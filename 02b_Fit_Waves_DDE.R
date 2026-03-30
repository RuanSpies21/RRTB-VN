## RR-TB Vietnam Model - 02b: Fit Waves
## Run HMER waves until 1000 unique plausible parameter sets
## DDE model with renamed parameters

library(parallel)
library(doParallel)
library(foreach)
library(rio)
library(here)
library(hmer)
library(lhs)
library(deSolve)
library(reshape2)
library(data.table)
library(beepr)
library(tictoc)

# =============================================================================
# 1. Setup
# =============================================================================

n_cores <- max(1, detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

targetsdb <- import(here("data", "fit", "targets.Rdata"))
targets <- setNames(
  lapply(1:nrow(targetsdb), function(i) c(targetsdb$lo[i], targetsdb$hi[i])),
  targetsdb$name
)

ranges <- list(
  beta = c(6, 20), kappa = c(0.62, 1.00), f = c(0.30, 1.00),
  eta = c(0.93, 3.30), p_n = c(0.04, 0.23), p_a = c(0.01, 0.10),
  nu = c(0.14, 0.23), xi = c(0.21, 0.28),
  zeta = c(1.24, 2.03), sigma = c(0.56, 0.94), gamma = c(0.46, 0.72),
  theta_s0 = c(0.00, 0.57), theta_s1 = c(0.57, 0.77),
  phi_s0 = c(0.11, 1.00), phi_s1 = c(0.07, 0.11),
  d_s0 = c(0.28, 0.38), d_s1 = c(0.00, 0.28),
  theta_r0 = c(0.00, 0.30), theta_r1 = c(0.50, 0.70),
  phi_r0 = c(0.29, 1.00), phi_r1 = c(0.07, 0.29),
  d_r0 = c(0.35, 0.50), d_r1 = c(0.15, 0.30),
  tau = c(0.50, 1.50), rho_s = c(0.01, 0.30), rho_r = c(0.01, 0.30),
  alpha_f = c(0.01, 0.20), alpha_r = c(0.01, 0.20)
)

# =============================================================================
# 2. DDE Model
# =============================================================================

Delta_s <- 0.5  # DS-TB treatment duration (years)

dde_2strain <- function(parms, end_time = 2024) {
  N <- 1e5; mu <- 0.013
  make_forcer <- function(p0, p1) {
    approxfun(c(1500, 1999, 2024), c(parms[p0], parms[p0], parms[p1]), method = "linear", rule = 2)
  }
  force_theta_s <- make_forcer('theta_s0', 'theta_s1')
  force_phi_s   <- make_forcer('phi_s0',   'phi_s1')
  force_d_s     <- make_forcer('d_s0',     'd_s1')
  force_theta_r <- make_forcer('theta_r0', 'theta_r1')
  force_phi_r   <- make_forcer('phi_r0',   'phi_r1')
  force_d_r     <- make_forcer('d_r0',     'd_r1')
  force_Delta_r <- approxfun(c(1500, 2020, 2024), c(1.5, 1.5, 0.5), method = "linear", rule = 2)

  des <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      infectious_s <- kappa * A_s + C_s
      infectious_r <- f * (kappa * A_r + C_r)
      lam_s <- (beta / N) * infectious_s
      lam_r <- (beta / N) * infectious_r
      lam <- lam_s + lam_r
      theta_s <- force_theta_s(time); phi_s <- force_phi_s(time); d_s <- force_d_s(time)
      theta_r <- force_theta_r(time); phi_r <- force_phi_r(time); d_r <- force_d_r(time)
      Delta_r <- force_Delta_r(time)

      # DDE: delayed treatment completion
      t0 <- 1500
      if (time - Delta_s > t0) {
        lag_s <- lagvalue(time - Delta_s); names(lag_s) <- names(state)
        completion_s <- max(0, theta_s * lag_s["C_s"] * exp(-mu * Delta_s))
      } else { completion_s <- (1 / Delta_s) * T_s }
      if (time - Delta_r > t0) {
        lag_r <- lagvalue(time - Delta_r); names(lag_r) <- names(state)
        completion_r <- max(0, theta_r * lag_r["C_r"] * exp(-mu * Delta_r))
      } else { completion_r <- (1 / Delta_r) * T_r }

      dS <- mu * N + d_s * C_s + d_r * C_r + eta * (I_s + I_r) + nu * (N_s + N_r) + tau * (P_s + P_r) - (lam + mu) * S
      dI_s <- lam_s * S - (eta + p_n + p_a + mu) * I_s
      dN_s <- p_n * I_s + zeta * A_s - (nu + xi + mu) * N_s
      dA_s <- p_a * I_s + xi * N_s + gamma * C_s + rho_s * (1 - alpha_r) * P_s - (zeta + sigma + mu) * A_s
      dC_s <- sigma * A_s + phi_s * (1 - alpha_f) * T_s - (gamma + theta_s + d_s + mu) * C_s
      dT_s <- theta_s * C_s - phi_s * T_s - completion_s - mu * T_s
      dP_s <- completion_s - (rho_s + tau + mu) * P_s
      Af <- phi_s * alpha_f * T_s; Ar <- rho_s * alpha_r * P_s
      dI_r <- lam_r * S - (eta + p_n + p_a + mu) * I_r
      dN_r <- p_n * I_r + zeta * A_r - (nu + xi + mu) * N_r
      dA_r <- p_a * I_r + xi * N_r + gamma * C_r + rho_r * P_r + Ar - (zeta + sigma + mu) * A_r
      dC_r <- sigma * A_r + phi_r * T_r + Af - (gamma + theta_r + d_r + mu) * C_r
      dT_r <- theta_r * C_r - phi_r * T_r - completion_r - mu * T_r
      dP_r <- completion_r - (rho_r + tau + mu) * P_r
      dA_r_acq <- Ar - (zeta + sigma + mu) * A_r_acq

      Pop <- S + I_s + N_s + A_s + C_s + T_s + P_s + I_r + N_r + A_r + C_r + T_r + P_r
      Inc_DS <- sigma * A_s; Inc_RR_acq <- Af + sigma * A_r_acq; Inc_RR <- sigma * A_r + Af
      return(list(c(dS, dI_s, dN_s, dA_s, dC_s, dT_s, dP_s,
                    dI_r, dN_r, dA_r, dC_r, dT_r, dP_r, dA_r_acq),
        Inc = ((Inc_DS + Inc_RR) / Pop) * 1e5, IncRR = (Inc_RR / Pop) * 1e5,
        pAc = Inc_RR_acq / (Inc_RR + 1e-10), pRR = Inc_RR / (Inc_DS + Inc_RR + 1e-10)))
    })
  }
  rr0 <- 0.025
  yini <- c(S = N - 500,
    I_s = 0, N_s = 0, A_s = 100*(1-rr0), C_s = 400*(1-rr0), T_s = 0, P_s = 0,
    I_r = 0, N_r = 0, A_r = 100*rr0, C_r = 400*rr0, T_r = 0, P_r = 0,
    A_r_acq = 0)
  deSolve::dede(yini, seq(1500, end_time, 1), des, parms, method = "lsoda")
}

hmer_res <- function(params, times, outputs) {
  all_res <- dde_2strain(params, end_time = max(times))
  actual_res <- all_res[all_res[, 'time'] %in% times, c('time', outputs)]
  shaped <- reshape2::melt(actual_res[, outputs])
  setNames(shaped$value, paste0(shaped$Var2, times))
}

run_model_parallel <- function(pts, calib_times, calib_outputs, targets) {
  clusterExport(cl, c("dde_2strain", "hmer_res", "calib_times", "calib_outputs", "targets", "Delta_s"),
                envir = environment())
  clusterEvalQ(cl, { library(deSolve); library(reshape2) })
  do.call(rbind, parLapply(cl, 1:nrow(pts), function(i) {
    tryCatch({
      p <- as.numeric(pts[i, ]); names(p) <- colnames(pts)
      data.frame(t(hmer_res(p, calib_times, calib_outputs)))[, names(targets)]
    }, error = function(e) { df <- data.frame(matrix(NA, 1, length(targets))); colnames(df) <- names(targets); df })
  }))
}

target_subset <- function(wave, sens = 0) {
  wave <- as.data.frame(wave)
  checks <- sapply(rownames(targetsdb), function(targ) {
    val <- wave[, grep(colnames(wave), pattern = targ)]
    lim <- targetsdb[targ, c('lo', 'hi')]
    val >= lim$lo & val <= lim$hi
  })
  wave[rowSums(checks) >= (ncol(checks) - sens), ]
}

# =============================================================================
# 3. Load Previous Waves & Run
# =============================================================================

get_wave_num <- function(x) as.integer(sub("^.*w0*([0-9]+)_.*$", "\\1", basename(x)))
wave_files <- sort(list.files(here("outputs", "waves"), "^w[0-9]+_wave\\.Rdata$", full.names = TRUE))
ems_files  <- sort(list.files(here("outputs", "ems"),   "^w[0-9]+_ems\\.Rdata$",  full.names = TRUE))
pts_files  <- sort(list.files(here("outputs", "pts"),   "^w[0-9]+_pts\\.Rdata$",  full.names = TRUE))

max_wave <- max(c(1, sapply(c(wave_files, ems_files, pts_files), get_wave_num)), na.rm = TRUE)
wave <- vector("list", max_wave); ems <- vector("list", max_wave); non_imp_pts <- vector("list", max_wave)
for (f_ in wave_files) wave[[get_wave_num(f_)]] <- import(f_)
for (f_ in ems_files) ems[[get_wave_num(f_)]] <- import(f_)
for (f_ in pts_files) non_imp_pts[[get_wave_num(f_)]] <- import(f_)

last_wave <- max(which(!vapply(non_imp_pts, is.null, logical(1))), 1)
cat("Continuing from wave:", last_wave + 1, "\n")

all_plausible <- do.call(rbind, lapply(seq_along(wave), function(i) {
  if (is.null(wave[[i]]) || nrow(wave[[i]]) == 0) return(NULL)
  wc <- target_subset(wave[[i]], sens = 0)
  if (nrow(wc) > 0) wc[, names(ranges)] else NULL
}))
if (!is.null(all_plausible)) all_plausible <- unique(all_plausible)
n_existing <- if (is.null(all_plausible)) 0L else nrow(all_plausible)
cat("Existing plausible points:", n_existing, "\n")

n_samples <- 20 * length(ranges)
target_plausible <- 1000
calib_times <- unique(targetsdb$time)
calib_outputs <- c("Inc", "IncRR", "pAc")

cat("\n=== Running waves until", target_plausible, "plausible points ===\n")

for (w in (last_wave + 1):(last_wave + 100)) {
  cat("\n=== Wave", w, "===\n"); tic()
  pts_to_run <- non_imp_pts[[w - 1]]
  if (is.null(pts_to_run)) stop("Missing points for wave ", w - 1)

  wave_res_w <- run_model_parallel(pts_to_run, calib_times, calib_outputs, targets)
  wave[[w]] <- cbind(pts_to_run, wave_res_w)[complete.cases(cbind(pts_to_run, wave_res_w)), ]
  cat("Successful:", nrow(wave[[w]]), "/", nrow(pts_to_run), "\n")
  export(wave[[w]], here("outputs", "waves", sprintf("w%03d_wave.Rdata", w)))

  wc <- target_subset(wave[[w]], sens = 0)
  cat("Hits all targets:", nrow(wc), "\n")
  if (nrow(wc) > 0) {
    new_pts <- wc[, names(ranges)]
    all_plausible <- if (is.null(all_plausible)) new_pts else unique(rbind(all_plausible, new_pts))
  }
  n_unique <- if (is.null(all_plausible)) 0L else nrow(all_plausible)
  cat("Cumulative plausible:", n_unique, "/", target_plausible, "\n")

  if (n_unique >= target_plausible) {
    cat("\n*** SUCCESS:", n_unique, "plausible sets found ***\n")
    export(all_plausible, here("outputs", "pts", "fitpts.Rdata"))
    toc(); beepr::beep(2); break
  }

  if (nrow(wave[[w]]) >= 20) {
    idx <- sample(nrow(wave[[w]]), nrow(wave[[w]]) %/% 2)
    wave_val <- wave[[w]][-idx, ]
    ems[[w]] <- emulator_from_data(wave[[w]][idx, ], names(targets), ranges, check.ranges = TRUE)
    for (j in seq_along(ems[[w]])) {
      for (it in 1:10) {
        mc <- tryCatch(nrow(classification_diag(ems[[w]][[j]], targets, wave_val, plt = FALSE)), error = function(e) 0)
        if (mc == 0) break; ems[[w]][[j]] <- ems[[w]][[j]]$mult_sigma(1.1)
      }
    }
    bad <- sapply(seq_along(ems[[w]]), function(j) {
      tryCatch(nrow(comparison_diag(ems[[w]][[j]], targets, wave_val, plt = FALSE)), error = function(e) 0) > nrow(wave_val) %/% 10
    })
    if (any(bad)) ems[[w]] <- ems[[w]][!bad]
    export(ems[[w]], here("outputs", "ems", sprintf("w%03d_ems.Rdata", w)))
    all_ems <- unlist(ems[!sapply(ems, is.null)], recursive = FALSE)
    non_imp_pts[[w]] <- generate_new_design(all_ems, n_samples * 2, targets, verbose = TRUE)
    if (nrow(wc) > 0) non_imp_pts[[w]] <- unique(rbind(non_imp_pts[[w]], wc[, names(ranges)]))
    export(non_imp_pts[[w]], here("outputs", "pts", sprintf("w%03d_pts.Rdata", w)))
  } else {
    non_imp_pts[[w]] <- non_imp_pts[[w-1]]
    for (col in names(ranges)) {
      rw <- diff(ranges[[col]])
      non_imp_pts[[w]][[col]] <- pmax(ranges[[col]][1], pmin(ranges[[col]][2],
        non_imp_pts[[w]][[col]] + rnorm(nrow(non_imp_pts[[w]]), 0, 0.05*rw)))
    }
    export(non_imp_pts[[w]], here("outputs", "pts", sprintf("w%03d_pts.Rdata", w)))
  }
  toc()
}

stopCluster(cl)
all_plausible <- unique(all_plausible)
export(all_plausible, here("outputs", "pts", "fitpts.Rdata"))
cat("\nFinal plausible sets:", nrow(all_plausible), "\n=== 02b_Fit_Waves.R completed ===\n")
beepr::beep(2)
