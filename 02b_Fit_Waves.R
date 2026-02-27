## RR-TB Vietnam Model - 02b: Fit Waves
## Run HMER waves until 1000 unique plausible parameter sets

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
  beta = c(6, 20), kappa = c(0.62, 1.00), fitness_rr = c(0.30, 1.00),
  omega_c = c(0.93, 3.30), omega_m = c(0.04, 0.23), omega_s = c(0.01, 0.10),
  mu_r = c(0.14, 0.23), mu_s = c(0.21, 0.28),
  sigma_m = c(1.24, 2.03), sigma_c = c(0.56, 0.94), gamma = c(0.46, 0.72),
  theta_ds_ini = c(0.00, 0.57), theta_ds_fin = c(0.57, 0.77),
  phi_ds_ini = c(0.11, 1.00), phi_ds_fin = c(0.07, 0.11),
  mutb_ds_ini = c(0.28, 0.38), mutb_ds_fin = c(0.00, 0.28),
  theta_rr_ini = c(0.00, 0.30), theta_rr_fin = c(0.50, 0.70),
  phi_rr_ini = c(0.29, 1.00), phi_rr_fin = c(0.07, 0.29),
  mutb_rr_ini = c(0.35, 0.50), mutb_rr_fin = c(0.15, 0.30),
  tau = c(0.50, 1.50), relapse_ds = c(0.01, 0.30), relapse_rr = c(0.01, 0.30),
  acq_fail = c(0.01, 0.20), acq_relapse = c(0.01, 0.20)
)

# =============================================================================
# 2. ODE Model (13 compartments, no CLE/REC)
# =============================================================================

ode_2strain <- function(parms, end_time = 2024) {
  N <- 1e5; mu <- 0.013; delta_ds <- 2
  make_forcer <- function(p_ini, p_fin) {
    approxfun(c(1900, 1999, 2024), c(parms[p_ini], parms[p_ini], parms[p_fin]), method = "linear", rule = 2)
  }
  force_theta_ds <- make_forcer('theta_ds_ini','theta_ds_fin')
  force_phi_ds   <- make_forcer('phi_ds_ini','phi_ds_fin')
  force_mutb_ds  <- make_forcer('mutb_ds_ini','mutb_ds_fin')
  force_theta_rr <- make_forcer('theta_rr_ini','theta_rr_fin')
  force_phi_rr   <- make_forcer('phi_rr_ini','phi_rr_fin')
  force_mutb_rr  <- make_forcer('mutb_rr_ini','mutb_rr_fin')
  force_delta_rr <- approxfun(c(1900, 2020, 2024), c(0.67, 0.67, 2.00), method = "linear", rule = 2)
  
  des <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      infectious_ds <- kappa*SUB_DS + CLN_DS
      infectious_rr <- fitness_rr*(kappa*SUB_RR + CLN_RR)
      lambda_ds <- beta/N*infectious_ds; lambda_rr <- beta/N*infectious_rr
      lambda_total <- lambda_ds + lambda_rr
      theta_ds_t <- force_theta_ds(time); phi_ds_t <- force_phi_ds(time); mutb_ds_t <- force_mutb_ds(time)
      theta_rr_t <- force_theta_rr(time); phi_rr_t <- force_phi_rr(time); mutb_rr_t <- force_mutb_rr(time)
      delta_rr_t <- force_delta_rr(time)
      
      dSUS <- mu*N + omega_c*(INF_DS+INF_RR) + mu_r*(MIN_DS+MIN_RR) + tau*(TRE_DS+TRE_RR) - (lambda_total+mu)*SUS
      dINF_DS <- lambda_ds*SUS - (omega_c+omega_m+omega_s+mu)*INF_DS
      dMIN_DS <- omega_m*INF_DS + sigma_m*SUB_DS - (mu_r+mu_s+mu)*MIN_DS
      dSUB_DS <- omega_s*INF_DS + mu_s*MIN_DS + gamma*CLN_DS + relapse_ds*(1-acq_relapse)*TRE_DS - (sigma_m+sigma_c+mu)*SUB_DS
      dCLN_DS <- sigma_c*SUB_DS - (gamma+theta_ds_t+mutb_ds_t+mu)*CLN_DS + phi_ds_t*(1-acq_fail)*TXT_DS
      dTXT_DS <- theta_ds_t*CLN_DS - (phi_ds_t+delta_ds+mu)*TXT_DS
      dTRE_DS <- delta_ds*TXT_DS - (relapse_ds+tau+mu)*TRE_DS
      acq_f <- phi_ds_t*acq_fail*TXT_DS; acq_r <- relapse_ds*acq_relapse*TRE_DS
      dINF_RR <- lambda_rr*SUS - (omega_c+omega_m+omega_s+mu)*INF_RR
      dMIN_RR <- omega_m*INF_RR + sigma_m*SUB_RR - (mu_r+mu_s+mu)*MIN_RR
      dSUB_RR <- omega_s*INF_RR + mu_s*MIN_RR + gamma*CLN_RR + relapse_rr*TRE_RR + acq_r - (sigma_m+sigma_c+mu)*SUB_RR
      dCLN_RR <- sigma_c*SUB_RR - (gamma+theta_rr_t+mutb_rr_t+mu)*CLN_RR + phi_rr_t*TXT_RR + acq_f
      dTXT_RR <- theta_rr_t*CLN_RR - (phi_rr_t+delta_rr_t+mu)*TXT_RR
      dTRE_RR <- delta_rr_t*TXT_RR - (relapse_rr+tau+mu)*TRE_RR
      dSUB_RR_acq <- acq_r - (sigma_m+sigma_c+mu)*SUB_RR_acq
      
      Pop <- SUS + INF_DS+MIN_DS+SUB_DS+CLN_DS+TXT_DS+TRE_DS + INF_RR+MIN_RR+SUB_RR+CLN_RR+TXT_RR+TRE_RR
      Inc_DS <- sigma_c*SUB_DS; Inc_RR_acq <- acq_f+sigma_c*SUB_RR_acq; Inc_RR <- sigma_c*SUB_RR+acq_f
      
      return(list(c(dSUS, dINF_DS,dMIN_DS,dSUB_DS,dCLN_DS,dTXT_DS,dTRE_DS,
                    dINF_RR,dMIN_RR,dSUB_RR,dCLN_RR,dTXT_RR,dTRE_RR,dSUB_RR_acq),
        Inc = ((Inc_DS+Inc_RR)/Pop)*1e5, IncRR = (Inc_RR/Pop)*1e5,
        pAc = Inc_RR_acq/(Inc_RR+1e-10), pRR = Inc_RR/(Inc_DS+Inc_RR+1e-10),
        Prev = ((SUB_DS+CLN_DS+SUB_RR+CLN_RR)/Pop)*1e5,
        Mor = ((mutb_ds_t*CLN_DS+mutb_rr_t*CLN_RR)/Pop)*1e5))
    })
  }
  rr0 <- 0.025
  yini <- c(SUS=N-500, INF_DS=0, MIN_DS=0, SUB_DS=100*(1-rr0), CLN_DS=400*(1-rr0), TXT_DS=0, TRE_DS=0,
            INF_RR=0, MIN_RR=0, SUB_RR=100*rr0, CLN_RR=400*rr0, TXT_RR=0, TRE_RR=0, SUB_RR_acq=0)
  deSolve::ode(yini, seq(1900, end_time, 1), des, parms, method = "lsoda")
}

hmer_res <- function(params, times, outputs) {
  all_res <- ode_2strain(params, end_time = max(times))
  actual_res <- all_res[all_res[, 'time'] %in% times, c('time', outputs)]
  shaped <- reshape2::melt(actual_res[, outputs])
  setNames(shaped$value, paste0(shaped$Var2, times))
}

run_model_parallel <- function(pts, calib_times, calib_outputs, targets) {
  clusterExport(cl, c("ode_2strain", "hmer_res", "calib_times", "calib_outputs", "targets"), envir = environment())
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
for (f in wave_files) wave[[get_wave_num(f)]] <- import(f)
for (f in ems_files) ems[[get_wave_num(f)]] <- import(f)
for (f in pts_files) non_imp_pts[[get_wave_num(f)]] <- import(f)

last_wave <- max(which(!vapply(non_imp_pts, is.null, logical(1))), 1)
cat("Continuing from wave:", last_wave + 1, "\n")

all_plausible <- do.call(rbind, lapply(seq_along(wave), function(i) {
  if (is.null(wave[[i]]) || nrow(wave[[i]]) == 0) return(NULL)
  wc <- target_subset(wave[[i]], sens = 0)
  if (nrow(wc) > 0) wc[, names(ranges)] else NULL
}))
all_plausible <- unique(all_plausible)
cat("Existing plausible points:", nrow(all_plausible), "\n")

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
  if (nrow(wc) > 0) all_plausible <- unique(rbind(all_plausible, wc[, names(ranges)]))
  n_unique <- nrow(all_plausible)
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
      non_imp_pts[[w]][[col]] <- pmax(ranges[[col]][1], pmin(ranges[[col]][2], non_imp_pts[[w]][[col]] + rnorm(nrow(non_imp_pts[[w]]), 0, 0.05*rw)))
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
