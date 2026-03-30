## RR-TB Vietnam Model - 02: Fit
## Two-strain DDE model with renamed parameters (Marc feedback)
## Treatment compartments use delay-differential equations

library(rio)
library(here)
library(hmer)
library(lhs)
library(deSolve)
library(reshape2)
library(data.table)
library(beepr)
library(tictoc)
library(progress)

# =============================================================================
# 1. Targets
# =============================================================================

targetsdb <- import(here("data", "fit", "targets.Rdata"))
targets_list <- import(here("data", "fit", "targets_list.Rdata"))

cat("=== Calibration targets loaded ===\n")
print(targetsdb[, c("var", "time", "lo", "hi")])

# =============================================================================
# 2. Parameters — renamed per Marc's feedback
# =============================================================================
# NAMING CONVENTION:
#   Variables: single letter + strain subscript (s=susceptible, r=resistant)
#     S       = Susceptible
#     I_s/I_r = Infection
#     N_s/N_r = Non-infectious TB
#     A_s/A_r = Asymptomatic TB
#     C_s/C_r = Clinical (symptomatic) TB
#     T_s/T_r = Treatment
#     P_s/P_r = Post-treatment (treated)
#     A_r_acq = Bookkeeping: acquired-resistance subset of A_r
#
#   Parameters: single Greek letter + subscript where needed
#     beta (β)        = transmission coefficient
#     kappa (κ)       = relative infectiousness of asymptomatic TB
#     f               = relative fitness of RR strain
#     eta (η)         = rate of clearance from infection → susceptible
#     p_n             = rate of progression: infection → non-infectious
#     p_a             = rate of progression: infection → asymptomatic
#     nu (ν)          = rate of recovery: non-infectious → susceptible
#     xi (ξ)          = rate of progression: non-infectious → asymptomatic
#     zeta (ζ)        = rate of regression: asymptomatic → non-infectious
#     sigma (σ)       = rate of progression: asymptomatic → symptomatic
#     gamma (γ)       = rate of regression: symptomatic → asymptomatic
#     theta (θ)       = rate of treatment initiation (time-varying)
#     phi (φ)         = rate of treatment failure (time-varying)
#     d               = TB-specific mortality rate (time-varying)
#     delta (δ)       = treatment completion rate (fixed delay for DDE)
#     tau (τ)         = rate of return: post-treatment → susceptible
#     rho (ρ)         = rate of relapse: post-treatment → asymptomatic
#     alpha_f (αf)    = prob. of acquiring RR during treatment failure
#     alpha_r (αr)    = prob. of acquiring RR during relapse
#     mu (μ)          = background mortality rate
#
# 13 compartments: S + 6 per strain (I, N, A, C, T, P)

parms <- c(
  beta = 12, kappa = 0.75, f = 0.65,
  eta = 1.83, p_n = 0.21, p_a = 0.07,
  nu = 0.16, xi = 0.25, zeta = 1.58, sigma = 0.77, gamma = 0.53,
  theta_s0 = 0.30, theta_s1 = 0.67,
  phi_s0 = 0.30, phi_s1 = 0.09,
  d_s0 = 0.33, d_s1 = 0.14,
  theta_r0 = 0.15, theta_r1 = 0.60,
  phi_r0 = 0.50, phi_r1 = 0.18,
  d_r0 = 0.42, d_r1 = 0.22,
  tau = 1.0, rho_s = 0.05, rho_r = 0.05,
  alpha_f = 0.05, alpha_r = 0.05
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

cat("Parameters:", length(parms), "| To calibrate:", length(ranges), "\n")

# =============================================================================
# 3. Two-Strain DDE Model (13 compartments)
#    Treatment compartments use fixed-delay DDEs
# =============================================================================

# Fixed treatment durations (years)
Delta_s <- 0.5   # DS-TB: 6 months

dde_2strain <- function(parms, end_time = 2024) {

  N <- 1e5; mu <- 0.013

  # Time-varying care cascade (piecewise linear 1999-2024)
  make_forcer <- function(p0, p1) {
    approxfun(c(1500, 1999, 2024), c(parms[p0], parms[p0], parms[p1]),
              method = "linear", rule = 2)
  }
  force_theta_s <- make_forcer('theta_s0', 'theta_s1')
  force_phi_s   <- make_forcer('phi_s0',   'phi_s1')
  force_d_s     <- make_forcer('d_s0',     'd_s1')
  force_theta_r <- make_forcer('theta_r0', 'theta_r1')
  force_phi_r   <- make_forcer('phi_r0',   'phi_r1')
  force_d_r     <- make_forcer('d_r0',     'd_r1')

  # RR-TB treatment duration: 18-month → 6-month (BPaL) over 2020-2024
  force_Delta_r <- approxfun(c(1500, 2020, 2024), c(1.5, 1.5, 0.5),
                              method = "linear", rule = 2)

  des <- function(time, state, parms) {
    with(as.list(c(state, parms)), {

      # Force of infection
      infectious_s <- kappa * A_s + C_s
      infectious_r <- f * (kappa * A_r + C_r)
      lam_s <- (beta / N) * infectious_s
      lam_r <- (beta / N) * infectious_r
      lam   <- lam_s + lam_r

      # Time-varying parameters
      theta_s <- force_theta_s(time)
      phi_s   <- force_phi_s(time)
      d_s     <- force_d_s(time)
      theta_r <- force_theta_r(time)
      phi_r   <- force_phi_r(time)
      d_r     <- force_d_r(time)
      Delta_r <- force_Delta_r(time)

      # --- DDE: delayed treatment completion ---
      # During early burn-in (before history exists), use rate-based approximation
      # Once history is available, use fixed-delay DDE
      t0 <- 1500  # model start time

      if (time - Delta_s > t0) {
        lag_s <- lagvalue(time - Delta_s)
        names(lag_s) <- names(state)
        completion_s <- max(0, theta_s * lag_s["C_s"] * exp(-mu * Delta_s))
      } else {
        # Rate-based fallback: 1/Delta_s per year
        completion_s <- (1 / Delta_s) * T_s
      }

      if (time - Delta_r > t0) {
        lag_r <- lagvalue(time - Delta_r)
        names(lag_r) <- names(state)
        completion_r <- max(0, theta_r * lag_r["C_r"] * exp(-mu * Delta_r))
      } else {
        completion_r <- (1 / Delta_r) * T_r
      }

      # --- Susceptible ---
      dS <- mu * N + d_s * C_s + d_r * C_r + eta * (I_s + I_r) + nu * (N_s + N_r) +
            tau * (P_s + P_r) - (lam + mu) * S

      # --- DS strain ---
      dI_s <- lam_s * S - (eta + p_n + p_a + mu) * I_s

      dN_s <- p_n * I_s + zeta * A_s - (nu + xi + mu) * N_s

      dA_s <- p_a * I_s + xi * N_s + gamma * C_s +
              rho_s * (1 - alpha_r) * P_s - (zeta + sigma + mu) * A_s

      dC_s <- sigma * A_s + phi_s * (1 - alpha_f) * T_s -
              (gamma + theta_s + d_s + mu) * C_s

      # Treatment: DDE — inflow from diagnosis, outflow from failure, delayed completion, death
      dT_s <- theta_s * C_s - phi_s * T_s - completion_s - mu * T_s

      dP_s <- completion_s - (rho_s + tau + mu) * P_s

      # Acquired resistance flows
      Af <- phi_s * alpha_f * T_s
      Ar <- rho_s * alpha_r * P_s

      # --- RR strain ---
      dI_r <- lam_r * S - (eta + p_n + p_a + mu) * I_r

      dN_r <- p_n * I_r + zeta * A_r - (nu + xi + mu) * N_r

      dA_r <- p_a * I_r + xi * N_r + gamma * C_r +
              rho_r * P_r + Ar - (zeta + sigma + mu) * A_r

      dC_r <- sigma * A_r + phi_r * T_r + Af -
              (gamma + theta_r + d_r + mu) * C_r

      # Treatment: DDE
      dT_r <- theta_r * C_r - phi_r * T_r - completion_r - mu * T_r

      dP_r <- completion_r - (rho_r + tau + mu) * P_r

      # --- Bookkeeping: acquired resistance subset of A_r ---
      dA_r_acq <- Ar - (zeta + sigma + mu) * A_r_acq

      # --- Outputs ---
      Pop <- S + I_s + N_s + A_s + C_s + T_s + P_s +
                 I_r + N_r + A_r + C_r + T_r + P_r

      Inc_DS  <- sigma * A_s
      Inc_RR_acq  <- Af + sigma * A_r_acq
      Inc_RR      <- sigma * A_r + Af
      Inc_total   <- Inc_DS + Inc_RR

      return(list(c(dS, dI_s, dN_s, dA_s, dC_s, dT_s, dP_s,
                        dI_r, dN_r, dA_r, dC_r, dT_r, dP_r,
                        dA_r_acq),
        Inc   = (Inc_total / Pop) * 1e5,
        IncRR = (Inc_RR / Pop) * 1e5,
        pAc   = Inc_RR_acq / (Inc_RR + 1e-10),
        pRR   = Inc_RR / (Inc_total + 1e-10),
        Prev  = ((A_s + C_s + A_r + C_r) / Pop) * 1e5,
        Mor   = ((d_s * C_s + d_r * C_r) / Pop) * 1e5,
        Notif = ((theta_s * C_s + theta_r * C_r) / Pop) * 1e5
      ))
    })
  }

  rr0 <- 0.025
  yini <- c(S = N - 500,
    I_s = 0, N_s = 0, A_s = 100*(1-rr0), C_s = 400*(1-rr0), T_s = 0, P_s = 0,
    I_r = 0, N_r = 0, A_r = 100*rr0,     C_r = 400*rr0,     T_r = 0, P_r = 0,
    A_r_acq = 0)

  times <- seq(1500, end_time, 1)

  deSolve::dede(yini, times, des, parms, method = "lsoda")
}

# Test
cat("\n=== Testing model ===\n")
test_df <- as.data.frame(dde_2strain(parms))
print(tail(test_df[, c("time", "Inc", "IncRR", "pRR", "pAc")], 5))

# =============================================================================
# 4. HMER Setup
# =============================================================================

hmer_res <- function(params, times, outputs) {
  all_res <- dde_2strain(params, end_time = max(times))
  actual_res <- all_res[all_res[, 'time'] %in% times, c('time', outputs)]
  shaped <- reshape2::melt(actual_res[, outputs])
  setNames(shaped$value, paste0(shaped$Var2, times))
}

targets <- setNames(
  lapply(1:nrow(targetsdb), function(i) c(targetsdb$lo[i], targetsdb$hi[i])),
  targetsdb$name
)

target_check <- function(wave, nwave, targets_df = targetsdb) {
  wave <- as.data.frame(wave)
  checks <- sapply(rownames(targets_df), function(targ) {
    val <- wave[, grep(colnames(wave), pattern = targ)]
    lim <- targets_df[targ, c('lo', 'hi')]
    val >= lim$lo & val <= lim$hi
  })
  list(runs = table(rowSums(checks)), targets = colSums(checks))
}

target_subset <- function(wave, sens = 0, targets_df = targetsdb) {
  wave <- as.data.frame(wave)
  checks <- sapply(rownames(targets_df), function(targ) {
    val <- wave[, grep(colnames(wave), pattern = targ)]
    lim <- targets_df[targ, c('lo', 'hi')]
    val >= lim$lo & val <= lim$hi
  })
  wave[rowSums(checks) >= (ncol(checks) - sens), ]
}

# =============================================================================
# 5. HMER Initial Wave
# =============================================================================

cat("\n=== Starting HMER Wave 1 ===\n")

n_samples <- 10 * length(ranges)
ini_LHS <- rbind(lhs::maximinLHS(n_samples, length(ranges)),
                  lhs::maximinLHS(n_samples, length(ranges)))

ini_pts <- setNames(
  data.frame(t(apply(ini_LHS, 1, function(x) {
    x * sapply(ranges, diff) + sapply(ranges, `[`, 1)
  }))),
  names(ranges)
)

cat("Generated", nrow(ini_pts), "initial parameter sets\n")

calib_times <- unique(targetsdb$time)
calib_outputs <- c("Inc", "IncRR", "pAc")

pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(ini_pts))
tic()
tmp <- lapply(seq_len(nrow(ini_pts)), function(i) {
  res <- tryCatch({
    r <- t(apply(ini_pts[i, ], 1, hmer_res, calib_times, calib_outputs))
    data.frame(r)[, names(targets)]
  }, error = function(e) {
    df <- data.frame(matrix(NA, 1, length(targets)))
    colnames(df) <- names(targets); df
  })
  pb$tick(); res
})
toc()

wave_res_1 <- do.call(rbind, tmp)
wave_1 <- cbind(ini_pts, wave_res_1)
wave_1 <- wave_1[complete.cases(wave_1), ]
wave_res_1 <- wave_res_1[complete.cases(wave_res_1), ]

cat("Successful runs:", nrow(wave_1), "/", nrow(ini_pts), "\n")

export(wave_res_1, here("outputs", "wave_res", "w001_waveres.Rdata"))
export(wave_1, here("outputs", "waves", "w001_wave.Rdata"))

checks_1 <- target_check(wave_res_1, 1)
cat("\nWave 1: Runs hitting all targets:", sum(checks_1$runs[as.character(length(targets))]), "\n")

# =============================================================================
# 6. Train Emulators
# =============================================================================

if (nrow(wave_1) >= 20) {
  cat("\n=== Training emulators ===\n")

  idx <- sample(nrow(wave_1), nrow(wave_1) %/% 2)
  wave_train <- wave_1[idx, ]
  wave_val <- wave_1[-idx, ]

  ems_1 <- emulator_from_data(wave_train, names(targets), ranges)

  for (j in seq_along(ems_1)) {
    for (iter in 1:10) {
      mc <- tryCatch(nrow(classification_diag(ems_1[[j]], targets, wave_val, plt = FALSE)),
                     error = function(e) 0)
      if (mc == 0) break
      ems_1[[j]] <- ems_1[[j]]$mult_sigma(1.1)
    }
  }

  bad <- sapply(seq_along(ems_1), function(j) {
    tryCatch(nrow(comparison_diag(ems_1[[j]], targets, wave_val, plt = FALSE)),
             error = function(e) 0) > nrow(wave_val) %/% 10
  })
  if (any(bad)) ems_1 <- ems_1[!bad]

  cat("Emulators after diagnostics:", length(ems_1), "\n")
  export(ems_1, here("outputs", "ems", "w001_ems.Rdata"))

  non_imp_pts_1 <- generate_new_design(ems_1, n_samples * 2, targets, verbose = TRUE)
  export(non_imp_pts_1, here("outputs", "pts", "w001_pts.Rdata"))
  cat("New points generated:", nrow(non_imp_pts_1), "\n")
} else {
  cat("Insufficient runs for emulator training.\n")
}

cat("\n=== 02_Fit.R completed ===\n")
beepr::beep(2)
