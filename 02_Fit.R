## RR-TB Vietnam Model - 02: Fit
## Define two-strain ODE model and run initial HMER calibration wave

library(rio)
library(here)
library(hmer)
library(lhs)
library(deSolve)
library(reshape2)
library(data.table)
library(beepr)
library(progress)
library(tictoc)

# =============================================================================
# 1. Load Data
# =============================================================================

targetsdb <- import(here("data", "fit", "targets.Rdata"))
targets_list <- import(here("data", "fit", "targets_list.Rdata"))

cat("=== Calibration targets loaded ===\n")
print(targetsdb[, c("var", "time", "lo", "hi")])

# =============================================================================
# 2. Parameters (aligned with model diagram)
# =============================================================================
# 13 compartments: SUS + 6 per strain (INF, MIN, SUB, CLN, TXT, TRE)
# No CLE/REC compartments — ωc and μr flow directly back to SUS

parms <- c(
  beta = 12, kappa = 0.75, fitness_rr = 0.65,
  omega_c = 1.83, omega_m = 0.21, omega_s = 0.07,
  mu_r = 0.16, mu_s = 0.25, sigma_m = 1.58, sigma_c = 0.77, gamma = 0.53,
  theta_ds_ini = 0.30, theta_ds_fin = 0.67,
  phi_ds_ini = 0.30, phi_ds_fin = 0.09,
  mutb_ds_ini = 0.33, mutb_ds_fin = 0.14,
  theta_rr_ini = 0.15, theta_rr_fin = 0.60,
  phi_rr_ini = 0.50, phi_rr_fin = 0.18,
  mutb_rr_ini = 0.42, mutb_rr_fin = 0.22,
  tau = 1.0, relapse_ds = 0.05, relapse_rr = 0.05,
  acq_fail = 0.05, acq_relapse = 0.05
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

cat("Parameters:", length(parms), "| To calibrate:", length(ranges), "\n")

# =============================================================================
# 3. Two-Strain ODE Model (13 compartments)
# =============================================================================

ode_2strain <- function(parms, end_time = 2024) {
  
  N <- 1e5; mu <- 0.013; delta_ds <- 2
  
  make_forcer <- function(p_ini, p_fin) {
    approxfun(c(1500, 1999, 2024), c(parms[p_ini], parms[p_ini], parms[p_fin]),
              method = "linear", rule = 2)
  }
  force_theta_ds <- make_forcer('theta_ds_ini', 'theta_ds_fin')
  force_phi_ds   <- make_forcer('phi_ds_ini',   'phi_ds_fin')
  force_mutb_ds  <- make_forcer('mutb_ds_ini',  'mutb_ds_fin')
  force_theta_rr <- make_forcer('theta_rr_ini', 'theta_rr_fin')
  force_phi_rr   <- make_forcer('phi_rr_ini',   'phi_rr_fin')
  force_mutb_rr  <- make_forcer('mutb_rr_ini',  'mutb_rr_fin')
  force_delta_rr <- approxfun(c(1500, 2020, 2024), c(0.67, 0.67, 2.00),
                               method = "linear", rule = 2)
  
  des <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      
      infectious_ds <- (kappa * SUB_DS) + CLN_DS
      infectious_rr <- fitness_rr * ((kappa * SUB_RR) + CLN_RR)
      lambda_ds <- (beta / N) * infectious_ds
      lambda_rr <- (beta / N) * infectious_rr
      lambda_total <- lambda_ds + lambda_rr
      
      theta_ds_t <- force_theta_ds(time)
      phi_ds_t   <- force_phi_ds(time)
      mutb_ds_t  <- force_mutb_ds(time)
      theta_rr_t <- force_theta_rr(time)
      phi_rr_t   <- force_phi_rr(time)
      mutb_rr_t  <- force_mutb_rr(time)
      delta_rr_t <- force_delta_rr(time)
      
      # Susceptible: inflows from births, TB deaths, clearance, recovery, treatment completion
      dSUS <- (mu * N) + mutb_ds_t*CLN_DS + mutb_rr_t*CLN_RR + omega_c*(INF_DS + INF_RR) + mu_r*(MIN_DS + MIN_RR) +
              tau*(TRE_DS + TRE_RR) - (lambda_total + mu) * SUS
      
      # DS strain
      dINF_DS <- lambda_ds * SUS -
                 (omega_c + omega_m + omega_s + mu) * INF_DS
      dMIN_DS <- (omega_m * INF_DS) + (sigma_m * SUB_DS) - (mu_r + mu_s + mu) * MIN_DS
      dSUB_DS <- (omega_s * INF_DS) + (mu_s * MIN_DS) + (gamma * CLN_DS) +
                 (relapse_ds * (1 - acq_relapse) * TRE_DS) -
                 (sigma_m + sigma_c + mu) * SUB_DS
      dCLN_DS <- (sigma_c * SUB_DS) - (gamma + theta_ds_t + mutb_ds_t + mu) * CLN_DS +
                 (phi_ds_t * (1 - acq_fail) * TXT_DS)
      dTXT_DS <- (theta_ds_t * CLN_DS) - (phi_ds_t + delta_ds + mu) * TXT_DS
      dTRE_DS <- (delta_ds * TXT_DS) - (relapse_ds + tau + mu) * TRE_DS
      
      # Acquired resistance flows
      acquired_from_failure <- phi_ds_t * acq_fail * TXT_DS
      acquired_from_relapse <- relapse_ds * acq_relapse * TRE_DS
      
      # RR strain
      dINF_RR <- lambda_rr * SUS -
                 (omega_c + omega_m + omega_s + mu) * INF_RR
      dMIN_RR <- (omega_m * INF_RR) + (sigma_m * SUB_RR) - (mu_r + mu_s + mu) * MIN_RR
      dSUB_RR <- (omega_s * INF_RR) + (mu_s * MIN_RR) + (gamma * CLN_RR) +
                 (relapse_rr * TRE_RR) + acquired_from_relapse -
                 (sigma_m + sigma_c + mu) * SUB_RR
      dCLN_RR <- (sigma_c * SUB_RR) - (gamma + theta_rr_t + mutb_rr_t + mu) * CLN_RR +
                 (phi_rr_t * TXT_RR) +
                 acquired_from_failure
      dTXT_RR <- (theta_rr_t * CLN_RR) - (phi_rr_t + delta_rr_t + mu) * TXT_RR
      dTRE_RR <- (delta_rr_t * TXT_RR) - (relapse_rr + tau + mu) * TRE_RR
      
      # Tracking state: acquired-via-relapse portion of SUB_RR
      # Same outflow rates as SUB_RR, only inflow from acq_r
      # Does NOT affect transmission — SUB_RR remains the full pool
      dSUB_RR_acq <- acquired_from_relapse - (sigma_m + sigma_c + mu) * SUB_RR_acq
      
      # Outputs
      Pop <- SUS + INF_DS + MIN_DS + SUB_DS + CLN_DS + TXT_DS + TRE_DS +
                   INF_RR + MIN_RR + SUB_RR + CLN_RR + TXT_RR + TRE_RR
      Inc_DS <- sigma_c * SUB_DS
      Inc_RR_acq <- acquired_from_failure + (sigma_c * SUB_RR_acq)
      Inc_RR_trans <- sigma_c * (SUB_RR - SUB_RR_acq)
      Inc_RR <- sigma_c * SUB_RR + acquired_from_failure
      Inc_total <- Inc_DS + Inc_RR
      
      return(list(c(dSUS, dINF_DS, dMIN_DS, dSUB_DS, dCLN_DS, dTXT_DS, dTRE_DS,
                    dINF_RR, dMIN_RR, dSUB_RR, dCLN_RR, dTXT_RR, dTRE_RR,
                    dSUB_RR_acq),
        Inc = (Inc_total / Pop) * 1e5,
        IncRR = (Inc_RR / Pop) * 1e5,
        pAc = Inc_RR_acq / (Inc_RR + 1e-10),
        pRR = Inc_RR / (Inc_total + 1e-10),
        Prev = ((SUB_DS + CLN_DS + SUB_RR + CLN_RR) / Pop) * 1e5,
        Mor = ((mutb_ds_t * CLN_DS + mutb_rr_t * CLN_RR) / Pop) * 1e5,
        Notif = ((theta_ds_t * CLN_DS + theta_rr_t * CLN_RR) / Pop) * 1e5
      ))
    })
  }
  
  rr0 <- 0.025
  yini <- c(SUS = N - 500,
    INF_DS = 0, MIN_DS = 0, SUB_DS = 100*(1-rr0), CLN_DS = 400*(1-rr0), TXT_DS = 0, TRE_DS = 0,
    INF_RR = 0, MIN_RR = 0, SUB_RR = 100*rr0, CLN_RR = 400*rr0, TXT_RR = 0, TRE_RR = 0,
    SUB_RR_acq = 0)
  
  deSolve::ode(yini, seq(1500, end_time, 1), des, parms, method = "lsoda")
}

# Test
cat("\n=== Testing model ===\n")
test_df <- as.data.frame(ode_2strain(parms))
print(tail(test_df[, c("time", "Inc", "IncRR", "pRR", "pAc")], 5))

# =============================================================================
# 4. HMER Setup
# =============================================================================

hmer_res <- function(params, times, outputs) {
  all_res <- ode_2strain(params, end_time = max(times))
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
