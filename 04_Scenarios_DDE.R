## RR-TB Vietnam Model - 04: Scenarios
## Evaluate impact of INH-R diagnosis/treatment on RR-TB burden
## DDE model with renamed parameters

library(rio)
library(here)
library(deSolve)
library(reshape2)
library(data.table)
library(progress)

fitpts <- import(here("outputs", "pts", "fitpts.Rdata"))
cat("Loaded", nrow(fitpts), "plausible parameter sets\n")

PAF <- list(estimate = 0.951, lo = 0.856, hi = 0.981)
set.seed(251195)

interventions <- list(
  baseline       = list(name = "Baseline",              coverage = 0,                              efficacy = 0),
  inhr_low       = list(name = "INH-R Dx Low (20%)",    coverage = c(lo=0.15, hi=0.25), efficacy = 0.975),
  inhr_med       = list(name = "INH-R Dx Medium (50%)", coverage = c(lo=0.40, hi=0.60), efficacy = 0.975),
  inhr_high      = list(name = "INH-R Dx High (80%)",   coverage = c(lo=0.70, hi=0.90), efficacy = 0.975),
  inhr_universal = list(name = "Universal DST (100%)",   coverage = c(lo=0.95, hi=1.00), efficacy = 0.975)
)

Delta_s <- 0.5  # DS-TB treatment duration (years)

# DDE with intervention
dde_2strain_intervention <- function(parms, intervention = NULL, interv_start = 2025, end_time = 2035) {
  N <- 1e5; mu <- 0.013
  if (is.null(intervention) || identical(intervention$coverage, 0)) {
    interv_eff <- 1
  } else {
    cov <- runif(1, intervention$coverage['lo'], intervention$coverage['hi'])
    eff <- if (length(intervention$efficacy) == 1) intervention$efficacy else runif(1, intervention$efficacy['lo'], intervention$efficacy['hi'])
    interv_eff <- 1 - cov * eff * runif(1, PAF$lo, PAF$hi)
  }
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

      ie <- ifelse(time >= interv_start, interv_eff, 1)
      alpha_f_t <- alpha_f * ie; alpha_r_t <- alpha_r * ie

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
      dA_s <- p_a * I_s + xi * N_s + gamma * C_s + rho_s * (1 - alpha_r_t) * P_s - (zeta + sigma + mu) * A_s
      dC_s <- sigma * A_s + phi_s * (1 - alpha_f_t) * T_s - (gamma + theta_s + d_s + mu) * C_s
      dT_s <- theta_s * C_s - phi_s * T_s - completion_s - mu * T_s
      dP_s <- completion_s - (rho_s + tau + mu) * P_s
      Af <- phi_s * alpha_f_t * T_s; Ar <- rho_s * alpha_r_t * P_s
      dI_r <- lam_r * S - (eta + p_n + p_a + mu) * I_r
      dN_r <- p_n * I_r + zeta * A_r - (nu + xi + mu) * N_r
      dA_r <- p_a * I_r + xi * N_r + gamma * C_r + rho_r * P_r + Ar - (zeta + sigma + mu) * A_r
      dC_r <- sigma * A_r + phi_r * T_r + Af - (gamma + theta_r + d_r + mu) * C_r
      dT_r <- theta_r * C_r - phi_r * T_r - completion_r - mu * T_r
      dP_r <- completion_r - (rho_r + tau + mu) * P_r
      dA_r_acq <- Ar - (zeta + sigma + mu) * A_r_acq

      Pop <- S + I_s + N_s + A_s + C_s + T_s + P_s + I_r + N_r + A_r + C_r + T_r + P_r
      Inc_DS <- sigma * A_s; Inc_RR_acq <- Af + sigma * A_r_acq
      Inc_RR_trans <- sigma * (A_r - A_r_acq); Inc_RR <- sigma * A_r + Af
      return(list(c(dS, dI_s, dN_s, dA_s, dC_s, dT_s, dP_s,
                    dI_r, dN_r, dA_r, dC_r, dT_r, dP_r, dA_r_acq),
        Inc = ((Inc_DS + Inc_RR) / Pop) * 1e5, IncDS = (Inc_DS / Pop) * 1e5,
        IncRR = (Inc_RR / Pop) * 1e5, IncRR_acq = (Inc_RR_acq / Pop) * 1e5,
        IncRR_trans = (Inc_RR_trans / Pop) * 1e5,
        pAc = Inc_RR_acq / (Inc_RR + 1e-10), pRR = Inc_RR / (Inc_DS + Inc_RR + 1e-10),
        Prev = ((A_s + C_s + A_r + C_r) / Pop) * 1e5, PrevRR = ((A_r + C_r) / Pop) * 1e5,
        Mor = ((d_s * C_s + d_r * C_r) / Pop) * 1e5, MorRR = (d_r * C_r / Pop) * 1e5,
        Notif = ((theta_s * C_s + theta_r * C_r) / Pop) * 1e5))
    })
  }
  rr0 <- 0.025
  yini <- c(S = N - 500,
    I_s = 0, N_s = 0, A_s = 100*(1-rr0), C_s = 400*(1-rr0), T_s = 0, P_s = 0,
    I_r = 0, N_r = 0, A_r = 100*rr0, C_r = 400*rr0, T_r = 0, P_r = 0,
    A_r_acq = 0)
  deSolve::dede(yini, seq(1500, end_time, 1), des, parms, method = "lsoda")
}

cat("\n=== Running Scenarios ===\n")
all_results <- list()
for (sname in names(interventions)) {
  cat("  ", interventions[[sname]]$name, "\n")
  pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(fitpts))
  res <- lapply(1:nrow(fitpts), function(i) {
    p <- as.numeric(fitpts[i,]); names(p) <- names(fitpts)
    out <- as.data.frame(dde_2strain_intervention(p, interventions[[sname]], 2025, 2035))
    out$run <- i; out$scenario <- sname; out$scenario_name <- interventions[[sname]]$name; pb$tick(); out
  })
  all_results[[sname]] <- rbindlist(res)
}
results_all <- rbindlist(all_results)[time >= 2020]
export(results_all, here("outputs", "results", "scenarios_raw.Rdata"))

output_vars <- c("Inc", "IncRR", "IncRR_acq", "IncRR_trans", "pAc", "pRR", "Prev", "PrevRR", "Mor", "MorRR")
results_summary <- rbindlist(lapply(output_vars, function(v) {
  results_all[, .(lo = quantile(get(v), 0.025), q25 = quantile(get(v), 0.25), median = quantile(get(v), 0.50),
                  q75 = quantile(get(v), 0.75), hi = quantile(get(v), 0.975), var = v), by = .(scenario, scenario_name, time)]
}))
export(results_summary, here("outputs", "results", "scenarios_summary.Rdata"))

bl <- results_all[scenario == "baseline" & time == 2035, .(run, IncRR_base = IncRR, MorRR_base = MorRR, IncRR_acq_base = IncRR_acq)]
impact_all <- rbindlist(lapply(setdiff(names(interventions), "baseline"), function(scen) {
  sc <- results_all[scenario == scen & time == 2035, .(run, IncRR_int = IncRR, MorRR_int = MorRR, IncRR_acq_int = IncRR_acq, scenario_name)]
  imp <- merge(bl, sc, by = "run")
  imp[, `:=`(IncRR_reduction = (IncRR_base - IncRR_int) / IncRR_base * 100,
             IncRR_acq_reduction = (IncRR_acq_base - IncRR_acq_int) / (IncRR_acq_base + 1e-10) * 100,
             MorRR_reduction = (MorRR_base - MorRR_int) / (MorRR_base + 1e-10) * 100,
             IncRR_averted = IncRR_base - IncRR_int, MorRR_averted = MorRR_base - MorRR_int, scenario = scen)]
}))
impact_summary <- impact_all[, .(scenario_name = first(scenario_name),
  IncRR_reduction_median = median(IncRR_reduction), IncRR_reduction_lo = quantile(IncRR_reduction, 0.025),
  IncRR_reduction_hi = quantile(IncRR_reduction, 0.975),
  IncRR_averted_median = median(IncRR_averted), IncRR_averted_lo = quantile(IncRR_averted, 0.025),
  IncRR_averted_hi = quantile(IncRR_averted, 0.975)), by = scenario]
export(impact_all, here("outputs", "results", "impact_raw.Rdata"))
export(impact_summary, here("outputs", "results", "impact_summary.Rdata"))

cat("\n===============================================\nIMPACT SUMMARY (by 2035)\n===============================================\n")
bl_sum <- results_all[scenario == "baseline" & time == 2035, .(m = median(IncRR), lo = quantile(IncRR, 0.025), hi = quantile(IncRR, 0.975))]
cat(sprintf("Baseline RR-TB incidence: %.2f (%.2f-%.2f) per 100k\n\n", bl_sum$m, bl_sum$lo, bl_sum$hi))
for (i in 1:nrow(impact_summary)) {
  r <- impact_summary[i]
  cat(sprintf("%s:\n  RR-TB reduction: %.1f%% (%.1f-%.1f%%)\n  Cases averted: %.2f (%.2f-%.2f) per 100k\n\n",
              r$scenario_name, r$IncRR_reduction_median, r$IncRR_reduction_lo, r$IncRR_reduction_hi,
              r$IncRR_averted_median, r$IncRR_averted_lo, r$IncRR_averted_hi))
}
cat("=== 04_Scenarios.R completed ===\n")
