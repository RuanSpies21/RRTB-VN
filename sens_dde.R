## RR-TB Vietnam Model - 04b: Global Sensitivity Analysis (DDE)
## Runs universal INH-R DST scenario across 9 PAF settings
## (3 INH-R prevalence levels × 3 RR levels)
## DDE model with renamed parameters

library(rio)
library(here)
library(deSolve)
library(data.table)
library(progress)
library(ggplot2)
library(patchwork)

fitpts <- as.data.table(import(here("outputs", "pts", "fitpts.Rdata")))
cat("Loaded", nrow(fitpts), "plausible parameter sets\n")

set.seed(251195)

# =============================================================================
# 1. Define PAF scenarios (Levin's formula with Monte Carlo CIs)
# =============================================================================

levin_paf <- function(pe, rr) {
  (pe * (rr - 1)) / (1 + pe * (rr - 1))
}

prevalence_scenarios <- list(
  low    = list(label = "5%",  estimate = 0.05, lo = 0.03, hi = 0.08),
  medium = list(label = "10%", estimate = 0.10, lo = 0.07, hi = 0.14),
  high   = list(label = "20%", estimate = 0.20, lo = 0.15, hi = 0.25)
)

# Updated RR: 105 (33-310) from weighted pseudo-cohort
rr_scenarios <- list(
  observed = list(label = "Observed",      estimate = 105,   lo = 33,    hi = 310),
  twofold  = list(label = "2-fold lower",  estimate = 52.5,  lo = 16.5,  hi = 155),
  fivefold = list(label = "5-fold lower",  estimate = 21,    lo = 6.6,   hi = 62)
)

# Pre-compute PAF ranges for each combination via Monte Carlo
n_mc <- 50000
paf_settings <- list()

for (p in names(prevalence_scenarios)) {
  for (r in names(rr_scenarios)) {
    pe_samp <- runif(n_mc, prevalence_scenarios[[p]]$lo, prevalence_scenarios[[p]]$hi)
    rr_samp <- runif(n_mc, rr_scenarios[[r]]$lo, rr_scenarios[[r]]$hi)
    paf_samp <- levin_paf(pe_samp, rr_samp)
    
    key <- paste0("prev_", p, "_rr_", r)
    paf_settings[[key]] <- list(
      prev = p, rr = r,
      prev_label = prevalence_scenarios[[p]]$label,
      rr_label = rr_scenarios[[r]]$label,
      name = sprintf("pe=%s, RR=%s", prevalence_scenarios[[p]]$label, rr_scenarios[[r]]$label),
      paf = list(estimate = median(paf_samp), lo = quantile(paf_samp, 0.025), hi = quantile(paf_samp, 0.975))
    )
    
    cat(sprintf("  %s: PAF = %.1f%% (%.1f-%.1f%%)\n",
                paf_settings[[key]]$name,
                median(paf_samp) * 100,
                quantile(paf_samp, 0.025) * 100,
                quantile(paf_samp, 0.975) * 100))
  }
}

# =============================================================================
# 2. DDE Model
# =============================================================================

Delta_s <- 0.5

dde_2strain_intervention <- function(parms, paf_range, coverage = 1.00,
                                      efficacy = 0.89, interv_start = 2025, end_time = 2035) {
  N <- 1e5; mu <- 0.013
  cov <- coverage
  paf_draw <- runif(1, paf_range$lo, paf_range$hi)
  interv_eff <- 1 - cov * efficacy * paf_draw

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
      lam_s <- (beta / N) * infectious_s; lam_r <- (beta / N) * infectious_r; lam <- lam_s + lam_r
      theta_s <- force_theta_s(time); phi_s <- force_phi_s(time); d_s <- force_d_s(time)
      theta_r <- force_theta_r(time); phi_r <- force_phi_r(time); d_r <- force_d_r(time)
      Delta_r <- force_Delta_r(time)
      ie <- ifelse(time >= interv_start, interv_eff, 1)
      alpha_f_t <- alpha_f * ie; alpha_r_t <- alpha_r * ie

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
      Inc_DS <- sigma * A_s; Inc_RR_acq <- Af + sigma * A_r_acq; Inc_RR <- sigma * A_r + Af
      return(list(c(dS, dI_s, dN_s, dA_s, dC_s, dT_s, dP_s,
                    dI_r, dN_r, dA_r, dC_r, dT_r, dP_r, dA_r_acq),
                  IncRR = (Inc_RR / Pop) * 1e5,
                  IncRR_acq = (Inc_RR_acq / Pop) * 1e5,
                  MorRR = (d_r * C_r / Pop) * 1e5))
    })
  }
  rr0 <- 0.025
  yini <- c(S = N - 500,
    I_s = 0, N_s = 0, A_s = 100*(1-rr0), C_s = 400*(1-rr0), T_s = 0, P_s = 0,
    I_r = 0, N_r = 0, A_r = 100*rr0, C_r = 400*rr0, T_r = 0, P_r = 0,
    A_r_acq = 0)
  deSolve::dede(yini, seq(1500, end_time, 1), des, parms, method = "lsoda")
}

# =============================================================================
# 3. Run baseline
# =============================================================================

cat("\n=== Running baseline ===\n")
pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(fitpts))
baseline_res <- rbindlist(lapply(1:nrow(fitpts), function(i) {
  p <- as.numeric(fitpts[i, ]); names(p) <- names(fitpts)
  paf_zero <- list(lo = 0, hi = 0)
  out <- tryCatch({
    as.data.frame(dde_2strain_intervention(p, paf_range = paf_zero,
                                            coverage = 0, efficacy = 0))
  }, error = function(e) NULL)
  if (is.null(out)) { pb$tick(); return(NULL) }
  out$run <- i; pb$tick(); out
}))[time >= 2020]
baseline_res$setting <- "baseline"

# =============================================================================
# 4. Run all 9 PAF settings at universal coverage
# =============================================================================

cat("\n=== Running 9 PAF sensitivity scenarios (universal coverage) ===\n")

sensitivity_results <- list()
sensitivity_results[["baseline"]] <- baseline_res

for (key in names(paf_settings)) {
  s <- paf_settings[[key]]
  cat(sprintf("\n  Running: %s (PAF = %.1f%%)\n", s$name, s$paf$estimate * 100))
  pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(fitpts))
  res <- rbindlist(lapply(1:nrow(fitpts), function(i) {
    p <- as.numeric(fitpts[i, ]); names(p) <- names(fitpts)
    out <- tryCatch({
      as.data.frame(dde_2strain_intervention(p, paf_range = s$paf))
    }, error = function(e) NULL)
    if (is.null(out)) { pb$tick(); return(NULL) }
    out$run <- i; pb$tick(); out
  }))[time >= 2020]
  res$setting <- key
  sensitivity_results[[key]] <- res
}

results_all <- rbindlist(sensitivity_results)

# =============================================================================
# 5. Summarise results
# =============================================================================

cat("\n=== Summarising ===\n")

bl_2035 <- results_all[setting == "baseline" & time == 2035, .(run, IncRR_base = IncRR, MorRR_base = MorRR)]

impact_sens <- rbindlist(lapply(names(paf_settings), function(key) {
  s <- paf_settings[[key]]
  sc_2035 <- results_all[setting == key & time == 2035, .(run, IncRR_int = IncRR, MorRR_int = MorRR)]
  imp <- merge(bl_2035, sc_2035, by = "run")
  imp[, `:=`(
    IncRR_reduction = (IncRR_base - IncRR_int) / IncRR_base * 100,
    IncRR_averted = IncRR_base - IncRR_int,
    MorRR_reduction = (MorRR_base - MorRR_int) / (MorRR_base + 1e-10) * 100
  )]
  imp[, .(
    setting = key, prev = s$prev, rr = s$rr,
    prev_label = s$prev_label, rr_label = s$rr_label,
    paf_median = s$paf$estimate,
    IncRR_reduction_median = median(IncRR_reduction),
    IncRR_reduction_lo = quantile(IncRR_reduction, 0.025),
    IncRR_reduction_hi = quantile(IncRR_reduction, 0.975),
    IncRR_averted_median = median(IncRR_averted),
    IncRR_averted_lo = quantile(IncRR_averted, 0.025),
    IncRR_averted_hi = quantile(IncRR_averted, 0.975),
    MorRR_reduction_median = median(MorRR_reduction),
    MorRR_reduction_lo = quantile(MorRR_reduction, 0.025),
    MorRR_reduction_hi = quantile(MorRR_reduction, 0.975)
  )]
}))

# Save
dir.create(here("outputs", "results", "sensitivity"), showWarnings = FALSE, recursive = TRUE)
export(results_all, here("outputs", "results", "sensitivity", "sens_raw.Rdata"))
export(impact_sens, here("outputs", "results", "sensitivity", "sens_impact.Rdata"))

# =============================================================================
# 6. Print results table
# =============================================================================

cat("\n===========================================================================\n")
cat("SENSITIVITY ANALYSIS: % Reduction in RR-TB Incidence by 2035\n")
cat("Universal INH-R DST (100% coverage, 89% efficacy)\n")
cat("===========================================================================\n\n")

cat(sprintf("%-10s %-15s %8s %25s\n", "INH-R prev", "RR", "PAF", "RR-TB reduction (95% UI)"))
cat(paste(rep("-", 62), collapse = ""), "\n")
for (i in 1:nrow(impact_sens)) {
  r <- impact_sens[i]
  cat(sprintf("%-10s %-15s %7.1f%% %7.1f%% (%.1f-%.1f%%)\n",
              r$prev_label, r$rr_label, r$paf_median * 100,
              r$IncRR_reduction_median, r$IncRR_reduction_lo, r$IncRR_reduction_hi))
}

# =============================================================================
# 7. Forest plot (inverted axes per Marc)
# =============================================================================

impact_sens[, prev_group := factor(
  ifelse(prev == "low", "Low (5%)",
         ifelse(prev == "medium", "Average (10%)", "High (20%)")),
  levels = c("Low (5%)", "Average (10%)", "High (20%)")
)]
impact_sens[, rr_display := factor(rr_label, levels = c("Observed", "2-fold lower", "5-fold lower"))]
impact_sens[, label := sprintf("%.1f\n(%.1f\u2013%.1f)",
                                IncRR_reduction_median, IncRR_reduction_lo, IncRR_reduction_hi)]

bg_colours <- c("Low (5%)" = "#C8E6C9", "Average (10%)" = "#FFF9C4", "High (20%)" = "#FFCDD2")

p <- ggplot(impact_sens, aes(x = rr_display, y = IncRR_reduction_median)) +
  geom_rect(aes(fill = prev_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            alpha = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = bg_colours) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3, linetype = "dashed") +
  geom_linerange(aes(ymin = IncRR_reduction_lo, ymax = IncRR_reduction_hi),
                 colour = "black", linewidth = 1) +
  geom_point(colour = "black", size = 3) +
  geom_text(aes(y = IncRR_reduction_hi + 3, label = label),
            size = 3, lineheight = 0.85, colour = "black") +
  facet_wrap(~ prev_group, ncol = 3, strip.position = "top") +
  scale_y_continuous(breaks = seq(0, 80, 20)) +
  coord_cartesian(ylim = c(0, 85)) +
  labs(title = "Prevalence of isoniazid resistance",
       y = "Reduction in RR-TB incidence (%)",
       x = "Relative risk for acquired\nrifampicin resistance") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey80"),
    axis.line.y = element_line(linewidth = 0.4, colour = "grey30"),
    axis.ticks.y = element_line(linewidth = 0.3, colour = "grey30"),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(10, 15, 10, 10)
  )

dir.create(here("plots"), showWarnings = FALSE)
pdf(here("plots", "sensitivity_forest.pdf"), width = 10, height = 6)
print(p)
dev.off()

cat("\nSaved: plots/sensitivity_forest.pdf\n")
