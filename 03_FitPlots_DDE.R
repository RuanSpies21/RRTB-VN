## RR-TB Vietnam Model - 03: Fit Plots
## Visualize calibration results — DDE model with renamed parameters

library(rio)
library(here)
library(ggplot2)
library(reshape2)
library(data.table)
library(patchwork)
library(deSolve)
library(progress)
library(scales)

targets <- as.data.table(import(here("data", "fit", "targets.Rdata")))
fitpts <- as.data.table(import(here("outputs", "pts", "fitpts.Rdata")))
cat("Loaded", nrow(fitpts), "plausible parameter sets\n")

Delta_s <- 0.5

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
      infectious_s <- kappa * A_s + C_s; infectious_r <- f * (kappa * A_r + C_r)
      lam_s <- (beta / N) * infectious_s; lam_r <- (beta / N) * infectious_r; lam <- lam_s + lam_r
      theta_s <- force_theta_s(time); phi_s <- force_phi_s(time); d_s <- force_d_s(time)
      theta_r <- force_theta_r(time); phi_r <- force_phi_r(time); d_r <- force_d_r(time)
      Delta_r <- force_Delta_r(time)
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
        pAc = Inc_RR_acq / (Inc_RR + 1e-10), pRR = Inc_RR / (Inc_DS + Inc_RR + 1e-10),
        Prev = ((A_s + C_s + A_r + C_r) / Pop) * 1e5,
        Mor = ((d_s * C_s + d_r * C_r) / Pop) * 1e5,
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

# Run model
cat("Running model...\n")
pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(fitpts))
all_results <- lapply(1:nrow(fitpts), function(i) {
  p <- as.numeric(fitpts[i,]); names(p) <- names(fitpts)
  out <- as.data.frame(dde_2strain(p)); out$run <- i; pb$tick(); out
})
results_raw <- rbindlist(all_results)[time >= 2000]

results_long <- data.table::melt(results_raw, id.vars = c("time", "run"),
  measure.vars = c("Inc", "IncRR", "pAc", "pRR", "Prev", "Mor", "Notif"),
  variable.name = "var", value.name = "value")
results <- results_long[, .(`2.5%` = quantile(value, 0.025), `25%` = quantile(value, 0.25),
  `50%` = quantile(value, 0.50), `75%` = quantile(value, 0.75),
  `97.5%` = quantile(value, 0.975)), by = .(time, var)]
export(results, here("outputs", "results.Rdata"))

# Plots
cl <- "#CE2931"; al <- 0.2; lw <- 0.6; brks <- seq(2000, 2024, 2); lims <- c(1999, 2025); yrs <- 2000:2024
theme_fit <- theme(panel.background = element_blank(), panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
  panel.grid.major.x = element_blank(), axis.line = element_line(linewidth = lw), axis.ticks = element_line(linewidth = lw),
  axis.title.x = element_text(vjust = -1, margin = margin(b = 10), size = 11),
  axis.title.y = element_text(vjust = 2, margin = margin(l = 10), size = 11),
  axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 10),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"), legend.position = "none")

targets_inc <- targets[grepl("^Inc[0-9]", name)][, var := "Inc"]
targets_incrr <- targets[grepl("^IncRR", name)][, var := "IncRR"]
targets_pac <- targets[grepl("^pAc", name)][, var := "pAc"]

rp <- function(d, tg, yb, yl, ttl, ylab, pct = FALSE) {
  p <- ggplot(d, aes(time)) + geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = cl, alpha = al) +
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`), fill = cl, alpha = al * 2) + geom_line(aes(y = `50%`), linewidth = lw, colour = cl) +
    geom_errorbar(data = tg, aes(x = time, ymin = lo, ymax = hi), linewidth = lw, width = 0.5) +
    scale_x_continuous(breaks = brks, limits = lims, expand = c(0, 0)) + coord_cartesian(ylim = yl) +
    labs(title = ttl, x = "Year", y = ylab) + theme_fit
  if (pct) p <- p + scale_y_continuous(breaks = yb, labels = percent, expand = c(0, 0)) else p <- p + scale_y_continuous(breaks = yb, expand = c(0, 0))
  p
}

max_y_A <- max(targets_inc$hi, results[var == "Inc" & time %in% yrs, `97.5%`], na.rm = TRUE)
y_upper_A <- ceiling(max_y_A / 50) * 50

a <- rp(results[var == "Inc" & time %in% yrs], targets_inc, seq(0, y_upper_A, 50), c(0, y_upper_A), "A. Total TB Incidence", "Per 100k")
b <- rp(results[var == "IncRR" & time %in% yrs], targets_incrr, seq(0, 20, 5), c(0, 20), "B. RR-TB Incidence", "Per 100k")
cc <- rp(results[var == "pAc" & time %in% yrs], targets_pac, seq(0, 0.5, 0.1), c(0, 0.5), "C. Proportion Acquired RR", "Proportion", TRUE)

pdf(here("plots", "01_calibration_fit.pdf"), width = 12, height = 4)
print(a | b | cc)
dev.off()

key <- c("beta", "f", "kappa", "alpha_f", "alpha_r", "theta_s1", "theta_r1", "rho_s", "sigma", "gamma")
fl <- data.table::melt(fitpts, measure.vars = names(fitpts), variable.name = "variable", value.name = "value")
pdf(here("plots", "03_parameters.pdf"), width = 12, height = 8)
print(ggplot(fl[variable %in% key], aes(value)) + geom_histogram(bins = 15, fill = cl, alpha = 0.7, colour = "white") +
  facet_wrap(~variable, scales = "free", ncol = 5) + labs(title = "Calibrated Parameters", x = "Value", y = "Count") +
  theme_fit + theme(strip.text = element_text(face = "bold")))
dev.off()

cat("\nCalibration (2024):\n")
for (v in c("Inc", "IncRR", "pAc", "pRR")) { r <- results[time == 2024 & var == v]; cat(sprintf("  %s: %.2f (%.2f-%.2f)\n", v, r$`50%`, r$`2.5%`, r$`97.5%`)) }
cat("=== 03_FitPlots.R completed ===\n")
