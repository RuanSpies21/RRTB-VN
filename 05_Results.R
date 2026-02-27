## RR-TB Vietnam Model - 05: Results
## Process and summarize intervention scenario results

library(rio)
library(here)
library(data.table)
options(scipen = 999)

# =============================================================================
# 1. Load Data
# =============================================================================

results_raw <- as.data.table(import(here("outputs", "results", "scenarios_raw.Rdata")))
results_summary <- as.data.table(import(here("outputs", "results", "scenarios_summary.Rdata")))
impact_raw <- as.data.table(import(here("outputs", "results", "impact_raw.Rdata")))
impact_summary <- as.data.table(import(here("outputs", "results", "impact_summary.Rdata")))

# =============================================================================
# 2. Cumulative Metrics (2025-2035)
# =============================================================================

cumulative <- results_raw[time >= 2025 & time == floor(time)][order(scenario, run, time)][, .(
  cumIncRR = sum(IncRR), cumIncRR_acq = sum(IncRR_acq), cumMorRR = sum(MorRR),
  cumInc = sum(Inc), IncRR_2035 = IncRR[time == 2035], pAc_2035 = pAc[time == 2035]
), by = .(scenario, scenario_name, run)]

# =============================================================================
# 3. Compare to Baseline
# =============================================================================

bl <- cumulative[scenario == "baseline", .(run, cumIncRR_base=cumIncRR, cumIncRR_acq_base=cumIncRR_acq,
                                            cumMorRR_base=cumMorRR, IncRR_2035_base=IncRR_2035)]

comparison_all <- rbindlist(lapply(setdiff(unique(cumulative$scenario), "baseline"), function(scen) {
  int <- cumulative[scenario == scen, .(run, scenario, scenario_name, cumIncRR_int=cumIncRR,
                                         cumIncRR_acq_int=cumIncRR_acq, cumMorRR_int=cumMorRR, IncRR_2035_int=IncRR_2035)]
  comp <- merge(bl, int, by = "run")
  comp[, `:=`(
    IncRR_averted_cum = cumIncRR_base - cumIncRR_int,
    IncRR_acq_averted_cum = cumIncRR_acq_base - cumIncRR_acq_int,
    MorRR_averted_cum = cumMorRR_base - cumMorRR_int,
    IncRR_reduction_pct = (cumIncRR_base - cumIncRR_int) / cumIncRR_base * 100,
    IncRR_2035_reduction_pct = (IncRR_2035_base - IncRR_2035_int) / IncRR_2035_base * 100
  )]
}))

# =============================================================================
# 4. Summary Table
# =============================================================================

qsum <- function(x, p) quantile(x, p, na.rm = TRUE)

results_table <- comparison_all[, .(
  scenario_name = first(scenario_name),
  IncRR_averted_cum_med = median(IncRR_averted_cum), IncRR_averted_cum_lo = qsum(IncRR_averted_cum, 0.025), IncRR_averted_cum_hi = qsum(IncRR_averted_cum, 0.975),
  MorRR_averted_cum_med = median(MorRR_averted_cum), MorRR_averted_cum_lo = qsum(MorRR_averted_cum, 0.025), MorRR_averted_cum_hi = qsum(MorRR_averted_cum, 0.975),
  IncRR_reduction_pct_med = median(IncRR_reduction_pct), IncRR_reduction_pct_lo = qsum(IncRR_reduction_pct, 0.025), IncRR_reduction_pct_hi = qsum(IncRR_reduction_pct, 0.975),
  IncRR_2035_reduction_med = median(IncRR_2035_reduction_pct), IncRR_2035_reduction_lo = qsum(IncRR_2035_reduction_pct, 0.025), IncRR_2035_reduction_hi = qsum(IncRR_2035_reduction_pct, 0.975)
), by = scenario]

baseline_summary <- cumulative[scenario == "baseline", .(
  cumIncRR_med = median(cumIncRR), cumIncRR_lo = qsum(cumIncRR, 0.025), cumIncRR_hi = qsum(cumIncRR, 0.975),
  cumMorRR_med = median(cumMorRR), cumMorRR_lo = qsum(cumMorRR, 0.025), cumMorRR_hi = qsum(cumMorRR, 0.975),
  IncRR_2035_med = median(IncRR_2035), IncRR_2035_lo = qsum(IncRR_2035, 0.025), IncRR_2035_hi = qsum(IncRR_2035, 0.975),
  pAc_2035_med = median(pAc_2035), pAc_2035_lo = qsum(pAc_2035, 0.025), pAc_2035_hi = qsum(pAc_2035, 0.975)
)]

export(cumulative, here("outputs", "outs", "cumulative.Rdata"))
export(comparison_all, here("outputs", "outs", "comparisons.Rdata"))
export(results_table, here("outputs", "outs", "results_table.Rdata"))
export(baseline_summary, here("outputs", "outs", "baseline_summary.Rdata"))

# =============================================================================
# 5. Print Results
# =============================================================================

cat("\n================================================================\n")
cat("INH-R INTERVENTION IMPACT ON RR-TB (2025-2035)\n")
cat("================================================================\n")

cat(sprintf("\nBaseline cumulative RR-TB: %.1f (%.1f-%.1f) per 100k\n",
            baseline_summary$cumIncRR_med, baseline_summary$cumIncRR_lo, baseline_summary$cumIncRR_hi))
cat(sprintf("Baseline RR-TB 2035: %.1f (%.1f-%.1f) per 100k\n",
            baseline_summary$IncRR_2035_med, baseline_summary$IncRR_2035_lo, baseline_summary$IncRR_2035_hi))

for (i in 1:nrow(results_table)) {
  r <- results_table[i]
  cat(sprintf("\n%s:\n  Cumulative cases averted: %.1f (%.1f-%.1f) per 100k\n  Cumulative reduction: %.1f%% (%.1f-%.1f%%)\n  2035 reduction: %.1f%% (%.1f-%.1f%%)\n",
    r$scenario_name, r$IncRR_averted_cum_med, r$IncRR_averted_cum_lo, r$IncRR_averted_cum_hi,
    r$IncRR_reduction_pct_med, r$IncRR_reduction_pct_lo, r$IncRR_reduction_pct_hi,
    r$IncRR_2035_reduction_med, r$IncRR_2035_reduction_lo, r$IncRR_2035_reduction_hi))
}

# Scale to Vietnam (~100M)
cat("\n--- Scaled to Vietnam (Pop ~100M) ---\n")
for (i in 1:nrow(results_table)) {
  r <- results_table[i]
  cat(sprintf("%s: %s (%s-%s) cases averted, %s (%s-%s) deaths averted\n",
    r$scenario_name,
    format(round(r$IncRR_averted_cum_med * 1000), big.mark=","),
    format(round(r$IncRR_averted_cum_lo * 1000), big.mark=","),
    format(round(r$IncRR_averted_cum_hi * 1000), big.mark=","),
    format(round(r$MorRR_averted_cum_med * 1000), big.mark=","),
    format(round(r$MorRR_averted_cum_lo * 1000), big.mark=","),
    format(round(r$MorRR_averted_cum_hi * 1000), big.mark=",")))
}

cat("\n=== 05_Results.R completed ===\n")
