###############################################################################
## Undetected isoniazid resistance and the emergence of rifampicin-resistant
## tuberculosis: a prospective observational study and transmission modelling
## analysis
##
## Statistical analysis of acquired rifampicin resistance
###############################################################################

library(tidyverse)
library(boot)
library(cobalt)
library(survey)
library(survival)
library(tableone)
library(patchwork)

dir.create("output", showWarnings = FALSE)
set.seed(2024)

# 1. Settings

# Population denominators from NTP notification data, Ho Chi Minh City,
# March 2020 - July 2024
pop_denominators <- c(
  subcohort1 = 42843 - 1775 - 550,  # RS-TB without recurrenc
  subcohort2 = 1775,                # RS-TB recurrence
  subcohort3 = 550                  # RR-TB recurrence
)

covars_primary <- c("age", "SEX", "DIABETES", "HIV", "fastlin_major_lineage")
trunc_primary  <- c(0.01, 0.99)
n_boot         <- 1000

# 2. Functions

# Sampling weight = population denominator / number analysed, per sub-cohort
add_sampling_weights <- function(data) {
  data %>%
    add_count(subcohort, name = "n_sub") %>%
    mutate(sampling_weight = pop_denominators[subcohort] / n_sub) %>%
    select(-n_sub)
}

# Stabilised inverse probability of exposure weights
# combined with sampling weights
add_ipew <- function(data, covars = covars_primary, trunc = trunc_primary) {
  ps_fit <- glm(reformulate(covars, "isoniazid_bin"), family = binomial, data = data)
  ps     <- predict(ps_fit, type = "response")
  p_exp  <- mean(data$isoniazid_bin)
  w      <- if_else(data$isoniazid_bin == 1, p_exp / ps, (1 - p_exp) / (1 - ps))
  w_raw  <- w
  if (!is.null(trunc)) {
    q <- quantile(w, trunc)
    w <- pmin(pmax(w, q[1]), q[2])
  }
  data %>% mutate(ps = ps, ipew_raw = w_raw, ipew = w,
                  combined_weight = sampling_weight * ipew)
}

# Weighted risks, risk difference, risk ratio and PAF
weighted_estimates <- function(data, outcome = "outcome", weight = "combined_weight") {
  y  <- data[[outcome]]
  w  <- data[[weight]]
  e  <- data$isoniazid_bin
  r1 <- weighted.mean(y[e == 1], w[e == 1])
  r0 <- weighted.mean(y[e == 0], w[e == 0])
  rr <- r1 / r0
  pe <- if (any(y == 1)) weighted.mean(e[y == 1], w[y == 1]) else NA
  c(risk_inh_r = r1, risk_inh_s = r0, RD = r1 - r0, RR = rr,
    PAF = pe * (rr - 1) / rr)
}

# Full estimation pipeline (propensity score model refitted each call)
estimate <- function(data, covars = covars_primary, trunc = trunc_primary,
                     outcome = "outcome") {
  weighted_estimates(add_ipew(data, covars, trunc), outcome)
}

# Bootstrap: sampling weights fixed, all propensity score steps repeated
boot_estimate <- function(data, R = n_boot, ...) {
  boot(data, function(d, i) estimate(d[i, ], ...), R = R)
}

# Percentile CI
boot_ci <- function(b, index) {
  t <- b$t[, index]
  ci <- boot.ci(b, type = "perc", index = index)$percent[4:5]
  tibble(estimate = b$t0[index], lower = ci[1], upper = ci[2],
         n_nonfinite = sum(!is.finite(t)))
}

# Kish effective sample size
kish_ess <- function(w) sum(w)^2 / sum(w^2)

# E-value for a risk ratio
e_value <- function(rr) rr + sqrt(rr * (rr - 1))

fmt_pct <- function(x, d = 2) sprintf(paste0("%.", d, "f"), 100 * x)

# 3. Weights

final_cohort <- final_cohort %>%
  add_sampling_weights() %>%
  add_ipew()

final_cohort %>%
  distinct(subcohort, sampling_weight) %>%
  arrange(subcohort) %>%
  mutate(N = pop_denominators[subcohort]) %>%
  print()

summary(final_cohort$ipew_raw)
summary(final_cohort$ipew)
summary(final_cohort$combined_weight)

# 4. Table 1: baseline characteristics by outcome

table1 <- CreateTableOne(
  vars       = c("age", "SEX", "HIV", "DIABETES", "fastlin_major_lineage",
                 "isoniazid_bin", "pyrazinamide", "ethambutol"),
  factorVars = c("isoniazid_bin"),
  strata     = "outcome",
  data       = final_cohort,
  test       = TRUE
)
print(table1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)
write.csv(print(table1, showAllLevels = TRUE, printToggle = FALSE),
          "output/table1.csv")

# 5. Primary analysis: weighted risks, risk ratio and PAF

primary <- weighted_estimates(final_cohort)

set.seed(2024)
boot_primary <- boot_estimate(final_cohort)

primary_results <- map_dfr(1:5, ~ boot_ci(boot_primary, .x)) %>%
  mutate(measure = names(primary), .before = 1)
print(primary_results)
write_csv(primary_results, "output/primary_results.csv")

# E-values
e_value(primary["RR"])
e_value(primary_results$lower[primary_results$measure == "RR"])

# Weighted prevalence of isoniazid resistance among RS-TB
weighted.mean(final_cohort$isoniazid_bin, final_cohort$sampling_weight)


# 6. Weighting diagnostics
# Propensity score overlap
p_ps <- ggplot(final_cohort, aes(ps, fill = factor(isoniazid_bin))) +
  geom_histogram(data = ~ filter(.x, isoniazid_bin == 1), aes(y = after_stat(count)),
                 bins = 30, colour = "white") +
  geom_histogram(data = ~ filter(.x, isoniazid_bin == 0), aes(y = -after_stat(count)),
                 bins = 30, colour = "white") +
  scale_fill_manual(values = c(`0` = "#3b82f6", `1` = "#ef4444"),
                    labels = c("INH-S", "INH-R"), name = NULL) +
  scale_y_continuous(labels = abs) +
  labs(title = "A", x = "Propensity score", y = "Count") +
  theme_minimal() + theme(legend.position = "top")

# Weight distribution after truncation
p_w <- ggplot(final_cohort, aes(ipew, fill = factor(isoniazid_bin))) +
  geom_histogram(bins = 30, colour = "white", position = "identity", alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c(`0` = "#3b82f6", `1` = "#ef4444"),
                    labels = c("INH-S", "INH-R"), name = NULL) +
  labs(title = "B", x = "Stabilised weight (truncated)", y = "Count") +
  theme_minimal() + theme(legend.position = "top")

# Covariate balance
bal_formula <- reformulate(covars_primary, "isoniazid_bin")
bal <- bal.tab(bal_formula, data = final_cohort, weights = final_cohort$combined_weight,
               un = TRUE, s.d.denom = "pooled", thresholds = c(m = 0.1))
print(bal)

p_bal <- love.plot(bal, abs = TRUE, var.order = "unadjusted",
                   colors = c("#ef4444", "#10b981"),
                   sample.names = c("Unweighted", "Weighted"), title = "C") +
  theme_minimal() + theme(legend.position = "top")

ggsave("output/figure_s1_weighting_diagnostics.pdf", p_ps + p_w + p_bal,
       width = 15, height = 5)

# Effective sample size
table_s7 <- final_cohort %>%
  group_by(isoniazid_bin) %>%
  summarise(N                = n(),
            ESS_ps_weights   = kish_ess(ipew),
            ESS_sampling     = kish_ess(sampling_weight),
            ESS_combined     = kish_ess(combined_weight)) %>%
  mutate(across(starts_with("ESS"), round))
print(table_s7)
write_csv(table_s7, "output/table_s7_ess.csv")

# 7. Sensitivity analyses: weighting and covariates

risk_row <- function(label, est) {
  tibble(analysis = label, risk_inh_r = est["risk_inh_r"],
         risk_inh_s = est["risk_inh_s"], RR = est["RR"])
}

table_s8 <- bind_rows(
  risk_row("Unweighted",
           weighted_estimates(mutate(final_cohort, w1 = 1), weight = "w1")),
  risk_row("Sampling weights only",
           weighted_estimates(final_cohort, weight = "sampling_weight")),
  risk_row("No truncation",                        estimate(final_cohort, trunc = NULL)),
  risk_row("Truncated 1st/99th (primary)",         estimate(final_cohort)),
  risk_row("Truncated 2.5th/97.5th",               estimate(final_cohort, trunc = c(0.025, 0.975))),
  risk_row("Truncated 5th/95th",                   estimate(final_cohort, trunc = c(0.05, 0.95))),
  risk_row("+ calendar year",
           estimate(final_cohort, covars = c(covars_primary, "index_year"))),
  risk_row("+ DTU",
           estimate(final_cohort, covars = c(covars_primary, "DTU"))),
  risk_row("+ calendar year and DTU",
           estimate(final_cohort, covars = c(covars_primary, "index_year", "DTU")))
)
print(table_s8)
write_csv(table_s8, "output/table_s8_sensitivity.csv")

# 8. Sensitivity analysis: <=5 SNP threshold for same-strain recurrence

final_cohort <- final_cohort %>%
  mutate(outcome_snp5 = as.integer(outcome == 1 & snp_distance <= 5))

sum(final_cohort$outcome != final_cohort$outcome_snp5)

set.seed(2024)
boot_snp5 <- boot_estimate(final_cohort, outcome = "outcome_snp5")
snp5_results <- map_dfr(1:5, ~ boot_ci(boot_snp5, .x)) %>%
  mutate(measure = names(primary), .before = 1)
print(snp5_results)
write_csv(snp5_results, "output/sensitivity_snp5.csv")

# 9. Time-to-event sensitivity analysis

run_tte <- function(data, end_date) {
  d <- data %>%
    mutate(
      has_event = !is.na(recur_date) & recur_date <= end_date,
      stop_date = if_else(has_event, recur_date, end_date),
      time_m    = as.numeric(stop_date - index_date) / 30.4375,
      event     = as.integer(has_event & outcome == 1)
    ) %>%
    filter(time_m > 0)
  
  person_time <- d %>%
    group_by(isoniazid_bin) %>%
    summarise(n = n(), events = sum(event),
              person_years     = sum(time_m) / 12,
              median_fu_months = median(time_m),
              w_rate_per_1000py = 1000 * sum(event * combined_weight) /
                sum(time_m / 12 * combined_weight),
              .groups = "drop")
  
  des  <- svydesign(ids = ~1, strata = ~subcohort, weights = ~combined_weight, data = d)
  cox  <- svycoxph(Surv(time_m, event) ~ isoniazid_bin, design = des)
  pois <- svyglm(event ~ isoniazid_bin + offset(log(time_m)), design = des,
                 family = quasipoisson())
  
  list(
    person_time = person_time,
    estimates = tibble(
      censoring = format(end_date),
      measure   = c("Hazard ratio", "Incidence rate ratio"),
      estimate  = c(exp(coef(cox)), exp(coef(pois))["isoniazid_bin"]),
      lower     = c(exp(confint(cox))[1], exp(confint(pois))["isoniazid_bin", 1]),
      upper     = c(exp(confint(cox))[2], exp(confint(pois))["isoniazid_bin", 2])
    )
  )
}

end_study <- as.Date("2024-07-31")
tte_main  <- run_tte(final_cohort, end_study)
tte_ext   <- run_tte(final_cohort, end_study + 365)

print(tte_main$person_time)
table_s9 <- bind_rows(tte_main$estimates, tte_ext$estimates)
print(table_s9)
write_csv(tte_main$person_time, "output/table_s9_person_time.csv")
write_csv(table_s9, "output/table_s9_estimates.csv")


# 10. Companion drug resistance


final_cohort <- final_cohort %>%
  mutate(
    pza_r = pyrazinamide %in% "R",
    emb_r = ethambutol %in% "R",
    profile = factor(case_when(
      isoniazid_bin == 0 ~ "Isoniazid-susceptible",
      !pza_r & !emb_r    ~ "Isoniazid only",
      pza_r & !emb_r     ~ "Isoniazid + pyrazinamide",
      !pza_r & emb_r     ~ "Isoniazid + ethambutol",
      TRUE               ~ "Isoniazid + pyrazinamide + ethambutol"),
      levels = c("Isoniazid-susceptible", "Isoniazid only", "Isoniazid + pyrazinamide",
                 "Isoniazid + ethambutol", "Isoniazid + pyrazinamide + ethambutol")),
    companion = factor(case_when(
      isoniazid_bin == 0 ~ "Isoniazid-susceptible",
      pza_r | emb_r      ~ "Isoniazid + any additional drug",
      TRUE               ~ "Isoniazid only"),
      levels = c("Isoniazid-susceptible", "Isoniazid only",
                 "Isoniazid + any additional drug"))
  )

# Companion drug resistance by isoniazid status
final_cohort %>%
  group_by(isoniazid_bin) %>%
  summarise(n = n(), companion_r = sum(pza_r | emb_r),
            weighted_prop = weighted.mean(pza_r | emb_r, sampling_weight))

# Share of HR-TB and of acquired RR among HR-TB with companion resistance
final_cohort %>%
  filter(isoniazid_bin == 1) %>%
  group_by(companion_r = pza_r | emb_r) %>%
  summarise(w_n = sum(sampling_weight), w_events = sum(outcome * sampling_weight)) %>%
  mutate(share_hr_tb = w_n / sum(w_n), share_events = w_events / sum(w_events))

# Counts and sampling-weighted risks with bootstrap CIs
profile_groups <- list(
  "Isoniazid-susceptible"                 = "Isoniazid-susceptible",
  "Isoniazid only"                        = "Isoniazid only",
  "Isoniazid + pyrazinamide"              = "Isoniazid + pyrazinamide",
  "Isoniazid + ethambutol"                = "Isoniazid + ethambutol",
  "Isoniazid + pyrazinamide + ethambutol" = "Isoniazid + pyrazinamide + ethambutol",
  "Isoniazid + any additional drug"       = c("Isoniazid + pyrazinamide",
                                              "Isoniazid + ethambutol",
                                              "Isoniazid + pyrazinamide + ethambutol")
)

profile_risks <- function(d, i) {
  b <- d[i, ]
  risks <- sapply(profile_groups, function(g) {
    x <- b[b$profile %in% g, ]
    weighted.mean(x$outcome, x$sampling_weight)
  })
  c(risks, risks[-1] / risks["Isoniazid only"])   # RR vs isoniazid only
}

set.seed(2024)
boot_profile <- boot(final_cohort, profile_risks, R = n_boot)
k <- length(profile_groups)

table_s4 <- tibble(
  profile   = names(profile_groups),
  n         = sapply(profile_groups, function(g) sum(final_cohort$profile %in% g)),
  events    = sapply(profile_groups, function(g) sum(final_cohort$outcome[final_cohort$profile %in% g])),
  risk_pct  = 100 * boot_profile$t0[1:k],
  risk_lo   = 100 * apply(boot_profile$t[, 1:k], 2, quantile, 0.025, na.rm = TRUE),
  risk_hi   = 100 * apply(boot_profile$t[, 1:k], 2, quantile, 0.975, na.rm = TRUE),
  rr_vs_mono    = c(NA, boot_profile$t0[(k + 1):(2 * k - 1)]),
  rr_vs_mono_lo = c(NA, apply(boot_profile$t[, (k + 1):(2 * k - 1)], 2, quantile, 0.025, na.rm = TRUE)),
  rr_vs_mono_hi = c(NA, apply(boot_profile$t[, (k + 1):(2 * k - 1)], 2, quantile, 0.975, na.rm = TRUE))
)

# Risk ratios vs isoniazid-susceptible
des_s <- svydesign(ids = ~1, strata = ~subcohort, weights = ~sampling_weight,
                   data = final_cohort)
rr_vs_s <- function(f) {
  m <- svyglm(f, design = des_s, family = quasipoisson())
  exp(cbind(coef(m), confint(m)))[-1, , drop = FALSE]
}
rr_s <- rbind(rr_vs_s(outcome ~ profile), rr_vs_s(outcome ~ companion)[2, , drop = FALSE])

table_s4 <- table_s4 %>%
  mutate(rr_vs_inh_s    = c(NA, rr_s[, 1]),
         rr_vs_inh_s_lo = c(NA, rr_s[, 2]),
         rr_vs_inh_s_hi = c(NA, rr_s[, 3]))
print(table_s4, width = Inf)
write_csv(table_s4, "output/table_s4_companion.csv")
