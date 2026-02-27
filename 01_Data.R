## RR-TB Vietnam Model - 01: Data
## Load WHO incidence data and set up calibration targets

library(data.table)
library(rio)
library(here)
library(readxl)

# Create directories ==========
for (d in c("data", "data/fit", "data/who", "outputs", "outputs/pts", "outputs/waves",
            "outputs/wave_res", "outputs/ems", "outputs/diagnostics", "outputs/results",
            "outputs/outs", "plots")) {
  dir.create(here(d), showWarnings = FALSE, recursive = TRUE)
}

# =============================================================================
# 1. Load & Process WHO Vietnam Data
# =============================================================================

who_data <- as.data.table(read_excel(here("data", "who", "vnm_data.xlsx")))

# Ensure numeric
num_cols <- c("year", "population", "dstb_incidence", "dstb_incidence_ll",
              "dstb_incidence_ul", "rr_tb_proportion", "rrtb_cases",
              "rrtb_cases_ll", "rrtb_cases_ul")
for (col in num_cols) who_data[, (col) := as.numeric(get(col))]

# Calculate RR-TB incidence rate (per 100k)
who_data[, rrtb_incidence := (rrtb_cases / population) * 1e5]
who_data[, rrtb_incidence_ll := (rrtb_cases_ll / population) * 1e5]
who_data[, rrtb_incidence_ul := (rrtb_cases_ul / population) * 1e5]

# Total TB incidence (rename for clarity)
who_data[, total_tb_incidence := dstb_incidence]
who_data[, total_tb_incidence_ll := dstb_incidence_ll]
who_data[, total_tb_incidence_ul := dstb_incidence_ul]

export(who_data, here("data", "who", "who_data_processed.Rdata"))

cat("=== WHO Data Loaded ===\n")
cat("Years:", min(who_data$year), "-", max(who_data$year), "\n")

# =============================================================================
# 2. Calibration Targets
# =============================================================================

calib_years_total <- c(2000, 2005, 2010, 2015, 2019, 2023, 2024)
calib_years_rr <- c(2015, 2017, 2019, 2023, 2024)

# Total TB incidence targets
targets_inc <- list()
for (yr in calib_years_total) {
  row <- who_data[year == yr]
  targets_inc[[paste0("Inc", yr)]] <- c(lo = row$total_tb_incidence_ll,
                                         hi = row$total_tb_incidence_ul)
}

# RR-TB incidence targets (per 100k)
targets_rr_inc <- list()
for (yr in calib_years_rr) {
  row <- who_data[year == yr]
  if (!is.na(row$rrtb_incidence)) {
    targets_rr_inc[[paste0("IncRR", yr)]] <- c(lo = row$rrtb_incidence_ll,
                                                 hi = row$rrtb_incidence_ul)
  }
}

# Proportion of RR-TB that is acquired: 13-28%
targets_acquired <- list(pAc2024 = c(lo = 0.13, hi = 0.28))

# Combine
targets_all <- c(targets_inc, targets_rr_inc, targets_acquired)

# Create targets dataframe for HMER
targetsdb <- data.frame(
  lo = sapply(targets_all, `[`, 1),
  hi = sapply(targets_all, `[`, 2)
)
targetsdb$name <- names(targets_all)
targetsdb$var <- substr(targetsdb$name, 1, 3)
targetsdb$time <- as.numeric(gsub("[^0-9]", "", targetsdb$name))
rownames(targetsdb) <- targetsdb$name

cat("\nCalibration targets:\n")
print(targetsdb[, c("var", "time", "lo", "hi")])

export(targetsdb, here("data", "fit", "targets.Rdata"))
export(targets_all, here("data", "fit", "targets_list.Rdata"))

# =============================================================================
# 3. Key Research Parameters
# =============================================================================

PAF_INHR <- list(estimate = 0.951, lo = 0.856, hi = 0.981)
prop_acquired <- list(estimate = 0.205, lo = 0.13, hi = 0.28)

research_params <- list(PAF_INHR = PAF_INHR, prop_acquired = prop_acquired)
export(research_params, here("data", "fit", "research_params.Rdata"))

cat("\nPAF for INH-R:", PAF_INHR$estimate, "(", PAF_INHR$lo, "-", PAF_INHR$hi, ")\n")
cat("Proportion acquired RR:", prop_acquired$estimate, "(", prop_acquired$lo, "-", prop_acquired$hi, ")\n")

cat("\n=== 01_Data.R completed ===\n")

rm(calib_years_total, calib_years_rr, yr, row, targets_inc, targets_rr_inc, targets_acquired, num_cols, col)
