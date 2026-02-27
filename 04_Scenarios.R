## RR-TB Vietnam Model - 04: Scenarios
## Evaluate impact of INH-R diagnosis/treatment on RR-TB burden

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
  inhr_low       = list(name = "INH-R Dx Low (20%)",    coverage = c(lo=0.15, hi=0.25), efficacy = c(lo=0.70, hi=0.90)),
  inhr_med       = list(name = "INH-R Dx Medium (50%)", coverage = c(lo=0.40, hi=0.60), efficacy = c(lo=0.70, hi=0.90)),
  inhr_high      = list(name = "INH-R Dx High (80%)",   coverage = c(lo=0.70, hi=0.90), efficacy = c(lo=0.70, hi=0.90)),
  inhr_universal = list(name = "Universal DST (100%)",   coverage = c(lo=0.95, hi=1.00), efficacy = c(lo=0.85, hi=0.95))
)

# ODE with intervention (13 compartments, no CLE/REC)
ode_2strain_intervention <- function(parms, intervention = NULL, interv_start = 2025, end_time = 2035) {
  N <- 1e5; mu <- 0.013; delta_ds <- 2
  if (is.null(intervention) || identical(intervention$coverage, 0)) {
    interv_eff <- 1
  } else {
    cov <- runif(1, intervention$coverage['lo'], intervention$coverage['hi'])
    eff <- runif(1, intervention$efficacy['lo'], intervention$efficacy['hi'])
    interv_eff <- 1 - cov * eff * runif(1, PAF$lo, PAF$hi)
  }
  make_forcer <- function(p_ini, p_fin) {
    approxfun(c(1900,1999,2024), c(parms[p_ini],parms[p_ini],parms[p_fin]), method="linear", rule=2)
  }
  force_theta_ds <- make_forcer('theta_ds_ini','theta_ds_fin')
  force_phi_ds   <- make_forcer('phi_ds_ini','phi_ds_fin')
  force_mutb_ds  <- make_forcer('mutb_ds_ini','mutb_ds_fin')
  force_theta_rr <- make_forcer('theta_rr_ini','theta_rr_fin')
  force_phi_rr   <- make_forcer('phi_rr_ini','phi_rr_fin')
  force_mutb_rr  <- make_forcer('mutb_rr_ini','mutb_rr_fin')
  force_delta_rr <- approxfun(c(1900,2020,2024), c(0.67,0.67,2.00), method="linear", rule=2)
  
  des <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      infectious_ds <- kappa*SUB_DS+CLN_DS; infectious_rr <- fitness_rr*(kappa*SUB_RR+CLN_RR)
      lambda_ds <- beta/N*infectious_ds; lambda_rr <- beta/N*infectious_rr; lambda_total <- lambda_ds+lambda_rr
      theta_ds_t <- force_theta_ds(time); phi_ds_t <- force_phi_ds(time); mutb_ds_t <- force_mutb_ds(time)
      theta_rr_t <- force_theta_rr(time); phi_rr_t <- force_phi_rr(time); mutb_rr_t <- force_mutb_rr(time)
      delta_rr_t <- force_delta_rr(time)
      ie <- ifelse(time >= interv_start, interv_eff, 1)
      acq_fail_t <- acq_fail*ie; acq_relapse_t <- acq_relapse*ie
      
      dSUS <- mu*N + omega_c*(INF_DS+INF_RR) + mu_r*(MIN_DS+MIN_RR) + tau*(TRE_DS+TRE_RR) - (lambda_total+mu)*SUS
      dINF_DS <- lambda_ds*SUS - (omega_c+omega_m+omega_s+mu)*INF_DS
      dMIN_DS <- omega_m*INF_DS + sigma_m*SUB_DS - (mu_r+mu_s+mu)*MIN_DS
      dSUB_DS <- omega_s*INF_DS + mu_s*MIN_DS + gamma*CLN_DS + relapse_ds*(1-acq_relapse_t)*TRE_DS - (sigma_m+sigma_c+mu)*SUB_DS
      dCLN_DS <- sigma_c*SUB_DS - (gamma+theta_ds_t+mutb_ds_t+mu)*CLN_DS + phi_ds_t*(1-acq_fail_t)*TXT_DS
      dTXT_DS <- theta_ds_t*CLN_DS - (phi_ds_t+delta_ds+mu)*TXT_DS
      dTRE_DS <- delta_ds*TXT_DS - (relapse_ds+tau+mu)*TRE_DS
      acq_f <- phi_ds_t*acq_fail_t*TXT_DS; acq_r <- relapse_ds*acq_relapse_t*TRE_DS
      dINF_RR <- lambda_rr*SUS - (omega_c+omega_m+omega_s+mu)*INF_RR
      dMIN_RR <- omega_m*INF_RR + sigma_m*SUB_RR - (mu_r+mu_s+mu)*MIN_RR
      dSUB_RR <- omega_s*INF_RR + mu_s*MIN_RR + gamma*CLN_RR + relapse_rr*TRE_RR + acq_r - (sigma_m+sigma_c+mu)*SUB_RR
      dCLN_RR <- sigma_c*SUB_RR - (gamma+theta_rr_t+mutb_rr_t+mu)*CLN_RR + phi_rr_t*TXT_RR + acq_f
      dTXT_RR <- theta_rr_t*CLN_RR - (phi_rr_t+delta_rr_t+mu)*TXT_RR
      dTRE_RR <- delta_rr_t*TXT_RR - (relapse_rr+tau+mu)*TRE_RR
      dSUB_RR_acq <- acq_r - (sigma_m+sigma_c+mu)*SUB_RR_acq
      
      Pop <- SUS+INF_DS+MIN_DS+SUB_DS+CLN_DS+TXT_DS+TRE_DS+INF_RR+MIN_RR+SUB_RR+CLN_RR+TXT_RR+TRE_RR
      Inc_DS <- sigma_c*SUB_DS; Inc_RR_acq <- acq_f+sigma_c*SUB_RR_acq; Inc_RR_trans <- sigma_c*(SUB_RR-SUB_RR_acq); Inc_RR <- sigma_c*SUB_RR+acq_f
      return(list(c(dSUS, dINF_DS,dMIN_DS,dSUB_DS,dCLN_DS,dTXT_DS,dTRE_DS,
                    dINF_RR,dMIN_RR,dSUB_RR,dCLN_RR,dTXT_RR,dTRE_RR,dSUB_RR_acq),
        Inc=((Inc_DS+Inc_RR)/Pop)*1e5, IncDS=(Inc_DS/Pop)*1e5, IncRR=(Inc_RR/Pop)*1e5,
        IncRR_acq=(Inc_RR_acq/Pop)*1e5, IncRR_trans=(Inc_RR_trans/Pop)*1e5,
        pAc=Inc_RR_acq/(Inc_RR+1e-10), pRR=Inc_RR/(Inc_DS+Inc_RR+1e-10),
        Prev=((SUB_DS+CLN_DS+SUB_RR+CLN_RR)/Pop)*1e5, PrevRR=((SUB_RR+CLN_RR)/Pop)*1e5,
        Mor=((mutb_ds_t*CLN_DS+mutb_rr_t*CLN_RR)/Pop)*1e5, MorRR=(mutb_rr_t*CLN_RR/Pop)*1e5,
        Notif=((theta_ds_t*CLN_DS+theta_rr_t*CLN_RR)/Pop)*1e5))
    })
  }
  rr0 <- 0.025
  yini <- c(SUS=N-500, INF_DS=0,MIN_DS=0,SUB_DS=100*(1-rr0),CLN_DS=400*(1-rr0),TXT_DS=0,TRE_DS=0,
            INF_RR=0,MIN_RR=0,SUB_RR=100*rr0,CLN_RR=400*rr0,TXT_RR=0,TRE_RR=0,SUB_RR_acq=0)
  deSolve::ode(yini, seq(1900, end_time, 1), des, parms, method="lsoda")
}

cat("\n=== Running Scenarios ===\n")
all_results <- list()
for (sname in names(interventions)) {
  cat("  ", interventions[[sname]]$name, "\n")
  pb <- progress_bar$new(format="[:bar] :percent :eta", total=nrow(fitpts))
  res <- lapply(1:nrow(fitpts), function(i) {
    p <- as.numeric(fitpts[i,]); names(p) <- names(fitpts)
    out <- as.data.frame(ode_2strain_intervention(p, interventions[[sname]], 2025, 2035))
    out$run <- i; out$scenario <- sname; out$scenario_name <- interventions[[sname]]$name; pb$tick(); out
  })
  all_results[[sname]] <- rbindlist(res)
}
results_all <- rbindlist(all_results)[time >= 2020]
export(results_all, here("outputs", "results", "scenarios_raw.Rdata"))

output_vars <- c("Inc","IncRR","IncRR_acq","IncRR_trans","pAc","pRR","Prev","PrevRR","Mor","MorRR")
results_summary <- rbindlist(lapply(output_vars, function(v) {
  results_all[, .(lo=quantile(get(v),0.025), q25=quantile(get(v),0.25), median=quantile(get(v),0.50),
                  q75=quantile(get(v),0.75), hi=quantile(get(v),0.975), var=v), by=.(scenario, scenario_name, time)]
}))
export(results_summary, here("outputs", "results", "scenarios_summary.Rdata"))

bl <- results_all[scenario=="baseline" & time==2035, .(run, IncRR_base=IncRR, MorRR_base=MorRR, IncRR_acq_base=IncRR_acq)]
impact_all <- rbindlist(lapply(setdiff(names(interventions), "baseline"), function(scen) {
  sc <- results_all[scenario==scen & time==2035, .(run, IncRR_int=IncRR, MorRR_int=MorRR, IncRR_acq_int=IncRR_acq, scenario_name)]
  imp <- merge(bl, sc, by="run")
  imp[, `:=`(IncRR_reduction=(IncRR_base-IncRR_int)/IncRR_base*100,
             IncRR_acq_reduction=(IncRR_acq_base-IncRR_acq_int)/(IncRR_acq_base+1e-10)*100,
             MorRR_reduction=(MorRR_base-MorRR_int)/(MorRR_base+1e-10)*100,
             IncRR_averted=IncRR_base-IncRR_int, MorRR_averted=MorRR_base-MorRR_int, scenario=scen)]
}))
impact_summary <- impact_all[, .(scenario_name=first(scenario_name),
  IncRR_reduction_median=median(IncRR_reduction), IncRR_reduction_lo=quantile(IncRR_reduction,0.025), IncRR_reduction_hi=quantile(IncRR_reduction,0.975),
  IncRR_averted_median=median(IncRR_averted), IncRR_averted_lo=quantile(IncRR_averted,0.025), IncRR_averted_hi=quantile(IncRR_averted,0.975)), by=scenario]
export(impact_all, here("outputs", "results", "impact_raw.Rdata"))
export(impact_summary, here("outputs", "results", "impact_summary.Rdata"))

cat("\n===============================================\nIMPACT SUMMARY (by 2035)\n===============================================\n")
bl_sum <- results_all[scenario=="baseline" & time==2035, .(m=median(IncRR), lo=quantile(IncRR,0.025), hi=quantile(IncRR,0.975))]
cat(sprintf("Baseline RR-TB incidence: %.2f (%.2f-%.2f) per 100k\n\n", bl_sum$m, bl_sum$lo, bl_sum$hi))
for (i in 1:nrow(impact_summary)) {
  r <- impact_summary[i]
  cat(sprintf("%s:\n  RR-TB reduction: %.1f%% (%.1f-%.1f%%)\n  Cases averted: %.2f (%.2f-%.2f) per 100k\n\n",
              r$scenario_name, r$IncRR_reduction_median, r$IncRR_reduction_lo, r$IncRR_reduction_hi,
              r$IncRR_averted_median, r$IncRR_averted_lo, r$IncRR_averted_hi))
}
cat("=== 04_Scenarios.R completed ===\n")
