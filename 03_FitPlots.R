## RR-TB Vietnam Model - 03: Fit Plots
## Visualize calibration results

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

# Model (13 compartments, no CLE/REC)
ode_2strain <- function(parms, end_time = 2024) {
  N <- 1e5; mu <- 0.013; delta_ds <- 2
  make_forcer <- function(p_ini, p_fin) {
    approxfun(c(1900,1999,2024), c(parms[p_ini],parms[p_ini],parms[p_fin]), method="linear", rule=2)
  }
  force_theta_ds <- make_forcer('theta_ds_ini','theta_ds_fin')
  force_phi_ds <- make_forcer('phi_ds_ini','phi_ds_fin')
  force_mutb_ds <- make_forcer('mutb_ds_ini','mutb_ds_fin')
  force_theta_rr <- make_forcer('theta_rr_ini','theta_rr_fin')
  force_phi_rr <- make_forcer('phi_rr_ini','phi_rr_fin')
  force_mutb_rr <- make_forcer('mutb_rr_ini','mutb_rr_fin')
  force_delta_rr <- approxfun(c(1900,2020,2024), c(0.67,0.67,2.00), method="linear", rule=2)
  des <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      infectious_ds <- kappa*SUB_DS+CLN_DS; infectious_rr <- fitness_rr*(kappa*SUB_RR+CLN_RR)
      lambda_ds <- beta/N*infectious_ds; lambda_rr <- beta/N*infectious_rr; lambda_total <- lambda_ds+lambda_rr
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
      Pop <- SUS+INF_DS+MIN_DS+SUB_DS+CLN_DS+TXT_DS+TRE_DS+INF_RR+MIN_RR+SUB_RR+CLN_RR+TXT_RR+TRE_RR
      Inc_DS <- sigma_c*SUB_DS; Inc_RR_acq <- acq_f+sigma_c*SUB_RR_acq; Inc_RR <- sigma_c*SUB_RR+acq_f
      return(list(c(dSUS,dINF_DS,dMIN_DS,dSUB_DS,dCLN_DS,dTXT_DS,dTRE_DS,
                    dINF_RR,dMIN_RR,dSUB_RR,dCLN_RR,dTXT_RR,dTRE_RR,dSUB_RR_acq),
        Inc=((Inc_DS+Inc_RR)/Pop)*1e5, IncRR=(Inc_RR/Pop)*1e5, pAc=Inc_RR_acq/(Inc_RR+1e-10),
        pRR=Inc_RR/(Inc_DS+Inc_RR+1e-10), Prev=((SUB_DS+CLN_DS+SUB_RR+CLN_RR)/Pop)*1e5,
        Mor=((mutb_ds_t*CLN_DS+mutb_rr_t*CLN_RR)/Pop)*1e5, Notif=((theta_ds_t*CLN_DS+theta_rr_t*CLN_RR)/Pop)*1e5))
    })
  }
  rr0 <- 0.025
  yini <- c(SUS=N-500, INF_DS=0,MIN_DS=0,SUB_DS=100*(1-rr0),CLN_DS=400*(1-rr0),TXT_DS=0,TRE_DS=0,
            INF_RR=0,MIN_RR=0,SUB_RR=100*rr0,CLN_RR=400*rr0,TXT_RR=0,TRE_RR=0,SUB_RR_acq=0)
  deSolve::ode(yini, seq(1900, end_time, 1), des, parms, method="lsoda")
}

# Run model
cat("Running model...\n")
pb <- progress_bar$new(format="[:bar] :percent :eta", total=nrow(fitpts))
all_results <- lapply(1:nrow(fitpts), function(i) {
  p <- as.numeric(fitpts[i,]); names(p) <- names(fitpts)
  out <- as.data.frame(ode_2strain(p)); out$run <- i; pb$tick(); out
})
results_raw <- rbindlist(all_results)[time >= 2000]

results_long <- melt(results_raw, id.vars=c("time","run"), measure.vars=c("Inc","IncRR","pAc","pRR","Prev","Mor","Notif"))
results <- as.data.table(results_long)[, .(`2.5%`=quantile(value,0.025), `25%`=quantile(value,0.25),
  `50%`=quantile(value,0.50), `75%`=quantile(value,0.75), `97.5%`=quantile(value,0.975)), by=.(time, variable)]
setnames(results, "variable", "var")
export(results, here("outputs", "results.Rdata"))

# Plots
cl <- "#CE2931"; al <- 0.2; lw <- 0.6; brks <- seq(2000,2024,2); lims <- c(1999,2025); yrs <- 2000:2024
theme_fit <- theme(panel.background=element_blank(), panel.grid.major.y=element_line(colour="grey90",linewidth=0.3),
  panel.grid.major.x=element_blank(), axis.line=element_line(linewidth=lw), axis.ticks=element_line(linewidth=lw),
  axis.title.x=element_text(vjust=-1,margin=margin(b=10),size=11), axis.title.y=element_text(vjust=2,margin=margin(l=10),size=11),
  axis.text.x=element_text(size=9,angle=45,hjust=1), axis.text.y=element_text(size=10),
  plot.title=element_text(size=13,hjust=0.5,face="bold"), legend.position="none")

targets_inc <- targets[grepl("^Inc[0-9]",name)][,var:="Inc"]
targets_incrr <- targets[grepl("^IncRR",name)][,var:="IncRR"]
targets_pac <- targets[grepl("^pAc",name)][,var:="pAc"]

rp <- function(d, tg, yb, yl, ttl, ylab, pct=FALSE) {
  p <- ggplot(d, aes(time)) + geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),fill=cl,alpha=al) +
    geom_ribbon(aes(ymin=`25%`,ymax=`75%`),fill=cl,alpha=al*2) + geom_line(aes(y=`50%`),linewidth=lw,colour=cl) +
    geom_errorbar(data=tg,aes(x=time,ymin=lo,ymax=hi),linewidth=lw,width=0.5) +
    scale_x_continuous(breaks=brks,limits=lims,expand=c(0,0)) + coord_cartesian(ylim=yl) +
    labs(title=ttl,x="Year",y=ylab) + theme_fit
  if (pct) p <- p + scale_y_continuous(breaks=yb,labels=percent,expand=c(0,0)) else p <- p + scale_y_continuous(breaks=yb,expand=c(0,0))
  p
}

a <- rp(results[var=="Inc" & time %in% yrs], targets_inc, seq(0,400,50), c(0,400), "A. Total TB Incidence", "Per 100k")
b <- rp(results[var=="IncRR" & time %in% yrs], targets_incrr, seq(0,20,5), c(0,20), "B. RR-TB Incidence", "Per 100k")
cc <- rp(results[var=="pAc" & time %in% yrs], targets_pac, seq(0,0.5,0.1), c(0,0.5), "C. Proportion Acquired RR", "Proportion", TRUE)
dd <- ggplot(results[var=="pRR" & time %in% yrs], aes(time)) +
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),fill=cl,alpha=al) + geom_ribbon(aes(ymin=`25%`,ymax=`75%`),fill=cl,alpha=al*2) +
  geom_line(aes(y=`50%`),linewidth=lw,colour=cl) + geom_hline(yintercept=0.04,linetype="dashed",colour="grey50") +
  scale_x_continuous(breaks=brks,limits=lims,expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0,0.10,0.02),labels=percent,expand=c(0,0)) + coord_cartesian(ylim=c(0,0.10)) +
  labs(title="D. RR Proportion",x="Year",y="Proportion") + theme_fit

png(here("plots","01_calibration_fit.png"), width=10,height=8,units='in',res=300); print((a|b)/(cc|dd)); dev.off()

key <- c("beta","fitness_rr","acq_fail","acq_relapse","theta_ds_fin","theta_rr_fin")
fl <- melt(as.data.table(fitpts), measure.vars=names(fitpts))
png(here("plots","03_parameters.png"), width=10,height=6,units='in',res=300)
print(ggplot(fl[variable %in% key], aes(value)) + geom_histogram(bins=15,fill=cl,alpha=0.7,colour="white") +
  facet_wrap(~variable,scales="free",ncol=3) + labs(title="Calibrated Parameters",x="Value",y="Count") +
  theme_fit + theme(strip.text=element_text(face="bold"))); dev.off()

cat("\nCalibration (2024):\n")
for (v in c("Inc","IncRR","pAc","pRR")) { r <- results[time==2024 & var==v]; cat(sprintf("  %s: %.2f (%.2f-%.2f)\n",v,r$`50%`,r$`2.5%`,r$`97.5%`)) }
cat("=== 03_FitPlots.R completed ===\n")
