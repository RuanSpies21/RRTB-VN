## RR-TB Vietnam Model - 06: Plots
## Publication figures for intervention scenarios

library(rio)
library(here)
library(ggplot2)
library(data.table)
library(patchwork)

results_summary <- as.data.table(import(here("outputs", "results", "scenarios_summary.Rdata")))
impact_summary <- as.data.table(import(here("outputs", "results", "impact_summary.Rdata")))

# Settings
sl <- c(baseline="Baseline", inhr_low="INH-R Dx Low (20%)", inhr_med="INH-R Dx Medium (50%)",
        inhr_high="INH-R Dx High (80%)", inhr_universal="Universal DST (100%)")
sc <- c(baseline="#828282", inhr_low="#FFBC47", inhr_med="#2984CE", inhr_high="#CE2931", inhr_universal="#1D2D5F")
so <- c("baseline", "inhr_low", "inhr_med", "inhr_high", "inhr_universal")
theme_rr <- theme_bw() + theme(legend.position="bottom", legend.direction="horizontal",
  plot.margin=margin(10,15,5,10), panel.spacing=unit(0.8,"lines"))

# A. RR-TB Incidence
d1 <- results_summary[var=="IncRR" & time>=2020]; d1$scenario <- factor(d1$scenario, so)
p1 <- ggplot(d1, aes(time, colour=scenario, fill=scenario)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, colour=NA) +
  geom_line(aes(y=median), linewidth=0.8) +
  geom_vline(xintercept=2025, linetype="dashed", colour="grey50") +
  scale_colour_manual(values=sc, labels=sl, name=NULL) +
  scale_fill_manual(values=sc, labels=sl, name=NULL) +
  scale_x_continuous(breaks=seq(2020,2035,5), limits=c(2020,2035), expand=c(0,0)) +
  labs(title="A. RR-TB Incidence", x="Year", y="Per 100k") + theme_rr +
  guides(colour=guide_legend(nrow=2), fill=guide_legend(nrow=2))

# B. Proportion Acquired
d2 <- results_summary[var=="pAc" & time>=2020]; d2$scenario <- factor(d2$scenario, so)
p2 <- ggplot(d2, aes(time, colour=scenario, fill=scenario)) +
  geom_ribbon(aes(ymin=lo*100, ymax=hi*100), alpha=0.15, colour=NA) +
  geom_line(aes(y=median*100), linewidth=0.8) +
  geom_vline(xintercept=2025, linetype="dashed", colour="grey50") +
  scale_colour_manual(values=sc, labels=sl, name=NULL) +
  scale_fill_manual(values=sc, labels=sl, name=NULL) +
  scale_x_continuous(breaks=seq(2020,2035,5), limits=c(2020,2035), expand=c(0,0)) +
  labs(title="B. Proportion RR-TB Acquired", x="Year", y="%") + theme_rr +
  guides(colour=guide_legend(nrow=2), fill=guide_legend(nrow=2))

# C. Impact Bar Chart
imp <- impact_summary[, .(scenario, scenario_name, y=IncRR_reduction_median,
                           ylo=IncRR_reduction_lo, yhi=IncRR_reduction_hi)]
imp$scenario <- factor(imp$scenario, so[-1])
p3 <- ggplot(imp, aes(scenario, y, fill=scenario)) +
  geom_col(width=0.7) + geom_errorbar(aes(ymin=ylo, ymax=yhi), width=0.2) +
  geom_text(aes(label=sprintf("%.1f%%", y)), vjust=-0.5, size=4, fontface="bold") +
  scale_fill_manual(values=sc[-1], guide="none") +
  scale_x_discrete(labels=c(inhr_low="Low\n(20%)", inhr_med="Medium\n(50%)", inhr_high="High\n(80%)", inhr_universal="Universal\n(100%)")) +
  scale_y_continuous(expand=c(0,0)) + coord_cartesian(ylim=c(0, max(imp$yhi)*1.3)) +
  labs(title="C. RR-TB Incidence Reduction by 2035", x="INH-R Dx Coverage", y="Reduction (%)") + theme_rr

# D. Remaining burden
bl_val <- results_summary[scenario=="baseline" & var=="IncRR" & time==2035, median]
wf <- data.table(
  cat = factor(c("Baseline", "Low (20%)", "Medium (50%)", "High (80%)", "Universal (100%)"),
               levels = c("Baseline", "Low (20%)", "Medium (50%)", "High (80%)", "Universal (100%)")),
  val = c(bl_val, bl_val - impact_summary$IncRR_averted_median)
)
p4 <- ggplot(wf, aes(cat, val, fill=ifelse(cat=="Baseline","a","b"))) +
  geom_col(width=0.7) + geom_text(aes(label=sprintf("%.1f",val)), vjust=-0.5, size=3.5, fontface="bold") +
  scale_fill_manual(values=c(a="#828282", b="#CE2931"), guide="none") +
  scale_y_continuous(expand=c(0,0)) + coord_cartesian(ylim=c(0, bl_val*1.3)) +
  labs(title="D. RR-TB Incidence 2035", x="Coverage", y="Per 100k") + theme_rr

# Combined
combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(title="Impact of INH-R Diagnosis and Treatment on RR-TB in Vietnam",
    subtitle="Model projections 2020-2035",
    theme=theme(plot.title=element_text(size=14, face="bold"), plot.subtitle=element_text(size=11, color="grey40")))

png(here("plots", "combined_figure.png"), width=14, height=10, units='in', res=300)
print(combined); dev.off()

# Individual panels
for (nm in c("p1","p2","p3","p4")) {
  png(here("plots", paste0("panel_", nm, ".png")), width=8, height=6, units='in', res=300)
  print(get(nm)); dev.off()
}

cat("Figures saved to plots/\n=== 06_Plots.R completed ===\n")
