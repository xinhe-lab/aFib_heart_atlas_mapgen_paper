# setting
pval.thr <- 1e-3 # p-value threshold for calling a SNP eQTL
z.thr <- qnorm(1-pval.thr/2)
n <- 500 # sample size

# power analysis
# Let z.hat be observed Z-score of a SNP, and z be its true effect (Z-score scale), the power analysis is based on the distribution: z.hat | z ~ N(z, 1). 
PVE <- seq(0, 0.1, 0.1/100) # true effect size of SNPs at the PVE scale
z <- sqrt(PVE * n) # true effect size of SNPs at the Z-score scale
power <- 1 - pnorm(z.thr, z, 1)
#power.reduced <- 1 - pnorm(z.thr, z*pi, 1)

# plotting
colors <- c("red", "blue", "green")
pdf(file="manuscript_figures/eQTL_power.pdf",  height=6, width=8)

plot(PVE, power, type="l", col=colors[1])
power.reduced <- 1 - pnorm(z.thr, z* 0.5, 1)
lines(PVE, power.reduced, type="l", col=colors[2])
power.reduced <- 1 - pnorm(z.thr, z* 0.3, 1)
lines(PVE, power.reduced, type="l", col=colors[3])
legend("topleft",c("pi = 1","pi = 0.5", "pi = 0.3"),col=colors,lty=1, cex=1.0)

dev.off()