library(abc)
rm(list = ls())
# Uncomment the target range
# Range 1
# rec_min <- 0.0
# rec_max <- 1e-06
# Range 2
# rec_min <- 3.16e-09
# rec_max <- 3.16e-07
# Range 3
# rec_min <- 3.16e-08
# rec_max <- 3.16e-07
# Range 4 - rec generally high
# rec_min <- 1e-08
# rec_max <- 1e-07
# Range 5 - rec generally high
rec_min <- 1e-07
rec_max <- 1e-06
sim <- read.table("./GC_stats.tsv", header=T)
unique(sim[,"rec"])
sim <- subset(sim, rec >= rec_min & rec <= rec_max)
param <- sim[, "sex"]+3.162e-06
sim <- subset(sim, select = -c(sex, rec)) # Remove sex and rec
statsets <- list(
  "H_i"         = "H_i",
  "F_is"        = "F_is",
  "r2"          = "r2",
  "kappa"       = "kappa",
  "delta_bias" = "delta_bias",
  "delta"       = "delta",
  "old"         = c("H_i", "F_is", "r2"),
  "new_1"       = c("kappa", "delta_bias"),
  "new_2"       = c("kappa", "delta")
)
statsets
sim<-sim[names(statsets)[1:6]]

# Select simulated and pseudo-observed data
o <- 10 # Number of pseudo-observed
oidx <- as.vector(sapply(seq(1, nrow(sim), 100), function(x)
  x:(x + o -1))) # First o per param combination
obs <- sim[oidx, ]
sim <- sim[-oidx, ]
obspar <- param[oidx]
param <- param[-oidx]

res<-abc(target=obs[1,],param=log10(param),sumstat=sim,tol=0.01,method="rejection")
as.vector(summary(res))[c(4,2,6)]

coeffs <- c()  
for (statset in statsets) {
  predict.full<-matrix(NA,nrow=nrow(obs),ncol=3)
  for (i in (1:nrow(obs))) {
    res <- abc(target=obs[i,statset],param=log10(param),sumstat=sim[,statset],tol=0.01,method="rejection")
    predict.full[i,]<-as.vector(summary(res))[c(4,2,6)]
  }
  # plot(log10(unique(obspar)),by(predict.full[,1],FUN=mean,obspar))
  # plot(predict.full[,1] ~ log10(obspar))
  m <- lm(predict.full[,1] ~ log10(obspar))
  coeff <- summary(m)$r.squared
  coeffs <- c(coeffs, summary(m)$r.squared)
}
names(coeffs) <- names(statsets)
coeffs <- as.data.frame(t(coeffs))

dir.create("abc_out", showWarnings = FALSE)
out_file <- paste0(
  "out/abc_results/abc_R2_rec_",
  format(rec_min, scientific = TRUE),
  "_to_",
  format(rec_max, scientific = TRUE),
  ".tsv"
)
write.table(coeffs, file=out_file, sep="\t", row.names=F, col.names=T, quote=F)
coeffs
