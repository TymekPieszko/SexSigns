rm(list = ls())

library(phangorn)
library(vegan)
library(abc)

setwd("/home/tymekpieszko/projects/RotSpec/abc")

n_inds <- 8 # 10 or 8, depending on whether it's run for Silwood or all samples
mode <- if (n_inds == 10) "total" else "silwood"
stats <- c("kappa","CL","AR","SX","delta")
tol <- 0.05
params <- paste0("inds~",n_inds,"_b~0.15_d~15")

# Simulated values
sim <- read.table(paste0("simulated/GC_stats_",params,".tsv"),header=T)

# Observed values
read.obs<-function(sp) {
    tmp<-read.table(paste0("observed/",sp,"_",mode,"_stats.tsv"),header=T)
    colMeans(tmp[,stats])
}
obs <- rbind(read.obs("Rmag"),read.obs("Rsoc"))
rownames(obs)<-c("Rmag","Rsoc")
print(obs)
write.table(obs, file=paste0("observed/abc_",mode,"_stats.tsv"), sep="\t")

# Run ABC using simple rejection method

params.mat<-sim[,1:2]
params.mat[,2]<-params.mat[,2]*1000000
print(params.mat)

res.mag <- abc(target=obs["Rmag",],
    param=params.mat, names=c("sex","rec"), sumstat=sim[,stats], tol=tol, method="rejection")
summary(res.mag )

res.soc <- abc(target=obs["Rsoc",],
    param=params.mat, names=c("sex","rec"), sumstat=sim[,stats], tol=tol, method="rejection")
summary(res.soc)

results<-cbind(summary(res.mag, print=F)[,1], summary(res.soc, print=F)[,1],summary(res.mag, print=F)[,2]/1000000, summary(res.soc, print=F)[,2]/1000000)
colnames(results)<-c("sex_rmag","sex_rsoc","rec_rmag","rec_rsoc")
results[, 3:4] <- results[, 3:4] * 5000
write.table(format(results, digits = 3, scientific = TRUE), file = paste0("Rotaria_inds~",n_inds,".tsv"), sep = "\t", quote = FALSE, col.names = NA)
classes <- function(res) table(cut(res$unadj.values[,"sex"], c(-Inf, 0, 0.001, Inf), labels = c("obligate", "rare", "frequent"))) / nrow(res$unadj.values)
write.table(rbind(Rmag = classes(res.mag), Rsoc = classes(res.soc)), file = paste0("Rotaria_inds~",n_inds,".classes.tsv"), sep = "\t", quote = FALSE, col.names = NA)

# Visualise PCA

pca.mat<-rbind(sim[,stats], obs)
pcas<-prcomp(pca.mat)
pve<-round(100*summary(pcas)$importance[2,1:2],1)

png(paste0("pca_",mode,".png"), width = 18, height = 13.7, units = "in", res = 150, pointsize = 36)
par(mai = c(2.5, 2.5, 0.3, 0.3))
palette(c("#0b4f8a","lightgreen","red"))
my.cols<-rep(1,nrow(pcas$x))
my.cols[sim[,1]>0]<-2
my.cols[sim[,1]>0.001]<-3
my.cols<-rev(c(my.cols,NA,NA))
plot(rev(pcas$x[,1]),
     rev(pcas$x[,2]),
     xlab=paste0("PC1 (",pve[1],"%)"),
     ylab=paste0("PC2 (",pve[2],"%)"),
     xlim=range(pcas$x[,1])+c(0, 0.5),
     ylim=range(pcas$x[,2])+c(-0.33, 0.05)*diff(range(pcas$x[,2]))
     ,col=my.cols,
     pch=20,
     cex.axis = 1.2,
     cex.lab = 1.2,
     cex=0.35,
     xaxt = "n",
     yaxt = "n")
axis(1, at = seq(floor(min(pcas$x[,1])), ceiling(max(pcas$x[,1])), by = 0.4),
     cex.axis = 1.3)
axis(2, at = seq(floor(min(pcas$x[,2])), ceiling(max(pcas$x[,2])), by = 0.4),
     cex.axis = 1.3)
points(pcas$x[(nrow(pcas$x)-1):nrow(pcas$x), 1],
       pcas$x[(nrow(pcas$x)-1):nrow(pcas$x), 2],
       pch = 21,                # circle
       bg = c(rgb(75, 41, 145, maxColorValue = 255), rgb(231, 90, 155, maxColorValue = 255)),  # fill colours
       col = "black",           # border colour
       cex = 1.7)
legend("topright",inset=0.01,cex=1.05,fill=1:3,legend=c("Obligate\nasex.","Rare\nsex","Freq.\nsex"),bty="n",y.intersp = 1.6)
tab <- results[c("2.5% Perc.:", "Median:", "97.5% Perc.:"), ] * 1000
usr <- par("usr")
lab <- formatC(tab, format = "f", digits = 3)
cw <- c(max(strwidth(rownames(tab), cex = 1.1)), apply(lab, 2, function(x) max(strwidth(x, cex = 1.1)))) + strwidth("00", cex = 1.1)
xs <- usr[2] - 0.01 * (usr[2] - usr[1]) - c(rev(cumsum(rev(cw[-1]))), 0)
ys <- usr[3] + 0.03 * (usr[4] - usr[3]) + (nrow(tab):0) * 1.8 * strheight("M", cex = 1.1)
rect(xs[1] - cw[1], usr[3], usr[2], ys[1] + strheight("M", cex = 1.1), col = "white", border = NA)
text(xs[-1], ys[1], expression(sigma[RM], sigma[RS], gamma[RM], gamma[RS]), adj = 1, cex = 1.1)
text(xs[1], ys[1], expression("(" * N^-1 * ")"), adj = 1, cex = 1.05)
segments(xs[1] - cw[1], mean(ys[1:2]), xs[5], mean(ys[1:2]))
text(xs[1], ys[-1], rownames(tab), adj = 1, cex = 1.1)
for (j in 1:4) text(xs[j + 1], ys[-1], lab[, j], adj = 1, cex = 1.1)
box()
dev.off()
# Visualise accepted scenarios

sex.lvls<-sort(unique(sim$sex))
rec.lvls<-sort(unique(sim$rec))
acc.tab<-function(res) table(factor(sim$sex[res$region],levels=sex.lvls),factor(sim$rec[res$region],levels=rec.lvls))
accs<-list(Rmag=acc.tab(res.mag),Rsoc=acc.tab(res.soc))
sp.names<-c(Rmag="R. magnacalcarata",Rsoc="R. socialis")
zmax<-max(unlist(accs))

png(paste0("accepted_",mode,".png"), width = 16, height = 8, units = "in", res = 150, pointsize = 20)
par(mfrow = c(1, 2), mai = c(1.5, 1.6, 0.7, 0.3))
for (sp in names(accs)) {
    acc<-accs[[sp]]
    image(1:length(sex.lvls),
          1:length(rec.lvls),
          acc,
          col=grey.colors(100, start=1, end=0),
          zlim=c(0,zmax),
          xlab="",
          ylab="",
          main=bquote(italic(.(sp.names[[sp]]))),
          xaxt="n",
          yaxt="n")
    axis(1, at=1:length(sex.lvls), labels=FALSE)
    text(1:length(sex.lvls), par("usr")[3]-0.05*diff(par("usr")[3:4]), signif(sex.lvls*1000,3), srt=45, adj=1, xpd=TRUE, cex=0.8)
    axis(2, at=1:length(rec.lvls), labels=signif(rec.lvls*5000*1000,3), las=1, cex.axis=0.8)
    mtext(expression(sigma * " (" * N^-1 * ")"), side=1, line=3.5)
    mtext(expression(gamma * " (" * N^-1 * ")"), side=2, line=3.5)
    nz<-which(acc>0, arr.ind=TRUE)
    text(nz[,1], nz[,2], acc[nz], cex=0.6, col=ifelse(acc[nz]>zmax/2,"white","black"))
    box()
}
dev.off()
