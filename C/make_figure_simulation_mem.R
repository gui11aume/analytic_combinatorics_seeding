d = 17
k = 1:210

L = list()

z1 = c(1.037418, 1.026633, 1.019079)
C1 = c(1.569487, 1.389335, 1.272535)


for (i in 1:3) {
   L[[i]] = 1-C1[i]/z1[i]^(k+1)
}

S = list(
   read.table("out-.08-MEM.txt"),
   read.table("out-.10-MEM.txt"),
   read.table("out-.12-MEM.txt")
)

pdf("simulp_mem.pdf", width=11, height=5.5, useDingbats=FALSE)
par(mfrow=c(1,2))
subs = seq(1,51,2)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e7, pch=19, cex=.5, ylim=c(0.5,1),
     plot.first=grid(),
     xlab="Read size", ylab="MEM seeding probability")
lines(k[1:154], L[[1]][1:154])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e7, pch=19, cex=.5)
lines(k[1:154], L[[2]][1:154])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e7, pch=19, cex=.5)
lines(k[1:154], L[[3]][1:154])

legend(x="bottomright", inset=0.05, legend=c("Simulation", "Estimate"),
       pch=c(19, NA), lty=c(NA, 1), pt.cex=.6, bg="white", box.col=NA)


subs = seq(51,76)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e7, pch=19, cex=.5, ylim=c(0.9,1),
     plot.first=grid(),
     xlab="Read size", ylab="MEM seeding probability")
lines(k[146:204], L[[1]][146:204])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e7, pch=19, cex=.5)
lines(k[146:204], L[[2]][146:204])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e7, pch=19, cex=.5)
lines(k[146:204], L[[3]][146:204])
dev.off()
