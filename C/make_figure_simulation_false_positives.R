d = 17
k = 1:160

L = list()

z2 = c(1.066122, 1.065963, 1.065649)
C2 = c(2.167762, 2.164468, 2.156257)

z1 = 1.065511
C1 = 2.152467

for (i in 1:3) {
   L[[i]] = C1/z1^(k+1) - C2[i]/z2[i]^(k+1)
}

S = list(
   read.table("out-.05-fp.txt"),
   read.table("out-.15-fp.txt"),
   read.table("out-.25-fp.txt")
)

pdf("simulp_false_positives.pdf", width=11, height=5.5, useDingbats=FALSE)
par(mfrow=c(1,2))
subs = seq(1,20)
plot(S[[1]][subs,1], S[[1]][subs,2]/1e7, pch=19, cex=.5,
     plot.first=grid(), ylim=c(0,.0035),
     xlab="Read size", ylab="Type I error rate")
lines(k[1:110], L[[1]][1:110])

points(S[[2]][subs,1], S[[2]][subs,2]/1e7, pch=19, cex=.5)
lines(k[1:110], L[[2]][1:110])

points(S[[3]][subs,1], S[[3]][subs,2]/1e7, pch=19, cex=.5)
lines(k[1:110], L[[3]][1:110])

legend(x="topright", inset=0.05, legend=c("Simulation", "Estimate"),
       pch=c(19, NA), lty=c(NA, 1), pt.cex=.6, bg="white", box.col=NA)

subs = seq(21,45)
plot(S[[1]][subs,1], S[[1]][subs,2]/1e8, pch=19, cex=.5,
     plot.first=grid(), ylim=c(0,.0002),
     xlab="Read size", ylab="Type I error rate")
lines(k[90:160], L[[1]][90:160])

points(S[[2]][subs,1], S[[2]][subs,2]/1e8, pch=19, cex=.5)
lines(k[90:160], L[[2]][90:160])

points(S[[3]][subs,1], S[[3]][subs,2]/1e8, pch=19, cex=.5)
lines(k[90:160], L[[3]][90:160])
dev.off()
