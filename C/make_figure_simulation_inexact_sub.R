
d = 17
k = 1:204

L = list()

# p = 0.08, 0.10, 0.12
z = c(1.10506131623070, 1.07924416591384, 1.06040791897610)
C = c(2.96172084803923, 2.32670984757655, 1.93944408541553)

for (i in 1:3) {
   L[[i]] = 1-C[i]/z[i]^(k+1)
}

S = list(
   read.table("out-.08-inexact.txt"),
   read.table("out-.10-inexact.txt"),
   read.table("out-.12-inexact.txt")
)

pdf("simulp-inexact.pdf", width=11, height=5.5, useDingbats=FALSE)
par(mfrow=c(1,2))
subs = seq(1,50,2)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e7, pch=19, cex=.5, ylim=c(0.98,1),
     plot.first=grid(),
     xlab="Read size", ylab="Seeding probability")
lines(k[1:154], L[[1]][1:154])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e7, pch=19, cex=.5)
lines(k[1:154], L[[2]][1:154])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e7, pch=19, cex=.5)
lines(k[1:154], L[[3]][1:154])

legend(x="bottomright", inset=0.05, legend=c("Simulation", "Estimate"),
       pch=c(19, NA), lty=c(NA, 1), pt.cex=.6, bg="white", box.col=NA)


subs = seq(51,76)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e8, pch=19, cex=.5, ylim=c(0.9996,1),
     plot.first=grid(),
     xlab="Read size", ylab="Seeding probability")
lines(k[146:204], L[[1]][146:204])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e8, pch=19, cex=.5)
lines(k[146:204], L[[2]][146:204])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e8, pch=19, cex=.5)
lines(k[146:204], L[[3]][146:204])
dev.off()

max(abs(1-S[[1]][,2]/1e7 - L[[1]][S[[1]][,1]]))
max(abs(1-S[[2]][,2]/1e7 - L[[2]][S[[2]][,1]]))
max(abs(1-S[[3]][,2]/1e7 - L[[3]][S[[3]][,1]]))
