f  = function(z, p, d) 1 - ((1-p)*z)^(d+1)
g  = function(z, p, d) 1 - p*z*(1-((1-p)*z)^d)/(1-(1-p)*z)
dg = function(z, p, d) (((d + 1)*p + (d*p^2 - d*p)*z)*(-(p - 1)*z)^d - p) / ((p^2 - 2*p + 1)*z^2 + 2*(p - 1)*z + 1)

Newton = function(p, d) {
   # Set initial value close to solution.
   oldz = 1
   while (TRUE) {
      newz = oldz - g(oldz, p, d)/ dg(oldz, p, d)
      dz = newz - oldz
      if (abs(dz) < 1e-15) { break }
      if (newz < 1) {
         # Backtrack.
         while (oldz + dz < 1) { dz = dz/2 }
         newz = oldz + dz
      }
      oldz = newz
   }
   return(newz)
}

bisect = function(p, d) {
   lo = 1.0
   hi = 1.2
   while (TRUE) {
      mid = (lo+hi) / 2
      if (g(mid, p, d) < 0) {
         hi = mid
      }
      else {
         lo = mid
      }
      if (hi-lo < 1e-12) break
   }
   return ((hi+lo)/2)
}

cst_term = function(z,p,d) {
   (1-(1-p)*z)^2 / p^2 / (1+d*((1-p)*z)^(d+1)-(d+1)*((1-p)*z)^d)
}

d = 17
k = 1:204

L = list()

for (p in c(0.08, 0.1, 0.12)) {
   z = Newton(p,d)
   C = cst_term(z,p,d)
   L[[as.character(p)]] = 1-C/z^(k+2)
}

S = list(
   read.table("out-.08.txt"),
   read.table("out-.10.txt"),
   read.table("out-.12.txt")
)

pdf("simulp.pdf", width=11, height=5.5, useDingbats=FALSE)
par(mfrow=c(1,2))
subs = seq(1,51,2)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e7, pch=19, cex=.5, ylim=c(0.5,1),
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
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e7, pch=19, cex=.5, ylim=c(0.9,1),
     plot.first=grid(),
     xlab="Read size", ylab="Seeding probability")
lines(k[146:204], L[[1]][146:204])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e7, pch=19, cex=.5)
lines(k[146:204], L[[2]][146:204])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e7, pch=19, cex=.5)
lines(k[146:204], L[[3]][146:204])
dev.off()

max(abs(1-S[[1]][,2]/1e7 - L[[1]][S[[1]][,1]]))
max(abs(1-S[[2]][,2]/1e7 - L[[2]][S[[2]][,1]]))
max(abs(1-S[[3]][,2]/1e7 - L[[3]][S[[3]][,1]]))
