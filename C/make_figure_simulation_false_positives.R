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
k = 1:108

L2 = list()

z = Newton(0.02,d)
C = cst_term(z,.02,d)
L1 = C/z^(k+2)

for (kappa in c(0.02, 0.05, 0.08)) {
   p = (1-0.02)*kappa + 0.02*2*kappa/3
   z = Newton(p,d)
   C = cst_term(z,p,d)
   L2[[as.character(p)]] = 1-C/z^(k+2)
}

S = list(
   read.table("out-.02-fp.txt"),
   read.table("out-.05-fp.txt"),
   read.table("out-.08-fp.txt")
)

pdf("simulp_false_positives.pdf", width=5.5, height=5.5, useDingbats=FALSE)
plot(S[[1]][,1], S[[1]][,2]/1e7, pch=19, cex=.5, ylim=c(0,0.002),
     plot.first=grid(),
     xlab="Read size", ylab="Seeding probability")
lines(k, L1*L2[[1]])

points(S[[2]][,1], S[[2]][,2]/1e7, pch=19, cex=.5)
lines(k, L1*L2[[2]])

points(S[[3]][,1], S[[3]][,2]/1e7, pch=19, cex=.5)
lines(k, L1*L2[[3]])

legend(x="topright", inset=0.05, legend=c("Simulation", "Estimate"),
       pch=c(19, NA), lty=c(NA, 1), pt.cex=.6, bg="white", box.col=NA)
dev.off()
