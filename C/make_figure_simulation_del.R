f = function(z, p, d, h) 1 + (1-h)*((-(p - 1)*z)^d*(-h + 1)^(d - 1) + (p - 1)*z)/((h - 1)*(p - 1)*z - 1)
g = function(z, p, d, h) -p*z + ((h - 1)*p*z - h)*((-(p - 1)*z)^d*(-h + 1)^(d - 1) + (p - 1)*z)/((h - 1)*(p - 1)*z - 1) + 1
dg = function(z, p, d, h) -((h - 1)*p*z - h)*((-(p - 1)*z)^d*(-h + 1)^(d - 1) + (p - 1)*z)*(h - 1)*(p - 1)/((h - 1)*(p - 1)*z - 1)^2 + ((-(p - 1)*z)^d*(-h + 1)^(d - 1) + (p - 1)*z)*(h - 1)*p/((h - 1)*(p - 1)*z - 1) - ((-(p - 1)*z)^(d - 1)*d*(-h + 1)^(d - 1)*(p - 1) - p + 1)*((h - 1)*p*z - h)/((h - 1)*(p - 1)*z - 1) - p

gamma = function(z,p,d,h) d*h - (1-h)*((d-1)*h-p*((d-1)*h+d+1))*z - d*(1-h)^2*p*(1-p)*z^2

alt_const = function(z,p,d,h) ((1 - h)*(1 - p)*z - 1)^2 / ((p+(1-p)*h)*z-gamma(z,p,d,h)*(1-h)^(d-1)*((1-p)*z)^d) / (h+p*z*(1-h))

Newton = function(p, d, h) {
   # Set initial value close to solution.
   oldz = 1
   while (TRUE) {
      newz = oldz - g(oldz, p, d, h) / dg(oldz, p, d, h)
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


bisect = function(p, d, h) {
   lo = 1.0
   hi = 1.2
   while (TRUE) {
      mid = (lo+hi) / 2
      if (g(mid, p, d, h) < 0) {
         hi = mid
      }
      else {
         lo = mid
      }
      if (hi-lo < 1e-12) break
   }
   return ((hi+lo)/2)
}


cst_term = function(z,p,d,h) {
   -f(z,p,d,h) / dg(z,p,d,h)
}

p = 0.05
d = 17
k = 1:504

L = list()

for (h in c(0.14, 0.15, 0.16)) {
   z = Newton(p,d,h)
   options(digits=22)
   C = alt_const(z,p,d,h)
   L[[as.character(h)]] = 1-C/z^k
}

S = list(
   read.table("out-.05-.14.txt"),
   read.table("out-.05-.15.txt"),
   read.table("out-.05-.16.txt")
)

pdf("simulpdel.pdf", width=11, height=5.5, useDingbats=FALSE)
par(mfrow=c(1,2))
subs = seq(1,64,3)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e7, pch=19, cex=.5, ylim=c(0.35,.95),
     plot.first=grid(),
     xlab="Read size", ylab="Seeding probability")
lines(k[90:365], L[[1]][90:365])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e7, pch=19, cex=.5)
lines(k[90:365], L[[2]][90:365])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e7, pch=19, cex=.5)
lines(k[90:365], L[[3]][90:365])

legend(x="bottomright", inset=0.05, legend=c("Simulation", "Estimate"),
       pch=c(19, NA), lty=c(NA, 1), pt.cex=.6, bg="white", box.col=NA)


subs = seq(64,138,3)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e7, pch=19, cex=.5, ylim=c(0.85,1),
     plot.first=grid(),
     xlab="Read size", ylab="Seeding probability")
lines(k[346:504], L[[1]][346:504])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e7, pch=19, cex=.5)
lines(k[346:504], L[[2]][346:504])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e7, pch=19, cex=.5)
lines(k[346:504], L[[3]][346:504])
dev.off()

max(abs(1-S[[1]][,2]/1e7 - L[[1]][S[[1]][,1]]))
max(abs(1-S[[2]][,2]/1e7 - L[[2]][S[[2]][,1]]))
max(abs(1-S[[3]][,2]/1e7 - L[[3]][S[[3]][,1]]))
