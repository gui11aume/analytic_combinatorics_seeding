f = function(z, p, d, h, r, a) -((a - r)*z - 1)*(((-(p + r - 1)*z)^d*(-h + 1)^(d - 1) + (p + r - 1)*z)*(h - 1)/((h - 1)*(p + r - 1)*z - 1) - 1)*(r - 1)
g = function(z, p, d, h, r, a) (((h - 1)*r^2 + a*h - ((a + 1)*h - a - 1)*r - a)*z^2 + (h*(-h + 1)^d*r + ((h*(-h + 1)^d - (-h + 1)^d)*p*r - (a*h*(-h + 1)^d - a*(-h + 1)^d)*p)*z^2 - h*(-h + 1)^d + (a*h*(-h + 1)^d + (h*(-h + 1)^d - (-h + 1)^d)*p - (a*h*(-h + 1)^d - (a - 1)*(-h + 1)^d + (h*(-h + 1)^d - (-h + 1)^d)*p)*r)*z)*(-(p + r - 1)*z)^d - (h - 1)*r - ((h - 1)*r^2 + (a + 1)*h - ((a + 2)*h - a - 2)*r - a - 1)*z + h - 1)/((h^2 - (h^2 - 2*h + 1)*p - (h^2 - 2*h + 1)*r - 2*h + 1)*z + h - 1)
dg = function(z, p, d, h, r, a) -(((h - 1)*p*r - (a*h - a)*p)*z^2 + h*r + (a*h + (h - 1)*p - (a*h + (h - 1)*p - a + 1)*r)*z - h)*((-(p + r - 1)*z)^d*(-h + 1)^(d - 1) + (p + r - 1)*z)*(h - 1)*(p + r - 1)/((h - 1)*(p + r - 1)*z - 1)^2 + (a + p)*r + 2*(a*p - p*r)*z - ((-(p + r - 1)*z)^(d - 1)*d*(-h + 1)^(d - 1)*(p + r - 1) - p - r + 1)*(((h - 1)*p*r - (a*h - a)*p)*z^2 + h*r + (a*h + (h - 1)*p - (a*h + (h - 1)*p - a + 1)*r)*z - h)/((h - 1)*(p + r - 1)*z - 1) + (a*h + (h - 1)*p - (a*h + (h - 1)*p - a + 1)*r + 2*((h - 1)*p*r - (a*h - a)*p)*z)*((-(p + r - 1)*z)^d*(-h + 1)^(d - 1) + (p + r - 1)*z)/((h - 1)*(p + r - 1)*z - 1) - a - p

Newton = function(p, d, h, r, a) {
   # Set initial value close to solution.
   oldz = 1
   while (TRUE) {
      newz = oldz - g(oldz, p, d, h, r, a) / dg(oldz, p, d, h, r, a)
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


bisect = function(p, d, h, r, a) {
   lo = 1.0
   hi = 1.2
   while (TRUE) {
      mid = (lo+hi) / 2
      if (g(mid, p, d, h, r, a) < 0) {
         hi = mid
      }
      else {
         lo = mid
      }
      if (hi-lo < 1e-12) break
   }
   return ((hi+lo)/2)
}


cst_term = function(z,p,d,h,r,a) {
   -f(z,p,d,h,r,a) / dg(z,p,d,h,r,a)
}

p = 0.05
d = 17
h = .15
a = 0.45
k = 1:1008

L = list()

for (r in c(0.04, 0.05, 0.06)) {
   options(digits=22)
   z = Newton(p,d,h,r,a)
   C = cst_term(z,p,d,h,r,a)
   L[[as.character(r)]] = 1-C/z^{k+1}
}

S = list(
   read.table("out-.05-.15-.04-.45.txt"),
   read.table("out-.05-.15-.05-.45.txt"),
   read.table("out-.05-.15-.06-.45.txt")
)

pdf("simulpins.pdf", width=11, height=5.5, useDingbats=FALSE)
par(mfrow=c(1,2))
subs1 = seq(1,44,2)
subs2 = seq(55,68,3)
plot(S[[1]][subs1,1], 1-S[[1]][subs1,2]/1e7,
     pch=19, cex=.5, ylim=c(0.35,.95),
     plot.first=grid(),
     xlab="Read size", ylab="Seeding probability")
points(S[[1]][subs2,1], 1-S[[1]][subs2,2]/1e7, pch=19, cex=.5)
lines(k[100:804], L[[1]][100:804])

points(S[[2]][subs1,1], 1-S[[2]][subs1,2]/1e7, pch=19, cex=.5)
points(S[[2]][subs2,1], 1-S[[2]][subs2,2]/1e7, pch=19, cex=.5)
lines(k[100:804], L[[2]][100:804])

points(S[[3]][subs1,1], 1-S[[3]][subs1,2]/1e7, pch=19, cex=.5)
points(S[[3]][subs2,1], 1-S[[3]][subs2,2]/1e7, pch=19, cex=.5)
lines(k[100:804], L[[3]][100:804])

legend(x="bottomright", inset=0.05, legend=c("Simulation", "Estimate"),
       pch=c(19, NA), lty=c(NA, 1), pt.cex=.6, bg="white", box.col=NA)


subs = seq(68,92)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e7, pch=19, cex=.5, ylim=c(0.85,1),
     plot.first=grid(),
     xlab="Read size", ylab="Seeding probability")
lines(k[792:1008], L[[1]][792:1008])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e7, pch=19, cex=.5)
lines(k[792:1008], L[[2]][792:1008])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e7, pch=19, cex=.5)
lines(k[792:1008], L[[3]][792:1008])
dev.off()

max(abs(1-S[[1]][,2]/1e7 - L[[1]][S[[1]][,1]]))
max(abs(1-S[[2]][,2]/1e7 - L[[2]][S[[2]][,1]]))
max(abs(1-S[[3]][,2]/1e7 - L[[3]][S[[3]][,1]]))
