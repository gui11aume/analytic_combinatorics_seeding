P1  = function(z, p, d) p*z*((1-((1-p)*z)^d) / (1-(1-p)*z))^2
dP1 = function(z, p, d) -2*((1-p)*z)^(d - 1)*(((1-p)*z)^d - 1)*d*((1-p) - 1)*(1-p)*z/((1-p)*z - 1)^2 + 2*(((1-p)*z)^d - 1)^2*(-p)*(1-p)*z/((1-p)*z - 1)^3 - (((1-p)*z)^d - 1)^2*(-p)/((1-p)*z - 1)^2
P2  = function(z, p, d) (1-((1-p)*z)^d) / (1-(1-p)*z)
#dP2 = function(z, p, d) (q*z)^(d - 1)*d*q/(q*z - 1) - ((q*z)^d - 1)*q/(q*z - 1)^2
Q   = function(z, p, d) 1 - p*z*(1-((1-p)*z)^d)/(1-(1-p)*z)
dQ  = function(z, p, d) (((d + 1)*p + (d*p^2 - d*p)*z)*(-(p - 1)*z)^d - p) / ((p^2 - 2*p + 1)*z^2 + 2*(p - 1)*z + 1)
d2Q = function(z, p, d) ((((d^2 - d)*(1-p)^3 - (d^2 - d)*(1-p)^2)*z^2 - d^2 + (d^2 + d)*(1-p) - 2*((d^2 - 1)*(1-p)^2 - (d^2 - 1)*(1-p))*z - d)*((1-p)*z)^d - 2*((1-p)^2 - (1-p))*z)/((1-p)^3*z^4 - 3*(1-p)^2*z^3 + 3*(1-p)*z^2 - z)


Newton = function(p, d) {
   # Set initial value close to solution.
   oldz = 1
   while (TRUE) {
      newz = oldz - Q(oldz, p, d)/ dQ(oldz, p, d)
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
      if (Q(mid, p, d) < 0) {
         hi = mid
      }
      else {
         lo = mid
      }
      if (hi-lo < 1e-12) break
   }
   return ((hi+lo)/2)
}

#cst_term_1 = function(z,p,d) -P1(z,p,d)/dQ(z,p,d)/z / P2(z,p,d)
cst_term_1 = function(z,p,d) {
   q = 1-p
   ((q*z - 1)*(q*z)^d - q*z + 1)/((d*q*z - d - 1)*(q*z)^d + 1)
}
#cst_term_2 = function(z,p,d) (dP1(z,p,d)/dQ(z,p,d) - P1(z,p,d)*d2Q(z,p,d)/dQ(z,p,d)^2) / P2(z,p,d)
cst_term_2 = function(z,p,d) {
   q = 1-p
   -((d^2*q^2*z^2 - (2*d^2 + 2*d + 1)*q*z + d^2 + 2*d + 1)*(q*z)^(2*d) + (d^2*q^2*z^2 - 2*(d^2 - d - 1)*q*z + d^2 - 2*d - 2)*(q*z)^d - q*z + 1)/((q*z)^d*d*q*z - (q*z)^d*d - (q*z)^d + 1)^2
}

d = 17
k = 1:160

L = list()

for (p in c(0.08, 0.1, 0.12)) {
   z = Newton(p,d)
   C1 = cst_term_1(z,p,d)
   C2 = cst_term_2(z,p,d)
   L[[as.character(p)]] = (k+1)*C1 + C2
}

S = list(
   read.table("out-.08-average-sub.txt"),
   read.table("out-.10-average-sub.txt"),
   read.table("out-.12-average-sub.txt")
)

pdf("simulp-average.pdf", width=5, height=5.5, useDingbats=FALSE)
subs = seq(1,51,2)
plot(S[[1]][subs,1], S[[1]][subs,2], pch=19, cex=.5, ylim=c(0,23),
     plot.first=grid(),
     xlab="Read size", ylab="Number of substitutions")
lines(k[1:156], L[[1]][1:156])

points(S[[2]][subs,1], S[[2]][subs,2], pch=19, cex=.5)
lines(k[1:156], L[[2]][1:156])

points(S[[3]][subs,1], S[[3]][subs,2], pch=19, cex=.5)
lines(k[1:156], L[[3]][1:156])

legend(x="bottomright", inset=0.05, legend=c("Simulation", "Estimate"),
       pch=c(19, NA), lty=c(NA, 1), pt.cex=.6, bg="white", box.col=NA)
dev.off()

max(abs(S[[1]][,2] - L[[1]][S[[1]][,1]]))
max(abs(S[[2]][,2] - L[[2]][S[[2]][,1]]))
max(abs(S[[3]][,2] - L[[3]][S[[3]][,1]]))
