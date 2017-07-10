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

d = 19
k = 75
p = 0.02

z = Newton(p,d)
C = cst_term(z,p,d)
print(C/z^(k+2))
