d = 17
p = .05
h = 0.15

#Q = function(z) 1-h-p*z*(1-h-(q*z)^d)/(1-q*z)
Q = function(z) 1-p*z-(p*z*(1-h)+h)*((1-p)*z-(1-h)^(d-1)*((1-p)*z)^d)/(1-(1-h)*(1-p)*z)
pdf("Q.pdf", width=5, height=5)
z = seq(-1.8, 1.8, 0.01)
plot(z, Q(z), type="l", ylim=c(-2,2), lwd=2, plot.first=grid())
abline(h=0, col="grey50")
dev.off()

mat = matrix(NA, nrow=512, ncol=512)
for (x in 1:512) {
for (y in 1:512) {
   z = 1.8 * (x-256) / 256 + (0+1i) * 1.8 * (y-256) / 256
   mat[x,y] = Mod(1/Q(z))
}
}

bw = colorRampPalette(c("white", "black"))(1024)
mat[mat > 7] = 7

png("singSdel.png", width=512, height=512)
par(mar=c(0,0,0,0))
image(mat, col=bw, bty="n", xaxt="n", yaxt="n")
dev.off()
