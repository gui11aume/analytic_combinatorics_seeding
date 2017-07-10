d = 17
p = .1
q = 1-p

Q = function(z) 1-p*z*(1-(q*z)^d)/(1-q*z)
pdf("Q.pdf", width=5, height=5)
z = seq(-1.5, 1.5, 0.01)
plot(z, Q(z), type="l", ylim=c(-2,2), lwd=2, plot.first=grid())
abline(h=0, col="grey50")
dev.off()

mat = matrix(NA, nrow=512, ncol=512)
for (x in 1:512) {
for (y in 1:512) {
   z = 1.5 * (x-256) / 256 + (0+1i) * 1.5 * (y-256) / 256
   mat[x,y] = Mod(1/Q(z))
}
}

bw = colorRampPalette(c("white", "black"))(1024)
mat[mat > 7] = 7

png("singS.png", width=512, height=512)
par(mar=c(0,0,0,0))
image(mat, col=bw, bty="n", xaxt="n", yaxt="n")
dev.off()
