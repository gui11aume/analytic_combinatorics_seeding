k = 1:160

L = list()

# Computed with a symbolic software.
C1 = c(0.1393347400, 0.1473025580, 0.1552921440)
C2 = c(-.4206813194, -.4265161201, -.4323752123)

for (i in 1:3) {
   L[[i]] = (k+1)*C1[i] + C2[i]
}

S = list(
   read.table("out-.05-.14-average-del.txt"),
   read.table("out-.05-.15-average-del.txt"),
   read.table("out-.05-.16-average-del.txt")
)

pdf("simuldel-average.pdf", width=5, height=5.5, useDingbats=FALSE)
subs = seq(1,51,2)
plot(S[[1]][subs,1], S[[1]][subs,2], pch=19, cex=.5,
     plot.first=grid(), ylim=c(0,23),
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
