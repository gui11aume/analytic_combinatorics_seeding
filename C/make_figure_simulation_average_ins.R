k = 1:160

L = list()

# Computed with a symbolic software.
C1 = c(0.07155318459, 0.08723441701, 0.1022467274)
C2 = c(-.1668078544, -.1930301413, -.2154411928)

for (i in 1:3) {
   L[[i]] = (k+1)*C1[i] + C2[i]
}

S = list(
   read.table("out-.05-.15-.04-.45-average.txt"),
   read.table("out-.05-.15-.05-.45-average.txt"),
   read.table("out-.05-.15-.06-.45-average.txt")
)

pdf("simulins-average.pdf", width=5, height=5.5, useDingbats=FALSE)
subs = seq(1,51,2)
plot(S[[1]][subs,1], S[[1]][subs,2], pch=19, cex=.5,
     plot.first=grid(), ylim=c(0,15),
     xlab="Read size", ylab="Number of insertions")
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
