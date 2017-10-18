d = 17
k = 1:160

L = list()
M = list()
N = list()

z2_1 = c(1.066122157, 1.065963331, 1.065648790)
C2_1 = c(2.167761949, 2.164467880, 2.156256928)

z2_2 = c(1.086538875+0.493986316i, 1.08674198+0.494062642i,
   1.086298145+0.493907639i)
C2_2 = c(1.134360181-0.013955522i, 1.133670719-0.012856294i,
  1.123364232-0.013466291i)

z2_3 = c(1.086538875-0.493986316i, 1.08674198-0.494062642i,
  1.086298145-0.493907639i)
C2_3 = c(1.134360181+0.013955522i, 1.133670719+0.012856294i,
  1.123364232+0.013466291i)

z3_1 = c(1.062359015, 1.065177621, 1.065473619)
C3_1 = c(1.894617782, 2.138059081, 2.151092322)

z3_2 = c(1.110242105+0.010966761i, 1.086523354+0.493793966i,
   1.08613424+0.49378375i)
C3_2 = c(0.1391783028+0.1400396942i, 1.12581633-0.017700243i,
   1.119625022-0.015123612i)

z3_3 = c(1.110242105-0.010966761i, 1.086523354-0.493793966i,
   1.08613424-0.49378375i)
C3_3 = c(0.1391783028-0.1400396942i, 1.12581633+0.017700243i,
   1.119625022+0.015123612i)

#.z3_2 = c(1.110242105+0.010966761i, 1, 1)
#.C3_2 = c(0.1391783028+0.1400396942i, 0, 0)
#
#.z3_3 = c(1.110242105-0.010966761i, 1, 1)
#.C3_3 = c(0.1391783028-0.1400396942i, 0, 0)


for (i in 1:3) {
   L[[i]] = C3_1[i]/z3_1[i]^(k+1) - C2_1[i]/z2_1[i]^(k+1)
   M[[i]] = C3_1[i]/z3_1[i]^(k+1) + C3_2[i]/z3_2[i]^(k+1) +
      C3_3[i]/z3_3[i]^(k+1) - C2_1[i]/z2_1[i]^(k+1)
#   N[[i]] = C3_1[i]/z3_1[i]^(k+1) + .C3_2[i]/.z3_2[i]^(k+1) +
#      .C3_3[i]/.z3_3[i]^(k+1) - C2_1[i]/z2_1[i]^(k+1)
}

S = list(
   read.table("out-.05-MEM_fp.txt"),
   read.table("out-.15-MEM_fp.txt"),
   read.table("out-.25-MEM_fp.txt")
)

options(scipen=22)
pdf("simulp_mem_fp.pdf", width=11, height=5.5, useDingbats=FALSE)
par(mfrow=c(1,2))
subs = seq(1,20)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e7, pch=19, cex=.5,
     plot.first=grid(), ylim=c(0,.007),
     xlab="Read size", ylab="Type I error rate")
lines(k[1:110], L[[1]][1:110])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e7, pch=19, cex=.5)
lines(k[1:110], L[[2]][1:110])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e7, pch=19, cex=.5)
lines(k[1:110], L[[3]][1:110])

legend(x="topright", inset=0.05, legend=c("Simulation", "Estimate"),
       pch=c(19, NA), lty=c(NA, 1), pt.cex=.6, bg="white", box.col=NA)

subs = seq(21,45)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e8, pch=19, cex=.5,
     plot.first=grid(), ylim=c(0,0.0008),
     xlab="Read size", ylab="Type I error rate")
lines(k[90:160], L[[1]][90:160])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e8, pch=19, cex=.5)
lines(k[90:160], L[[2]][90:160])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e8, pch=19, cex=.5)
lines(k[90:160], L[[3]][90:160])
dev.off()

pdf("simulp_mem_fp2.pdf", width=5.5, height=5.5, useDingbats=FALSE)
#par(mfrow=c(1,2))
subs = seq(1,20)
plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e7, pch=19, cex=.5,
     plot.first=grid(), ylim=c(-0.003,0.016),
     xlab="Read size", ylab="Type I error rate")
lines(k[1:110], M[[1]][1:110])

points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e7, pch=19, cex=.5)
lines(k[1:110], M[[2]][1:110])

points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e7, pch=19, cex=.5)
lines(k[1:110], M[[3]][1:110])

legend(x="topright", inset=0.05, legend=c("Simulation", "Estimate"),
       pch=c(19, NA), lty=c(NA, 1), pt.cex=.6, bg="white", box.col=NA)

#subs = seq(21,45)
#plot(S[[1]][subs,1], 1-S[[1]][subs,2]/1e8, pch=19, cex=.5,
#     plot.first=grid(), ylim=c(0,0.0008),
#     xlab="Read size", ylab="Type I error rate")
#lines(k[90:160], M[[1]][90:160])
#
#points(S[[2]][subs,1], 1-S[[2]][subs,2]/1e8, pch=19, cex=.5)
#lines(k[90:160], M[[2]][90:160])
#
#points(S[[3]][subs,1], 1-S[[3]][subs,2]/1e8, pch=19, cex=.5)
#lines(k[90:160], M[[3]][90:160])
dev.off()
