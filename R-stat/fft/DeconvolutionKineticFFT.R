
args <- commandArgs(trailingOnly = T)
args
fname <- args[1]
td <- read.table(paste0(fname, ".dat"), sep=",", head=TRUE)

t <- td$Time/3600
f <- td$Normalized_heat_flow
ff <- fft(f)
#ff

pdf(paste0(fname, ".pdf"), family="NimbusSan", encoding="KOI8-R.enc")

plot(t, f, type="l", xlab = "t, h", ylab = "Normalized heat flow", col="black", lwd=2)

spectrum(f)
plot(ff, type="h", xlab = "Freq", ylab = "Int", col="black", lwd=2)

de <- 1.0 # sampling interval
x.spec <- spectrum(f,log="no",span=10,plot=FALSE)
x.spec
spx <- x.spec$freq/de
spy <- 2*x.spec$spec
plot(spy~spx, xlab="frequency", ylab="spectral density", type="l")

dev.off()
