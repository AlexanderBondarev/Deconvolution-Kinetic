
t <- seq(0,200,by=0.1)
x <- cos(2*pi*t/16) + 0.75*sin(2*pi*t/5)
ftx <- fft(x)

f <- runif(10, -10.0, 10.0)
f
"FFT:"
ft <- fft(f)
ft

pdf("fft.pdf", family="NimbusSan", encoding="KOI8-R.enc")

plot(t, x, type="l", xlab = "", ylab = "Int", col="black", lwd=2)
spectrum(x)
#plot(ftx, type="h", xlab = "Freq", ylab = "Int", col="black", lwd=2)

plot(f, type="h", xlab = "", ylab = "Int", col="black", lwd=2)
spectrum(f)
plot(ft, type="h", xlab = "Freq", ylab = "Int", col="black", lwd=2)

dev.off()
