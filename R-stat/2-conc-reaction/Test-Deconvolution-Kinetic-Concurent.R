
args <- commandArgs(trailingOnly = T)
args
fname <- args[1]
p_A0  <- as.numeric(args[2])
p_dH1 <- as.numeric(args[3])
p_dH2 <- as.numeric(args[4])
p_k1  <- as.numeric(args[5])
p_k2  <- as.numeric(args[6])
p_C01 <- as.numeric(args[7])
p_C02 <- as.numeric(args[8])
td = read.table(paste0(fname, ".dat"), sep=",", head=TRUE)

pdf(paste0(fname, "_Concurent.pdf"), family="NimbusSan", encoding="KOI8-R.enc")

test_myf <- function(x, caption) {
 m1 <- 1.0
 m2 <- 1.0

 k1 <- x[1]
 k2 <- x[2]
 A01 <- p_A0
 A02 <- p_A0
 C01 <- x[3]
 C02 <- x[4]
 H1 <- x[5]
 H2 <- x[6]

 tm <- as.array(td$Time)
 flow <- as.array(td$Normalized_heat_flow)

 C1 <- C01
 C2 <- C02
 tmax <- max(tm)
 dt <- tmax / length(tm)
 Q <- 0
 mQ <- c()
 mQ1 <- c()
 mQ2 <- c()
 n <- 1
 t <- 0
 mt <- c()
    print( sprintf(" N: %d", n) )
 while (t < tmax) {
    A1 <- A01-m1*(C1-C01)
    dC1 <- k1*(A01-m1*(C1-C01))*C1
    P1 <- -1000.0*H1*dC1/dt
    A2 <- A02-m2*(C2-C02)
    dC2 <- k2*(A02-m2*(C2-C02))*C2
    P2 <- -1000.0*H2*dC2/dt
    Psum <- P1+P2;
    C1 <- C1 + dC1
    C2 <- C2 + dC2
    if(t<tm[1]) { t <- t + dt; next }
    dQ <- flow[n]-Psum;
    mQ1[n] <- P1
    mQ2[n] <- P2
    mQ[n] <- Psum
    mt[n] <- t
    Q <- Q + dQ*dQ;
    n <- n + 1
    t <- t + dt
  }
 print( sprintf("Quality: %d %f", n, Q/n ) )
 cap <- sprintf("%s (Quality: %.3e)", caption, Q/n )
 yrange <- range(td$Normalized_heat_flow, mQ1, mQ2, mQ)
 plot(td$Time/3600, td$Normalized_heat_flow, type="l", xlab = "t, h", ylab = "Normalized heat flow", col="black", main=cap, lwd=3, ylim = yrange)
 lines(mt/3600, mQ1, col="blue", lwd=2)
 lines(mt/3600, mQ2, col="red", lwd=2)
 lines(mt/3600, mQ, col="green", lwd=2)
 main = "Heat flow"
 location = "topright"
 labels = c("Experimental", "P1", "P2", "Psum")
 colors = c("black", "blue", "red", "green")
 legend(location, labels, title = main, fill=colors)
 return(Q/n)
}

par0 <- c(p_k1, p_k2, p_C01, p_C02, p_dH1, p_dH2)

test_myf(par0, "Test Opt")

dev.off()
