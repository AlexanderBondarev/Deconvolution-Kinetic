
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

pdf(paste0(fname, "_Sequential.pdf"), family="NimbusSan", encoding="KOI8-R.enc")

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
#    C1 <- C1 + dC1
#    C2 <- C2 + dC2
    A1 <- A1 - dC1;
    C1 <- C1 + (dC1-dC2);
    C2 <- C2 + dC2;
    A2 <- C1;

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
 plot(td$Time/3600, td$Normalized_heat_flow, type="l", xlab = "t, h", ylab = "Normalized heat flow", col="black", main=caption, lwd=3)
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


myf <- function(x) {
 m1 <- 1.0
 m2 <- 1.0

 k1 <- x[1]
 k2 <- x[2]
 A01 <- p_A0
 A02 <- p_A0
#A01 <- 0.0035185
#A02 <- 0.0035185
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
 n <- 1
 t <- 0
 while (t < tmax) {
    A1 <- A01-m1*(C1-C01)
    dC1 <- k1*(A01-m1*(C1-C01))*C1
    P1 <- -1000.0*H1*dC1/dt
    A2 <- A02-m2*(C2-C02)
    dC2 <- k2*(A02-m2*(C2-C02))*C2
    P2 <- -1000.0*H2*dC2/dt
    Psum <- P1+P2;
#    C1 <- C1 + dC1
#    C2 <- C2 + dC2
    A1 <- A1 - dC1;
    C1 <- C1 + (dC1-dC2);
    C2 <- C2 + dC2;
    A2 <- C1;
    if(t<tm[1]) { t <- t + dt; next }
    dQ <- flow[n]-Psum;
    Q <- Q + dQ*dQ;
    n <- n + 1
    t <- t + dt
  }
 return(Q/n)
}

k <- 0.1

krnd <- function(x) { runif(1, x-k*abs(x), x+k*abs(x)) }

par0 <- c(p_k1, p_k2, p_C01, p_C02, p_dH1, p_dH2)

###"STANDART OPT"; r <- optim(par, myf, control=list(trace=1, maxit=100) )
#"BFGS OPT"; optim(par, myf, NULL, method = "BFGS", hessian = TRUE, control=list(trace=1))
#"SANN OPT"; optim(par, myf, method = "SANN", control=list(trace=1))
###"Result: "; r$par
# -k1 0.65 -k2 0.070 -C01 0.0009 -C02 0.0027 -H1 40 -H2 -430
#test_myf(c(0.65, 0.070, 0.0009, 0.0027, 40.0, -430.0))
###test_myf(r$par, "red")

ControlPar <- list(trace=0, maxit=30000, eltol=1e-11)
"Nelder-Mead Opt 30000"; r1 <- optim(sapply(par0, krnd), myf, control=ControlPar )
"Result: "; r1$par
test_myf(r1$par, "Nelder-Mead Opt 30000")

ControlPar <- list(trace=0, maxit=50000)
"SANN Opt 50000"; r2 <- optim(sapply(par0, krnd), myf,  method = "SANN", control=ControlPar )
"Result: "; r2$par
test_myf(r2$par, "SANN Opt 50000")

dev.off()
