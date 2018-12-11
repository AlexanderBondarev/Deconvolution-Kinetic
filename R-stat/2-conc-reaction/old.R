
args <- commandArgs(trailingOnly = T)
args
fname <- args[1]
td = read.table(paste0(fname, ".dat"), sep=",", head=TRUE)

pdf(paste0(fname, ".pdf"), family="NimbusSan", encoding="KOI8-R.enc")
plot(td$Time/3600, td$Normalized_heat_flow, type="l", xlab = "t, h", ylab = "Normalized heat flow", col="black", lwd=3)

test_myf <- function(x, color) {
 m1 <- 1.0
 m2 <- 1.0

 k1 <- x[1]
 k2 <- x[2]
 A01 <- 0.0035185
 A02 <- 0.0035185
 C01 <- x[3]
 C02 <- x[4]
 H1 <- x[5]
 H2 <- x[6]
 tm <- as.array(td$Time)
 flow <- as.array(td$Normalized_heat_flow)

# printf("t, k1=%f k2=%f A01=%f A02=%f C01=%f C02=%f, 0, 0, Hf1, Hf2, Hf\n",k1,k2,A01,A02,C01,C02);
 C1 <- C01
 C2 <- C02
 tmax <- max(tm)
 dt <- tmax / length(tm)
 Q <- 0
 mQ <- c()
 n <- 1
 t <- 0
 mt <- c()
    print( sprintf(" N: %d", n) )
# printf("dt=%f\n",dt);
# for(t=0;t<tmax;t+=dt)
 while (t < tmax) {
    A1 <- A01-m1*(C1-C01)
    dC1 <- k1*(A01-m1*(C1-C01))*C1
    P1 <- -1000.0*H1*dC1/dt
    A2 <- A02-m2*(C2-C02)
    dC2 <- k2*(A02-m2*(C2-C02))*C2
    P2 <- -1000.0*H2*dC2/dt
    Psum <- P1+P2;
#    printf("%.1f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",t,dC1/dt,dC2/dt,A1,A2,C1,C2,Hf1,Hf2,Hf);
    C1 <- C1 + dC1
    C2 <- C2 + dC2
    if(t<tm[1]) { t <- t + dt; next }
    dQ <- flow[n]-Psum;
    mQ[n] <- Psum
    mt[n] <- t
    Q <- Q + dQ*dQ;
#    printf("t=%f %f\n",t,tm[n]);
    n <- n + 1
    t <- t + dt
#    print( sprintf(" Q: %d %f %f", n, t, Q) )
  }
 print( sprintf("Quality: %d %f", n, Q/n ) )
 lines(mt/3600, mQ, col=color, lwd=2)
 return(Q/n)
}

myf <- function(x) {
 m1 <- 1.0
 m2 <- 1.0

 k1 <- x[1]
 k2 <- x[2]
 A01 <- 0.0035185
 A02 <- 0.0035185
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
    C1 <- C1 + dC1
    C2 <- C2 + dC2
    if(t<tm[1]) { t <- t + dt; next }
    dQ <- flow[n]-Psum;
    Q <- Q + dQ*dQ;
    n <- n + 1
    t <- t + dt
  }
 return(Q/n)
}

rnd30 <- function(x) { runif(1, x-0.3*abs(x), x+0.3*abs(x)) }

par0 <- c(0.65, 0.070, 0.0009, 0.0027, 40.0, -430.0)
ControlPar <- list(trace=0, maxit=1000)

"Nelder-Mead OPT 1"; optim(sapply(par0, rnd30), myf, control=ControlPar );
"Nelder-Mead OPT 2"; optim(sapply(par0, rnd30), myf, control=ControlPar );
"Nelder-Mead OPT 3"; optim(sapply(par0, rnd30), myf, control=ControlPar );
"Nelder-Mead OPT 4"; optim(sapply(par0, rnd30), myf, control=ControlPar );
"Nelder-Mead OPT 5"; optim(sapply(par0, rnd30), myf, control=ControlPar );

###"STANDART OPT"; r <- optim(par, myf, control=list(trace=1, maxit=100) )
#"BFGS OPT"; optim(par, myf, NULL, method = "BFGS", hessian = TRUE, control=list(trace=1))
#"SANN OPT"; optim(par, myf, method = "SANN", control=list(trace=1))
###"Result: "; r$par
# -k1 0.65 -k2 0.070 -C01 0.0009 -C02 0.0027 -H1 40 -H2 -430
#test_myf(c(0.65, 0.070, 0.0009, 0.0027, 40.0, -430.0))
###test_myf(r$par, "red")

###"STANDART OPT"; r <- optim(par, myf, control=list(trace=1, maxit=1000) )
###"Result: "; r$par
# -k1 0.65 -k2 0.070 -C01 0.0009 -C02 0.0027 -H1 40 -H2 -430
#test_myf(c(0.65, 0.070, 0.0009, 0.0027, 40.0, -430.0))
###test_myf(r$par, "green")

###"STANDART OPT"; r <- optim(par, myf, control=list(trace=1, maxit=10000) )
###"Result: "; r$par
# -k1 0.65 -k2 0.070 -C01 0.0009 -C02 0.0027 -H1 40 -H2 -430
#test_myf(c(0.65, 0.070, 0.0009, 0.0027, 40.0, -430.0))
###test_myf(r$par, "blue")

main = "Море"
location = "topright"
labels = c("Experimental", "Opt 100", "Opt 1000", "Opt 10000")
colors = c("black", "red", "green", "blue")
legend(location, labels, title = main, fill=colors)

###"SANN OPT"; optim(par, myf, method = "SANN", control=list(trace=1))

dev.off()
