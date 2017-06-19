require(opLabR)
rm(list=ls())
dev.off()

# In this script, we
# a) compare the difference between the Carr-Madan approach and
# the standard approach to price puts with the Characteristic Function to
# compute the CDF
# b) compare various option pricing models to compute the CDF

# Diffusion and costs parameters (same as in jump_cdfs)
r <- 0
v <- sqrt(log(1+0.8^2))
p <- 100
sv <- 20
cs <- 44

# A grid of lead-times
tau <- seq(0,1,1e-2)

# Log-normal CDF obtained using the analytic expression
myc <- logNormalCDF(sv = sv, cs = cs, sigma = v)
cdf1 <- myc$cdf

# Jump parameters (as in jump_cdfs)
lambda <- 0.1
muY <- 1.4
vY <- 1
# Transformations for the standard parametrization
muJ <- muY-1
vJ <- log(1 + vY^2)

# Stochastic volatility parameters
v0 <- v^2
vT <- v^2
k <- 1
sigma <- 0.4
rho <- 0.4

# In the case of the pure diffusion, the Carr-Madan approach is about 10-20 faster than the standard method
# (albeit obviously much slower than the analytic formula)
system.time(cdf1b <- CDF(cfBSM, p = p, cs = cs, sv = sv, v = v^2))
system.time(cdf1c <- CDF(cfBSM, carrMadan = FALSE, p = p, cs = cs, sv = sv, v = v^2))
plot(tau, abs(cdf1-cdf1b), type = "l", ylab = "Absolute error", xlab = "Lead-time", col = "red")
lines(tau, abs(cdf1-cdf1c)) # Notice the y-scale (i.e., negligible error) !!!
legend("topright", legend = c("Standard option price", "Carr-Madan approach"), lty = 1, col = c(1,2))

# With the Heston model, the Carr-Madan is about 20 times faster
system.time(cdf2 <- CDF(cfHeston, p = p, cs = cs, sv = sv,
                          v0 = v0, vT = vT, k = k, rho = rho, sigma = sigma))
system.time(cdf2b <- CDF(cfHeston, carrMadan = FALSE, p = p, cs = cs, sv = sv,
                          v0 = v0, vT = vT, k = k, rho = rho, sigma = sigma))
plot(tau, abs(cdf2-cdf2b), type = "l", ylab = "Difference between standard option price and carr-madan approach", xlab = "Lead-time")
# This time, there is not analytical formula to compute the absolute error,
# but notice that the difference between the standard and Carr-Madan approach is negligible

# With the Merton model, the speed-up is between three and four
system.time(cdf3 <- CDF(cfMerton, p = p, cs = cs, sv = sv, v = v^2,
              lambda = lambda, muJ = muJ, vJ = vJ))
system.time(cdf3b <- CDF(cfMerton, carrMadan = FALSE, p = p, cs = cs, sv = sv, v = v^2,
                          lambda = lambda, muJ = muJ, vJ = vJ))
plot(tau, abs(cdf3-cdf3b), type = "l", ylab = "Difference between standard option price and carr-madan approach", xlab = "Lead-time")
# Same comment as for Heston concerning the differences

# With the Bates model, the speed-up is about 20 again
system.time(cdf4 <- CDF(cfBates, p = p, cs = cs, sv = sv,
                          lambda = lambda, muJ = muJ, vJ = vJ,
                          v0 = v0, vT = vT, k = k, rho = rho, sigma = sigma))
system.time(cdf4b <- CDF(cfBates, carrMadan = FALSE, p = p, cs = cs, sv = sv,
                           lambda = lambda, muJ = muJ, vJ = vJ,
                           v0 = v0, vT = vT, k = k, rho = rho, sigma = sigma))
plot(tau, abs(cdf4-cdf4b), type = "l", ylab = "Difference between standard option price and carr-madan approach", xlab = "Lead-time")
# Same comment as for Heston concerning the differences

# We can finally plot the four cost-differential frontiers
plot(tau,cdf1*100, type = "l", xlab = "Lead Time", ylab = "CDF",
     main = "Cost differential frontier", ylim = c(0,50))
lines(tau,cdf2*100, col = "red")
lines(tau,cdf3*100, col = "green")
lines(tau,cdf4*100, col = "blue")
legend("topleft", c("LN", "SV", "Jumps", "SV+Jumps"),
       lty = 1, col = 1:4, cex = 0.7)

