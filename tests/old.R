# rm(list=ls())
#
# putCF <- function(cf, S, X, tau, r, ...) {
#   if (tau == 0) {
#     max(X - S, 0)
#   } else {
#     callCF(cf, S, X, tau, r, q = 0, ...) + X * exp(-r*tau) - S
#   }
# }
#
# profit <- function(cf, Q, tau, r, p, cs, sv, ...) {
#   (p-cs) * Q - (p-sv) * putCF(cf, 1, Q, tau, r, ...)
# }
#
# # profit2 <- function(cf, Q, tau, r, crf, ...) {
# #   crf * Q - putCF(cf, 1, Q, tau, r, ...)
# # }
#
# dev.off()
# r <- 0
# v <- sqrt(log(1+0.8^2))
# p <- 100
# sv <- 20
# cs <- 44
# tau <- seq(0,1,1e-2)
# myc <- logNormalCDF(sv = sv, cs = cs, sigma = v)
# cdf1 <- myc$cdf
#
# lambda <- 0.1
# muY <- 1.4
# vY <- 1
# muJ <- muY-1
# vJ <- log(1 + vY^2)
#
# v0 <- 0.64
# vT <- 0.64
# k <- 1
# sigma <- 0.1
# rho <- 0.4
#
# myCDF <- function(cf,  carrMadan = TRUE,
#                   r = 0, tau = seq(0,1,1e-2), p, cl, sv, ...) {
#
#   if (carrMadan == TRUE) {
#     puts <- carrMadanPut(cf, r = r, tau = tau, ...)
#     Q <- puts[,1]
#     P <- puts[,-1]
#     betaS <- sapply(1:dim(P)[2], function(i)
#       max((p-cl[i]) * Q - (p-sv) * P[,i]))/(p-sv)
#   } else {
#     betaS <- -sapply(1:length(tau), function(i)
#       optim(1, function(Q) -profit(cf, Q, tau[i], r, p, cl[i], sv, ...),
#             method = "Brent", lower = 0, upper = 10)$value)/(p-sv)
#   }
#
#   cs <- p-betaS*(p-sv)
#   cdf <- (cs-cl)/cs
#   if (any(tau == 0)) {
#     cdf[tau == 0] <- 0
#   }
#   return(cdf)
# }
#
# #system.time(cdf1b <- myCDF(cfBSM, p = p, cl = myc$cl, sv = sv, v = v^2))
# # system.time(cdf1c <- myCDF(cfBSM, carrMadan = FALSE, p = p, cl = myc$cl, sv = sv, v = v^2))
# # plot(cdf1-cdfb, type = "l")
# # lines(cdf1-cdf1c, col = "red")
#
# system.time(cdf2 <- myCDF(cfHeston, p = p, cl = myc$cl, sv = sv,
#                           v0 = v0, vT = vT, k = k, rho = rho, sigma = sigma))
# # system.time(cdf2b <- myCDF(cfHeston, carrMadan = FALSE, p = p, cl = myc$cl, sv = sv,
# #                           v0 = v0, vT = vT, k = k, rho = rho, sigma = sigma))
# # plot(cdf2-cdf2b, type = "l")
#
# system.time(cdf3 <- myCDF(cfMerton, p = p, cl = myc$cl, sv = sv, v = v^2,
#               lambda = lambda, muJ = muJ, vJ = vJ))
# # system.time(cdf3b <- myCDF(cfMerton, carrMadan = FALSE, p = p, cl = myc$cl, sv = sv, v = v^2,
# #                           lambda = lambda, muJ = muJ, vJ = vJ))
# # plot(cdf3-cdf3b, type = "l")
#
# system.time(cdf4 <- myCDF(cfBates, p = p, cl = myc$cl, sv = sv,
#                           lambda = lambda, muJ = muJ, vJ = vJ,
#                           v0 = v0, vT = vT, k = k, rho = rho, sigma = sigma))
# # system.time(cdf4b <- myCDF(cfBates, carrMadan = FALSE, p = p, cl = myc$cl, sv = sv,
# #                            lambda = lambda, muJ = muJ, vJ = vJ,
# #                            v0 = v0, vT = vT, k = k, rho = rho, sigma = sigma))
# # plot(cdf4-cdf4b, type = "l")
#
# plot(tau,cdf1*100, type = "l", xlab = "Lead Time", ylab = "CDF",
#      main = "Cost differential frontier", ylim = c(0,50))
# lines(tau,cdf2*100, col = "red")
# lines(tau,cdf3*100, col = "green")
# lines(tau,cdf4*100, col = "blue")
# legend("topleft", c("LN", "SV", "Jumps", "SV+Jumps"),
#        lty = 1, col = 1:4, cex = 0.7)
#
# head(cbind(cdf1, cdf2, cdf3, cdf4))
# tail(cbind(cdf1, cdf2, cdf3, cdf4))
#
#
