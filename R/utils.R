# Implementation of the log-normal model (fully analytic)
logNormalCDF <- function(price = 100, cs = 50, sv = 30, 
                         cf = (price-cs)/(price-sv),
                         sigma = 0.8, tau = seq(0,1,1e-2)) {

  assert_that(is.number(price) && price > 0)
  assert_that(is.number(cs) && cs > 0 && cs < price)
  assert_that(is.number(sv) && sv > 0 && sv < price && sv < cs)
  assert_that(is.numeric(cf) && cf >= 0 && cf < 1)
  assert_that(is.number(sigma) && sigma > 0)
  assert_that(is.numeric(tau) && all(tau >= 0) && all(tau <= 1))

  ci <- ((price-sv) * cf - (price-cs)) /(1 - cf)
  zl <- qnorm(cf) + sigma*sqrt(tau)
  cl <- price + ci - pnorm(zl) * (price - sv + ci)

  out <- data.frame(tau = tau, cl = cl, cdf = (cs-cl)/cs)
  return(out)
}

# Carr-Madan approach
carrMadanPut <- function(cf, r = 0, tau = seq(0,1,1e-2), ...) {

  # Parameters for the Carr Madan approach
  N <- 2^12 # Number of strikes
  alpha <- 3/4 # Damping factor
  eta <- 1/4 # 1/step in units of log-strike

  # Needed for Carr Madan
  n <- seq(0,N-1)
  u <- n*eta
  S <- 100
  s0 <- log(S)
  lambda <- 2*pi/(N*eta)

  # Computation of the characteristic functions
  chi <- sapply(tau, function(t)
    cf(u-(1+alpha)*1i, S = S, tau = t, r = r, q = 0, ...))
  psi <- t(exp(-r*tau)*t(chi/(alpha^2+alpha-u^2+1i*(1+2*alpha)*u)))

  # FFT step
  x <- exp(1i*u*(1/2*N*lambda-s0))*psi*eta
  x[1,] <- x[1,]*1/2
  X <- mvfft(x)

  # Compute the put prices
  k <- -1/2*N*lambda+lambda*n+s0
  P <- exp(-alpha*k)*Re(X/pi) + exp(k) %o% exp(-r*tau) - S

  return(cbind(exp(k)/S,P/S))
}

# Compute the calls by standard integration of the characteristic function
putCF <- function(cf, S, X, tau, r, ...) {
  if (tau == 0) {
    max(X - S, 0)
  } else {
    callCF(cf, S, X, tau, r, q = 0, ...) + X * exp(-r*tau) - S
  }
}

# Compute the profit
profit <- function(cf, Q, tau, r, p, cs, sv, ...) {
  (p-cs) * Q - (p-sv) * putCF(cf, 1, Q, tau, r, ...)
}

# Compute the CDF by using either the Carr-Madan approach or the standard way of obtaining option prices
CDF <- function(cf,  carrMadan = TRUE,
                  r = 0, tau = seq(0,1,1e-2), p, cl, sv, ...) {

  if (carrMadan == TRUE) {
    puts <- carrMadanPut(cf, r = r, tau = tau, ...) # Compute the grid of options prices and quantities at each time-step
    Q <- puts[,1]
    P <- puts[,-1]
    profitlong <- sapply(1:dim(P)[2], function(i) # Solve for the profit maximizing quantity at each time-step by grid-search
      max((p-cl[i]) * Q - (p-sv) * P[,i]))
  } else {
    profitlong <- -sapply(1:length(tau), function(i) # Solve for the profit maximizing quantity by brute-force
      optim(1, function(Q) -profit(cf, Q, tau[i], r, p, cl[i], sv, ...),
            method = "Brent", lower = 0, upper = 10)$value)
  }

  # Recover the cost short and the cdf
  cs <- p-profitlong
  cdf <- (cs-cl)/cs
  if (any(tau == 0)) {
    cdf[tau == 0] <- 0
  }
  return(cdf)
}


# carrMadanOptions <- function(cf, S, r, tau, q) {
#   N <- 2^12
#   alpha <- 3/4
#   eta <- 2
#
#   n <- seq(0,N-1)
#   u <- n*eta
#   ucf <- u-(1+alpha)*1i
#   s0 <- log(S)
#   lambda <- 2*pi/(N*eta);
#   k <- -1/2*N*lambda+lambda*n+s0
#
#   psi <- cf(ucf)
#   chi <- exp(-r*tau)*psi/(alpha^2+alpha-u^2+1i*(1+2*alpha)*u)
#   x <- exp(1i*u*(1/2*N*lambda-s0))*chi*eta
#   x[1] <- x[1]*1/2
#   X <- fft(x)
#   callPrice <- exp(-alpha*k)*Re(X/pi)
#   putPrice <- callPrice + exp(k)*exp(-r*tau)-S*exp(-q*tau)
#
#   return(cbind(exp(k),callPrice,putPrice))
# }
