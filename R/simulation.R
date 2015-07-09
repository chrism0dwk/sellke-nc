# Sellke simulation function

# Implements homogeneous epidemic
# with exponential removal times.

sellkeSim <- function(beta, gamma, nu)
{
  S <- length(nu)-1
  I <- 1
  R <- 0
  t <- 0
  
  # Results
  sim <- data.frame(nu=nu, i=rep(Inf, length(nu)), r=rep(Inf, length(nu)), row.names=1:length(nu))
  names(nu) <- 1:length(nu)
  
  # Initial conditions
  I1 <- sample(length(nu),size=1)
  seeds <- sort(nu[-I1])
  iTimes <- 0
  rTimes <- rexp(1, gamma) # Removal times here
  sim$i[I1] <- iTimes
  sim$r[I1] <- rTimes
  
  cumPressure <- 0

  while(I > 0)
  {
    minSeed <- seeds[1]
    nextInfec <- (minSeed-cumPressure)/(beta*I) + t
    nextRemove <- min(rTimes[rTimes>t])
    tOld <- t
    if (nextInfec < nextRemove)
    {
      toInfect <- names(minSeed)
      S <- S-1
      I <- I+1
      t <- nextInfec
      iTimes <- c(iTimes,t)
      seeds <- seeds[-which.min(seeds)]
      rt <- rexp(1,gamma) + t
      rTimes <- c(rTimes,rt)
      sim[toInfect,2:3] <- c(t,rt)
    }
    else
    {
      I <- I-1
      R <- R+1
      t <- nextRemove
    }
    
    cumPressure <- cumPressure + beta*I*(t-tOld)
  }
  
  message("S: ", S, ", I: ", I, ", R: ", R)
  message("Max IIP: ", cumPressure)
  class(sim) <- "epidemic"
  sim
}

