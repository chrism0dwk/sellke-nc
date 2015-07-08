# Sellke simulation function

# Implements homogeneous epidemic
# with exponential removal times.

sellkeSim <- function(beta, gamma, nu)
{
  S <- length(nu)-1
  I <- 1
  R <- 0
  t <- 0
  
  seeds <- sort(nu[-sample(length(nu),size=1)])
  iTimes <- 0
  rTimes <- rexp(1, gamma) # Removal times here
  
  cumPressure <- 0

  while(I > 0)
  {
    minSeed <- seeds[1]
    nextInfec <- (minSeed-cumPressure)/(beta*I) + t
    nextRemove <- min(rTimes[rTimes>t])
    tOld <- t
    if (nextInfec < nextRemove)
    {
      S <- S-1
      I <- I+1
      t <- nextInfec
      iTimes <- c(iTimes,t)
      seeds <- seeds[-which.min(seeds)]
      rTimes <- c(rTimes,rexp(1,gamma)+t)
    }
    else
    {
      I <- I-1
      R <- R+1
      t <- nextRemove
    }
    
    cumPressure <- cumPressure + beta*I*(t-tOld)
  }
  
  rv <- data.frame(time=c(iTimes,rTimes),event=c(rep("i",length(iTimes)),rep("r",length(rTimes))))
  message("S: ", S, ", I: ", I, ", R: ", R)
  rv[order(rv$time),]
}

plotsim <- function(sim)
{
  events <- split(sim, sim$event)
  I <- sapply(sim$time, function(t) {
    sum(events$i$time < t) - sum(events$r$time < t)
  } )  
  plot(sim$time, I, type='l')
}