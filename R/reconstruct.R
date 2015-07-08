# Implements naive Sellke reconstruction of an epidemic

reconstruct <- function(beta, sim)
{
  oldI <- sim$i
  
  S <- nrow(sim) - 1
  I <- 1
  R <- 0
  t <- 0
  
  I1 <- which.min(sim$i)
  nu <- sim$nu
  names(nu) <- rownames(sim)
  seeds <- sort(nu[-I1])
  
  cumPressure <- 0
  maxR <- max(sim$r[sim$r < Inf])
  
  while (t < maxR)
  {
    if(I == 0) { # Check to make sure we have infectious pressure.
      warning("Discontinuous epidemic. Abort.",immediate. = T)
      browser()
      sim$i <- oldI
      break
    }
    
    minSeed <- seeds[1]
    nextInfec <- ifelse(length(seeds)>0, (minSeed - cumPressure) / (beta * I) + t, Inf)
    toRemove <- which(sim$r == min(sim$r[sim$r > t]))
    nextRemove <- sim$r[toRemove]
    tOld <- t
  
    if (nextRemove > maxR & nextInfec > maxR) break
    
    if (nextInfec < nextRemove)
    {
      S <- S - 1
      I <- I + 1
      t <- nextInfec
      seeds <- seeds[-1]
      sim[names(minSeed), 'i'] <- t
    }
    else
    {
      I <- I - 1
      R <- R + 1
      t <- nextRemove
      # Check if I > R
      if(sim$r[toRemove] < sim$i[toRemove]) {
        warning("Invalidated epidemic. Abort.",immediate. = T)
        sim$i <- oldI
        break
      }
    }
    
    
      
    cumPressure <- cumPressure + beta * I * (t - tOld)
  }
  sim
}
