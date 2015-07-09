# Implements naive Sellke reconstruction of an epidemic

#' @title Recalculate infection times using Sellke construction.
#' 
#' @description
#' \code{reconstruct} returns a rescaled \code{epidemic}, conditioned on
#' removal times, Sellke thresholds, and a new \code{beta}.
#' 
#' @details 
#' This function takes a Sellke-defined epidemic, and rescales the infection times based on a 
#' new transmission rate \code{beta}.  If the rescaling results in either a discontinuous
#' epidemic, or a violation of I < R, then the original epidemic is returned with a warning.
#' 
#' @param beta The new infection rate
#' @param epidemic an epidemic simulated using \code{sellkeSim}
#' @param condition.on.I Should we condition on only known infections?
#' 
#' @return A new epidemic, with rescaled infection times.
reconstruct <- function(beta, epidemic, condition.on.I=F)
{
  epidemic <- data.frame(nu=epidemic$nu, i=epidemic$i, r=epidemic$r, row.names=rownames(epidemic))  
  if(condition.on.I) {
    # If we're not allowing 
    epidemic$nu[epidemic$i == Inf] <- Inf
  }
  
  oldI <- epidemic$i
  
  S <- nrow(epidemic) - 1
  I <- 1
  R <- 0
  t <- 0
  
  I1 <- which.min(epidemic$i)
  nu <- epidemic$nu
  names(nu) <- rownames(epidemic)
  seeds <- sort(nu[-I1])
  
  cumPressure <- 0
  maxR <- max(epidemic$r[epidemic$r < Inf])
  
  while (t < maxR)
  {
    if(I == 0) { # Check to make sure we have infectious pressure.
      warning("Discontinuous epidemic. Abort.",immediate. = T)
      epidemic$i <- oldI
      break
    }
    
    minSeed <- seeds[1]
    nextInfec <- ifelse(length(seeds)>0, (minSeed - cumPressure) / (beta * I) + t, Inf)
    toRemove <- which(epidemic$r == min(epidemic$r[epidemic$r > t]))
    nextRemove <- epidemic$r[toRemove]
    tOld <- t
  
    if (nextRemove > maxR & nextInfec > maxR) break
    
    if (nextInfec < nextRemove)
    {
      S <- S - 1
      I <- I + 1
      t <- nextInfec
      seeds <- seeds[-1]
      epidemic[names(minSeed),'i'] <- t
    }
    else
    {
      I <- I - 1
      R <- R + 1
      t <- nextRemove
      # Check if I > R
      if(epidemic$r[toRemove] < epidemic$i[toRemove]) {
        warning("Invalidated epidemic. Abort.",immediate. = T)
        epidemic$i <- oldI
        break
      }
    }
    cumPressure <- cumPressure + beta * I * (t - tOld)
  }
  message("Max IIP: ", cumPressure)
  epidemic$nu <- epidemic$nu # Reset nu's
  class(epidemic) <- "epidemic"
  epidemic
}
