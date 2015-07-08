# Plot functions for epidemics

plotsim <- function(sim,...)
{
  sim <- sim[sim$i < Inf,]
  eventTimes <- sort(c(sim$i,sim$r))
  I <- sapply(eventTimes, function(t) {
    sum(sim$i <= t & t < sim$r)
  } )  

  # Plot curve
  plot(eventTimes, I, type='s',...)
  
  # Plot events
  points(sim$r, rep(0,nrow(sim)), pch=19, cex=1, col='blue')
  points(sim$i, rep(0,nrow(sim)), pch=19, cex=1, col='red')
  
}
