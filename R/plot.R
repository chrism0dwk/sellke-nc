# Plot functions for epidemics

plot.epidemic <- function(epidemic,type='i',...)
{
  epidemic <- lapply(epidemic, function(x) x[epidemic$i < Inf])
  N <- length(epidemic$i)
  eventTimes <- sort(c(epidemic$i,epidemic$r))
  
  n <- NULL
  n <- switch (type,
    i=sapply(eventTimes, function(t) {
      sum(epidemic$i <= t & t < epidemic$r)
    } ),
    r=sapply(eventTimes, function(t) {
      sum(epidemic$r <= t)
    })
  )

  # Plot curve
  plot(eventTimes, n, type='s', ylab=switch(type, i="Infected", r="Removed"),...)
  
  # Plot events
  points(epidemic$r, rep(0,N), pch=19, cex=1, col='blue')
  points(epidemic$i[epidemic$r<Inf], rep(0,sum(epidemic$r<Inf)), pch=19, cex=1, col='red')
  points(epidemic$i[epidemic$r==Inf],rep(0,sum(epidemic$r==Inf)), pch=19, cex=1, col='green')
  
}
