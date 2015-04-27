addlog <- function(a,b) {
    if(b > a + 200) return(b)
    if(a > b + 200) return(a)
    a + log(1+exp(b-a))
}

nct <-
    function(x, dof, ncp, n_iterations=600)
{
      tol <- 1e-14

      prob = pnorm(-ncp, 0, 1, log.p=TRUE)
      value <- tol+1
      lastvalue <- value
      y <- x*x/(x*x + dof)
      result <- matrix(nrow=n_iterations, ncol=6)
      p <- -ncp*ncp/2
      q = -ncp*ncp/2 - 0.5 * log(2) - lgamma(3/2) + log(ncp)
      j <- 0
      result <- prob
      while(j < n_iterations) {
          lastvalue = value
          if(j>0) {
              p = p + 2*log(ncp) - log(2) - log(j)
              q = q + 2*log(ncp) - log(2) - log(j+0.5)
          }
          value = addlog(p + pbeta(y, j+0.5, dof/2, log.p=TRUE),
                         q + pbeta(y, j+1, dof/2, log.p=TRUE))
          prob = addlog(prob, value - log(2))
          j <- j + 1
          result <- c(result, value)
      }

      message(exp(prob))
      result
  }
