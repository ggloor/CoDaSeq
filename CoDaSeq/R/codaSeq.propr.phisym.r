codaSeq.propr.phisym <- function (X)
{
  Cov    <- stats::var(X)
  tmp    <- 2 * Cov / outer(diag(Cov), diag(Cov), "+")
  return((1-tmp)/(1+tmp))
}

