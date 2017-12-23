#' @title Print n prime numbers
#'
#' @description Function by CloudStat posted on R-bloggers on November 28, 2011. Link: https://www.r-bloggers.com/prime-number-in-r-language-cloudstat/
#'
#' @param n number of prime numbers to be printed
#'
#' @export

prime <- function(n) {
  n <- as.integer(n)
  if(n > 1e8) stop("n too large")
  primes = rep(TRUE, n)
  primes[1] = FALSE
  last.prime = 2L
  fsqr = floor(sqrt(n))
  while (last.prime <= fsqr)
  {
    primes[seq.int(2L*last.prime, n, last.prime)] = FALSE
    sel = which(primes[(last.prime+1):(fsqr+1)])
    if(any(sel)){
      last.prime = last.prime + min(sel)
    }else last.prime = fsqr+1
  }
  which(primes)
}
