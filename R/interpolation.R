## Newton’s divided differences
## vector x
## vector y
ndd <- function(x,y){
  if(!is.vector(x) || !is.vector(y) || length(x) != length(y) || length(x) < 2){
    stop("invalid args")
  }
  n <- length(x)
  cof <- c(y[1])
  for(i in 2:n){
    y <- (y[-1] - y[-length(y)]) / (x[i:n] - x[1:(n-i+1)])
    cof[i] = y[1]
  }
  cof
}

## Chebyshev interpolation nodes
## interval [a,b]
## no of points n
chnd <- function(a = -1, b = 1, n=10){
  if(!is.numeric(a) || !is.numeric(b) || !is.numeric(n)){
    stop("invalid args")
  }
  ((b + a)/2) + ((b - a)/2) * cos(((2*(1:n) - 1)*pi)/(2*n))
}
