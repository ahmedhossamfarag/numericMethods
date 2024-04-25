
## linear congruential generator
## x(i) = a*x(i−1) + b (mod m)
## ui = x(i) / m
## stretch (s(1) − s(0))r + s(0)
lcg <- function(a = 7^5, b = 0, m = 2^31-1, x0 = 6754, n = 10, s = c(0, 1)){
  g <- c(x0 / m)
  length(g) <- n
  x <- x0;
  for(i in 2:n){
    x <- (a * x + b) %% m
    g[i] <- x / m
  }
  g * s[2] - s[1] + s[1]
}

## plotting linear congruential generator

lcgPlot <- function(x01 = 6565, x02 = 1458){
  plot(lcg(x0 = x01, n = 1000), lcg(x0 = x02, n = 1000))
}


## random numbers generator
## x(i) = x(i-2)*(a*x(i-1) mod m) mod m
## ui = x(i) / m
## stretch (s(1) − s(0))r + s(0)
rng <- function(a = 7^5, b = 0, m = 2^31-1, x0 = 6754, x1 = 4857, n = 10, s = c(0, 1)){
  g <- c(x0 / m, x1 / m)
  length(g) <- n
  xi2 <- x0
  xi1 <- x1
  for(i in 2:n){
    x <- (xi2*((a * xi1 + b) %% m)) %% m
    xi2 <- xi1
    xi1 <- x
    g[i] <- x / m
  }
  g * s[2] - s[1] + s[1]
}

## plotting random numbers generator

rngPlot <- function(x10 = 6565, x11 = 1458, x20 = 4258, x21 = 1634){
  plot(rng(x0 = x10, x1 = x11, n = 1000), rng(x0 = x20, x1 = x21, n = 1000))
}
