## Bisection method
## function f
## interval [a, b]
## tolerance t
bis <- function(f, a, b, t = 0.0001){
  if(!is.function(f) || !is.numeric(a) || !is.numeric(b) || !is.numeric(t) || b <= a){
    stop("invalid arguments")
  }
  fa <- f(a)
  fb <- f(b)
  if(fa * fb > 0){
    stop("invalid interval")
  }
  if(fa == 0){
    return(a)
  }
  if(fb == 0){
    return(b)
  }
  print(format(c("a", "b", "c", "fc"), width = 10), quote = FALSE)
  while(b-a > t){
    c <- (a+b) / 2
    fc <- f(c)
    print(format(c(a, b, c, fc), width = 10, digits = 3,justify="left"), quote = FALSE)
    if(fc == 0){
      return(c)
    }
    else if(fc * fa < 0){
      b <- c
      fb <- fc
    }
    else{
      a <- c
      fa <- fc
    }
  }
  return((a+b) / 2)
}

## fixed point iteration
## function g
## initial guess x0
## tolerance t
fpi <- function(g, x0, t = 0.0001){
  if(!is.function(g) || !is.numeric(x0) || !is.numeric(t)){
    stop("invalid arguments")
  }
  cat("i", "x", "\n")
  cat(0, x0, "\n")
  xp <- NA
  x <- x0
  i <- 1
  repeat{
    xp <- x
    x <- g(x)
    if(abs(x - xp) <= t){
      break
    }
    cat(i, x, "\n")
    i <- i + 1
    if(is.nan(x) || is.infinite(x))
    {
      stop("invalid iteration")
    }
  }
  x
}

## newton method iteration
## function f
## function derivative fd
## initial guess x0
## tolerance t
nmi <- function(f, fd, x0, t = 0.0001){
  if(!is.function(f) || !is.function(fd) || !is.numeric(x0) || !is.numeric(t)){
    stop("invalid arguments")
  }
  cat("i", "x", "\n")
  cat(0, x0, "\n")
  x <- x0
  i <- 1
  repeat{
    fx <- f(x)
    if(abs(fx) <= t){
      break
    }
    x <- x - (fx / fd(x))
    cat(i, x, "\n")
    i <- i + 1
    if(is.nan(x) || is.infinite(x))
    {
      stop("invalid iteration")
    }
  }
  x
}

## secant method iteration
## function f
## initial guess x0, x1
## tolerance t
smi <- function(f, x0, x1, t = 0.0001){
  if(!is.function(f) || !is.numeric(x0) || !is.numeric(x1) || !is.numeric(t)){
    stop("invalid arguments")
  }
  cat("i", "x", "\n")
  cat(0, x0, "\n")
  cat(1, x1, "\n")
  xi1 <- x0
  xi <- x1
  i <- 2

  repeat{
    fxi1 <- f(xi1)
    fxi <- f(xi)
    if(abs(fxi) <= t){
      break
    }
    xt <- xi
    xi <- xi - (fxi * (xi - xi1))/(fxi - fxi1)
    xi1 <- xt
    cat(i, xi, "\n")
    i <- i + 1
    if(is.nan(xi) || is.infinite(xi))
    {
      stop("invalid iteration")
    }
  }
  xi
}


