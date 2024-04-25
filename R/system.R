## Gaussian elimination
## matrix m
gausseli <- function(m){
  if(!is.matrix(m)){
    stop("invalid args")
  }
  for(j in 1:(min(dim(m))-1)){
   for (i in (j+1):nrow(m)) {
     cat(sprintf("subtract %.2f/%.2f * row %i from row %i \n", m[i,j] , m[j,j], j, i))
     m[i,] <- m[i,] - (m[j,] * m[i,j] / m[j,j])
     print(m)
   }
  }
  m
}

## LU factorization
## matrix m
luf <- function(m){
  if(!is.matrix(m)){
    stop("invalid args")
  }
  l <- diag(1, nrow(m), ncol(m))
  for(j in 1:(min(dim(m))-1)){
    for (i in (j+1):nrow(m)) {
      l[i, j] <- m[i,j] / m[j,j]
      m[i,] <- m[i,] - (m[j,] * l[i, j])
    }
  }
  list(L = l, U = m)
}

## PA = LU factorization
## matrix m
pluf <- function(m){
  if(!is.matrix(m) || nrow(m) != ncol(m)){
    stop("invalid args")
  }
  n <- nrow(m)
  p <- 1:n

  for(j in 1:(n-1)){
    mx <- max.col(t(abs(m[j:n,j])))[1] + j - 1
    if(mx != j){

      tmp <- m[mx,]
      m[mx,] <- m[j,]
      m[j,] <- tmp

      tmp <- p[mx]
      p[mx] <- p[j]
      p[j] <- tmp

    }
    for (i in (j+1):n) {
      d <- m[i,j] / m[j,j]
      rng <- j:n
      m[i,rng] <- m[i,rng] - (m[j,rng] * d)
      m[i, j] <- d
    }
  }

  mp <- matrix(0, n, n)
  ind <- array(c(1:n, p), dim = c(n,2))
  mp[ind] <- 1

  l <- diag(1, n, ncol(m))
  l[lower.tri(l)] <- m[lower.tri(m)]

  m[lower.tri(m)] <- 0

  list(P = mp, L = l, U = m)
}

## Jacobi Method
## matrix A
## column b
## initial guess x0
## no of iterations n
jacobi <- function(A, b, x0, n = 10){
  if(!is.matrix(A) || !is.matrix(b) || !is.matrix(x0) || !all(dim(b) == dim(x0), dim(b) == c(nrow(A), 1)) || !is.numeric(n)){
    stop("invalid args")
  }
  x = x0
  LU = A
  diag(LU) <- 0
  D = diag(diag(A)^-1, nrow(A), ncol(A))
  for(i in 1:n){
    x = D %*% (b - (LU %*% x))
    cat(sprintf("[%i] %f", i, x), "\n")
  }
  x
}

## Gauss–Seidel Method
## matrix A
## column b
## initial guess x0
## no of iterations n
gsm <- function(A, b, x0, n = 10){
  if(!is.matrix(A) || !is.matrix(b) || !is.matrix(x0) || !all(dim(b) == dim(x0), dim(b) == c(nrow(A), 1)) || !is.numeric(n)){
    stop("invalid args")
  }
  x = x0
  LU = A
  diag(LU) <- 0
  L = LU
  L[upper.tri(L)] = 0
  U = LU
  U[lower.tri(U)] = 0
  D = diag(A)^-1
  for(i in 1:n){
    for(j in 1:nrow(A)){
      x[j,] = D[j] * (b[j,] - (L[j,] %*% x) - (U[j,] %*% x))
    }
    cat(sprintf("[%i] %f", i, x), "\n")
  }
  x
}

## Cholesky factorization
## matrix A
chf <- function(A){
  if(!is.matrix(A) || !isSymmetric(A) || nrow(A) < 2 || !all(eigen(A, only.values = TRUE)$values >= 0)){
    stop("invalid args")
  }
  fac <- function(m){
    n = nrow(m)
    r = matrix(0, n, n)
    r[1,1] = sqrt(m[1,1])
    r[1, 2:n] = m[2:n, 1] / r[1,1]
    if(n == 2){
      r[2,2] = m[2,2] - m[1,2]^2/m[1,1]
    }
    else{
      u = r[1, 2:n]
      r[2:n, 2:n] = fac(m[2:n, 2:n] - (u %*% t(u)))
    }
    r
  }
  fac(A)
}
