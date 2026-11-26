# Real-data design: Boston housing predictors with a simulated response.
# Base R plus MASS (ships with R), run with Rscript.

library(MASS)

set.seed(2026)
Xfull <- as.matrix(Boston[, names(Boston) != "medv"])
N <- nrow(Xfull); p <- ncol(Xfull)
n <- 100; k <- 3; nrep <- 500
beta0 <- c(1.5, 1, 0.5)

sym_pow <- function(A, pw) {
  e <- eigen(A, symmetric = TRUE)
  e$vectors %*% diag(pmax(e$values, 0)^pw) %*% t(e$vectors)
}

lasso_supports <- function(A, y, nlam = 30) {
  p <- ncol(A)
  lam <- max(abs(crossprod(A, y))) * exp(seq(0, log(0.01), length.out = nlam))
  cn <- colSums(A^2)
  b <- numeric(p); r <- y
  out <- matrix(FALSE, nlam, p)
  for (k in 1:nlam) {
    for (it in 1:200) {
      old <- b
      for (j in 1:p) {
        rj <- r + A[, j] * b[j]
        z <- sum(A[, j] * rj)
        b[j] <- sign(z) * max(abs(z) - lam[k], 0) / cn[j]
        r <- rj - A[, j] * b[j]
      }
      if (max(abs(b - old)) < 1e-6) break
    }
    out[k, ] <- b != 0
  }
  out
}

recovers <- function(A, y, truth) {
  S <- lasso_supports(A, y)
  any(apply(S, 1, function(s) identical(which(s), as.integer(truth))))
}

chol_fit <- function(X, y, ord, truth) {
  qrX <- qr(X[, ord])
  Q <- qr.Q(qrX) %*% diag(sign(diag(qr.R(qrX))))
  recovers(Q, y, sort(match(truth, ord)))
}

methods <- c("Lasso", "Puffer", "ZCA", "PCA", "Cholesky (dataset order)",
             "Cholesky (true variables first)", "Cholesky (ranked by |OLS|)")
hits <- matrix(NA, nrep, length(methods), dimnames = list(NULL, methods))

for (r in 1:nrep) {
  repeat {
    rows <- sample(N, n)
    if (all(apply(Xfull[rows, ], 2, sd) > 0)) break
  }
  X <- scale(Xfull[rows, ])
  truth <- sort(sample(p, k))
  beta <- numeric(p)
  beta[truth] <- beta0
  y <- drop(X %*% beta + rnorm(n))
  y <- y - mean(y)
  sv <- svd(X)
  first <- c(truth, setdiff(1:p, truth))
  hits[r, ] <- c(
    recovers(X, y, truth),
    recovers(sv$u %*% t(sv$v), drop(sv$u %*% diag(1 / sv$d) %*% t(sv$u) %*% y), truth),
    recovers(X %*% sym_pow(crossprod(X), -1/2), y, truth),
    recovers(X %*% sv$v %*% diag(1 / sv$d), y, truth),
    chol_fit(X, y, 1:p, truth),
    chol_fit(X, y, first, truth),
    chol_fit(X, y, order(-abs(solve(crossprod(X), crossprod(X, y)))), truth)
  )
}

m <- colMeans(hits)
res <- data.frame(method = methods,
                  recovery = sprintf("%.2f (%.2f)", m, sqrt(m * (1 - m) / nrep)))
cat("Boston predictors, n =", n, ", p =", p, ",", k, "true variables drawn at random\n")
cat("Share of", nrep, "data sets in which the Lasso path contains the true support (standard error)\n\n")
print(res)