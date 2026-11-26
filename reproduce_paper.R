# Code for "Decorrelation Before Variable Selection: Left and Right
# Transforms and Their Trade-Off". Base R only, run with Rscript.

tol <- 1e-8

sym_pow <- function(A, pw) {
  e <- eigen(A, symmetric = TRUE)
  e$vectors %*% diag(pmax(e$values, 0)^pw) %*% t(e$vectors)
}

## 1. the (T1, N, T2) example from Section 3
Sigma <- matrix(c(1.0, 0.6, 0.5,
                  0.6, 1.0, 0.0,
                  0.5, 0.0, 1.0), 3, 3, byrow = TRUE,
                dimnames = list(c("T1", "N", "T2"), c("T1", "N", "T2")))
beta  <- c(T1 = 1, N = 0, T2 = 1)
gamma <- drop(chol(Sigma) %*% beta)    # chol() gives R with R'R = Sigma

cat("Example (T1, N, T2)\n")
cat("  gamma =", round(gamma, 3), "\n")
cat("  N is active:", abs(gamma["N"]) > tol, "\n\n")

## 2. checks behind Table 3 and Proposition 3
set.seed(2026)
n <- 60
p <- 5
X <- matrix(rnorm(n * p), n, p) %*% chol(0.7^abs(outer(1:p, 1:p, "-")))
X <- scale(X, scale = FALSE)
X <- X / rep(sqrt(colSums(X^2)), each = n)   # unit-norm columns
beta <- c(2, 1.5, 0, 0, 0)                   # true variables first
eps  <- rnorm(n)
y    <- drop(X %*% beta + eps)
S    <- crossprod(X)

sv  <- svd(X)
U   <- sv$u
D   <- sv$d
V   <- sv$v
qrX <- qr(X)
sg  <- sign(diag(qr.R(qrX)))
Q   <- qr.Q(qrX) %*% diag(sg)
Rm  <- diag(sg) %*% qr.R(qrX)

# W for the projection transform, with Lambda = lambda_min * I
lam <- min(eigen(S, symmetric = TRUE)$values)
Z   <- qr.Q(qrX, complete = TRUE)[, (p + 1):(2 * p)] %*% sym_pow(diag(p) - lam * solve(S), 1/2)
W   <- X %*% solve(S) * sqrt(lam) + Z

uncorrelated <- function(M) {
  G <- crossprod(M)
  off <- G
  diag(off) <- 0
  max(abs(off)) < tol * max(diag(G))
}
white <- function(C) max(abs(C - diag(nrow(C)))) < tol

left <- function(F) {
  list(design = F %*% X, noise = F %*% t(F),
       keeps = max(abs(F %*% y - F %*% X %*% beta - F %*% eps)) < tol)
}
right <- function(G) {
  g <- drop(solve(G, beta))
  keeps <- if (max(abs(g - beta)) < tol) "yes"
           else if (max(abs(G - diag(diag(G)))) < tol) "rescaled"
           else "no"
  list(design = X %*% G, noise = diag(n), gamma = g, keeps = keeps)
}

L <- list("Puffer" = left(U %*% diag(1 / D) %*% t(U)),
          "Cholesky Q'" = left(t(Q)),
          "Projection W'" = left(t(W)))
R <- list("ZCA" = right(sym_pow(S, -1/2)),
          "PCA" = right(V %*% diag(1 / D)),
          "Cholesky (Gram-Schmidt)" = right(solve(Rm)),
          "Weighted Lasso" = right(diag(1 / (1:5))))

make_row <- function(nm, side, m) {
  data.frame(method = nm, side = side, uncorrelated = uncorrelated(m$design),
             keeps_beta = if (is.logical(m$keeps)) ifelse(m$keeps, "yes", "no") else m$keeps,
             white_noise = white(m$noise))
}
tab <- rbind(do.call(rbind, Map(make_row, names(L), "left", L)),
             do.call(rbind, Map(make_row, names(R), "right", R)))
rownames(tab) <- NULL
cat("Table 3 properties (n = 60, p = 5)\n")
print(tab)

# left Cholesky keeps the least-squares estimate
cat("\nleast squares unchanged by Q':",
    max(abs(solve(S, crossprod(X, y)) - solve(Rm, crossprod(Q, y)))) < tol, "\n")

# true variables first: gamma should be 0 on the noise variables for Cholesky only
cat("\nmax |gamma| on the noise variables\n")
for (nm in c("ZCA", "PCA", "Cholesky (Gram-Schmidt)"))
  cat(sprintf("  %-24s %.3f\n", nm, max(abs(R[[nm]]$gamma[3:5]))))

# Proposition 3
cat("\nProposition 3\n")
cat("  W'W = I:", max(abs(crossprod(W) - diag(p))) < tol, "\n")
cat("  X'W W'X = lambda_min * I:", max(abs(t(X) %*% W %*% t(W) %*% X - diag(lam, p))) < tol, "\n")
cat("  X'X - Lambda is PSD:", min(eigen(S - diag(lam, p), symmetric = TRUE)$values) > -tol, "\n")
v <- rbind(reduced = rep(1 / lam, p), full = diag(solve(S)))
colnames(v) <- paste0("beta", 1:p)
cat("  variances (sigma^2 = 1):\n")
print(round(v, 3))
cat("  reduced >= full for all j:", all(v["reduced", ] >= v["full", ] - tol), "\n")
