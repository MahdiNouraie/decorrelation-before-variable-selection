# Simulation for Section 6: grouped designs with p < n, as in the DVS paper.
# Base R only, run with Rscript. Writes fig_simulation.pdf to the working directory.

set.seed(2026)
n <- 100; p <- 50; nrep <- 500

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

methods <- c("Lasso", "Puffer", "ZCA", "PCA",
             "Cholesky (given order)", "Cholesky (ranked by |OLS|)")

run_design <- function(gsize, rho, beta0) {
  ng <- p / gsize
  Sig <- matrix(0, p, p)
  for (g in 1:ng) {
    i <- (g - 1) * gsize + 1:gsize
    Sig[i, i] <- rho[g]
  }
  diag(Sig) <- 1
  Sig <- chol(Sig)
  starts <- (0:(ng - 1)) * gsize
  
  res <- data.frame(method = methods)
  for (s in c("first", "last")) {
    truth <- starts + if (s == "first") 1 else gsize
    beta <- numeric(p)
    beta[truth] <- beta0
    hits <- matrix(NA, nrep, 6, dimnames = list(NULL, methods))
    for (r in 1:nrep) {
      X <- scale(matrix(rnorm(n * p), n, p) %*% Sig)
      y <- drop(X %*% beta + rnorm(n))
      y <- y - mean(y)
      sv <- svd(X)
      hits[r, ] <- c(
        recovers(X, y, truth),
        recovers(sv$u %*% t(sv$v), drop(sv$u %*% diag(1 / sv$d) %*% t(sv$u) %*% y), truth),
        recovers(X %*% sym_pow(crossprod(X), -1/2), y, truth),
        recovers(X %*% sv$v %*% diag(1 / sv$d), y, truth),
        chol_fit(X, y, 1:p, truth),
        chol_fit(X, y, order(-abs(solve(crossprod(X), crossprod(X, y)))), truth)
      )
    }
    m <- colMeans(hits)
    res[[s]] <- m
    res[[paste0(s, "_se")]] <- sqrt(m * (1 - m) / nrep)
  }
  res
}

show <- function(res) {
  data.frame(method = res$method,
             first = sprintf("%.2f (%.2f)", res$first, res$first_se),
             last = sprintf("%.2f (%.2f)", res$last, res$last_se))
}

A <- run_design(10, c(0.2, 0.4, 0.6, 0.8, 0.9), c(4, 3.5, 3, 2.5, 2))
B <- run_design(5, seq(0.5, 0.9, length.out = 10), seq(4, 1, length.out = 10))

cat("Share of", nrep, "simulations in which the Lasso path contains the true support (standard error)\n")
cat("n =", n, ", p =", p, "\n\n")
cat("Design A: five groups of ten, one true variable per group\n")
print(show(A))
cat("\nDesign B: ten groups of five, one true variable per group, decaying signal\n")
print(show(B))

panel <- function(res, title, labels) {
  y <- length(methods):1
  plot(NA, xlim = c(0, 1), ylim = c(0.5, length(methods) + 0.5), yaxt = "n",
       xlab = "Share with exact recovery", ylab = "", main = title, bty = "l")
  abline(h = y, col = "gray85")
  for (s in c("first", "last")) {
    yy <- y + if (s == "first") 0.15 else -0.15
    lo <- pmax(res[[s]] - 1.96 * res[[paste0(s, "_se")]], 0)
    hi <- pmin(res[[s]] + 1.96 * res[[paste0(s, "_se")]], 1)
    segments(lo, yy, hi, yy)
    points(res[[s]], yy, pch = if (s == "first") 16 else 2, cex = 1.1)
  }
  if (labels) axis(2, at = y, labels = methods, las = 1)
}

pdf("fig_simulation.pdf", width = 7.5, height = 3.9)
layout(matrix(1:2, 1), widths = c(1.75, 1))
par(mar = c(4, 12, 2, 1), oma = c(2, 0, 0, 0))
panel(A, "Design A", TRUE)
par(mar = c(4, 1, 2, 1))
panel(B, "Design B", FALSE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("bottom", pch = c(16, 2), horiz = TRUE, bty = "n",
       legend = c("True variable first", "True variable last"))
dev.off()