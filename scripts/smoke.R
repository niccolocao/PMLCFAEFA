library(PMLCFAEFA)
stopifnot(PMLR_status()$bundle_present)
cat("Compiled Julia bundle found.\n")
stopifnot(!any(grepl("pmlr_run_isolated", deparse(body(PMLcfaefa)), fixed = TRUE)))

set.seed(511)
n <- 120L
z <- rnorm(n)
X <- replicate(5, {
  as.integer(cut(z + rnorm(n, sd = 0.75),
                 breaks = c(-Inf, -0.8, 0, 0.8, Inf)))
})

cat("Testing an actual in-process native Julia model fit.\n")
fit <- PMLcfaefa(
  X, q = 1L, K = 1L, c = 4L,
  ics = matrix(1, n, 1),
  tol = 1e-2, maxit = 50L, maxNumRestart = 2L, replix = 2L,
  verbose = TRUE, fast = TRUE
)
stopifnot(inherits(fit, "pmlcfaefa_fit"), length(fit$raw) > 10)
cat("Julia returned a model result to R.\n")
