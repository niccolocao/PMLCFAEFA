.pmlr_initialize <- local({
  initialized <- FALSE
  function() {
    if (initialized) return(invisible(TRUE))
    root <- system.file("julia", package = "PMLCFAEFA")
    if (!nzchar(root) || !dir.exists(root)) {
      stop("Compiled Julia libraries are missing. Install a built PMLCFAEFA package containing the Julia bundle.", call. = FALSE)
    }
    .Call(C_pmlr_init, root)
    initialized <<- TRUE
    invisible(TRUE)
  }
})

PMLR_status <- function() {
  root <- system.file("julia", package = "PMLCFAEFA")
  list(bundle_present = nzchar(root) && dir.exists(root),
       runner_present = nzchar(root) &&
         file.exists(file.path(root, "bin", "pmlr-runner")),
       bundle_path = root)
}

PMLcfaefa <- function(X, p = ncol(X), q, K, c = max(X),
                      ics = matrix(1, nrow(X), 1), tol = 1e-6,
                      maxit = 1000L, maxNumRestart = 20L, replix = 10L,
                      fixthrehsolds = 0, L_str = FALSE,
                      latentmeans = 1, parametr = FALSE) {
  X <- as.matrix(X)
  if (!is.numeric(X) || anyNA(X) || any(!is.finite(X)) ||
      any(X != round(X))) stop("X must contain finite integer-valued categories.", call. = FALSE)
  if (length(dim(X)) != 2L || nrow(X) < 2L || ncol(X) < 2L)
    stop("X must be a matrix with at least two rows and two columns.", call. = FALSE)
  if (length(p) != 1L || is.na(p) || p != ncol(X))
    stop("p must equal ncol(X).", call. = FALSE)
  if (any(X < 1) || any(X > c))
    stop("X must contain categories between 1 and c.", call. = FALSE)
  if (!identical(parametr, FALSE))
    stop("Custom 'parametr' initialization is not yet supported by the native R interface.", call. = FALSE)
  if (length(q) != 1L || length(K) != 1L || q < 1 || K < 1 ||
      q != as.integer(q) || K != as.integer(K) || K > p)
    stop("q and K must be valid positive integer factor counts.", call. = FALSE)
  if (length(c) != 1L || is.na(c) || c != as.integer(c) || c < 2L)
    stop("c must be an integer >= 2.", call. = FALSE)
  if (isFALSE(L_str) && p %% q != 0L)
    stop("Without L_str, p must be divisible by q in the original Julia implementation.", call. = FALSE)
  if (c < 4L && !(c == 2L && latentmeans == 1 && fixthrehsolds == -1))
    stop("For this model configuration, the original Julia code requires c >= 4.", call. = FALSE)

  ics <- as.matrix(ics)
  storage.mode(ics) <- "double"
  if (nrow(ics) != nrow(X) || ncol(ics) < 1L || anyNA(ics) || any(!is.finite(ics)))
    stop("ics must be a finite numeric matrix with nrow(X) rows.", call. = FALSE)
  if (isFALSE(L_str)) {
    loading <- FALSE
  } else {
    loading <- as.matrix(L_str)
    storage.mode(loading) <- "double"
    if (!identical(dim(loading), c(ncol(X), as.integer(q))) || anyNA(loading) ||
        any(!is.finite(loading))) stop("L_str must be a finite p x q matrix.", call. = FALSE)
  }
  if (!is.finite(tol) || tol <= 0 || maxit < 1 || maxNumRestart < 1 || replix < 1)
    stop("tol must be positive; maxit, maxNumRestart and replix must be >= 1.", call. = FALSE)
  if (latentmeans == 0 && !(fixthrehsolds %in% 0:7) ||
      latentmeans == 1 && !(fixthrehsolds %in% c(-1, 0:5)) ||
      !(latentmeans %in% c(0, 1, 2)))
    stop("Unsupported latentmeans/fixthrehsolds combination.", call. = FALSE)

  raw <- .pmlr_run_isolated(
    X, ics, as.integer(q), as.integer(K), as.integer(c),
    as.double(tol), as.integer(maxit), as.integer(maxNumRestart),
    as.integer(replix), as.double(fixthrehsolds), loading,
    as.double(latentmeans)
  )
  structure(list(raw = raw, status_code = raw[length(raw) - 1L],
                 restart_count = raw[length(raw)],
                 n = nrow(X), p = ncol(X), q = q, K = K, c = c),
            class = "pmlcfaefa_fit")
}

print.pmlcfaefa_fit <- function(x, ...) {
  cat("PML CFA-EFA fit\n")
  cat("Observations:", x$n, " Items:", x$p, " CFA factors:", x$q,
      " EFA factors:", x$K, "\n")
  cat("Status code:", x$status_code, " Restart count:", x$restart_count, "\n")
  cat("Raw result length:", length(x$raw), "\n")
  invisible(x)
}
