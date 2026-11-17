.pmlr_run_isolated <- function(X, ics, q, K, c, tol, maxit,
                               maxNumRestart, replix, fixthrehsolds,
                               loading, latentmeans) {
  root <- system.file("julia", package = "PMLCFAEFA")
  if (!nzchar(root)) stop("Compiled Julia bundle is missing.", call. = FALSE)
  runner <- file.path(root, "bin", "pmlr-runner")
  if (!file.exists(runner))
    stop("The isolated Julia runner is missing. Install version 0.1.1 or newer.", call. = FALSE)

  input <- tempfile("pmlr_input_")
  output <- tempfile("pmlr_output_")
  on.exit(unlink(c(input, output)), add = TRUE)

  fmt <- function(x) formatC(as.numeric(x), digits = 17, format = "g",
                              decimal.mark = ".")
  has_l <- is.matrix(loading)

  header <- c(
    nrow(X), ncol(X), q, K, c, ncol(ics),
    maxit, maxNumRestart, replix, as.integer(has_l),
    fmt(tol), fmt(fixthrehsolds), fmt(latentmeans)
  )

  writeLines(
    c(
      paste(header, collapse = " "),
      as.character(as.integer(X)),
      fmt(ics),
      if (has_l) fmt(loading)
    ),
    con = input,
    useBytes = TRUE
  )

  message("Running PMLCFAEFA in a separate process. RStudio will remain open if Julia crashes.")
  debug <- identical(Sys.getenv("PMLR_DEBUG_GDB"), "1")
  if (debug) {
    gdb <- Sys.which("gdb")
    if (!nzchar(gdb)) stop("gdb is required for PMLR_DEBUG_GDB=1", call. = FALSE)
    executable <- gdb
    arguments <- c("-q", "-batch", "-ex", shQuote("set pagination off"),
                   "-ex", "run", "-ex", shQuote("bt full"), "--args",
                   shQuote(c(runner, root, input, output)))
  } else {
    executable <- runner
    arguments <- shQuote(c(root, input, output))
  }
  log <- suppressWarnings(
    system2(executable, args = arguments, stdout = TRUE, stderr = TRUE)
  )
  status <- attr(log, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L || !file.exists(output)) {
    diagnostic <- paste(tail(log, 30L), collapse = "\n")
    stop(sprintf("Compiled Julia model failed (exit code %s):\n%s",
                 status, diagnostic), call. = FALSE)
  }
  result <- scan(output, what = double(), quiet = TRUE)
  if (!length(result) || result[1] < 1L || length(result) != result[1] + 1)
    stop("The Julia runner returned an invalid result.", call. = FALSE)
  result[-1L]
}
