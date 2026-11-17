if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("CTN-0094/public.ctn0094data")
remotes::install_url("https://github.com/niccolocao/PMLCFAEFA/releases/download/v0.1.0/PMLCFAEFA_0.1.0.tar.gz")

library(PMLCFAEFA)
data("everybody", package = "public.ctn0094data")
data("randomization", package = "public.ctn0094data")
data("rbs", package = "public.ctn0094data")

ids <- intersect(everybody$who[as.character(everybody$project) == "51"],
                 randomization$who)
d <- rbs[rbs$who %in% ids, ]
d$days[is.na(d$days) & d$did_use == "No"] <- 0L
wide <- reshape(as.data.frame(d[c("who", "what", "days")]),
                idvar = "who", timevar = "what", direction = "wide")

fit_rbs <- function(items, K) {
  selected <- paste0("days.", items)
  stopifnot(all(selected %in% names(wide)))
  sub <- wide[complete.cases(wide[, selected, drop = FALSE]), , drop = FALSE]
  days <- as.matrix(sub[, selected, drop = FALSE])
  stopifnot(all(days %in% c(0L, 4L, 14L, 30L)))
  X <- matrix(match(days, c(0L, 4L, 14L, 30L)), nrow = nrow(days),
              dimnames = dimnames(days))
  freq <- t(apply(X, 2, tabulate, nbins = 4))
  print(freq)
  if (any(freq == 0L)) warning("At least one item has an empty response category.")
  set.seed(123)
  fit <- PMLcfaefa(X, q = 1L, K = K, c = 4L,
                   ics = matrix(1, nrow(X), 1), tol = 1e-2,
                   maxit = 200L, maxNumRestart = 5L, replix = 3L)
  list(fit = fit, who = sub$who, X = X, frequencies = freq)
}

all_rbs <- fit_rbs(c("cocaine", "heroin", "opioid", "speed", "speedball"), K = 2L)
opioid_rbs <- fit_rbs(c("heroin", "opioid", "speedball"), K = 1L)

print(all_rbs$fit)
print(opioid_rbs$fit)
