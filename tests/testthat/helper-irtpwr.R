# Access to non-exported functions
ns <- function(x) get(x, envir = asNamespace("irtpwr"))

# Small 1PLvs2PL hypothesis (3 items, fast analytical path)
hyp_1pl3 <- function() {
  setup.hypothesis("1PLvs2PL", list(a = c(0.8, 1.3, 1.6), d = c(-0.5, 0.2, 0.9)))
}
