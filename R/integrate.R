# Vectorised pattern-space integration (internal)
#
# g_all(model, patterns, pars)[p] is bit-identical to
# load.functions(model)$g(patterns[p, ], pars), and row p of
# ldot_all(model, patterns, pars) to load.functions(model)$ldot(patterns[p, ], pars).
# What keeps them identical:
# - f(th, a_j, d_j, x) (and its derivatives) is elementwise, so its value for
#   pattern p and item j depends only on (th, a_j, d_j, x_pj). It is evaluated
#   once per node, item and response category, with exactly the same
#   floating-point operations as in the per-pattern call, and looked up for
#   every pattern;
# - the product over items of each pattern is accumulated like prod() (long
#   double accumulator, items in the same order; src/rowprods.c), as is
#   prod(res[-i]) in fdot();
# - the quadrature nodes are accumulated in R, in the same order.
# Models without a vectorised version (3PL, multidimensional 2PL) fall back
# to the per-pattern functions.

gh_order <- function(model) {
    if (model == "gpcm") 50 else 30
}

vec_model <- function(model, pars) {
    model %in% c("2PL", "gpcm") && !is.matrix(pars$a)
}

use_long_double <- function() {
    ld <- .gh_cache[["long.double"]]
    if (is.null(ld)) {
        ld <- isTRUE(capabilities("long.double"))
        assign("long.double", ld, envir = .gh_cache)
    }
    ld
}

# response categories of a pattern matrix as 0-based codes into `cats`
pattern_codes <- function(patterns) {
    cats <- sort(unique(as.vector(patterns)))
    idx <- match(patterns, cats) - 1L
    dim(idx) <- dim(patterns)
    list(idx = idx, cats = as.numeric(cats))
}

# f (fun = funs$f) or one of its derivatives for every response category
# x every item j and node th: array [category, item, node]
f_table <- function(model, fun, cats, pars, xx) {
    nc <- length(cats)
    k <- length(pars$a)
    nq <- length(xx)
    TH <- rep(xx, each = nc * k)
    X <- rep(cats, times = k * nq)
    A <- rep(rep(pars$a, each = nc), times = nq)
    if (model == "2PL") {
        D <- rep(rep(pars$d, each = nc), times = nq)
        v <- fun(TH, A, D, X)
    } else {
        D1 <- rep(rep(pars$d[, 2], each = nc), times = nq)
        D2 <- rep(rep(pars$d[, 3], each = nc), times = nq)
        v <- fun(TH, A, D1, D2, X)
    }
    array(as.numeric(v), dim = c(nc, k, nq))
}

g_all <- function(model, patterns, pars) {
    patterns <- as.matrix(patterns)
    if (!vec_model(model, pars)) {
        funs <- load.functions(model)
        return(vapply(seq_len(nrow(patterns)), function(i) funs$g(patterns[i,
            ], pars), 0))
    }
    pc <- pattern_codes(patterns)
    nd <- gh_nodes(gh_order(model))
    xx <- nd$xx
    ww <- nd$ww
    tab <- f_table(model, load.functions(model)$f, pc$cats, pars, xx)
    # RP[p, q] = prod(f(xx[q], pars$a, pars$d, patterns[p, ]))
    RP <- .Call(irtpwr_pattern_prods, pc$idx, tab, use_long_double())
    ans <- 0
    for (i in seq_along(xx)) ans <- ans + ww[i] * RP[, i]
    ans
}

gdot_all <- function(model, patterns, pars) {
    # row p is load.functions(model)$gdot(patterns[p, ], pars)
    pc <- pattern_codes(patterns)
    nd <- gh_nodes(gh_order(model))
    xx <- nd$xx
    ww <- nd$ww
    funs <- load.functions(model)
    # the derivative functions used by funs$fdot, in its order
    fenv <- environment(funs$fdot)
    derivs <- if (model == "2PL")
        list(fenv$f1_2pl, fenv$f2_2pl) else list(fenv$f1_gpcm, fenv$f2_gpcm,
        fenv$f3_gpcm)
    tab <- f_table(model, funs$f, pc$cats, pars, xx)
    dtabs <- lapply(derivs, function(fd) f_table(model, fd, pc$cats,
        pars, xx))
    # G[p, , q] = fdot(xx[q], ..., patterns[p, ]), i.e.
    # prod(res[-i]) * c(temp1[i], temp2[i], ...) for each item i
    G <- .Call(irtpwr_pattern_fdot, pc$idx, tab, dtabs, use_long_double())
    ans <- 0
    dG <- dim(G)[1:2]
    for (i in seq_along(xx)) ans <- ans + ww[i] * array(G[, , i, drop = FALSE],
        dG)
    ans
}

ldot_all <- function(model, patterns, pars) {
    # row p is load.functions(model)$ldot(patterns[p, ], pars)
    patterns <- as.matrix(patterns)
    if (!vec_model(model, pars)) {
        funs <- load.functions(model)
        return(do.call(rbind, lapply(seq_len(nrow(patterns)), function(i) funs$ldot(patterns[i,
            ], pars))))
    }
    gdot_all(model, patterns, pars)/g_all(model, patterns, pars)
}

#' @useDynLib irtpwr, .registration = TRUE
NULL
