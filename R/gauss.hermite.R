# Gauss-Hermite quadrature for expectations under the standard normal
#
# Nodes and weights come from the Golub-Welsch algorithm (Golub & Welsch, 1969,
# Math. Comp. 23, 221-230): for the probabilists' Hermite polynomials, which are
# orthogonal with respect to the N(0, 1) density, the three-term recurrence
# He_{k+1}(x) = x He_k(x) - k He_{k-1}(x) gives a symmetric tridiagonal Jacobi
# matrix with zero diagonal and off-diagonal sqrt(1), ..., sqrt(n - 1). Its
# eigenvalues are the nodes, and the squared first components of the normalised
# eigenvectors are the weights (they sum to 1 because the density integrates
# to 1).
#
# Nodes and weights depend only on the order, so they are computed once per
# order and cached.

.gh_cache <- new.env(parent = emptyenv())

gh_nodes <- function(order) {
    # nodes xx and weights ww for E f(Z), Z ~ N(0, 1)
    key <- as.character(order)
    nd <- .gh_cache[[key]]
    if (is.null(nd)) {
        if (order == 1) {
            nd <- list(xx = 0, ww = 1)
        } else {
            J <- matrix(0, order, order)
            off <- sqrt(seq_len(order - 1))
            J[cbind(seq_len(order - 1), 2:order)] <- off
            J[cbind(2:order, seq_len(order - 1))] <- off
            e <- eigen(J, symmetric = TRUE)
            o <- order(e$values)
            nd <- list(xx = e$values[o], ww = e$vectors[1, o]^2)
        }
        assign(key, nd, envir = .gh_cache)
    }
    nd
}

gauss_hermite <- function(f, order = 5) {
    # E f(Z), Z ~ N(0, 1)
    nd <- gh_nodes(order)
    xx <- nd$xx
    ww <- nd$ww
    ans <- 0
    for (i in seq_along(xx)) ans <- ans + ww[i] * f(xx[i])
    return(ans)
}
