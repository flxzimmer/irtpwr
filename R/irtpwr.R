

#' Perform Power Analysis
#'
#' Perform analytical or sampling-based power analysis for the Wald, LR, score, or gradient statistic.
#'
#' @param hyp Hypothesis Object created by the setup.hypothesis function
#' @param stat character vector containing the statistics to be calculated. Options are 'Wald','LR','Score', and 'Gradient'. By default, all statistics are included
#' @param sampling.npers integer, sample size for the sampling-based approach. An artificial data set of this size is generated to fit a model and later estimate the noncentrality parameter from.
#' @param approx.npers integer, sample size for approximating the Fisher expected information matrix in the sampling-based approach. An artificial data set is calculated of this size to calculate the Fisher expected information matrix from. In contrast to the data created with the sampling.npers sample size, this sample is not used to fit a model.
#' @param method character, indicating the method used. The options are 'analytical'(default) for the analytical power analysis method or 'sampling' for the sampling-based method. The sampling-based method is generally recommended for higher numbers of items.
#' @param SE.type Method for calculation of the observed information matrix used for calculating the statistics in the sampling-based approach ('Oakes' by default). Another option is 'Fisher'.
#' @param sampling.mat Approach to calculate the information matrix used for calculating the statistics in the sampling-based approach. By default ('ApproxFisher'), an sampling-based approximation of the expected Fisher matrix is calculated using an observed information matrix of the type SE.type
#' @param power numeric, statistical power for which the necessary sample size is calculated
#' @param N integer, sample size for which the statistical power is calculated.
#' @param alpha numeric, alpha level
#'
#' @return function returns an object of class irtpwrresult
#' @export
#'
#' @examples \donttest{
#'
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1,verbose = FALSE)
#' hyp <- setup.hypothesis(type = '1PLvs2PL', altpars = mirtfit)
#' res <- irtpwr(hyp=hyp,alpha=.05,power =.8)
#' summary(res)
#' }
irtpwr <- function(hyp, stat = c("Wald", "LR", "Score",
    "Gradient"), method = "analytical", sampling.npers = 10^5,
    approx.npers = 10^5, SE.type = "Oakes", sampling.mat = "ApproxFisher",
    power = NULL, N = NULL, alpha = NULL) {


    if (method == "sampling") {

        data <- setup.data(hyp, n = sampling.npers)
        fitted <- mml.fit(hyp, data = data, infmat.unres = sampling.mat,
            infmat.res = sampling.mat, approx.npers = approx.npers,
            SE.type = SE.type)
        obs <- stat_obs(fitted, stat)
        ncps <- sapply(obs, function(x) get_ncp(x,
            df = nrow(hyp$resmod$Amat))$ncp)
        ncps <- ncps/sampling.npers
        names(ncps) <- stat
    }

    if (method == "analytical") {

        a <- rep(NA, 4)  # Prepare Result Vector

        lx <- NULL  # Prepare empty lx Vector for Gradient stat

        # Calculate ML Restricted parameters
        # centrally for use in Score and LR and
        # Gradient
        if ("LR" %in% stat | "Score" %in% stat | "Gradient" %in%
            stat) {
            parsr <- hyp$type$maximizeL(hyp)
        }

        if ("Wald" %in% stat) {
            a[1] <- ncp.wald(hyp = hyp)
        }
        if ("LR" %in% stat) {
            a[2] <- ncp.lr(hyp = hyp, parsr = parsr)
        }
        if ("Score" %in% stat) {
            sctemp <- ncp.score(hyp = hyp, parsr = parsr)
            a[3] <- sctemp$ncp
            lx <- sctemp$lx
        }
        if ("Gradient" %in% stat) {
            a[4] <- ncp.grad(hyp = hyp, parsr = parsr,
                lx = lx)
        }
        ncps <- a[!is.na(a)]
        names(ncps) <- stat
    }

    re <- list(call = match.call(), hyp = hyp, ncps = ncps,
        power = power, N = N, alpha = alpha, method = method)
    class(re) <- "irtpwrresult"
    return(re)
}


ncp.wald <- function(hyp, n = 1) {

    # Analytical noncentrality parameter for the
    # Wald statistic

    resmod <- hyp$resmod
    unresmod <- hyp$unresmod

    multigroup <- isTRUE(unresmod$multigroup)

    method <- "Fisher"

    A <- resmod$Amat
    dif <- A %*% unresmod$longpars - resmod$cvec
    sigma <- infmat(unresmod$parsets, method = method,
        model = hyp$unresmod$model, multigroup = multigroup,
        itemtype = hyp$unresmod$itemtype) |>
        solve()
    re <- t(dif) %*% solve(A %*% sigma %*% t(A)) %*%
        dif |>
        c()
    return(re * n)
}



ncp.lr <- function(hyp, n = 1, parsr = NULL) {
    # Analytical noncentrality parameter for the
    # LR statistic

    resmod <- hyp$resmod
    unresmod <- hyp$unresmod

    if (is.null(parsr)) {
        parsr <- hyp$type$maximizeL(hyp)
    }
    res <- c()
    unres <- c()
    funs <- load.functions(resmod$itemtype)


    if (isTRUE(unresmod$multigroup)) {

        pars1 <- unresmod$parsets[[1]]
        pars2 <- unresmod$parsets[[2]]
        n.kat <- max(ncol(pars1$d), 2)
        patterns <- as.matrix(expand.grid(lapply(1:resmod$n.items,
            function(x) 0:(n.kat - 1))))

        for (i in seq_len(nrow(patterns))) {
            gr <- funs$g(patterns[i, ], parsr)
            g1 <- funs$g(patterns[i, ], pars1)
            g2 <- funs$g(patterns[i, ], pars2)
            res[i] <- log(gr) * (0.5 * g1 + 0.5 * g2)
            unres[i] <- 0.5 * log(g1) * g1 + 0.5 *
                log(g2) * g2
        }

    } else {
        # unresmod$multigroup NA or FALSE

        pars <- unresmod$parsets
        n.kat <- max(ncol(pars$d), 2)
        patterns <- as.matrix(expand.grid(lapply(1:resmod$n.items,
            function(x) 0:(n.kat - 1))))

        for (i in seq_len(nrow(patterns))) {
            gr <- funs$g(patterns[i, ], parsr)
            gu <- funs$g(patterns[i, ], pars)
            res[i] <- log(gr) * (gu)
            unres[i] <- log(gu) * (gu)
        }

    }

    res <- sum(res)
    unres <- sum(unres)

    re <- 2 * (unres - res) |>
        c()

    return(re * n)

}



ncp.score <- function(hyp, n = 1, parsr = NULL) {
    # Analytical noncentrality parameter for the
    # Score statistic


    resmod <- hyp$resmod
    unresmod <- hyp$unresmod

    if (is.null(parsr)) {
        parsr <- hyp$type$maximizeL(hyp)
    }

    funs <- load.functions(resmod$itemtype)
    multigroup <- isTRUE(unresmod$multigroup)

    if (multigroup) {
        # Multigroup Model

        relpars <- resmod$relpars

        pars1 <- unresmod$parsets[[1]]
        pars2 <- unresmod$parsets[[2]]
        n.kat <- max(ncol(pars1$d), 2)
        patterns <- lapply(1:resmod$n.items, function(x) 0:(n.kat -
            1)) |>
            expand.grid() |>
            as.matrix() |>
            (function(x) split(x, seq(nrow(x))))()

        freq12 <- lapply(patterns, function(x) c(funs$g(x,
            pars1)/2, funs$g(x, pars2)/2)) |>
            (function(x) do.call(rbind, x))()
        freq <- cbind(rowSums(freq12), freq12)
        ly <- lapply(seq_len(length(patterns)), function(i) {
            l <- funs$ldot(patterns[[i]], parsr)
            list(freq[i, 1] * l, freq[i, 2] * l, freq[i,
                3] * l)
        })
        lx <- lapply(ly, function(x) x[[1]]) |>
            (function(x) do.call(rbind, x))() |>
            colSums() |>
            (function(x) array(x, dim = c(length(x),
                1)))()
        lx1 <- lapply(ly, function(x) x[[2]]) |>
            (function(x) do.call(rbind, x))() |>
            colSums() |>
            (function(x) array(x, dim = c(length(x),
                1)))()
        lx2 <- lapply(ly, function(x) x[[3]]) |>
            (function(x) do.call(rbind, x))() |>
            colSums() |>
            (function(x) array(x, dim = c(length(x),
                1)))()

        lx[relpars, ] <- lx1[relpars, ]
        lx <- c(lx, lx2[relpars, ]) |>
            (function(x) array(x, dim = c(length(x),
                1)))()

        # transform restricted pars for infmat
        # calculation
        parsr <- list(parsr, parsr)

    } else {
        # Single Group Model

        pars <- unresmod$parsets
        n.kat <- max(ncol(pars$d), 2)
        patterns <- lapply(1:resmod$n.items, function(x) 0:(n.kat -
            1)) |>
            expand.grid() |>
            as.matrix() |>
            (function(x) split(x, seq(nrow(x))))()
        i <- 0
        ly <- lapply(patterns, function(x) {
            funs$ldot(x, parsr) * funs$g(x, pars)
        })
        lx <- do.call(rbind, ly) |>
            colSums() |>
            (function(x) array(x, dim = c(length(x),
                1)))()
    }

    sigma <- infmat(parsr, method = "Fisher", model = hyp$unresmod$model,
        multigroup = multigroup, itemtype = hyp$unresmod$itemtype) |>
        solve()

    ncp <- t(lx) %*% sigma %*% lx |>
        c()

    re <- list(ncp = ncp * n, lx = lx)

    return(re)
}


ncp.grad <- function(hyp, n = 1, parsr = NULL, lx = NULL) {
    # Analytical noncentrality parameter for the
    # Gradient statistic

    resmod <- hyp$resmod
    unresmod <- hyp$unresmod

    if (is.null(parsr)) {
        parsr <- hyp$type$maximizeL(hyp)
    }

    funs <- load.functions(resmod$itemtype)
    multigroup <- isTRUE(unresmod$multigroup)

    if (multigroup & is.null(lx)) {
        # Multigroup Model

        relpars <- resmod$relpars

        pars1 <- unresmod$parsets[[1]]
        pars2 <- unresmod$parsets[[2]]
        n.kat <- max(ncol(pars1$d), 2)
        patterns <- lapply(1:resmod$n.items, function(x) 0:(n.kat -
            1)) |>
            expand.grid() |>
            as.matrix() |>
            (function(x) split(x, seq(nrow(x))))()

        freq12 <- lapply(patterns, function(x) c(funs$g(x,
            pars1)/2, funs$g(x, pars2)/2)) |>
            (function(x) do.call(rbind, x))()
        freq <- cbind(rowSums(freq12), freq12)
        ly <- lapply(seq_len(length(patterns)), function(i) {

            l <- funs$ldot(patterns[[i]], parsr)
            list(freq[i, 1] * l, freq[i, 2] * l, freq[i,
                3] * l)
        })
        lx <- lapply(ly, function(x) x[[1]]) |>
            (function(x) do.call(rbind, x))() |>
            colSums() |>
            (function(x) array(x, dim = c(length(x),
                1)))()
        lx1 <- lapply(ly, function(x) x[[2]]) |>
            (function(x) do.call(rbind, x))() |>
            colSums() |>
            (function(x) array(x, dim = c(length(x),
                1)))()
        lx2 <- lapply(ly, function(x) x[[3]]) |>
            (function(x) do.call(rbind, x))() |>
            colSums() |>
            (function(x) array(x, dim = c(length(x),
                1)))()

        lx[relpars, ] <- lx1[relpars, ]
        lx <- c(lx, lx2[relpars, ]) |>
            (function(x) array(x, dim = c(length(x),
                1)))()

        # transform restricted pars for infmat
        # calculation
        parsr <- list(parsr, parsr)

    }

    if (!multigroup & is.null(lx)) {
        # Single Group Model

        pars <- unresmod$parsets
        n.kat <- max(ncol(pars$d), 2)
        patterns <- lapply(1:resmod$n.items, function(x) 0:(n.kat -
            1)) |>
            expand.grid() |>
            as.matrix() |>
            (function(x) split(x, seq(nrow(x))))()
        i <- 0
        ly <- lapply(patterns, function(x) {
            funs$ldot(x, parsr) * funs$g(x, pars)
        })
        lx <- do.call(rbind, ly) |>
            colSums() |>
            (function(x) array(x, dim = c(length(x),
                1)))()
    }

    A <- resmod$Amat
    dif <- A %*% unresmod$longpars - resmod$cvec

    lambda <- findlambda(lx, A)

    re <- t(lambda) %*% dif |>
        as.numeric() |>
        abs()

    return(re * n)
}



