


#' Setup null and alternative hypothesis
#'
#' @param type hypothesis type. Either a string for one of the included presets or a list with the required functions for a custom hypothesis. The presets are '1PLvs2PL' (Rasch model against 2PL model), 'DIF2PL' (DIF in a 2PL model, two groups), 'PCMvsGPCM' (PCM against GPCM, items with 3 response categories), 'basic' (2PL model; tests whether the intercept d of the first item is 0), '3PL_basic' (3PL model; tests whether a, d and g of the first item are 1, 0 and 0.2; sampling-based method only) and 'multi_basic' (two-dimensional 2PL model; tests whether the first two items have the same intercept d; sampling-based method only). For a tutorial on custom hypotheses, see the vignette "adding_hypotheses" (\code{vignette("adding_hypotheses", package = "irtpwr")}); the vignette "hypothesis_templates" contains templates based on the presets.
#' @param altpars List of model parameters following the alternative hypothesis. The format depends on the preset. When parameters are derived from observed data, a model fitted with mirt can be passed directly; it is converted to the right format.
#' @param nullpars List of model parameters following the null hypothesis. The format depends on the preset. Null parameters are not necessary for some hypothesis presets, e.g. 1PLvs2PL.
#'
#' @return a list specifying the hypothesis for usage in further functions, e.g. estimation of the noncentrality parameters.
#' @export
#'
#' @examples \donttest{
#'
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1,verbose = FALSE)
#' hyp <- setup.hypothesis(type = '1PLvs2PL', altpars = mirtfit)
#' }
setup.hypothesis <- function(type, altpars = NULL,
    nullpars = NULL) {
    # setup a null and alternative hypothesis
    # according to a given type

    if (!is.null(attr(class(altpars), "package"))) {
        altpars <- coef_short(altpars)
    }
    # altpars is in a coef_short format

    if (!is.list(type)) {

        if (type == "1PLvs2PL") {
            type <- h_1PLvs2PL
        } else if (type == "DIF2PL") {
            type <- h_DIF2PL
        } else if (type == "basic") {
            type <- h_2PL_basic
        } else if (type == "PCMvsGPCM") {
            type <- h_PCMvsGPCM
        } else if (type == "multi_basic") {
            type <- h_multi_basic
        } else if (type == "3PL_basic") {
            type <- h_3PL_basic
        } else {
            stop("Unknown hypothesis type '", type, "'. Available types are ",
                "'1PLvs2PL', 'DIF2PL', 'basic', 'PCMvsGPCM', 'multi_basic' ",
                "and '3PL_basic', or a list with the functions of a custom ",
                "hypothesis.", call. = FALSE)
        }
    }

    resmod <- type$res(altpars = altpars, nullpars = nullpars)
    unresmod <- type$unres(altpars = altpars)

    re <- list(resmod = resmod, unresmod = unresmod,
        type = type)
    return(re)
}

