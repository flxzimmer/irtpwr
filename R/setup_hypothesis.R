


#' Setup null and alternative hypothesis
#'
#' @param type preset, e.g. "1PLvs2PL". Either a string for an existing preset or a list with the required functions for a custom  preset. See the section extending the package for a tutorial on how to create a custom preset.
#' @param altpars List of model parameters following the alternative hypothesis. The format depends on the preset. When parameters are derived from observed data by mirt, the coef_short function can be used to convert the parameters to the right format.
#' @param nullpars List of model parameters following the null hypothesis. The format depends on the preset. Null parameters are not necessary for some hypothesis presets, e.g. 1PLvs2PL.
#'
#' @return a list specifying the hypothesis for usage in further functions, e.g. estimation of the noncentrality parameters.
#' @export
#'
#' @examples
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1)
#' pars = coef_short(mirtfit)
#' hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = pars)
#'
setup_hypothesis = function(type, altpars = NULL,nullpars = NULL) {
# setup a null and alternative hypothesis according to a given type

  # altpars is in a coef_short format

  if (!is.list(type)) {

    if (type == "1PLvs2PL") {
      type = h_1PLvs2PL
    }
    else if (type == "DIF2PL") {
      type = h_DIF2PL
    }
    else {print("error..")}
  }


  resmod = type$res(altpars = altpars,nullpars = nullpars)
  unresmod = type$unres(altpars = altpars)

  re = list(resmod=resmod,unresmod = unresmod,type=type)
  return(re)
}

