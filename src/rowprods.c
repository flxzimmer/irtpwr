/*
 * Pattern products for the pattern-space integration (see R/integrate.R).
 *
 * idx is the n x k integer matrix of 0-based response category codes,
 * tab the [category, item, node] array of item response functions f.
 * Each product over items is accumulated exactly like base R's prod()
 * (rprod() in R's src/main/summary.c): an LDOUBLE accumulator starting at
 * 1, multiplied by the elements in item order, clamped to +-Inf beyond
 * DBL_MAX and then rounded to double. That keeps the results bit-identical
 * to prod(f(...)) and prod(res[-i]) in load.functions() on the same platform.
 * LDOUBLE is long double unless R was built without long double support;
 * that is not visible in the public headers, so the R side passes
 * capabilities("long.double") as `ld`.
 * Only multiplications happen here (no additions), so there is nothing a
 * compiler could contract into fused multiply-adds; the quadrature sums
 * stay in R.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <float.h>

static double ld_to_double(long double s)
{
    if (s > DBL_MAX) return R_PosInf;
    else if (s < -DBL_MAX) return R_NegInf;
    else return (double) s;
}

static double d_to_double(double s)
{
    if (s > DBL_MAX) return R_PosInf;
    else if (s < -DBL_MAX) return R_NegInf;
    else return s;
}

/* product over items j != skip of tq[idx[r, j] + nc * j] */
#define PAT_PROD(TYPE, CONV, res, skip) do {                          \
        TYPE s = 1.0;                                                 \
        for (int j = 0; j < k; j++)                                   \
            if (j != (skip))                                          \
                s *= tq[ix[r + (R_xlen_t) n * j] + (R_xlen_t) nc * j]; \
        (res) = CONV(s);                                              \
    } while (0)

static void check_args(SEXP idx, SEXP tab, int *n, int *k, int *nc, int *nq)
{
    if (!Rf_isInteger(idx) || !Rf_isMatrix(idx))
        Rf_error("idx must be an integer matrix");
    SEXP dim = Rf_getAttrib(tab, R_DimSymbol);
    if (!Rf_isReal(tab) || LENGTH(dim) != 3)
        Rf_error("tab must be a 3-dimensional double array");
    *n = Rf_nrows(idx);
    *k = Rf_ncols(idx);
    *nc = INTEGER(dim)[0];
    *nq = INTEGER(dim)[2];
    if (INTEGER(dim)[1] != *k) Rf_error("tab does not match idx");
    const int *ix = INTEGER(idx);
    for (R_xlen_t e = 0; e < XLENGTH(idx); e++)
        if (ix[e] == NA_INTEGER || ix[e] < 0 || ix[e] >= *nc)
            Rf_error("invalid category code");
}

/* RP[r, q] = prod over items j of tab[idx[r, j], j, q] */
SEXP irtpwr_pattern_prods(SEXP idx, SEXP tab, SEXP ld)
{
    int n, k, nc, nq;
    check_args(idx, tab, &n, &k, &nc, &nq);
    int use_ld = Rf_asLogical(ld);
    const int *ix = INTEGER(idx);
    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, n, nq));
    double *out = REAL(ans);
    for (int q = 0; q < nq; q++) {
        const double *tq = REAL(tab) + (R_xlen_t) nc * k * q;
        double *oq = out + (R_xlen_t) n * q;
        for (int r = 0; r < n; r++) {
            if (use_ld) PAT_PROD(long double, ld_to_double, oq[r], -1);
            else PAT_PROD(double, d_to_double, oq[r], -1);
        }
    }
    UNPROTECT(1);
    return ans;
}

/*
 * fdot() for every pattern and node:
 * G[r, i * m + l, q] = prod_{j != i} tab[idx[r, j], j, q] * dtabs[[l]][idx[r, i], i, q]
 * (the order of fdot()'s unlist(lapply(items, prod(res[-i]) * c(temp1[i], ...))))
 */
SEXP irtpwr_pattern_fdot(SEXP idx, SEXP tab, SEXP dtabs, SEXP ld)
{
    int n, k, nc, nq;
    check_args(idx, tab, &n, &k, &nc, &nq);
    if (!Rf_isNewList(dtabs)) Rf_error("dtabs must be a list");
    int m = LENGTH(dtabs);
    for (int l = 0; l < m; l++) {
        SEXP dt = VECTOR_ELT(dtabs, l);
        if (!Rf_isReal(dt) || XLENGTH(dt) != XLENGTH(tab))
            Rf_error("dtabs do not match tab");
    }
    int use_ld = Rf_asLogical(ld);
    const int *ix = INTEGER(idx);
    SEXP dim = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(dim)[0] = n;
    INTEGER(dim)[1] = m * k;
    INTEGER(dim)[2] = nq;
    SEXP ans = PROTECT(Rf_allocArray(REALSXP, dim));
    double *out = REAL(ans);
    R_xlen_t ncol = (R_xlen_t) m * k;
    for (int q = 0; q < nq; q++) {
        const double *tq = REAL(tab) + (R_xlen_t) nc * k * q;
        double *oq = out + (R_xlen_t) n * ncol * q;
        for (int r = 0; r < n; r++) {
            for (int i = 0; i < k; i++) {
                double e;
                if (use_ld) PAT_PROD(long double, ld_to_double, e, i);
                else PAT_PROD(double, d_to_double, e, i);
                R_xlen_t cell = ix[r + (R_xlen_t) n * i] + (R_xlen_t) nc * i +
                    (R_xlen_t) nc * k * q;
                for (int l = 0; l < m; l++)
                    oq[r + (R_xlen_t) n * ((R_xlen_t) i * m + l)] =
                        e * REAL(VECTOR_ELT(dtabs, l))[cell];
            }
        }
    }
    UNPROTECT(2);
    return ans;
}

static const R_CallMethodDef CallEntries[] = {
    {"irtpwr_pattern_prods", (DL_FUNC) &irtpwr_pattern_prods, 3},
    {"irtpwr_pattern_fdot", (DL_FUNC) &irtpwr_pattern_fdot, 4},
    {NULL, NULL, 0}
};

void R_init_irtpwr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
