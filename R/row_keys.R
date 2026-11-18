row_keys <- function(patterns) {
    # One string key per row of a response matrix: rows get the same key
    # if and only if they are identical. Replaces
    # apply(patterns, 1, digest::digest), which serialised and hashed
    # every row; only key equality and first-appearance order are used.
    # Response data are integer category codes (or NA), so paste() is exact.
    patterns <- as.matrix(patterns)
    cols <- lapply(seq_len(ncol(patterns)), function(j) patterns[, j])
    do.call(paste, c(cols, sep = "\r"))
}
