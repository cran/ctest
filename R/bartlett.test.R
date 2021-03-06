bartlett.test <- function(x, g) {
    LM <- FALSE
    if (is.list(x)) {
        if (length(x) < 2)
            stop("x must be a list with at least 2 elements")
        DNAME <- deparse(substitute(x))
        if (all(sapply(x, function(obj) inherits(obj, "lm"))))
            LM <- TRUE
        else
            x <- lapply(x, function(x) x <- x[is.finite(x)])
        k <- length(x)
    }
    else {
        if (length(x) != length(g))
            stop("x and g must have the same length")
        DNAME <- paste(deparse(substitute(x)), "and",
                       deparse(substitute(g)))
        OK <- complete.cases(x, g)
        x <- x[OK]
        g <- as.factor(g[OK])
        k <- nlevels(g)
        if (k < 2)
            stop("all observations are in the same group")
        x <- split(x, g)
    }

    if (LM) {
        n <- sapply(x, function(obj) obj$df.resid)
        v <- sapply(x, function(obj) sum(obj$residuals^2))
    } else {
        n <- sapply(x, "length") - 1
        if (any(n <= 0))
            stop("there must be at least 2 observations in each group")
        v <- sapply(x, "var")
    }

    n.total <- sum(n)
    v.total <- sum(n * v) / n.total
    STATISTIC <- ((n.total * log(v.total) - sum(n * log(v))) /
                  (1 + (sum(1 / n) - 1 / n.total) / (3 * (k - 1))))
    names(STATISTIC) <- "Bartlett's K-square"
    PARAMETER <- k - 1
    names(PARAMETER) <- "df"
  
    RVAL <- list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = 1 - pchisq(STATISTIC, PARAMETER),
                 data.name = DNAME,
                 method = "Bartlett test for homogeneity of variances")
    class(RVAL) <- "htest"
    return(RVAL)
}
