ks.test <- function(x, y, ..., alternative = "two.sided")
{
    CHOICES <- c("two.sided", "less", "greater")
    alternative <- CHOICES[pmatch(alternative, CHOICES)]
    if (length(alternative) > 1 || is.na(alternative)) 
        stop("alternative must be \"two.sided\", \"less\" or \"greater\"")

    DNAME <- deparse(substitute(x))      
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 1)
        stop("Not enough x data")

    if (is.numeric(y)) {
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        y <- y[!is.na(y)]
        n.x <- n
        n.y <- length(y)
        if (n.y < 1)
            stop("Not enough y data")
        METHOD <- "Two-sample Kolmogorov-Smirnov test"
        n <- n.x * n.y / (n.x + n.y)
        w <- c(x, y)
        z <- cumsum(ifelse(order(w) <= n.x, 1 / n.x, - 1 / n.y))
        if (length(unique(w)) < (n.x + n.y)) {
            warning("cannot compute correct p-values with ties")
            z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
            print(z)
        }
        STATISTIC <- switch(alternative,
                            "two.sided" = max(abs(z)),
                            "greater" = max(z),
                            "less" = - min(z))
    }
    else {
        if (is.character(y))
            y <- get(y, mode="function")
        if (mode(y) != "function")
            stop("y must be numeric or a string naming a valid function")
        METHOD <- "One-sample Kolmogorov-Smirnov test"
        n <- length(x)
        x <- y(sort(x), ...) - (0 : (n-1)) / n
        STATISTIC <- switch(alternative,
                            "two.sided" = max(abs(c(x, x-1/n))),
                            "greater" = max(c(x, x-1/n)),
                            "less" = - min(c(x, x-1/n)))
    }

    names(STATISTIC) <- switch(alternative,
                               "two.sided" = "D",
                               "greater" = "D^+",
                               "less" = "D^-")
    PVAL <- ifelse(alternative == "two.sided",
                   1 - pks(sqrt(n) * STATISTIC),
                   exp(- 2 * n * STATISTIC^2))
    
    RVAL <- list(statistic = STATISTIC,
                 p.value = PVAL,
                 alternative = alternative,
                 method = METHOD,
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}
