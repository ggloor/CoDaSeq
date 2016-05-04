biplot2Color <- function (x, choices = 1L:2L, scale = 1, pc.biplot = FALSE, ...)
{
    if (length(choices) != 2L)
        stop("length of choices must be 2")
    if (!length(scores <- x$x))
        stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
            domain = NA)
    if (is.complex(scores))
        stop("biplots are not defined for complex PCA")
    lam <- x$sdev[choices]
    n <- NROW(scores)
    lam <- lam * sqrt(n)
    if (scale < 0 || scale > 1)
        warning("'scale' is outside [0, 1]")
    if (scale != 0)
        lam <- lam^scale
    else lam <- 1
    if (pc.biplot)
        lam <- lam/sqrt(n)
    biplot.default2Color(t(t(scores[, choices])/lam), t(t(x$rotation[,
        choices]) * lam), ...)
    invisible()
}
