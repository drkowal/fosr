library(refund) # fosr.vs function
library(grpreg) # for cv.grpreg functions used by refund::fosr.vs
library(splines) # for bs

modified_fosr.vs_for_lasso = function(
    formula,
    data,
    lambda,
    nbasis = 10,
    method = "grLasso",
    epsilon = 0.00001,
    max.iter_num = 100)
{
    stopifnot(method == "grLasso")
    cl = match.call()
    mf = match.call(expand.dots = FALSE)
    m = match(c("formula", "data"), names(mf), 0)
    mf = mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] = as.name("model.frame")
    mf = eval.parent(mf)
    Y = model.response(mf)
    W.des = model.matrix(formula, data = mf)
    N = dim(Y)[1]
    D = dim(Y)[2]
    K = dim(W.des)[2]
    Theta = bs(1:D, df = nbasis, intercept = TRUE, degree = 3)
    Z = kronecker(W.des, Theta)
    Y.vec = as.vector(t(Y))
    ls <- lm(Y.vec ~ Z + 0)
    a <- matrix(ls$residuals, nrow = N, ncol = D, byrow = T)
    w <- fpca.sc(a, var = T)
    b <- Reduce("+", lapply(1:length(w$evalues), function(x) {
        w$evalues[x] * w$efunctions[, x] %*% t(w$efunctions[,
            x])
    })) + w$sigma2 * diag(1, D, D)
    f <- coef(ls)
    d <- rep(0, nbasis * K)
    num.iter = 0
    cat("Beginning iterative algorithm \n")
    while (sum((f - d)^2) > epsilon & num.iter <= max.iter_num) {
        num.iter = num.iter + 1
        d <- f
        c <- chol(solve(b))
        Y_new = Y %*% t(c)
        Theta_new = c %*% Theta
        Z_new = kronecker(W.des, Theta_new)
        Y.vecnew = as.vector(t(Y_new))
        if (method == "ls") {
            ls <- lm(Y.vecnew ~ Z_new)
            f <- as.vector(coef(ls))[-1]
        }
        else {
            cvlam <- grpreg(Z_new, Y.vecnew, group = switch(colnames(W.des)[1],
                `(Intercept)` = c(rep(1:K, each = nbasis)) -
                  1, c(rep(1:K, each = nbasis))), penalty = method,
                lambda = lambda)
            f <- as.vector(coef(cvlam))[-1]
        }
        a <- matrix(Y.vec - Z %*% f, nrow = N, ncol = D, byrow = T)
        w <- fpca.sc(a, var = T)
        b <- Reduce("+", lapply(1:length(w$evalues), function(x) {
            w$evalues[x] * w$efunctions[, x] %*% t(w$efunctions[,
                x])
        })) + w$sigma2 * diag(1, D, D)
        if (num.iter%%10 == 1)
            cat(".")
    }
    B <- matrix(f, nrow = K, byrow = T)
    est <- B %*% t(Theta)
    rownames(est) <- colnames(W.des)
    fitted <- W.des %*% est
    res <- Y - fitted
    ret <- list(call = cl, formula = formula, coefficients = est,
        fitted.values = fitted, residuals = res, vcov = b, method = method)
    class(ret) <- "fosr.vs"
    return(ret)
}

