library(RTMB)


compiler::enableJIT(0)

autodiff = function(fn, arg0=rep(0, 4)) {
    fn = rlang::as_function(fn)
    fn_ad = MakeTape(fn, arg0)
    fn_ad$simplify()
    fn_ad
}

regr = function(S, v1, v2) {
    solve(S[v1, v1], S[v1, v2])
}

regr_se = function(S, v1, v2, n = 1) {
    num = S[v2, v2] - S[v2, v1] %*% solve(S[v1, v1], S[v1, v2])
    p = length(v1)
    denom = S[v1[1], v1[1]] - S[v1[1], v1[-1]] %*% solve(S[v1[-1], v1[-1]], S[v1[-1], v1[1]])
    sqrt(c(num / denom) / (n - p))
}

sigma_ext = function(S0, pars = rep(0, 3)) {
  `[<-` <- RTMB::ADoverload("[<-")
    c <- RTMB::ADoverload("c")
    L = diag(4)
    L[2:4, 2:4] = t(chol(S0))

    L[1, 2:4] = pars * c(1, cumprod(sqrt(1 - pars[1:2]^2)))

    L[1, 1] = sqrt(1 - sum(L[1, 2:4]^2))
    rownames(L) = c("u", rownames(S0))

    # Su = tcrossprod(L)
    Su = L %*% t(L)

    coef_zhat = regr(Su, 2, 3)
    zhat = coef_zhat * L[2, ]

    L = rbind(L, zhat = c(zhat))
    tcrossprod(L)
}

ests = function(S) {
    c(
        soo = unname(regr(S, c(3,2), 4)[1]),
        iv = unname(regr(S, 5, 4)[1])
    )
}

ses = function(S, n) {
    c(
        soo = unname(regr_se(S, 3:2, 4, n)),
        iv = unname(regr_se(S, 5, 4, n))
    )
}

tau <- function(S) {
    tau = regr(S, 3:1, 4)[1]
}

biases = function(S) {
    tau = regr(S, 3:1, 4)[1]
    ests(S) - tau
}

part_r2 = function(S, idxs) {
    P = solve(S[idxs, idxs])
    -P[1, 2] /  sqrt(P[1, 1] * P[2, 2])
}

arrows = function(S) {
    c(
        soo_z = part_r2(S, c(1, 3, 2)),
        soo_y = part_r2(S, c(1, 4, 2:3)),
        iv_z  = part_r2(S, c(1, 2)),
        iv_y  = part_r2(S, c(4, 2, 1, 3))
    )
}

arrows_idx = list(
    soo = c("soo_z", "soo_y"),
    iv  = c("iv_z", "iv_y")
)


calc_rv <- function(
    S0,
    b,
    est = ests(sigma_ext(S0)),
    lb = rep(-1, 3),
    ub = rep(1, 3),
    method = c("soo", "iv"),
    tol = 1e-10
) {
    method = match.arg(method)

    stopifnot(
        nrow(S0) == 3,
        ncol(S0) == 3,
        length(lb) == 3,
        length(ub) == 3
    )

    max_rho = 1 - 1e-5
    ub = pmax(pmin(ub, max_rho), -max_rho)
    lb = pmax(pmin(lb, max_rho), -max_rho)
    stopifnot(all(lb <= ub))

    learn = which(lb < ub)
    midpt = (lb + ub) / 2

    ctou = function(rho) atanh(rho[learn])


    utoc = function(x) {
      `[<-` <- RTMB::ADoverload("[<-")
        out = midpt
        out[learn] = tanh(x)
        out
    }

    # Residual variance of Y given WZ and Z
    resvar_y =
        S0[3, 3] -
        S0[3, 1:2] %*%
        solve(S0[1:2, 1:2], S0[1:2, 3])

    # Residual variance of Z given WZ
    resvar_z =
        S0[2, 2] -
        S0[2, 1]^2 / S0[1, 1]

    bias_fac = sqrt(c(resvar_y / resvar_z))

    # Find an initial value having zero bias for the selected estimator
    b0 = unname(est[method] - est["soo"])

    rho0 = uniroot(
        \(x) b0^2 * (1 - x^2) - x^4 * bias_fac^2,
        c(0, max_rho),
        tol = 1e-16
    )$root

    # rho = (U-WZ, U-Z|WZ, U-Y|WZ,Z)
    init1 = c(0, -rho0, rho0)
    init2 = c(0,  rho0, rho0)

    # True coefficient on Z in Y ~ Z + WZ + U
    tau0_1 = regr(sigma_ext(S0, init1), 3:1, 4)[1]
    tau0_2 = regr(sigma_ext(S0, init2), 3:1, 4)[1]

    if (abs(tau0_1 - est[method]) <
        abs(tau0_2 - est[method])) {
        init = init1
        tau0 = tau0_1
    } else {
        init = init2
        tau0 = tau0_2
    }

    lb = pmin(lb, init)
    ub = pmax(ub, init)

    arr_i = arrows_idx[[method]]

    run_optim = function(x0, pen) {
        fn = autodiff(
            function(x) {
                rho = utoc(x)
                Sr = sigma_ext(S0, rho)

                # True effect controlling for WZ and U
                tau = regr(Sr, 3:1, 4)[1]

                max(arrows(Sr)[arr_i]^2) +
                    pen * ((tau0 - tau) - b)^2 +
                    1e-2 * sum(rho^2)
            },
            x0
        )

        x_sc = abs(fn$jacobian(x0))

        while (any(x_sc < 1e-8)) {
            x0[x_sc < 1e-8] =
                x0[x_sc < 1e-8] + 1e-4
            x_sc = abs(fn$jacobian(x0))
        }

        optim(
            x0,
            fn,
            fn$jacobian,
            method = "L-BFGS-B",
            lower = ctou(lb),
            upper = ctou(ub),
            control = list(
                parscale = x_sc,
                factr = 1e4,
                pgtol = 1e-9,
                lmm = 15,
                maxit = 100000
            )
        )
    }

    x = ctou(init)
    pen = 1e2

    while (pen <= 1e15) {
        res = run_optim(x, pen)
        rho = utoc(res$par)
        Sr = sigma_ext(S0, rho)
        gap = abs(biases(Sr)[method] - b)

        if (gap < tol &&
            max(abs(res$par - x)) < tol) {
            break
        }

        pen = pen * 10
        x = res$par
    }

    arr = arrows(Sr)
    rv = max(arr[arr_i]^2)

    if (res$convergence != 0 || gap > 10 * tol) {
        rv = NA_real_
    }

    list(
        rho = zapsmall(rho),
        rv = rv,
        bias = biases(Sr)[method],
        arrows = zapsmall(arr),
        n_iter = log10(pen / 1e2) + 1
    )
}

calc_rvs <- function(S0, b) {
    methods = c("soo", "iv")
    ub = c(1, 0.99, 1)
    est = ests(sigma_ext(S0))

    if (is.null(names(b))) {
        if (length(b) != length(methods)) {
            stop(
                "Unnamed `b` must have length ",
                length(methods),
                "."
            )
        }

        names(b) = methods
    }

    missing_b = setdiff(methods, names(b))

    if (length(missing_b) > 0L) {
        stop(
            "`b` is missing: ",
            paste(missing_b, collapse = ", ")
        )
    }

    purrr::map_dfr(methods, function(method) {
        target_b = unname(b[method])

        withCallingHandlers(
            {
                res = calc_rv(
                    S0 = S0,
                    b = target_b,
                    est = est,
                    lb = -ub,
                    ub = ub,
                    method = method,
                    tol = 1e-4
                )

                tibble::tibble(
                    method = method,
                    target_b = target_b,
                    achieved_b = unname(res$bias),
                    gap = res$gap,
                    rv = res$rv,
                    converged = res$converged,
                    optim_code = res$optim_code,
                    optim_message = res$optim_message,
                    rho_1 = res$rho[1],
                    rho_2 = res$rho[2],
                    rho_3 = res$rho[3],
                    error = NA_character_
                )
            },
            error = function(e) {
                tibble::tibble(
                    method = method,
                    target_b = target_b,
                    achieved_b = NA_real_,
                    gap = NA_real_,
                    rv = NA_real_,
                    converged = FALSE,
                    optim_code = NA_integer_,
                    optim_message = NA_character_,
                    rho_1 = NA_real_,
                    rho_2 = NA_real_,
                    rho_3 = NA_real_,
                    error = conditionMessage(e)
                )
            }
        )
    })
}

