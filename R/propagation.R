BinLigandExpr <- function(lidardata, bin = NULL, n_bins = 10L, 
    imp = TRUE) {
    .CheckObject(lidardata)
    # Include the terminal zero-weight bin in the internal bin count.
    n_bins <- .Integer(n_bins, "n_bins", min = 1L) + 1L
    expr <- .Expression(lidardata, imp)
    pos <- .AlignRows(lidardata@Image@pos, lidardata@cellname, 
        "Coordinates")
    if (!all(c("x", "y") %in% names(pos))) 
        stop("Coordinates must have x/y columns; recreate object")
    distance <- as.matrix(stats::dist(pos[, c("x", "y"), drop = FALSE]))
    if (any(!is.finite(distance))) 
        stop("Invalid distances")
    if (is.null(bin)) {
        positive <- distance[distance > 0]
        if (!length(positive)) 
            stop("No positive distances; specify bin explicitly")
        bin <- min(positive)
    }
    if (!is.numeric(bin) || length(bin) != 1L || !is.finite(bin) || 
        bin <= 0) 
        stop("bin must be positive")
    raw_bin <- floor(distance/bin) + 1
    raw_bin[raw_bin > n_bins] <- 0
    dimnames(raw_bin) <- list(lidardata@cellname, lidardata@cellname)
    valid <- character()
    for (l in rownames(expr$ligand)) {
        err <- tryCatch({
            .Nonnegative(expr$ligand[l, ], "Ligand expression")
            .Receptors(lidardata, l, expr$receptor)
            NULL
        }, error = function(e) conditionMessage(e))
        if (is.null(err)) 
            valid <- c(valid, l)
        else lidardata <- .Failure(lidardata, l, "distance", 
            err)
    }
    arr <- array(0, c(length(valid), length(lidardata@cellname), 
        n_bins), dimnames = list(Ligand = valid, Cell = lidardata@cellname, 
        Distance = as.character((seq_len(n_bins) - 1) * bin)))
    for (j in seq_len(n_bins)) {
        mask <- Matrix::Matrix(raw_bin == j, sparse = TRUE)
        arr[, , j] <- as.matrix(expr$ligand[valid, , drop = FALSE] %*% 
            mask)
    }
    lidardata@Image@dist <- distance
    lidardata@Image@dist_bin_all <- Matrix::Matrix(raw_bin, sparse = TRUE)
    attr(arr, "imp") <- imp
    lidardata@Diffusion@ligand_dist <- arr
    lidardata@Diffusion@dist_threshold <- as.numeric(n_bins)
    lidardata@Diffusion@bin <- bin
    lidardata@Diffusion@ligand_diff <- .EmptyFit()
    lidardata@Diffusion@ligand_conc <- .EmptyMatrix(lidardata@cellname)
    lidardata@Lidarscore <- .ScoresEmpty(lidardata@cellname)
    lidardata
}
.DecayWeights <- function(par, distance, B) {
    if (length(par) != 2L || any(!is.finite(par)) || any(par <= 
        0) || length(B) != 1L || !is.finite(B) || B <= 0 || any(!is.finite(distance)) || 
        any(distance < 0)) 
        stop("Invalid decay parameters/distances")
    z <- pmin(distance/B, 1)
    result <- -expm1(par[2] * log(z)) * stats::plogis(-par[2] * 
        (log(z) - log(par[1])))
    result[z == 0] <- 1
    result[z >= 1] <- 0
    if (any(!is.finite(result)) || any(result < 0) || any(result > 
        1)) 
        stop("Invalid decay curve")
    result
}
CalcRanges <- function(ki, n, B) {
    if (length(ki) != 1L || length(n) != 1L || length(B) != 1L || 
        any(!is.finite(c(ki, n, B))) || any(c(ki, n, B) <= 0)) 
        stop("ki, n and B must be finite positive scalars")
    get_range <- function(q) {
        a <- log(q) - n * log(ki)
        logden <- max(a, 0) + log1p(exp(-abs(a)))
        exp(log(B) + (log1p(-q) - logden)/n)
    }
    data.frame(r50 = get_range(0.5), r01 = get_range(0.01))
}
CalcExposure <- function(ligand_dist, decay_weights) {
    X <- as.matrix(ligand_dist)
    if (ncol(X) != length(decay_weights)) 
        stop("Distance bins and weights differ")
    .Nonnegative(X, "Binned expression")
    .Nonnegative(decay_weights, "Decay weights")
    out <- X %*% as.numeric(decay_weights)
    .Nonnegative(out, "Exposure")
    Matrix::Matrix(out, sparse = TRUE)
}
.PrepCoupling <- function(X, y) {
    X <- as.matrix(X)
    y <- as.numeric(y)
    if (nrow(X) != length(y) || length(y) < 3L) 
        stop("Insufficient matched receivers")
    .Nonnegative(X, "Binned expression")
    .Nonnegative(y, "Response")
    yc <- y - mean(y)
    ss <- sum(yc^2)
    if (!is.finite(ss) || ss <= 0) 
        stop("Response is constant or has invalid variance")
    Xc <- sweep(X, 2, colMeans(X), "-")
    list(Xc = Xc, A = crossprod(Xc), yc = yc, ss = ss)
}
.FitDecayFast <- function(prepared, yc, distance, B) {
    if (length(distance) != ncol(prepared$Xc) || length(yc) != 
        nrow(prepared$Xc)) 
        stop("Distance/response dimensions do not match prepared statistics")
    g <- as.vector(crossprod(prepared$Xc, yc))
    objective <- function(par) {
        w <- .DecayWeights(par, distance, B)
        v <- as.numeric(crossprod(w, prepared$A %*% w))
        size <- sum(abs(outer(w, w) * prepared$A))
        if (!is.finite(v) || v <= 1e-12 * max(size, .Machine$double.xmin)) 
            v <- sum(as.vector(prepared$Xc %*% w)^2)
        if (!is.finite(v) || v <= 0) 
            stop("Exposure has zero or invalid variance")
        coupling <- sum(w * g)/sqrt(v * prepared$ss)
        if (!is.finite(coupling)) 
            stop("Non-finite coupling")
        -max(-1, min(1, coupling))
    }
    fit <- stats::optim(c(0.01, 0.01), objective, method = "L-BFGS-B", 
        lower = c(0.01, 0.01), upper = c(1, 10), control = list(trace = 0))
    fit$value <- objective(fit$par)
    fit
}
.Timed <- function(fun, seconds) {
    if (length(seconds) != 1L || !is.finite(seconds) || seconds <= 
        0) 
        stop("Timeout must be positive")
    setTimeLimit(cpu = Inf, elapsed = seconds, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), 
        add = TRUE)
    fun()
}
FitDecay <- function(response, ligand_dist, distance, B, fit_timeout = 300) {
    p <- .PrepCoupling(ligand_dist, response)
    fit <- .Timed(function() .FitDecayFast(p, p$yc, distance, 
        B), fit_timeout)
    if (fit$convergence != 0L) 
        warning("Decay optimizer did not converge; valid finite fit accepted", 
            call. = FALSE)
    fit
}
SummarizeDecay <- function(ligand, response, exposure, fit, B) {
    y <- as.numeric(response)
    exposure <- as.numeric(exposure)
    .Nonnegative(y, "Response")
    .Nonnegative(exposure, "Exposure")
    if (length(y) != length(exposure) || length(y) < 3 || stats::sd(y) == 
        0 || stats::sd(exposure) == 0) 
        stop("Coupling needs at least three receivers with nonconstant response/exposure")
    test <- stats::cor.test(y, exposure, method = "pearson")
    ranges <- CalcRanges(fit$par[1], fit$par[2], B)
    data.frame(Ligand = ligand, ki = fit$par[1], n = fit$par[2], 
        coupling = unname(test$estimate), p.value = test$p.value, r50 = ranges$r50, 
        r01 = ranges$r01, convergence = as.integer(fit$convergence), 
        n_receivers = length(y))
}
.PermuteCoupling <- function(prepared, observed, distance, B, 
    n_perm, fit_timeout) {
    null <- rep(NA_real_, n_perm)
    conv <- rep(NA_integer_, n_perm)
    errors <- rep(NA_character_, n_perm)
    rng_start <- get(".Random.seed", envir = .GlobalEnv)
    for (m in seq_len(n_perm)) {
        yc <- prepared$yc[sample.int(length(prepared$yc))]
        fit <- tryCatch(.Timed(function() .FitDecayFast(prepared, 
            yc, distance, B), fit_timeout), error = function(e) e)
        if (inherits(fit, "error")) 
            errors[m] <- conditionMessage(fit)
        else {
            null[m] <- -fit$value
            conv[m] <- fit$convergence
        }
    }
    good <- is.finite(null)
    complete <- all(good)
    nonconverged <- sum(conv != 0L, na.rm = TRUE)
    if (nonconverged) 
        warning(nonconverged, " permutation fits did not converge; valid finite results accepted", 
            call. = FALSE)
    if (!complete) 
        warning("Incomplete permutation test: empirical P value is NA", 
            call. = FALSE)
    list(p.empirical = if (complete) (1 + sum(null >= observed))/(n_perm + 
        1) else NA_real_, n_perm_requested = n_perm, n_perm_valid = sum(good), 
        n_perm_failed = sum(!good), n_perm_nonconverged = nonconverged, 
        status = if (complete) "complete" else "incomplete", 
        observed = observed, null_cor = null, convergence = conv, 
        errors = errors, rng_start = rng_start, rng_kind = RNGkind())
}
.LigandBins <- function(arr, ligand, cells) {
    if (!ligand %in% dimnames(arr)[[1]] || !all(cells %in% dimnames(arr)[[2]])) 
        stop("Missing binned ligand/cells")
    matrix(arr[ligand, cells, , drop = FALSE], nrow = length(cells), 
        ncol = dim(arr)[3], dimnames = list(cells, dimnames(arr)[[3]]))
}
FitPropagation <- function(lidardata, CalcLigandConc = TRUE, 
    imp = TRUE, ligand = NULL, permutation_p = FALSE, n_perm = 999L, 
    perm_seed = 123L, fit_timeout = 300) {
    .CheckObject(lidardata)
    .Flag(CalcLigandConc, "CalcLigandConc")
    .Flag(permutation_p, "permutation_p")
    if (length(fit_timeout) != 1L || !is.finite(fit_timeout) || 
        fit_timeout <= 0) 
        stop("fit_timeout must be positive")
    if (permutation_p) {
        n_perm <- .Integer(n_perm, "n_perm")
        perm_seed <- .Integer(perm_seed, "perm_seed", min = 0L)
        had_seed <- exists(".Random.seed", envir = .GlobalEnv, 
            inherits = FALSE)
        if (had_seed) 
            old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit({
            if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", 
                envir = .GlobalEnv, inherits = FALSE)) rm(list = ".Random.seed", 
                envir = .GlobalEnv)
        }, add = TRUE)
        set.seed(perm_seed)
    }
    expr <- .Expression(lidardata, imp)
    cells <- lidardata@cellname
    response <- .AlignRows(lidardata@Diffusion@ligand_act, cells, 
        "Ligand response")
    arr <- lidardata@Diffusion@ligand_dist
    if (length(dim(arr)) != 3L || !dim(arr)[2] || dim(arr)[3] < 
        2L) 
        stop("Run BinLigandExpr first")
    if (!identical(attr(arr, "imp"), imp)) 
        stop("imp differs from BinLigandExpr; recompute distance expression")
    n_bins <- dim(arr)[3]
    bin <- lidardata@Diffusion@bin
    distance <- (seq_len(n_bins) - 1) * bin
    B <- (n_bins - 1) * bin
    ligands <- if (is.null(ligand)) 
        colnames(response)
    else unique(as.character(ligand))
    fits <- list()
    exposures <- list()
    permutations <- list()
    for (l in ligands) {
        result <- tryCatch({
            if (!l %in% colnames(response) || !l %in% rownames(expr$ligand)) 
                stop("Ligand lacks response/expression")
            .Nonnegative(response[, l], "Response")
            .Nonnegative(expr$ligand[l, ], "Ligand expression")
            receptor <- .Receptors(lidardata, l, expr$receptor)
            receivers <- cells[receptor > 0]
            y <- as.numeric(response[receivers, l])
            X <- .LigandBins(arr, l, receivers)
            if (sum(y > 0) <= 10L) 
                stop("Not enough positive responses (requires >10)")
            fit <- FitDecay(y, X, distance, B, fit_timeout)
            w <- .DecayWeights(fit$par, distance, B)
            exposure <- CalcExposure(X, w)
            if (sum(exposure > 0) <= 10L) 
                stop("Not enough positive exposures (requires >10)")
            row <- SummarizeDecay(l, y, exposure, fit, B)
            full <- if (CalcLigandConc) 
                CalcExposure(.LigandBins(arr, l, cells), w)
            else NULL
            perm <- NULL
            if (permutation_p) {
                message("Permutation test: ", l, " (", n_perm, 
                  " draws)")
                perm <- .PermuteCoupling(.PrepCoupling(X, y), 
                  -fit$value, distance, B, n_perm, fit_timeout)
                perm$receivers <- receivers
                perm$seed <- perm_seed
                row$p.empirical <- perm$p.empirical
                row$n_perm_requested <- perm$n_perm_requested
                row$n_perm_valid <- perm$n_perm_valid
                row$n_perm_failed <- perm$n_perm_failed
                row$n_perm_nonconverged <- perm$n_perm_nonconverged
                row$perm_status <- perm$status
                row$perm_seed <- perm_seed
            }
            list(row = row, exposure = full, permutation = perm)
        }, error = function(e) e)
        if (inherits(result, "error")) {
            lidardata <- .Failure(lidardata, l, "fit", conditionMessage(result))
            next
        }
        fits[[l]] <- result$row
        if (CalcLigandConc) 
            exposures[[l]] <- result$exposure
        if (permutation_p) 
            permutations[[l]] <- result$permutation
    }
    df <- if (length(fits)) 
        do.call(rbind, fits)
    else .EmptyFit()
    if (permutation_p) {
        if (!nrow(df)) {
            df$p.empirical <- numeric()
            df$n_perm_requested <- integer()
            df$n_perm_valid <- integer()
            df$n_perm_failed <- integer()
            df$n_perm_nonconverged <- integer()
            df$perm_status <- character()
            df$perm_seed <- integer()
        }
        attr(df, "permutation_results") <- permutations
    }
    attr(df, "fit_config") <- list(imp = imp, bin = bin, n_bins = n_bins, 
        B = B, distance = distance, initial = c(0.01, 0.01), 
        lower = c(0.01, 0.01), upper = c(1, 10), method = "L-BFGS-B", 
        fit_timeout = fit_timeout, range_definition = "continuous curve", 
        coordinate_units = "input units")
    lidardata@Diffusion@ligand_diff <- df
    conc <- if (length(exposures)) 
        do.call(cbind, exposures)
    else .EmptyMatrix(cells)
    if (length(exposures)) {
        rownames(conc) <- cells
        colnames(conc) <- names(exposures)
    }
    lidardata@Diffusion@ligand_conc <- conc
    lidardata@Lidarscore <- .ScoresEmpty(cells)
    lidardata
}
