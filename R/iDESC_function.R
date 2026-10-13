#' iDESC: differential expression in single-cell data with a negative binomial mixed model
#'
#' For every gene a negative binomial mixed model with the subject as random intercept is fitted
#' with `glmmTMB`: `y ~ group + covariates + offset(log(norm_sf)) + (1 | sub)`. The group effect is
#' tested with a likelihood-ratio test (full model against the same model without `group`), and
#' every group coefficient (and, optionally, every remaining pairwise contrast) with a Wald test.
#'
#' @section Model:
#' `zi_model = "none"` (default) fits the negative binomial mixed model without zero inflation
#' (`ziformula = ~ 0`). `zi_model = "constant"` adds one zero-inflation probability per gene
#' (`ziformula = ~ 1`). The zero-inflation formula of iDESC <= 0.2.0
#' (`~ offset(predict_pi_offset) - 1 + ind_zero + (1 | sub)`) was removed in 0.3.0, see `NEWS.md`.
#'
#' @section Gene status:
#' There is one status vocabulary, used in the `Status` column of `model_results`, in the `status`
#' column of `gene_status` and in the `Reason` column of `problematic_genes_cells_expressed`:
#' \describe{
#'   \item{`ok`}{the model was fitted and every group coefficient has a Wald p-value.}
#'   \item{`no_standard_errors`}{the model was fitted, but at least one group coefficient has no
#'   standard error / Wald p-value (typically a Hessian that is not positive definite). The
#'   estimates and the likelihood-ratio test are still reported.}
#'   \item{`fit_error`}{the model could not be fitted or its results could not be extracted; the
#'   error message is in `gene_status$message`.}
#'   \item{`one_group`}{the gene is expressed in cells of fewer than two groups; not modelled.}
#'   \item{`filtered`}{the gene was removed by the gene filter of [preprocessing()]; it appears
#'   only in `gene_status`.}
#' }
#' The optimizer convergence code and the positive-definite-Hessian flag are not folded into the
#' status: they are the columns `convergence` (0 = converged), `pdHess` and `convergence_null`
#' (null model of the likelihood-ratio test) of `gene_status`.
#'
#' @param mat Count matrix (genes x cells) with gene names as row names and cell names as column
#' names
#' @param meta Data frame including information for cells (row names = cell names, same order as
#' the columns of `mat`)
#' @param subject_var The name of subject information in meta
#' @param group_var The name of group/disease information for DE analysis in meta. The variable is
#' used as a factor; its first level is the reference.
#' @param norm_opt Option for normalizing factors: `"SeqDepth"` (default), `"SizeFactor"`, `"User"`
#' or `"None"`
#' @param user_sf Option for user-specific normalizing factors
#' @param sub_cell_filtering Filtering on cells within each subject
#' @param gene_sub_filtering Filtering on genes based on expression in each subject: a gene is kept
#' if, in at least one group, it is expressed in at least `gene_cell_filtering` of the cells in at
#' least `gene_sub_filtering` x (number of subjects of that group) subjects. The default 0 keeps
#' every gene.
#' @param gene_cell_filtering Minimum proportion of cells of a subject expressing the gene (used
#' together with `gene_sub_filtering`)
#' @param ncell_filtering Filtering on cells based on the number of genes expressed
#' @param covariates Optional. A vector of covariates (columns of `meta`) to include in the model.
#' @param cores Number of cores for parallel processing (default is 1)
#' @param zi_model Zero-inflation component: `"none"` (default, `ziformula = ~ 0`) or `"constant"`
#' (`ziformula = ~ 1`, the intercept is reported in `Theta` on the logit scale).
#' @param pairwise If `TRUE` (default) and the group has three or more levels, adds for every pair
#' of non-reference levels the contrast `Beta_<B>_vs_<A>` (= coefficient of B minus coefficient of
#' A, A being the earlier factor level), its Wald p-value `Pval_Beta_<B>_vs_<A>` (from the
#' covariance matrix of the fixed effects) and `FDR_Beta_<B>_vs_<A>`. Ignored with two levels.
#' @param fit_stats If `TRUE`, adds the columns `LogLik`, `Npar` and `AIC` of the fitted model.
#'
#' @return A list containing:
#' \itemize{
#'   \item{\strong{model_results}: A data frame with one row per gene that passed the gene filter
#'   (in the order of `mat`). Columns: `Status` (see the section Gene status); `Alpha` (intercept);
#'   for every coefficient of the group (`group<level>`) and of the covariates `Beta_<coef>`,
#'   `Pval_Beta_<coef>` (Wald) and `FDR_Beta_<coef>`; the pairwise contrasts (see `pairwise`);
#'   `LRT_stat`, `Pval_LRT_group` and `FDR_LRT_group` (likelihood-ratio test of the whole group
#'   effect, chi-square with one degree of freedom per group coefficient); `Dispersion`
#'   (1 / negative binomial size parameter); `Sigma2` (variance of the subject random intercept);
#'   `Theta` (zero-inflation intercept, `NA` for `zi_model = "none"`); `Minus2LogLik`
#'   (-2 x maximised log-likelihood of the full model) and, with `fit_stats = TRUE`, `LogLik`,
#'   `Npar` and `AIC`. All FDR columns are Benjamini-Hochberg adjustments across the genes that
#'   have the corresponding p-value. Genes with status `one_group` or `fit_error` have `NA` in
#'   every numeric column.}
#'   \item{\strong{problematic_genes_cells_expressed}: For the genes of `model_results` whose
#'   status is not `ok`: number of cells with non-zero expression in each group and `Reason`
#'   (the status). `NULL` when there is no such gene.}
#'   \item{\strong{total_cells}: A one-row data frame with the number of cells per group after
#'   filtering.}
#'   \item{\strong{gene_status}: A data frame with one row per input gene (row of `mat`) and the
#'   columns `gene`, `status`, `convergence`, `pdHess`, `convergence_null` and `message`
#'   (see the section Gene status).}
#' }
#' The log-likelihoods are the maximised values reported by the optimizer
#' (`-fit$fit$objective`), which are available even when the Hessian is not positive definite
#' (`stats::logLik()` of a `glmmTMB` fit returns `NA` in that case).
#' @export
#'
#' @examples
#' \donttest{
#' library(iDESC)
#' data(IPF_example)
#' mat <- IPF_example$mat
#' meta <- IPF_example$meta
#' sequencing_depth <- IPF_example$meta$sequencing_depth
#' result <- iDESC(mat, meta, subject_var = "subject", group_var = "disease",
#'                 norm_opt = "User", user_sf = sequencing_depth, cores = 2)
#' head(result$model_results)
#' table(result$gene_status$status)
#' }
iDESC <- function(mat, meta, subject_var, group_var,
                  norm_opt = c("SeqDepth", "SizeFactor", "User", "None"), user_sf = NULL,
                  sub_cell_filtering = 5, gene_sub_filtering = 0, gene_cell_filtering = 0.05,
                  ncell_filtering = 1, covariates = NULL, cores = 1,
                  zi_model = c("none", "constant"), pairwise = TRUE, fit_stats = FALSE) {
  norm_opt <- match.arg(norm_opt)
  zi_model <- match.arg(zi_model)

  # ---- input checks -------------------------------------------------------------------------
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("`mat` needs gene names (row names) and cell names (column names).")
  }
  missing_vars <- setdiff(c(subject_var, group_var, covariates), colnames(meta))
  if (length(missing_vars) > 0) {
    stop("Not found in `meta`: ", paste(missing_vars, collapse = ", "))
  }
  reserved <- c("y", "norm_sf", "group", "sub")
  if (any(covariates %in% reserved)) {
    stop("Covariates cannot be called ", paste(reserved, collapse = ", "), " (names used internally).")
  }
  input_genes <- rownames(mat)
  if (!is.null(user_sf) && is.null(names(user_sf))) {
    if (length(user_sf) != ncol(mat)) stop("`user_sf` must have one value per cell (column of `mat`).")
    names(user_sf) <- colnames(mat)
  }

  # ---- filtering and normalisation ----------------------------------------------------------
  dat <- preprocessing(mat, meta, subject_var, group_var, sub_cell_filtering = sub_cell_filtering,
                       gene_sub_filtering = gene_sub_filtering, gene_cell_filtering = gene_cell_filtering,
                       ncell_filtering = ncell_filtering, cores = cores)
  mat <- dat$mat
  meta <- dat$meta
  if (ncol(mat) == 0) stop("No cell left after filtering.")
  if (!is.null(user_sf)) user_sf <- user_sf[colnames(mat)]
  norm_factor <- normalization_factors(mat, norm_opt, user_sf)
  if (any(!is.finite(norm_factor)) || any(norm_factor <= 0)) {
    stop("Normalization factors must be finite and positive for every cell that is kept.")
  }

  if (any(is.na(meta[[group_var]])) || any(is.na(meta[[subject_var]]))) {
    stop("Missing values in `", group_var, "` or `", subject_var, "`.")
  }
  group <- droplevels(factor(meta[[group_var]]))
  group_levels <- levels(group)
  if (length(group_levels) < 2) stop("`", group_var, "` needs at least two groups after filtering.")
  subject <- as.character(meta[[subject_var]])

  # Data shared by all genes; the count (`y`) is added per gene
  base_df <- data.frame(norm_sf = as.numeric(norm_factor), group = group, sub = factor(subject))
  for (cv in covariates) {
    x <- meta[[cv]]
    if (any(is.na(x))) stop("Missing values in covariate `", cv, "`.")
    base_df[[cv]] <- if (is.character(x) || is.logical(x)) factor(x) else if (is.factor(x)) droplevels(x) else x
  }

  # ---- cells expressing each gene, per group -------------------------------------------------
  gene_group_expressed <- do.call(cbind, lapply(group_levels, function(l) {
    rowSums(mat[, group == l, drop = FALSE] > 0)
  }))
  dimnames(gene_group_expressed) <- list(rownames(mat), group_levels)
  gene_group_expressed_df <- as.data.frame(gene_group_expressed, check.names = FALSE)

  total_cells_per_group <- as.data.frame(as.list(table(group)))
  rownames(total_cells_per_group) <- "Total_cells"

  # Genes expressed in fewer than two groups are not modelled
  genes_to_exclude <- rownames(mat)[rowSums(gene_group_expressed > 0) < 2]
  if (length(genes_to_exclude) > 0) {
    message("The following genes were excluded because they were expressed in only one group: ",
            paste(genes_to_exclude, collapse = ", "))
  }
  model_genes <- setdiff(rownames(mat), genes_to_exclude)

  # ---- formulas, coefficients and contrasts --------------------------------------------------
  rhs_null <- paste(c(if (length(covariates) > 0) paste0("`", covariates, "`") else "1"), collapse = " + ")
  rhs_full <- paste(c("group", if (length(covariates) > 0) paste0("`", covariates, "`")), collapse = " + ")
  full_formula <- stats::as.formula(paste("y ~", rhs_full, "+ offset(log(norm_sf)) + (1 | sub)"))
  null_formula <- stats::as.formula(paste("y ~", rhs_null, "+ offset(log(norm_sf)) + (1 | sub)"))
  zi_formula <- switch(zi_model, none = ~ 0, constant = ~ 1)

  # Fixed-effect coefficients of the full model, in model order: group first, then covariates
  coef_model <- colnames(stats::model.matrix(stats::as.formula(paste("~", rhs_full)), base_df))[-1]
  group_coefs <- paste0("group", group_levels[-1])
  # Column names are made syntactically valid (as data.frame() did in earlier versions)
  coef_cols <- make.names(paste0("Beta_", coef_model))

  # Contrasts of the group effect: every coefficient against the reference level and, with
  # `pairwise`, the differences between non-reference levels
  n_gc <- length(group_coefs)
  contrast_mat <- diag(n_gc)
  contrast_names <- make.names(paste0("Beta_", group_coefs))
  pair_names <- character(0)
  if (isTRUE(pairwise) && n_gc >= 2) {
    non_ref <- group_levels[-1]
    for (a in seq_len(n_gc - 1)) {
      for (b in (a + 1):n_gc) {
        l <- numeric(n_gc)
        l[b] <- 1
        l[a] <- -1
        contrast_mat <- rbind(contrast_mat, l)
        pair_names <- c(pair_names, make.names(paste0("Beta_", non_ref[b], "_vs_", non_ref[a])))
      }
    }
  }
  contrast_names <- c(contrast_names, pair_names)
  dimnames(contrast_mat) <- list(contrast_names, group_coefs)

  # ---- model fitting --------------------------------------------------------------------------
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  lrt_tol <- 1e-6

  fit_model <- function(y, grp, full = TRUE, se = TRUE, start = NULL) {
    df <- base_df
    df$y <- y
    df$group <- grp
    glmmTMB::glmmTMB(formula = if (full) full_formula else null_formula, data = df,
                     family = glmmTMB::nbinom2, ziformula = zi_formula, se = se, start = start)
  }
  loglik_of <- function(fit) -as.numeric(fit$fit$objective)[1]

  # Null model (no group term): log-likelihood and starting values for a restart of the full model
  fit_null <- function(y) {
    f0 <- tryCatch(quiet(fit_model(y, group, full = FALSE, se = FALSE)), error = function(e) e)
    if (inherits(f0, "error")) return(NULL)
    ll0 <- loglik_of(f0)
    if (!is.finite(ll0)) return(NULL)
    start <- tryCatch({
      pl <- f0$obj$env$parList(f0$fit$par, f0$fit$parfull)
      b0 <- glmmTMB::fixef(f0)$cond
      beta <- stats::setNames(numeric(length(coef_model) + 1), c("(Intercept)", coef_model))
      if (!all(names(b0) %in% names(beta)) || any(!is.finite(b0))) stop("unexpected coefficients")
      beta[names(b0)] <- b0
      # the dispersion parameter is called betadisp (betad in older versions of glmmTMB)
      st <- c(list(beta = unname(beta)), pl[intersect(c("betadisp", "betad"), names(pl))], list(theta = pl$theta))
      if (zi_model == "constant") st$betazi <- pl$betazi
      st
    }, error = function(e) NULL)
    list(ll = ll0, convergence = as.integer(f0$fit$convergence)[1], start = start)
  }

  # Full model for a given vector of group labels. If the fit fails or ends below the null model
  # (the optimizer did not reach the maximum), it is fitted again starting from the null estimates
  # and the better of the two fits is kept. The same rule is used for observed and permuted labels.
  fit_full <- function(y, grp, null, se = TRUE) {
    f <- tryCatch(quiet(fit_model(y, grp, full = TRUE, se = se)), error = function(e) e)
    ll <- if (inherits(f, "error")) NA_real_ else loglik_of(f)
    if (!is.null(null) && !is.null(null$start) && (!is.finite(ll) || 2 * (ll - null$ll) < -lrt_tol)) {
      f2 <- tryCatch(quiet(fit_model(y, grp, full = TRUE, se = se, start = null$start)), error = function(e) e)
      if (!inherits(f2, "error")) {
        ll2 <- loglik_of(f2)
        if (is.finite(ll2) && (!is.finite(ll) || ll2 > ll)) f <- f2
      }
    }
    f
  }

  # Likelihood-ratio statistic of the group effect; NA when it cannot be computed
  lrt_statistic <- function(ll1, null) {
    if (is.null(null) || !is.finite(ll1)) return(NA_real_)
    stat <- 2 * (ll1 - null$ll)
    if (!is.finite(stat) || stat < -lrt_tol) return(NA_real_)
    max(stat, 0)
  }

  value_cols <- c("Alpha", coef_cols, paste0("Pval_", coef_cols),
                  unlist(lapply(pair_names, function(pn) c(pn, paste0("Pval_", pn)))),
                  "LRT_stat", "Pval_LRT_group", "Dispersion", "Sigma2", "Theta", "Minus2LogLik",
                  "LogLik", "Npar", "AIC")

  fit_gene <- function(gene_name) {
    out <- list(gene = gene_name, status = "fit_error", convergence = NA_integer_, pdHess = NA,
                convergence_null = NA_integer_, message = NA_character_, values = NULL,
                null = NULL, contrasts = NULL)
    y <- mat[gene_name, ]
    null <- fit_null(y)
    if (!is.null(null)) out$convergence_null <- null$convergence
    f1 <- fit_full(y, group, null, se = TRUE)
    # The failure is returned to the parent process (a forked worker cannot modify its variables)
    if (inherits(f1, "error")) {
      out$message <- conditionMessage(f1)
      return(out)
    }
    extracted <- tryCatch({
      values <- stats::setNames(rep(NA_real_, length(value_cols)), value_cols)
      ctab <- summary(f1)$coefficients$cond
      est <- ctab[, "Estimate"]
      values["Alpha"] <- est["(Intercept)"]
      has <- coef_model %in% rownames(ctab)
      values[coef_cols[has]] <- est[coef_model[has]]
      values[paste0("Pval_", coef_cols[has])] <- ctab[coef_model[has], "Pr(>|z|)"]

      # Group contrasts (also used as permutation statistics)
      beta_group <- stats::setNames(rep(NA_real_, n_gc), group_coefs)
      beta_group[group_coefs[group_coefs %in% names(est)]] <- est[group_coefs[group_coefs %in% names(est)]]
      contrasts <- drop(contrast_mat %*% beta_group)
      names(contrasts) <- contrast_names
      if (length(pair_names) > 0) {
        V <- tryCatch(quiet(stats::vcov(f1)$cond), error = function(e) NULL)
        for (pn in pair_names) {
          values[pn] <- contrasts[pn]
          l <- contrast_mat[pn, ]
          if (!is.null(V) && all(group_coefs %in% rownames(V))) {
            v <- drop(t(l) %*% V[group_coefs, group_coefs, drop = FALSE] %*% l)
            if (is.finite(v) && v > 0 && is.finite(contrasts[pn])) {
              values[paste0("Pval_", pn)] <- 2 * stats::pnorm(-abs(contrasts[pn] / sqrt(v)))
            }
          }
        }
      }

      # Likelihood-ratio test of the group effect
      ll1 <- loglik_of(f1)
      lrt_df <- sum(!is.na(beta_group))
      lrt_stat <- lrt_statistic(ll1, null)
      if (!is.na(lrt_stat) && lrt_df > 0) {
        values["LRT_stat"] <- lrt_stat
        values["Pval_LRT_group"] <- stats::pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)
      }

      values["Dispersion"] <- 1 / stats::sigma(f1)
      values["Sigma2"] <- exp(f1$fit$par[names(f1$fit$par) == "theta"][1])^2
      if (zi_model == "constant") values["Theta"] <- glmmTMB::fixef(f1)$zi[1]
      npar <- length(f1$fit$par)
      values["Minus2LogLik"] <- -2 * ll1
      values["LogLik"] <- ll1
      values["Npar"] <- npar
      values["AIC"] <- -2 * ll1 + 2 * npar

      group_p <- values[paste0("Pval_", make.names(paste0("Beta_", group_coefs)))]
      list(values = values, contrasts = contrasts,
           status = if (any(is.na(group_p))) "no_standard_errors" else "ok",
           convergence = as.integer(f1$fit$convergence)[1],
           pdHess = isTRUE(f1$sdr$pdHess))
    }, error = function(e) e)

    if (inherits(extracted, "error")) {
      out$message <- paste0("result extraction failed: ", conditionMessage(extracted))
      return(out)
    }
    out[names(extracted)] <- extracted
    out$null <- null
    out
  }

  res_list <- parallel::mclapply(model_genes, fit_gene, mc.cores = cores)
  # mclapply can return try-error (or NULL) elements, e.g. when a worker is killed
  for (i in seq_along(res_list)) {
    r <- res_list[[i]]
    if (!is.list(r) || is.null(r$status)) {
      msg <- if (inherits(r, "try-error")) trimws(as.character(r)[1]) else "worker returned no result"
      res_list[[i]] <- list(gene = model_genes[i], status = "fit_error", convergence = NA_integer_,
                            pdHess = NA, convergence_null = NA_integer_, message = msg, values = NULL,
                            null = NULL, contrasts = NULL)
    }
  }
  names(res_list) <- model_genes

  # ---- per-gene status: one row per input gene ---------------------------------------------
  n_in <- length(input_genes)
  gene_status <- data.frame(gene = input_genes, status = rep("filtered", n_in),
                            convergence = rep(NA_integer_, n_in), pdHess = rep(NA, n_in),
                            convergence_null = rep(NA_integer_, n_in),
                            message = rep(NA_character_, n_in), stringsAsFactors = FALSE)
  gene_status$status[input_genes %in% genes_to_exclude] <- "one_group"
  if (length(res_list) > 0) {
    idx <- match(model_genes, input_genes)
    gene_status$status[idx] <- vapply(res_list, function(r) r$status, character(1))
    gene_status$convergence[idx] <- vapply(res_list, function(r) as.integer(r$convergence), integer(1))
    gene_status$pdHess[idx] <- vapply(res_list, function(r) as.logical(r$pdHess), logical(1))
    gene_status$convergence_null[idx] <- vapply(res_list, function(r) as.integer(r$convergence_null), integer(1))
    gene_status$message[idx] <- vapply(res_list, function(r) as.character(r$message), character(1))
  }
  failed <- gene_status$gene[gene_status$status == "fit_error"]
  if (length(failed) > 0) {
    message("The model could not be fitted for: ", paste(failed, collapse = ", "))
  }

  # ---- results table: one row per gene that passed the gene filter ---------------------------
  res_genes <- rownames(mat)
  values <- matrix(NA_real_, length(res_genes), length(value_cols), dimnames = list(res_genes, value_cols))
  for (r in res_list) {
    if (!is.null(r$values)) values[r$gene, ] <- r$values
  }
  bh <- function(p) stats::p.adjust(p, method = "BH")
  res_tb <- data.frame(Status = gene_status$status[match(res_genes, gene_status$gene)],
                       Alpha = values[, "Alpha"], row.names = res_genes, stringsAsFactors = FALSE)
  for (cn in c(coef_cols, pair_names)) {
    res_tb[[cn]] <- values[, cn]
    res_tb[[paste0("Pval_", cn)]] <- values[, paste0("Pval_", cn)]
    res_tb[[paste0("FDR_", cn)]] <- bh(values[, paste0("Pval_", cn)])
  }
  res_tb[["LRT_stat"]] <- values[, "LRT_stat"]
  res_tb[["Pval_LRT_group"]] <- values[, "Pval_LRT_group"]
  res_tb[["FDR_LRT_group"]] <- bh(values[, "Pval_LRT_group"])
  for (cn in c("Dispersion", "Sigma2", "Theta", "Minus2LogLik", if (isTRUE(fit_stats)) c("LogLik", "Npar", "AIC"))) {
    res_tb[[cn]] <- values[, cn]
  }

  # Cells expressing the genes that are not "ok"
  problematic_genes <- res_genes[res_tb$Status != "ok"]
  if (length(problematic_genes) > 0) {
    problematic_genes_expressed <- gene_group_expressed_df[problematic_genes, , drop = FALSE]
    problematic_genes_expressed$Reason <- res_tb[problematic_genes, "Status"]
  } else {
    problematic_genes_expressed <- NULL
  }

  list(
    model_results = res_tb,
    problematic_genes_cells_expressed = problematic_genes_expressed,
    total_cells = total_cells_per_group,
    gene_status = gene_status
  )
}
