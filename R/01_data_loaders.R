## ============================================================
## 01_data_loaders.R
##
## Standardized data loaders for the three benchmark datasets
## referenced in Section IV-A of the paper (METABRIC, SUPPORT,
## TCGA-BRCA).
##
## Each loader returns a list with the following slots:
##
##   $time              numeric vector, length n, time in days
##                      (T_tilde_i in Section III-A notation)
##   $event             integer vector in {0,1}, length n
##                      (delta_i in Section III-A notation)
##   $X                 numeric matrix, n x p_raw, covariates
##                      BEFORE preprocessing (no scaling, no
##                      imputation, no encoding); preprocessing
##                      is the responsibility of 02_preprocessing.R
##                      and is performed fold-wise.
##   $feature_names     character vector of length p_raw
##   $feature_types     character vector of length p_raw, each
##                      element in {"continuous","categorical"}
##   $dataset_name      character scalar
##   $baseline_censoring numeric scalar in (0,1), the empirical
##                      censoring proportion 1 - mean(event)
##   $time_unit         character scalar, always "days" for
##                      consistency with the existing METABRIC code
##
## NO train/test splitting happens here. Splitting is the job of
## the 5x2 cross-validation scaffolding in 04_cv_scaffolding.R,
## as required by Section IV-D.
##
## NO censoring augmentation happens here either. Augmentation is
## applied per training fold in 03_censoring_augmentation.R, as
## required by Section IV-C.
## ============================================================

suppressPackageStartupMessages({
  library(survival)
})

## ------------------------------------------------------------
## Helper: validate the returned object so downstream modules
## can rely on the schema without defensive checks everywhere.
## ------------------------------------------------------------
.validate_dataset_object <- function(obj) {
  required <- c("time", "event", "X", "feature_names", "feature_types",
                "dataset_name", "baseline_censoring", "time_unit")
  missing_slots <- setdiff(required, names(obj))
  if (length(missing_slots) > 0) {
    stop(sprintf("Dataset object missing slots: %s",
                 paste(missing_slots, collapse = ", ")))
  }
  n <- length(obj$time)
  if (length(obj$event) != n) stop("time and event length mismatch")
  if (nrow(obj$X) != n) stop("X row count does not match time length")
  if (ncol(obj$X) != length(obj$feature_names)) {
    stop("X column count does not match feature_names length")
  }
  if (length(obj$feature_types) != length(obj$feature_names)) {
    stop("feature_types and feature_names length mismatch")
  }
  if (!all(obj$feature_types %in% c("continuous", "categorical"))) {
    stop("feature_types must contain only 'continuous' or 'categorical'")
  }
  if (!all(obj$event %in% c(0L, 1L))) {
    stop("event must be 0/1 integer")
  }
  if (any(obj$time <= 0)) {
    stop("time must be strictly positive")
  }
  invisible(obj)
}


## ============================================================
## METABRIC loader
##
## Adapted from the user's existing METABRIC preprocessing code.
## Differences from the original:
##   - Returns the RAW gene matrix (no scaling, no train/test
##     split). Scaling moves to 02_preprocessing.R per fold.
##   - Time is converted to days (months * 30.44) to match the
##     existing convention.
##   - The 10th-percentile variance filter is RETAINED here
##     because it is an unsupervised, dataset-level filter that
##     does not use outcome information and is identical across
##     folds. This matches the spirit of Section IV-B.
## ============================================================
load_metabric <- function(data_file = "METABRIC_RNA_Mutation.csv") {
  if (!file.exists(data_file)) {
    stop(sprintf("METABRIC file not found: %s", data_file))
  }
  cat(sprintf("Loading METABRIC from %s ...\n", data_file))
  
  raw <- read.csv(data_file, stringsAsFactors = FALSE)
  
  ## --- Survival columns ---
  time_days <- raw$overall_survival_months * 30.44
  event     <- as.integer(raw$overall_survival)
  
  keep <- !is.na(time_days) & !is.na(event) & time_days > 0
  raw       <- raw[keep, , drop = FALSE]
  time_days <- time_days[keep]
  event     <- event[keep]
  
  ## --- Locate the gene expression block ---
  ## In the standard METABRIC CSV from cBioPortal the expression
  ## block runs from "brca1" to the column just before the first
  ## "_mut" mutation column. Lower-case fallback is checked first
  ## because that is what the existing user code assumes.
  gene_start <- which(names(raw) == "brca1")
  if (length(gene_start) == 0) gene_start <- which(names(raw) == "BRCA1")
  if (length(gene_start) == 0) {
    stop("Could not locate gene expression block in METABRIC file")
  }
  mut_cols <- which(grepl("_mut$", names(raw)))
  gene_end <- if (length(mut_cols) > 0) mut_cols[1] - 1 else ncol(raw)
  
  gene_matrix <- as.matrix(raw[, gene_start:gene_end, drop = FALSE])
  
  ## --- Drop genes with any missing values across the cohort ---
  complete_genes <- complete.cases(t(gene_matrix))
  gene_matrix <- gene_matrix[, complete_genes, drop = FALSE]
  
  ## --- Unsupervised low-variance filter (10th percentile) ---
  ## This is the same filter as the user's existing code. It is
  ## unsupervised and outcome-free so it is safe to apply at
  ## load time rather than per fold.
  vars <- apply(gene_matrix, 2, var, na.rm = TRUE)
  keep_var <- vars >= quantile(vars, 0.10, na.rm = TRUE)
  gene_matrix <- gene_matrix[, keep_var, drop = FALSE]
  
  feature_names <- colnames(gene_matrix)
  feature_types <- rep("continuous", length(feature_names))
  
  obj <- list(
    time               = time_days,
    event              = event,
    X                  = gene_matrix,
    feature_names      = feature_names,
    feature_types      = feature_types,
    dataset_name       = "METABRIC",
    baseline_censoring = 1 - mean(event),
    time_unit          = "days"
  )
  .validate_dataset_object(obj)
  
  cat(sprintf("  METABRIC loaded: n = %d, p = %d, baseline censoring = %.3f\n",
              length(obj$time), ncol(obj$X), obj$baseline_censoring))
  obj
}


## ============================================================
## SUPPORT loader
##
## Reads the publicly distributed SUPPORT TSV. Columns observed
## in the user's file:
##   age, death, sex, hospdead, slos, d.time, dzgroup, dzclass,
##   num.co, edu, income, scoma, charges, totcst, totmcst,
##   avtisst, race, meanbp, wblc, hrt, resp, temp, pafi, alb,
##   bili, crea, sod, ph, glucose, bun, urine, adlp, adls,
##   sfdm2, adlsc
##
## Survival endpoint per Section IV-A:
##   time  = d.time  (days)
##   event = death   (0/1)
##
## NOTE FOR PAPER: The user's SUPPORT file has n = 1000, while
## Section IV-A states n ~ 9000. This is the publicly distributed
## subset; either obtain the full cohort or update Section IV-A
## before running the benchmark. A warning is printed at load
## time so this does not get forgotten.
##
## Covariates excluded by design (rationale in comments):
##   hospdead -- in-hospital death, leaks the outcome
##   slos     -- length of hospital stay, leaks the outcome
##   sfdm2    -- functional disability at follow-up, leaks the
##               outcome
##   charges, totcst, totmcst -- post-hoc cost measures
##   adlp, adls, adlsc -- post-hoc activities of daily living
## ============================================================
load_support <- function(data_file = "support.tsv") {
  if (!file.exists(data_file)) {
    stop(sprintf("SUPPORT file not found: %s", data_file))
  }
  cat(sprintf("Loading SUPPORT from %s ...\n", data_file))
  
  raw <- read.delim(data_file, stringsAsFactors = FALSE)
  
  ## --- Survival columns ---
  if (!all(c("d.time", "death") %in% names(raw))) {
    stop("SUPPORT file must contain 'd.time' and 'death' columns")
  }
  time_days <- as.numeric(raw$d.time)
  event     <- as.integer(raw$death)
  
  keep <- !is.na(time_days) & !is.na(event) & time_days > 0
  raw       <- raw[keep, , drop = FALSE]
  time_days <- time_days[keep]
  event     <- event[keep]
  
  ## --- Covariate selection ---
  ## Exclude the survival columns themselves and the outcome-leaking
  ## or post-hoc variables listed in the file header comment.
  exclude <- c("d.time", "death", "hospdead", "slos", "sfdm2",
               "charges", "totcst", "totmcst",
               "adlp", "adls", "adlsc")
  cov_names <- setdiff(names(raw), exclude)
  
  ## --- Classify each covariate as continuous or categorical ---
  ## Numeric columns -> continuous; character columns -> categorical.
  ## Missingness is preserved here and handled in 02_preprocessing.R.
  feature_types <- vapply(cov_names, function(nm) {
    if (is.numeric(raw[[nm]])) "continuous" else "categorical"
  }, character(1))
  
  ## --- Build the raw covariate matrix ---
  ## Categorical columns are kept as character at this stage and
  ## one-hot encoded later in 02_preprocessing.R, which is the
  ## point at which fold-wise level handling can be done correctly.
  ## To keep the schema uniform (X is a numeric matrix), categorical
  ## columns are temporarily encoded as integer factor codes here
  ## and re-encoded properly in the preprocessor. The original
  ## character labels are preserved in an attribute.
  X_list <- lapply(cov_names, function(nm) {
    col <- raw[[nm]]
    if (is.character(col)) {
      f <- factor(col)
      attr_levels <- levels(f)
      out <- as.integer(f)
      attr(out, "levels") <- attr_levels
      out
    } else {
      as.numeric(col)
    }
  })
  X <- do.call(cbind, X_list)
  colnames(X) <- cov_names
  
  ## Preserve factor levels for categorical columns so the
  ## preprocessor can reconstruct the original labels.
  cat_levels <- list()
  for (j in seq_along(cov_names)) {
    if (feature_types[j] == "categorical") {
      cat_levels[[ cov_names[j] ]] <- attr(X_list[[j]], "levels")
    }
  }
  attr(X, "categorical_levels") <- cat_levels
  
  obj <- list(
    time               = time_days,
    event              = event,
    X                  = X,
    feature_names      = cov_names,
    feature_types      = unname(feature_types),
    dataset_name       = "SUPPORT",
    baseline_censoring = 1 - mean(event),
    time_unit          = "days"
  )
  .validate_dataset_object(obj)
  
  cat(sprintf("  SUPPORT loaded: n = %d, p = %d, baseline censoring = %.3f\n",
              length(obj$time), ncol(obj$X), obj$baseline_censoring))
  
  if (length(obj$time) < 5000) {
    warning(sprintf(
      "SUPPORT cohort has n = %d, but Section IV-A states n ~ 9000. ",
      length(obj$time)),
      "This is the publicly distributed subset. Either obtain the ",
      "full cohort or update Section IV-A before running the benchmark.",
      call. = FALSE)
  }
  
  obj
}


## ============================================================
## TCGA-BRCA loader
##
## Reads the cBioPortal-style clinical TSV the user provided.
## The file contains 63 variables, all clinical (no expression
## matrix). Survival endpoint per Section IV-A:
##   time  = Overall.Survival..Months.  (converted to days)
##   event = Overall.Survival.Status    (decoded to 0/1)
##
## NOTE FOR PAPER: Section IV-B describes TCGA-BRCA as a
## high-dimensional MOLECULAR cohort reduced to p = 500 most
## variable probes. The clinical-only file the user provided
## has p ~ 30 after preprocessing. Two options:
##
##   (a) Provide an expression matrix file (e.g., from
##       cBioPortal mRNA Z-scores or GDC). The loader accepts
##       an optional `expression_file` argument that, when
##       supplied, is merged on Patient.ID and concatenated to
##       the clinical block.
##   (b) Update Section IV-B to describe TCGA-BRCA as the
##       clinical cohort.
##
## A warning is printed at load time when no expression file is
## supplied so this decision does not get forgotten.
##
## Covariate selection: identifier columns, free-text columns,
## and outcome-related columns are excluded by name. Disease-Free
## and Progression-Free fields are also excluded since they refer
## to alternative endpoints, not predictors.
## ============================================================
load_tcga_brca <- function(data_file = "tcga_brca_data.tsv",
                           expression_file = NULL) {
  if (!file.exists(data_file)) {
    stop(sprintf("TCGA-BRCA clinical file not found: %s", data_file))
  }
  cat(sprintf("Loading TCGA-BRCA clinical from %s ...\n", data_file))
  
  raw <- read.delim(data_file, stringsAsFactors = FALSE,
                    na.strings = c("", "NA", "[Not Available]",
                                   "[Unknown]", "[Not Evaluated]"))
  
  ## --- Survival columns ---
  surv_time_col  <- "Overall.Survival..Months."
  surv_event_col <- "Overall.Survival.Status"
  if (!all(c(surv_time_col, surv_event_col) %in% names(raw))) {
    stop(sprintf("TCGA-BRCA file must contain '%s' and '%s'",
                 surv_time_col, surv_event_col))
  }
  
  time_days <- as.numeric(raw[[surv_time_col]]) * 30.44
  
  ## OS status in cBioPortal exports is encoded as
  ## "0:LIVING" / "1:DECEASED". Decode robustly.
  os_status <- raw[[surv_event_col]]
  event <- rep(NA_integer_, length(os_status))
  event[grepl("^0", os_status) | os_status == "LIVING"]   <- 0L
  event[grepl("^1", os_status) | os_status == "DECEASED"] <- 1L
  
  keep <- !is.na(time_days) & !is.na(event) & time_days > 0
  raw       <- raw[keep, , drop = FALSE]
  time_days <- time_days[keep]
  event     <- event[keep]
  
  ## --- Covariate selection ---
  ## Exclude identifiers, free-text fields, alternative endpoints,
  ## and date columns that are not directly usable as predictors.
  exclude <- c(
    surv_time_col, surv_event_col,
    "Study.ID", "Patient.ID", "Sample.ID", "Other.Patient.ID",
    "Form.completion.date",
    "Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date",
    "Birth.from.Initial.Pathologic.Diagnosis.Date",
    "Last.Alive.Less.Initial.Pathologic.Diagnosis.Date.Calculated.Day.Value",
    "Disease.Free..Months.", "Disease.Free.Status",
    "Months.of.disease.specific.survival",
    "Disease.specific.Survival.status",
    "Progress.Free.Survival..Months.", "Progression.Free.Status",
    "Informed.consent.verified",
    "Cancer.Type", "TCGA.PanCanAtlas.Cancer.Type.Acronym",
    "Cancer.Type.Detailed",
    "American.Joint.Committee.on.Cancer.Publication.Version.Type",
    "ICD.10.Classification",
    "International.Classification.of.Diseases.for.Oncology..Third.Edition.ICD.O.3.Histology.Code",
    "International.Classification.of.Diseases.for.Oncology..Third.Edition.ICD.O.3.Site.Code",
    "Oncotree.Code",
    "In.PanCan.Pathway.Analysis",
    "Tissue.Source.Site", "Tissue.Source.Site.Code",
    "Number.of.Samples.Per.Patient",
    "Sample.Type",
    "Somatic.Status",
    "New.Neoplasm.Event.Post.Initial.Therapy.Indicator",
    "Neoadjuvant.Therapy.Type.Administered.Prior.To.Resection.Text"
  )
  cov_names <- setdiff(names(raw), exclude)
  
  ## Drop covariates that are constant or all-NA after filtering.
  not_useless <- vapply(cov_names, function(nm) {
    col <- raw[[nm]]
    n_unique <- length(unique(col[!is.na(col)]))
    n_unique >= 2
  }, logical(1))
  cov_names <- cov_names[not_useless]
  
  ## --- Classify ---
  feature_types <- vapply(cov_names, function(nm) {
    if (is.numeric(raw[[nm]])) "continuous" else "categorical"
  }, character(1))
  
  ## --- Build the raw covariate matrix (same convention as SUPPORT) ---
  X_list <- lapply(cov_names, function(nm) {
    col <- raw[[nm]]
    if (is.character(col) || is.factor(col)) {
      f <- factor(col)
      out <- as.integer(f)
      attr(out, "levels") <- levels(f)
      out
    } else {
      as.numeric(col)
    }
  })
  X <- do.call(cbind, X_list)
  colnames(X) <- cov_names
  
  cat_levels <- list()
  for (j in seq_along(cov_names)) {
    if (feature_types[j] == "categorical") {
      cat_levels[[ cov_names[j] ]] <- attr(X_list[[j]], "levels")
    }
  }
  attr(X, "categorical_levels") <- cat_levels
  
  ## --- Optional: merge an expression matrix ---
  ## When supplied, expression_file should be a TSV with one row
  ## per patient and the first column equal to Patient.ID. All
  ## other columns are treated as continuous molecular features.
  ## This is the path that would restore the p ~ 500 setting of
  ## Section IV-B once the user provides an expression file.
  if (!is.null(expression_file)) {
    if (!file.exists(expression_file)) {
      stop(sprintf("Expression file not found: %s", expression_file))
    }
    cat(sprintf("  Merging expression matrix from %s ...\n",
                expression_file))
    expr <- read.delim(expression_file, stringsAsFactors = FALSE,
                       check.names = FALSE)
    if (!"Patient.ID" %in% names(expr)) {
      stop("Expression file must have a 'Patient.ID' column")
    }
    pid <- raw$Patient.ID
    expr_match <- expr[match(pid, expr$Patient.ID), , drop = FALSE]
    expr_match$Patient.ID <- NULL
    expr_mat <- as.matrix(expr_match)
    storage.mode(expr_mat) <- "numeric"
    X <- cbind(X, expr_mat)
    expr_names <- colnames(expr_mat)
    cov_names <- c(cov_names, expr_names)
    feature_types <- c(feature_types, rep("continuous", length(expr_names)))
  }
  
  obj <- list(
    time               = time_days,
    event              = event,
    X                  = X,
    feature_names      = cov_names,
    feature_types      = unname(feature_types),
    dataset_name       = "TCGA-BRCA",
    baseline_censoring = 1 - mean(event),
    time_unit          = "days"
  )
  .validate_dataset_object(obj)
  
  cat(sprintf("  TCGA-BRCA loaded: n = %d, p = %d, baseline censoring = %.3f\n",
              length(obj$time), ncol(obj$X), obj$baseline_censoring))
  
  if (is.null(expression_file)) {
    warning(
      "TCGA-BRCA loaded WITHOUT an expression matrix. ",
      "Section IV-B describes a high-dimensional cohort with p = 500. ",
      "Either supply expression_file= or update Section IV-B to ",
      "describe the clinical cohort before running the benchmark.",
      call. = FALSE)
  }
  
  obj
}


## ============================================================
## Convenience: load all three datasets at once.
## ============================================================
load_all_datasets <- function(metabric_file = "METABRIC_RNA_Mutation.csv",
                              support_file  = "support.tsv",
                              tcga_file     = "tcga_brca_data.tsv",
                              tcga_expression_file = NULL) {
  list(
    METABRIC  = load_metabric(metabric_file),
    SUPPORT   = load_support(support_file),
    `TCGA-BRCA` = load_tcga_brca(tcga_file, tcga_expression_file)
  )
}

#---------------------------------------------------------------------------------

#Example run
#datasets <- load_all_datasets(
#  metabric_file = "METABRIC_RNA_Mutation.csv",
#  support_file  = "support.tsv",
#  tcga_file     = "tcga_brca_data.tsv"
#)
#str(datasets, max.level = 2)

#str(datasets$METABRIC, max.level = 1)
#str(datasets$SUPPORT, max.level = 1)
#str(datasets$`TCGA-BRCA`, max.level = 1)
