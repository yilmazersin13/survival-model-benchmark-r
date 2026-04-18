
## 01_data_loaders.R
##

suppressPackageStartupMessages({
  library(survival)
})


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
  

  feature_types <- vapply(cov_names, function(nm) {
    if (is.numeric(raw[[nm]])) "continuous" else "categorical"
  }, character(1))
  

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
