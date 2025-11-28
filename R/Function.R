#' alignment_wrapper: count 'good' and 'total' cases per TR
#'
#' @description
#' Compare input region indices against the reference TR database to count
#' the number of 'good' and total informative cases for each TR.
#'
#' @param input_vec A numeric vector of indices from \code{import_input_regions()}.
#' @param bin_width Width of bin in base pairs.
#' @param genome The genome of TR ChIP-seq data, either \code{"hg38"} or \code{"mm10"}.
#'
#' @return A data frame with columns \code{TR}, \code{GOOD}, and \code{TOTAL}.
#' @export
alignment_wrapper <- function(input_vec, filter_vec=NULL, bin_width, genome = c("hg38", "mm10")) {

  # Validate genome input
  genome <- match.arg(genome)

  # Validate input vector
  if (!is.numeric(input_vec) || length(input_vec) == 0) {
    stop("input_vec must be a non-empty numeric vector")
  }

  # Load the meta table
  message("Loading meta table...")
  meta_table_path <- system.file("meta_table.rds", package = "TRex")

  if (!file.exists(meta_table_path)) {
    stop("Meta table not found in package directory. Please ensure the package is installed correctly.")
  }

  meta_table <- readRDS(meta_table_path)

  # Get the appropriate ChIP-seq reference data
  meta_key <- paste0("meta_", genome, "_", bin_width)
  file_table <- meta_table[[meta_key]]

  if (is.null(file_table)) {
    stop("ChIP-seq files not found for genome '", genome, "' and bin width '", bin_width, "'. ",
         "Please download and load the ChIP-seq data first using load_chip_data().\n",
         "You may follow the tutorial on: https://github.com/ZeyuL01/TRex")
  }

  # Validate file table structure
  if (!all(c("TR", "File_Path") %in% colnames(file_table))) {
    stop("Invalid meta table structure. Expected columns 'TR' and 'File_Path'.")
  }

  if (nrow(file_table) == 0) {
    stop("No ChIP-seq files found in the meta table.")
  }

  # Initialize the results table
  chip_table <- data.frame(
    TR = file_table$TR,
    GOOD = numeric(nrow(file_table)),
    TOTAL = numeric(nrow(file_table)),
    stringsAsFactors = FALSE
  )

  # Progress bar setup
  message("Starting alignment process for ", nrow(chip_table), " TR ChIP-seq datasets...")
  pb <- txtProgressBar(min = 0, max = nrow(chip_table), style = 3, width = 50, char = "=")

  # Preallocate vectors for storing results
  good_vec <- numeric(nrow(chip_table))
  total_vec <- numeric(nrow(chip_table))

  # Loop over each ChIP-seq file and perform alignment
  for (i in seq_len(nrow(chip_table))) {
    # Read reference vector from file
    ref_vec <- data.table::fread(file_table$File_Path[i])[[1]]
    if(!is.null(filter_vec)){
      ref_vec <- intersect(ref_vec,filter_vec)
    }
    # Perform alignment using C++ function
    alignment_result <- Alignment(input_vec, ref_vec)

    # Store results in pre-allocated vectors
    good_vec[i] <- alignment_result$Xi_GOOD
    total_vec[i] <- alignment_result$Ni_TOTAL

    # Update progress bar
    setTxtProgressBar(pb, i)
  }

  # Close the progress bar
  close(pb)

  # Update chip_table with results
  chip_table$GOOD <- good_vec
  chip_table$TOTAL <- total_vec

  message("Alignment completed successfully. Processed ", nrow(chip_table), " TRs.")

  return(chip_table)
}


#' Transform the input file to vector
#'
#' @param file file path to the user-input.
#' @param format format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
#' @param bin_width desired width of bin, default: 1000.
#' @param genome the genome of TR ChIP-seq data, either as "hg38" or "mm10".
#'
#' @return A numeric vector contains the index of peaks with pre-specified number of bins in each chromosome.
#' @export
import_input_regions <- function(file, format = NULL, bin_width = 1000, genome=c("hg38", "mm10")) {
  # Determine the format if not provided
  if (is.null(format)) {
    extensions <- tools::file_ext(file)
    format <- switch(extensions,
                     bed = "bed",
                     bb = "bigNarrowPeak",
                     narrowPeak = "narrowPeak",
                     broadPeak = "broadPeak",
                     csv = "csv",
                     stop("File type not in (bed, bigNarrowPeak, narrowPeak, broadPeak, csv), please check the file type.")
    )
  }

  # Import data based on the format
  peak_dat <- switch(format,
                     bed = rtracklayer::import(file, format = "BED"),
                     narrowPeak = rtracklayer::import(file, format = "BED", extraCols = c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")),
                     broadPeak = rtracklayer::import(file, format = "BED", extraCols = c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")),
                     bigNarrowPeak = rtracklayer::import(file),
                     csv = read.csv(file),
                     stop("Unsupported format")
  )

  # Define fixed window numbers
  if (genome == "hg38") {
    N <- 3031030 * (1000 / bin_width)
    chr_windows <- c(0, 248956, 242193, 198295, 190214, 181538, 170805, 159345, 145138, 138394, 133797, 135086, 133275, 114364, 107043, 101991, 90338, 83257, 80373, 58617, 64444, 46709, 50818, 156040)
    chr_labels <- paste0("chr", c(1:22, "X"))
  } else if (genome == "mm10") {
    N <- 2631636 * (1000 / bin_width)
    chr_windows <- c(0, 195372, 182014, 159940, 156359, 151735, 149587, 145342, 129302, 124496, 130595, 121983, 120029, 120322, 124803, 103944, 98108, 94888, 90603, 61332, 170882)
    chr_labels <- paste0("chr", c(1:19, "X"))
  } else {
    stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
  }

  chr_windows <- chr_windows * (1000 / bin_width)
  chr_windows_cs <- cumsum(chr_windows)

  bin_inds <- c()

  for (i in seq_along(chr_labels)) {
    chr_lab <- chr_labels[i]
    chr_win_num <- chr_windows[i + 1]

    inds <- switch(format,
                   bed = {
                     gr <- peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab]
                     if (length(gr)) (GenomicRanges::start(gr) + GenomicRanges::end(gr)) %/% (2 * bin_width) + 1 else integer()
                   },
                   narrowPeak = {
                     gr <- peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab]
                     if (length(gr) && "peak" %in% names(mcols(gr))) {
                       (GenomicRanges::start(gr) + mcols(gr)$peak) %/% bin_width + 1
                     } else integer()
                   },
                   broadPeak = {
                     gr <- peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab]
                     if (length(gr)) (GenomicRanges::start(gr) + GenomicRanges::end(gr)) %/% (2 * bin_width) + 1 else integer()
                   },
                   bigNarrowPeak = {
                     gr <- peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab]
                     if (length(gr) && "abs_summit" %in% names(mcols(gr))) {
                       mcols(gr)$abs_summit %/% bin_width + 1
                     } else integer()
                   },
                   csv = {
                     idx <- which(peak_dat$Chrom == chr_lab)
                     if (length(idx)) (peak_dat$End[idx] + peak_dat$Start[idx]) %/% (2 * bin_width) + 1 else integer()
                   }
    )

    # keep in-range, unique, finite
    inds <- inds[is.finite(inds) & inds >= 1 & inds <= chr_win_num]
    inds <- unique(inds)
    if (length(inds)) bin_inds <- c(bin_inds, chr_windows_cs[i] + inds)
  }

  return(sort(bin_inds))
}


#' Load pre-compiled ChIP-seq data
#'
#' @description
#' Load and organize pre-compiled ChIP-seq data for TRex analysis. This function
#' scans the specified directory for ChIP-seq files, extracts transcription regulator
#' (TR) labels, and creates a meta table for efficient data access. Please follow
#' the tutorial on: https://github.com/ZeyuL01/TRex.
#'
#' @param data_path Character string. Path to the ChIP-seq data folder, can be
#'   absolute or relative path.
#' @param bin_width Integer. Width of bin in base pairs, must be one of: 100, 200, 500, 1000.
#'   This should match your ChIP-seq data resolution.
#' @param genome Character string. The reference genome, either "hg38" or "mm10".
#' @param overwrite Logical. Whether to overwrite existing meta table for the same
#'   genome and bin width combination. Default: FALSE.
#'
#' @return Invisibly returns the created meta table data frame.
#' @export
load_chip_data <- function(data_path, bin_width, genome = c("hg38", "mm10"), overwrite = FALSE) {

  # Validate genome input
  genome <- match.arg(genome)

  # Validate bin_width
  supported_bin_widths <- c(100, 200, 500, 1000)
  if (!bin_width %in% supported_bin_widths) {
    stop("bin_width must be one of: ", paste(supported_bin_widths, collapse = ", "))
  }

  # Convert to absolute path and validate directory
  data_path <- R.utils::getAbsolutePath(data_path)
  if (!dir.exists(data_path)) {
    stop("ChIP-seq data directory does not exist: ", data_path)
  }

  # Get meta table file path
  meta_table_path <- file.path(system.file(package = "TRex"), "meta_table.rds")

  # Check if meta table exists
  meta_table_exists <- file.exists(meta_table_path)

  if (meta_table_exists) {
    message("Loading existing meta table...")
    data_list <- readRDS(meta_table_path)

    # Check if entry already exists for this genome and bin width
    meta_key <- paste0("meta_", genome, "_", bin_width)
    if (!is.null(data_list[[meta_key]]) && !overwrite) {
      stop("Meta table for genome '", genome, "' and bin width '", bin_width, "' already exists. ",
           "Use overwrite = TRUE to replace it.")
    }
  } else {
    message("Creating new meta table...")
    data_list <- list()
  }

  # Scan directory for ChIP-seq files
  message("Scanning directory for ChIP-seq files...")
  chip_seq_files <- list.files(data_path, full.names = FALSE)

  if (length(chip_seq_files) == 0) {
    stop("No files found in the ChIP-seq data directory: ", data_path)
  }

  # Extract TR labels from filenames
  # Expected format: TR_label_*.txt or similar
  tr_labels <- sapply(strsplit(chip_seq_files, "_", fixed = TRUE), function(x) {
    if (length(x) >= 1) return(x[1]) else return(NA)
  })

  # Remove any NA values
  valid_indices <- !is.na(tr_labels)
  chip_seq_files <- chip_seq_files[valid_indices]
  tr_labels <- tr_labels[valid_indices]

  if (length(tr_labels) == 0) {
    stop("No valid TR labels could be extracted from filenames. ",
         "Expected format: TR_label_*.txt")
  }

  # Create meta table
  meta_table <- data.frame(
    TR = tr_labels,
    File_Path = file.path(data_path, chip_seq_files),
    stringsAsFactors = FALSE
  )

  # Validate file existence
  file_exists <- file.exists(meta_table$File_Path)
  if (!all(file_exists)) {
    missing_files <- meta_table$File_Path[!file_exists]
    warning("Some files do not exist: ", paste(basename(missing_files), collapse = ", "))
  }

  # Store in data list
  meta_key <- paste0("meta_", genome, "_", bin_width)
  data_list[[meta_key]] <- meta_table

  # Save meta table
  message("Saving meta table...")
  saveRDS(data_list, meta_table_path)

    # Summary
  message("ChIP-seq data successfully loaded!")
  message("  - Genome: ", genome)
  message("  - Bin width: ", bin_width)
  message("  - Number of TRs: ", nrow(meta_table))
  message("  - Files processed: ", sum(file_exists), "/", length(file_exists))
  message("")
  message("You can use check_loaded_chip_data() to view all loaded ChIP-seq data.")

  return(invisible(meta_table))
}

#' Check loaded ChIP-seq data
#'
#' @description
#' Retrieve and display information about loaded ChIP-seq data meta tables.
#' This function provides an overview of all available genome and bin width
#' combinations that have been loaded.
#'
#' @param genome Character string. Optional filter for specific genome ("hg38" or "mm10").
#'   If NULL, returns data for all genomes.
#' @param bin_width Integer. Optional filter for specific bin width (100, 200, 500, 1000).
#'   If NULL, returns data for all bin widths.
#'
#' @return A data frame containing information about loaded ChIP-seq data, or NULL
#'   if no data is loaded.
#' @export
check_loaded_chip_data <- function(genome = NULL, bin_width = NULL) {

  # Validate optional filters
  if (!is.null(genome)) {
    genome <- match.arg(genome, choices = c("hg38", "mm10"))
  }

  if (!is.null(bin_width)) {
    supported_bin_widths <- c(100, 200, 500, 1000)
    if (!bin_width %in% supported_bin_widths) {
      stop("bin_width must be one of: ", paste(supported_bin_widths, collapse = ", "))
    }
  }

  # Get meta table file path
  meta_table_path <- file.path(system.file(package = "TRex"), "meta_table.rds")

  if (!file.exists(meta_table_path)) {
    message("No ChIP-seq data has been loaded yet.")
    return(NULL)
  }

  # Load meta table
  data_list <- readRDS(meta_table_path)

  # Extract available keys
  available_keys <- names(data_list)
  meta_keys <- available_keys[grepl("^meta_", available_keys)]

  if (length(meta_keys) == 0) {
    message("No ChIP-seq meta tables found.")
    return(NULL)
  }

  # Parse keys to extract genome and bin width information
  key_info <- data.frame(
    Key = meta_keys,
    Genome = sapply(strsplit(meta_keys, "_"), function(x) x[2]),
    Bin_Width = as.integer(sapply(strsplit(meta_keys, "_"), function(x) x[3])),
    stringsAsFactors = FALSE
  )

  # Apply filters if specified
  if (!is.null(genome)) {
    key_info <- key_info[key_info$Genome == genome, , drop = FALSE]
  }

  if (!is.null(bin_width)) {
    key_info <- key_info[key_info$Bin_Width == bin_width, , drop = FALSE]
  }

  if (nrow(key_info) == 0) {
    message("No ChIP-seq data found matching the specified criteria.")
    return(NULL)
  }

  # Get detailed information for each key
  result_list <- list()
  for (i in seq_len(nrow(key_info))) {
    key <- key_info$Key[i]
    meta_table <- data_list[[key]]

    if (!is.null(meta_table) && nrow(meta_table) > 0) {
      # Check file existence
      file_exists <- file.exists(meta_table$File_Path)

      result_list[[i]] <- data.frame(
        Genome = key_info$Genome[i],
        Bin_Width = key_info$Bin_Width[i],
        Num_TRs = nrow(meta_table),
        Files_Available = sum(file_exists),
        Files_Total = length(file_exists),
        Data_Path = dirname(meta_table$File_Path[1]),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(result_list) == 0) {
    message("No valid ChIP-seq data found.")
    return(NULL)
  }

  # Combine results
  result <- do.call(rbind, result_list)
  rownames(result) <- NULL

  return(result)
}


#' display_tables
#' @description To show the ranking table by inspecting the results of variational inference.
#' @param file_path path to the saved TRex variational inference results.
#' @param output_path path to save the rank table.
#' @param burnin number of samples used for burn-in. If not specify, TRex will use the half of the iterations as burn in.
#' @return a data.frame object contains TR names, theta_i, TRex scores for each TR.
#' @export
display_tables<-function(file_path, output_path){
  dat<-readRDS(file_path)

  TR_names <- dat[["TR_names"]]
  TR_index<-!duplicated(TR_names)

  theta_i <- dat$theta_i[TR_index]
  TR_names <- TR_names[TR_index]

  TRex_score<-pnorm(theta_i)

  data("TR_meta_table")

  results_theta_i=data.frame(TR=TR_names,
                             Theta_i=theta_i,
                             TRex_score=TRex_score,
                             Is.TF = TR_meta_table$is.TF[match(TR_names,TR_meta_table$TR)],
                             Type = TR_meta_table$Type[match(TR_names,TR_meta_table$TR)]
                             )

  results_theta_i=results_theta_i[order(-results_theta_i$Theta_i),]
  results_theta_i$Rank=rank(-results_theta_i$Theta_i)
  row.names(results_theta_i) <- NULL

  write.csv(results_theta_i,paste0(tools::file_path_as_absolute(output_path),"/",
                                   tools::file_path_sans_ext(basename(file_path)),"_rank_table.csv"))

  return(results_theta_i)
}



#' RankAggre_MRRF: Aggregate Multiple TR Profiles Using MRRF
#'
#' @description
#' Computes an aggregated transcription regulator (TR) profile from a collection
#' of individual TR ranking tables using the Mean Reciprocal Rank with Smoothing
#' (MRRF) method. This produces a high-level summary of functional TRs by
#' integrating evidence across many datasets.
#'
#' @param profiles_dir Character string specifying the directory containing
#'   individual TR profile files. These files must follow the naming convention
#'   `*_rank_table.csv`.
#' @param output_dir Optional character string specifying the directory where the
#'   aggregated MRRF table will be saved. If `NULL` (default), the aggregated
#'   results are returned as a data.frame instead of being written to disk.
#' @param k Numeric smoothing parameter for MRRF calculation. Default is `50`.
#'
#' @return
#' If `output_dir` is `NULL`, returns a data.frame containing:
#' \itemize{
#'   \item \code{TR}: TR names
#'   \item \code{MRRF}: aggregated MRRF scores
#'   \item \code{Rank}: aggregated rank (1 = highest)
#'   \item \code{Is.TF}: logical indicator of whether the TR is a TF
#'   \item \code{Type}: functional category of the TR
#' }
#' If `output_dir` is not `NULL`, the function writes `MRRF_Table.csv` to the
#' specified directory and returns `NULL`.
#'
#' @export
RankAggre_MRRF <- function(profiles_dir, output_dir=NULL, k = 50) {

  # Get all TR profile files
  profile_files <- list.files(
    profiles_dir,
    pattern = "_rank_table\\.csv$",
    full.names = TRUE
  )

  if (length(profile_files) == 0) {
    stop("No TR profile files found in the provided directory.")
  }

  n_files <- length(profile_files)

  # ---- Read the first file to initialize structure ----
  first_profile <- read.csv(profile_files[1])
  tr_names <- sort(first_profile$TR)
  n_tr <- length(tr_names)

  rank_mat <- matrix(NA_real_, nrow = n_files, ncol = n_tr)
  colnames(rank_mat) <- tr_names

  # ---- Fill ranking matrix ----
  for (i in seq_along(profile_files)) {
    tr_df <- read.csv(profile_files[i])

    # Match TR order
    matched_ranks <- tr_df$Rank[match(tr_names, tr_df$TR)]
    rank_mat[i, ] <- matched_ranks
  }

  data("TR_meta_table")

  # ---- Compute MRRF ----
  output_table <- data.frame(
    TR   = tr_names,
    MRRF = colMeans(1 / (rank_mat + k), na.rm = TRUE),
    Rank = seq_len(n_tr),
    Is.TF = TR_meta_table$is.TF[match(tr_names,TR_meta_table$TR)],
    Type = TR_meta_table$Type[match(tr_names,TR_meta_table$TR)]
  )

  # ---- Sort by MRRF ----
  output_table <- output_table[order(output_table$MRRF, decreasing = TRUE), ]
  rownames(output_table) <- NULL

  if(!is.null(output_dir)){
    write.csv(output_table,paste0(output_dir,"MRRF_Table.csv"))
    return(NULL)
  }else{
    return(output_table)
  }
}

