#' TRex: variational Bayesian Transcriptional Regulators explorer
#'
#' @description
#' Performs variational Bayesian inference on ChIP-seq data by comparing input regions
#' with pre-compiled reference data and generating comprehensive analysis results.
#'
#' @param file Character. Path to the input file containing genomic regions.
#' @param output_path Character. Directory path where results will be saved.
#' @param filter_path Character, optional. Path to filter regions file. Default: NULL.
#' @param show Logical. Whether to display results tables. Default: TRUE.
#' @param plot.bar Logical. Whether to generate bar plots (currently disabled). Default: TRUE.
#' @param format Character, optional. Input file format. Default: NULL.
#' @param N Integer. Number of iterations for variational inference. Default: 5000.
#' @param bin_width Integer. Width of genomic bins in base pairs. Default: 1000.
#' @param genome Character vector. Supported genomes: "hg38", "mm10". Default: c("hg38", "mm10").
#'
#' @return Invisible NULL. Results are saved to the specified output path.
#'
#' @examples
#' \dontrun{
#' TRex(
#'   file = "peaks.bed",
#'   output_path = "./results",
#'   N = 3000,
#'   bin_width = 500
#' )
#' }
#'
#' @export
TRex <- function(
  file,
  output_path,
  filter_path = NULL,
  show = TRUE,
  plot.bar = TRUE,
  format = NULL,
  N = 5000,
  bin_width = 1000,
  genome = c("hg38", "mm10")
) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  # Validate required parameters
  if (missing(file) || missing(output_path)) {
    stop("Both 'file' and 'output_path' are required parameters")
  }

  if (!file.exists(file)) {
    stop("Input file does not exist: ", file)
  }

  # Validate genome parameter
  genome <- match.arg(genome, several.ok = TRUE)

  # Validate numeric parameters
  if (!is.numeric(N) || N <= 0) {
    stop("N must be a positive integer")
  }

  if (!is.numeric(bin_width) || bin_width <= 0) {
    stop("bin_width must be a positive integer")
  }

  # ============================================================================
  # INITIALIZATION
  # ============================================================================

  cat("Starting TRex Analysis\n")
  cat("=", rep("=", 50), "\n", sep = "")

  # Convert to absolute path
  output_path <- R.utils::getAbsolutePath(output_path)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
    cat("Created output directory:", output_path, "\n")
  }

  # ============================================================================
  # STEP 1: LOAD AND MAP PEAKS TO BINS
  # ============================================================================

  cat("\n Step 1: Loading and mapping peaks to bins...\n")

  peak_inds <- import_input_regions(
    file = file,
    format = format,
    bin_width = bin_width,
    genome = genome
  )

  # Apply filtering if filter_path is provided
  if (!is.null(filter_path)) {
    cat("Applying region filtering...\n")
    filter_inds <- import_input_regions(
      filter_path,
      format = format,
      bin_width = bin_width,
      genome = genome
    )
    cat("  - Original peaks: ", length(peak_inds), " regions\n")
    cat("  - Filter regions: ", length(filter_inds), " regions\n")
    filtered_peak_inds <- filter_peaks(peak_inds, filter_inds)
    cat("  - Filtered peaks: ", length(filtered_peak_inds), " regions\n")
    cat("  - Regions removed: ", length(peak_inds) - length(filtered_peak_inds), " regions\n")
    cat("Filtering completed\n")
  } else {
    filtered_peak_inds <- peak_inds
    cat("No filtering applied - using all ", length(peak_inds), " regions\n")
  }

  cat("Peak mapping completed\n")

  # ============================================================================
  # STEP 2: ALIGNMENT WITH REFERENCE DATA
  # ============================================================================

  cat("\n Step 2: Comparing with reference ChIP-seq data...\n")
  cat("   Bin width:", bin_width, "bps\n")

  alignment_results <- alignment_wrapper(
    filtered_peak_inds,
    bin_width = bin_width,
    genome = genome
  )

  cat("Alignment completed\n")

  # ============================================================================
  # STEP 3: PREPARE DATA FOR INFERENCE
  # ============================================================================

  cat("\n Step 3: Preparing data for inference...\n")

  # Extract alignment results
  xct <- alignment_results$GOOD
  nct <- alignment_results$TOTAL
  tr_labels <- as.numeric(factor(alignment_results$TR))

  cat("Data preparation completed\n")

  # ============================================================================
  # STEP 4: VARIATIONAL BAYESIAN INFERENCE
  # ============================================================================

  cat("\n Step 4: Running variational Bayesian inference...\n")
  cat("   Iterations:", N, "\n")

  VI_chain_results <- Main_Function(N, xct, nct, tr_labels)

  # Add transcription factor names to results
  VI_chain_results[["TR_names"]] <- alignment_results$TR

  cat("Inference completed\n")

  # ============================================================================
  # STEP 5: SAVE RESULTS
  # ============================================================================

  cat("\n Step 5: Saving results...\n")

  # Generate output filename
  base_name <- tools::file_path_sans_ext(basename(file))
  file_name <- file.path(output_path, paste0(base_name, ".rds"))

  # Save results
  saveRDS(VI_chain_results, file_name)

  cat("Results saved to:", file_name, "\n")

  # ============================================================================
  # STEP 6: DISPLAY RESULTS (OPTIONAL)
  # ============================================================================

  if (show) {
    cat("\n Step 6: Displaying results tables...\n")
    display_tables(file_path = file_name, output_path = output_path)
    cat("Tables displayed\n")
  }

  # ============================================================================
  # STEP 7: GENERATE PLOTS (OPTIONAL - CURRENTLY DISABLED)
  # ============================================================================

  if (plot.bar) {
    cat("\n Step 7: Plot generation (currently disabled)...\n")
    # TODO: Implement plot generation when burnin parameter is available
    # rank_plot(file_path = file_name, output_path = output_path, burnin = burnin)
    cat("Plot generation is currently disabled\n")
  }

  # ============================================================================
  # COMPLETION
  # ============================================================================

  cat("\n", "=", rep("=", 50), "\n", sep = "")
  cat("TRex Analysis Completed Successfully!\n")
  cat("Results location:", output_path, "\n")
  cat("Main TR profile:", paste0(output_path, "/", base_name, "_rank_table.csv"),"\n")
  cat("VI file:", file_name, "\n")
  cat("=", rep("=", 50), "\n", sep = "")

  # Return invisibly
  invisible(NULL)
}







