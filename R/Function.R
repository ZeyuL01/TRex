##functions to transform regions to bins.

#' Count the 'good' and 'informative' cases by comparing input with the reference database.
#'
#' @param input_vec A numeric vector containing indices of transformed regions
#'   obtained by applying \code{import_input_regions}.
#' @param bin_width Width of bin in base pairs.
#' @param genome The genome of TR ChIP-seq data, either "hg38" or "mm10".
#'
#' @return A data frame with three columns: TR labels, number of 'good' cases,
#'   and number of 'total' informative cases.
#' @export
alignment_wrapper <- function(input_vec, bin_width, genome = c("hg38", "mm10")) {

  # Validate genome input
  genome <- match.arg(genome)

  # Validate input vector
  if (!is.numeric(input_vec) || length(input_vec) == 0) {
    stop("input_vec must be a non-empty numeric vector")
  }

  # Load the meta table
  message("Loading meta table...")
  meta_table_path <- system.file("meta_table.rds", package = "vBIT")

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
         "You may follow the tutorial on: https://github.com/ZeyuL01/BIT")
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
                     start <- rtracklayer::start(peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab])
                     end <- rtracklayer::end(peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab])
                     (end + start) %/% (2 * bin_width) + 1
                   },
                   narrowPeak = {
                     start <- rtracklayer::start(peak_dat[GenomicRanges::seqnames(peak_dat)==chr_lab])
                     (start+peak_dat[GenomicRanges::seqnames(peak_dat)==chr_lab]$peak)%/% bin_width + 1
                   },
                   broadPeak = {
                     start <- rtracklayer::start(peak_dat[GenomicRanges::seqnames(peak_dat)==chr_lab])
                     (start+peak_dat[GenomicRanges::seqnames(peak_dat)==chr_lab]$abs_summit)%/% bin_width + 1
                   },
                   bigNarrowPeak = peak_dat[GenomicRanges::seqnames(peak_dat) == chr_lab]$abs_summit %/% bin_width + 1,
                   csv = {
                     start <- peak_dat$Start[peak_dat$Chrom == chr_lab]
                     end <- peak_dat$End[peak_dat$Chrom == chr_lab]
                     (end + start) %/% (2 * bin_width) + 1
                   }
    )

    inds <- inds[inds <= chr_win_num]
    inds <- inds[!duplicated(inds)]
    bin_inds <- c(bin_inds, (chr_windows_cs[i] + inds))
  }

  return(sort(bin_inds))
}



##functions to load chip-seq datasets
##just need to run once.

#' Load pre-compiled ChIP-seq data
#'
#' @description
#' Load and organize pre-compiled ChIP-seq data for vBIT analysis. This function
#' scans the specified directory for ChIP-seq files, extracts transcription regulator
#' (TR) labels, and creates a meta table for efficient data access. Please follow
#' the tutorial on: https://github.com/ZeyuL01/BIT.
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
  meta_table_path <- file.path(system.file(package = "vBIT"), "meta_table.rds")

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
  meta_table_path <- file.path(system.file(package = "vBIT"), "meta_table.rds")

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
#' @description To show the ranking table by inspecting the results of Gibbs sampler.
#' @param file_path path to the saved BIT Gibbs sampling results.
#' @param output_path path to save the rank table.
#' @param burnin number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.
#' @return a data.frame object contains TR names, theta_i, BIT scores and 95 CIs for each TR.
#' @export
display_tables<-function(file_path, output_path){
  dat<-readRDS(file_path)

  TR_names <- dat[["TR_names"]]
  TR_index<-!duplicated(TR_names)

  theta_i <- dat$theta_i[TR_index]
  TR_names <- TR_names[TR_index]

  BIT_score<-logistic(theta_i)
  results_theta_i=data.frame(TR=TR_names,Theta_i=theta_i,BIT_score=BIT_score)
  results_theta_i=results_theta_i[order(-results_theta_i$Theta_i),]
  results_theta_i$Rank=rank(-results_theta_i$Theta_i)
  row.names(results_theta_i) <- NULL

  write.csv(results_theta_i,paste0(tools::file_path_as_absolute(output_path),"/",
                                   tools::file_path_sans_ext(basename(file_path)),"_rank_table.csv"))

  return(results_theta_i)
}

#' rank_plot
#' @description To draw a barplot for the top n TRs.
#' @param file_path path to the saved BIT Gibbs sampling results.
#' @param output_path path to save the barplot.
#' @param burnin number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.
#' @param n top n TRs will show in the barplot, default: 10.
#' @param colors colors for each bar, default "NPG" for n<=10, has to be manually specified if n>10.
#' @param main main title for the barplot, default: NULL.
#' @param xlab x axis label, default: BIT score.
#' @param ylab y axis label, default: TR symbols.
#' @return a data.frame object contains TR names, theta_i, BIT scores and 95 CIs for each TR.
#' @export
rank_plot<-function(file_path=NULL, output_path, burnin=NULL, n=10, colors="NPG", main=NULL, xlab="BIT score", ylab="TR symbols"){
  dat<-readRDS(file_path)

  TR_names <- dat[["TR_names"]]
  theta_i_mat<-dat$theta_i

  if(is.null(burnin)){
    burnin=dim(theta_i_mat)[2]%/%2
  }

  theta_i_mat<-theta_i_mat[!duplicated(TR_names),(dim(theta_i_mat)[2]-burnin):dim(theta_i_mat)[2]]
  TR_names<-TR_names[!duplicated(TR_names)]
  tr_results_i<-rowMeans(theta_i_mat)
  BIT_score<-logistic(tr_results_i)
  CI_intervals<-t(apply(theta_i_mat, 1, function(x) quantile(x, probs = c(0.025, 0.975))))

  results_theta_i=data.frame(TR=TR_names,Theta_i=tr_results_i,lower=CI_intervals[,1],upper=CI_intervals[,2],
                             BIT_score=BIT_score,BIT_score_lower=logistic(CI_intervals[,1]),BIT_score_upper=logistic(CI_intervals[,2]))
  results_theta_i=results_theta_i[order(-results_theta_i$Theta_i),]
  results_theta_i$Rank=rank(-results_theta_i$Theta_i)
  row.names(results_theta_i) <- NULL

  output_file_path_complete<-paste0(tools::file_path_as_absolute(output_path),"/",
                                    tools::file_path_sans_ext(basename(file_path)),".pdf")

  TR_names_used<-results_theta_i[n:1,"TR"]
  BIT_score_used<-results_theta_i[n:1,"BIT_score"]
  BIT_score_lower_used<-results_theta_i[n:1,"BIT_score_lower"]
  BIT_score_upper_used<-results_theta_i[n:1,"BIT_score_upper"]

  if(colors=="NPG"){
    colors=ggsci::pal_npg()(n)
  }

  pdf(output_file_path_complete)
  par(mfrow=c(1,1),oma = c(1,1,1,1) + 0.1,mar = c(4,6.5,2,1) + 0.1,mgp=c(3, 0.5, 0),font.axis=2)
  bp<-barplot(BIT_score_used,horiz=TRUE,yaxt="n",xlim=c(0,max(BIT_score_upper_used)*1.05),cex.axis=1.5,col=colors,tcl=-0.2,main=main,cex.main=1.3,font.axis=1)
  text(BIT_score_used-max(BIT_score_used)/10,bp,round(BIT_score_used,3),font=1,cex=1.2)
  points(BIT_score_used,bp,pch=16)
  segments(BIT_score_used,bp,BIT_score_upper_used,bp,lwd=2)
  segwidth<-(bp[2]-bp[1])/4
  segments(BIT_score_upper_used,bp+(bp[2]-bp[1])/4,BIT_score_upper_used,bp-(bp[2]-bp[1])/4,lwd=2)

  axis(side=2, las=1, at=bp, labels=TR_names_used,tcl=-0.2,cex.axis=1.2,font.axis=1)
  title(xlab=xlab,line = 3, cex.lab=1.4,font.lab=2)
  title(ylab=ylab,line = 5.5, cex.lab=1.4,font.lab=2)
  title(main=main,cex.main=1.4,font.main=2)
  box()
  dev.off()

  return()
}

#' compare_scatter_plot
#' @description To draw a barplot for the top n TRs.
#' @param file1_path path to the saved BIT Gibbs sampling results of input 1.
#' @param file2_path path to the saved BIT Gibbs sampling results of input 2.
#' @param output_path path to save the barplot.
#' @param burnin number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.
#' @return NULL
#' @export
compare_scatter_plot<-function(file1_path, file2_path, output_path, burnin=NULL){
  dat1<-readRDS(file1_path)
  dat2<-readRDS(file2_path)

  TR_names_dat1 <- dat1[["TR_names"]]
  dat1_theta_i_mat<-dat1$theta_i

  TR_names_dat2 <- dat2[["TR_names"]]
  dat2_theta_i_mat<-dat2$theta_i

  if(is.null(burnin)){
    burnin=dim(dat1_theta_i_mat)[2]%/%2
  }

  output_file_path_complete<-paste0(tools::file_path_as_absolute(output_path),"/",
                                    tools::file_path_sans_ext(basename(file1_path)),"_",
                                    tools::file_path_sans_ext(basename(file2_path)),"_compare.pdf")

  theta_i_mat_dat1<-dat1_theta_i_mat[!duplicated(TR_names_dat1),(dim(dat1_theta_i_mat)[2]-burnin):dim(dat1_theta_i_mat)[2]]
  TR_names_dat1<-TR_names_dat1[!duplicated(TR_names_dat1)]
  tr_results_i_dat1<-rowMeans(theta_i_mat_dat1)
  BIT_score<-logistic(tr_results_i_dat1)
  CI_intervals<-t(apply(theta_i_mat_dat1, 1, function(x) quantile(x, probs = c(0.025, 0.975))))

  results_theta_i_dat1=data.frame(TR=TR_names_dat1,Theta_i=tr_results_i_dat1,BIT_score=BIT_score,BIT_score_lower=CI_intervals[,1],BIT_score_upper=CI_intervals[,2])
  results_theta_i_dat1=results_theta_i_dat1[order(-results_theta_i_dat1$Theta_i),]
  results_theta_i_dat1$Rank=rank(-results_theta_i_dat1$Theta_i)
  row.names(results_theta_i_dat1) <- NULL

  theta_i_mat_dat2<-dat2_theta_i_mat[!duplicated(TR_names_dat2),(dim(dat2_theta_i_mat)[2]-burnin):dim(dat2_theta_i_mat)[2]]
  TR_names_dat2<-TR_names_dat2[!duplicated(TR_names_dat2)]
  tr_results_i_dat2<-rowMeans(theta_i_mat_dat2)
  BIT_score<-logistic(tr_results_i_dat2)
  CI_intervals<-t(apply(theta_i_mat_dat2, 1, function(x) quantile(x, probs = c(0.025, 0.975))))

  results_theta_i_dat2=data.frame(TR=TR_names_dat2,Theta_i=tr_results_i_dat2,BIT_score=BIT_score,BIT_score_lower=CI_intervals[,1],BIT_score_upper=CI_intervals[,2])
  results_theta_i_dat2=results_theta_i_dat2[order(-results_theta_i_dat2$Theta_i),]
  results_theta_i_dat2$Rank=rank(-results_theta_i_dat2$Theta_i)
  row.names(results_theta_i_dat2) <- NULL

  file1_values<-results_theta_i_dat1[,"BIT_score"]
  file2_values<-results_theta_i_dat2[match(results_theta_i_dat1[,"TR"],results_theta_i_dat2[,"TR"]),"BIT_score"]

  aligned_tables<-data.frame(TR=results_theta_i_dat1[,"TR"],file1_values=file1_values,file2_values=file2_values)

  TR_names1<-results_theta_i_dat1[1:10,"TR"]
  TR_names2<-results_theta_i_dat2[1:10,"TR"]
  TR_union<-union(TR_names1,TR_names2)

  file1_coords<-aligned_tables[which(aligned_tables[,"TR"]%in%TR_union),"file1_values"]
  file2_coords<-aligned_tables[which(aligned_tables[,"TR"]%in%TR_union),"file2_values"]
  TR_plot_names<-aligned_tables[which(aligned_tables[,"TR"]%in%TR_union),"TR"]

  scale<-c()
  scale1_vec<-(file1_values/max(file1_values))
  scale2_vec<-(file2_values/max(file2_values))
  for(i in 1:length(file1_values)){
    scale<-c(scale,max(c(scale1_vec[i],scale2_vec[i])))
  }

  size<-1-median(scale)+scale

  colors<-rgb(file1_values/max(file1_values),0.4,file2_values/max(file2_values),scale)

  pdf(output_file_path_complete)
  par(mfrow=c(1,1),oma = c(1,1,1,1) + 0.1,mar = c(4,4,2,1) + 0.1,mgp=c(3, 0.5, 0),font.axis=2)
  plot(file1_values,file2_values,type="p",pch=19,col=colors,cex=size, xlab="Region set 1 BIT score",ylab="Region set 2 BIT score",font.lab=2,cex.lab=1.4)
  basicPlotteR::addTextLabels(file1_coords, file2_coords, TR_plot_names, cex.label=1.4,lwd=2)
  box()
  dev.off()

  return()
}


#' Add non-overlapping text labels to plot
#' The function addTextLabels is copied and used under GPL 3.0 license,
#' we thank the contribution by original author: Joseph Crispell https://github.com/JosephCrispell/basicPlotteR
#' @param xCoords A vector containing the X coordinates for labels
#' @param yCoords A vector containing the Y coordinates for labels
#' @param labels A vector containing the labels to be plotted
#' @param cex.label A number to scale the size of the plotted labels. Defaults to 1
#' @param col.label The colour of the plotted labels. Defaults to "red". Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param col.line The colour of the line to plot from relocated labels to original location. Defaults to "black". Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param col.background An optional colour for a background polygon plotted behind labels. Defaults to NULL - won't be plotted. Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param lty A number detailing the type of line to plot from relocated labels to original location. 0: blank, 1: solid, 2: dashed, 3: dotted, 4: dotdash, 5: longdash, and 6: twodash. Defaults to 1. Multiple line types can be provided. If more options than labels provided types will be recycled.
#' @param lwd A number to scale the size of line from relocated labels to original location. Defaults to 1. Multiple line widths can be provided. If more options than labels provided widths will be recycled.
#' @param border The colour of the border to be plotted around the polygon. Defaults to NA - won't be plotted. Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param avoidPoints A logical variable indicating whether labels shouldn't be plotted on top of points. Defaults to TRUE
#' @param keepLabelsInside A logical variable indicating whether the labels shouldn't be plotted outside of plotting region. Defaults to TRUE
#' @param cex.pt A number used to scale the points plotted on the graph that labels are to be added to. Defaults to 1
#' @keywords text label plot
#' @examples
#' # Create some random points
#' n <- 50
#' coords <- data.frame(X=runif(n), Y=runif(n), Name="Test Label")
#'
#' # Plot them without labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#'
#' # With potentially overlapping labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#' text(coords$X, coords$Y, labels=coords$Name, xpd=TRUE)
#'
#' # Plot them with non-overlapping labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#' addTextLabels(coords$X, coords$Y, coords$Name, cex.label=1, col.label="black")
#'
#' # Plot them with non-overlapping labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#' addTextLabels(coords$X, coords$Y, coords$Name, cex.label=1, col.background=rgb(0,0,0, 0.75), col.label="white")
addTextLabels <- function(xCoords, yCoords, labels, cex.label=1, col.label="red", col.line="black", col.background=NULL,
                          lty=1, lwd=1, border=NA, avoidPoints=TRUE, keepLabelsInside=TRUE, cex.pt=1){

  #######################################################
  # Check that the input data are in the correct format #
  #######################################################

  # Are each of coordinate vectors the same length?
  if(length(xCoords) != length(yCoords)){
    stop("addTextLabels() The vectors containing the X and Y coodinates must be the same length.")
  }
  if(length(xCoords) != length(labels)){
    stop("addTextLabels() The vector of labels must be the same length as the coordinate vectors.")
  }

  #######################
  # Get the axis limits #
  #######################

  # Get the axis limits
  axisLimits <- graphics::par("usr")

  ############################
  # Check for NA coordinates #
  ############################

  # Check if any NA coordinates present
  indicesOfNAs <- which(is.na(xCoords) | is.na(yCoords))
  if(length(indicesOfNAs) > 0){

    # Send warning
    warning("NA values present in coordinates provided. These are ignored.")

    # Check for each of the parameters that can have multiple parameters
    if(length(col.line) == length(xCoords)){
      col.line = col.line[-indicesOfNAs]
    }
    if(length(col.background) == length(xCoords)){
      col.background = col.background[-indicesOfNAs]
    }
    if(length(lty) == length(xCoords)){
      lty = lty[-indicesOfNAs]
    }
    if(length(lwd) == length(xCoords)){
      lwd = lwd[-indicesOfNAs]
    }
    if(length(border) == length(xCoords)){
      border = border[-indicesOfNAs]
    }

    # Remove the NA coordinates
    xCoords <- xCoords[-indicesOfNAs]
    yCoords <- yCoords[-indicesOfNAs]

    # Remove the respective labels
    labels <- labels[-indicesOfNAs]
  }

  ############################
  # Check if axes are logged #
  ############################

  # Check X axis
  xAxisLogged <- FALSE
  if(graphics::par("xlog")){

    # Note that X axis was logged
    xAxisLogged <- TRUE

    # Log the X coordinates
    xCoords <- log10(xCoords)

    # Reset the X axis logged flag - fools points and polygon commands below
    graphics::par(xlog=FALSE)
  }

  # Check Y axis
  yAxisLogged <- FALSE
  if(graphics::par("ylog")){

    # Note that Y axis was logged
    yAxisLogged <- TRUE

    # Log the Y coordinates
    yCoords <- log10(yCoords)

    # Reset the Y axis logged flag - fools points and polygon commands below
    graphics::par(ylog=FALSE)
  }

  ###############################
  # Store the point information #
  ###############################

  # Store the input coordinates and labels
  pointInfo <- list("X"=xCoords, "Y"=yCoords, "Labels"=labels, "N"=length(xCoords), "cex"=cex.pt)

  # Set the amount to pad onto height and width
  heightPad <- 0.5
  widthPad <- 0.04
  if(is.null(col.background)){
    heightPad <- 0
    widthPad <- 0
  }

  # Calculate the label heights and widths
  pointInfo <- calculateLabelHeightsAndWidths(pointInfo=pointInfo, cex=cex.label,
                                              heightPad=heightPad, widthPad=widthPad)

  ###########################################
  # Produce a list of alternative locations #
  ###########################################

  # Generate the alternative locations
  alternativeLocations <- generateAlternativeLocations(axisLimits)

  # Calculate the distance between the actual and alternative points - rescale X axis remove axis range bias
  distances <- euclideanDistancesWithRescaledXAxis(pointInfo, alternativeLocations, axisLimits)

  ###############################################################
  # Create a list to store the information about plotted points #
  ###############################################################

  # Initialise the list to store the information about plotted labels
  plottedLabelInfo <- list("X"=c(), "Y"=c(), "Height"=c(), "Width"=c(), "N"=0)

  ##############################################################
  # Add labels to plot assigning new locations where necessary #
  ##############################################################

  # Plot the point label
  for(i in seq_len(pointInfo$N)){

    # Set the colours for plotting the label - allows multiple colours and cycling through colours
    labelColour <- setOption(options=col.label, index=i)
    backgroundColour <- setOption(options=col.background, index=i)
    borderColour <- setOption(options=border, index=i)

    # Set the line characteristics
    lineColour <- setOption(options=col.line, index=i)
    lineType <- setOption(options=lty, index=i)
    lineWidth <- setOption(options=lwd, index=i)

    # Get the information for the current point
    x <- pointInfo$X[i]
    y <- pointInfo$Y[i]
    label <- pointInfo$Labels[i]
    height <- pointInfo$Heights[i]
    width <- pointInfo$Widths[i]

    # Get a new location
    newLocationIndex <- chooseNewLocation(pointInfo, i, alternativeLocations, distances, plottedLabelInfo, axisLimits, keepLabelsInside)

    # Is the current point too close to others?
    if(alternativeLocations$N != 0 && newLocationIndex != -1 &&
       (avoidPoints == TRUE || tooClose(x, y, height, width, plottedLabelInfo) || outsidePlot(x, y, height, width, axisLimits))){

      # Get the coordinates for the chosen alternate location
      altX <- alternativeLocations$X[newLocationIndex]
      altY <- alternativeLocations$Y[newLocationIndex]

      # Add line back to previous location
      addLineBackToOriginalLocation(altX=altX, altY=altY, x=x, y=y, label=label,
                                    cex=cex.label, col=lineColour, lty=lineType, lwd=lineWidth, heightPad=heightPad,
                                    widthPad=widthPad)

      # Add label
      addLabel(x=altX, y=altY, label=label,
               cex=cex.label, col=labelColour, bg=backgroundColour, border=borderColour, heightPad=heightPad, widthPad=widthPad)

      # Append the plotted label information
      plottedLabelInfo <- addPlottedLabel(x=altX, y=altY, height=height, width=width,
                                          plottedLabelInfo=plottedLabelInfo)

      # Remove the alternative plotting location used
      alternativeLocations$X <- alternativeLocations$X[-newLocationIndex]
      alternativeLocations$Y <- alternativeLocations$Y[-newLocationIndex]
      alternativeLocations$N <- alternativeLocations$N - 1
      distances <- distances[, -newLocationIndex]

    }else{

      # Add label
      addLabel(x=x, y=y, label=label,
               cex=cex.label, col=labelColour, bg=backgroundColour, border=borderColour,
               heightPad=heightPad, widthPad=widthPad)

      # Append the plotted label information
      plottedLabelInfo <- addPlottedLabel(x=x, y=y, height=height, width=width,
                                          plottedLabelInfo=plottedLabelInfo)
    }
  }

  #####################################################################################
  # Return axes logged flags to original state - for if person makes any future plots #
  #####################################################################################

  graphics::par(xlog=xAxisLogged)
  graphics::par(ylog=yAxisLogged)

}

#' function used to compute the logistic transformation
#' @param x vectorized object
#'
#' @return None
#'
logistic <- function(x) {
  1 / (1 + exp(-x))
}


#' function used to filter converted peaks based on distal elements
#' @param input_vec vectorized object
#'
#' @return index of peak index after the filtering
#'
filter_peaks <- function(input_vec, filter_vec){

  filtered_input_vec<-intersect(input_vec,filter_vec)

  return(filtered_input_vec)
}


#' Plot Gene Regulatory Network
#'
#' @description
#' This function reads peak data and transcription regulator (TR) information
#' to construct and visualize a gene regulatory network, integrating data
#' from STRINGdb and TFLink databases.
#'
#' @param input_file_path Character string. Path to the peak file (e.g., BED format).
#' @param input_table_path Character string. Path to the CSV file containing TRs,
#'   expected to have a column named 'TR'.
#' @param thres_TRs Integer. The number of top TRs to include from `input_table_path`. Default is 10.
#' @param genome Character string. The reference genome, either "hg38" or "mm10".
#' @param link Character string. Method to select target genes:
#'   "pathway" - selects genes from the top enriched GO BP pathway.
#'   "number" - selects genes based on peak frequency (top proportion).
#'   Default is "number".
#' @param tss_region Numeric vector of length 2. Region around TSS for peak annotation (e.g., c(-3000, 3000)). Default is c(-3000, 3000).
#' @param number_proportion Numeric. Proportion of top genes to select when `link = "number"`. Default is 0.005 (top 0.5\%).
#' @param string_score_threshold Integer. Minimum interaction score for STRINGdb connections (0-1000). Default is 400.
#' @param node_size Numeric. Base size for nodes in the plot. Default is 5.
#' @param label_size Numeric. Size for node labels in the plot. Default is 3.
#' @param layout_algorithm Character string. Layout algorithm for ggraph (e.g., 'fr', 'kk', 'nicely'). Default is 'fr' (Fruchterman-Reingold).
#' @param vbit_package_name Character string. The name of the package where TFLink data resides, needed for `system.file`. Default is "vBIT".
#'
#' @return A ggraph plot object representing the gene regulatory network.
#'
#' @importFrom ChIPseeker readPeakFile annotatePeak
#' @importFrom clusterProfiler enrichGO
#' @importFrom biomaRt useMart getBM
#' @importFrom STRINGdb STRINGdb
#' @importFrom igraph graph_from_data_frame V E
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text theme_graph create_layout
#' @importFrom dplyr filter select rename left_join bind_rows distinct arrange desc %>% slice_head
#' @importFrom readr read_csv cols
#' @importFrom data.table fread
#' @importFrom methods as
#' @importFrom stats quantile
#' @importFrom AnnotationDbi select
#' @importFrom utils head installed.packages
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage (requires valid input files and package installation)
#' peak_file <- "path/to/your/peaks.bed"
#' tr_file <- "path/to/your/tr_table.csv"
#'
#' # Create dummy files for demonstration if needed
#' # writeLines(c("chr1\t1000\t1500", "chr1\t2000\t2500"), peak_file)
#' # write.csv(data.frame(TR = c("TF1", "TF2", "TF3")), tr_file, row.names = FALSE)
#'
#' grn_plot <- plot_GRN(
#'   input_file_path = peak_file,
#'   input_table_path = tr_file,
#'   thres_TRs = 3,
#'   genome = "hg38", # or "mm10"
#'   link = "number"
#' )
#'
#' print(grn_plot)
#' }
plot_GRN <- function(input_file_path,
                     input_table_path,
                     thres_TRs = 10,
                     genome = c("hg38", "mm10"),
                     link = c("number", "pathway"), # Removed "score" as it wasn't implemented
                     tss_region = c(-3000, 3000),
                     number_proportion = 0.005,
                     string_score_threshold = 400,
                     node_size = 5,
                     label_size = 3,
                     layout_algorithm = 'nicely',
                     vbit_package_name = "vBIT") {

  # --- Input Validation and Argument Matching ---
  genome <- match.arg(genome)
  link <- match.arg(link)

  if (!file.exists(input_file_path)) {
    stop("Input peak file not found: ", input_file_path)
  }
  if (!file.exists(input_table_path)) {
    stop("Input TR table file not found: ", input_table_path)
  }
  if (! (vbit_package_name %in% rownames(installed.packages())) ) {
    warning("Package '", vbit_package_name, "' not found. Cannot load TFLink data from package path.")
    # Consider adding alternative TFLink loading mechanism or stopping execution
    # For now, we'll let it fail later if system.file returns empty string.
  }


  # --- Genome-Specific Setup ---
  message("Setting up resources for genome: ", genome)
  if (genome == "hg38") {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    anno_db <- "org.Hs.eg.db"
    species_id <- 9606
    ensembl_dataset <- "hsapiens_gene_ensembl"
    symbol_col <- "hgnc_symbol" # Column name for gene symbols in biomaRt result
    tflink_ss_file <- paste0(system.file(package = vbit_package_name), "/data/TFLink_Human_SS_v1.0.tsv")
    tflink_ls_file <- paste0(system.file(package = vbit_package_name), "/data/TFLink_Human_LS_v1.0.tsv")
  } else if (genome == "mm10") {
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene # Corrected TxDb for mm10
    anno_db <- "org.Mm.eg.db"
    species_id <- 10090
    ensembl_dataset <- "mmusculus_gene_ensembl"
    symbol_col <- "mgi_symbol" # Column name for gene symbols in biomaRt result
    tflink_ss_file <- paste0(system.file(package = vbit_package_name), "/data/TFLink_Mice_SS_v1.0.tsv")
    tflink_ls_file <- paste0(system.file(package = vbit_package_name), "/data/TFLink_Mice_LS_v1.0.tsv")
  } else {
    # This case is redundant due to match.arg but kept for clarity
    stop("Unsupported genome specified.")
  }

  # --- Load and Annotate Peaks ---
  message("Loading and annotating peaks...")
  peaks <- readPeakFile(input_file_path)
  peak_anno <- annotatePeak(peaks,
                            tssRegion = tss_region,
                            TxDb = txdb,
                            annoDb = anno_db)
  peak_anno_df <- as.data.frame(peak_anno)

  # --- Load TRs ---
  message("Loading Transcription Regulators (TRs)...")
  tr_table <- readr::read_csv(input_table_path, col_types = cols()) # More robust reading
  if (!"TR" %in% colnames(tr_table)) {
    stop("Input TR table must contain a column named 'TR'.")
  }
  # Ensure we don't select more TRs than available
  n_trs_available <- nrow(tr_table)
  if (thres_TRs > n_trs_available) {
    warning("Requested ", thres_TRs, " TRs, but only ", n_trs_available, " are available in the table. Using all available TRs.")
    thres_TRs <- n_trs_available
  }
  trs_selected <- tr_table$TR[1:thres_TRs]


  # --- Get Protein Coding Genes List ---
  message("Fetching protein-coding gene list from Ensembl...")
  ensembl <- useMart("ensembl", dataset = ensembl_dataset)
  protein_coding_genes <- getBM(
    attributes = c('ensembl_gene_id', symbol_col, 'gene_biotype'),
    filters = 'biotype',
    values = 'protein_coding',
    mart = ensembl
  )
  # Filter for genes that have a symbol
  protein_coding_genes <- protein_coding_genes %>%
    filter(!!sym(symbol_col) != "")

  # --- Filter Annotated Peaks for Protein Coding Genes ---
  # Ensure column names match ('SYMBOL' from ChIPseeker, symbol_col from biomaRt)
  if (!"SYMBOL" %in% names(peak_anno_df)) {
    stop("ChIPseeker annotation result does not contain 'SYMBOL' column.")
  }
  peak_anno_df_coding <- peak_anno_df %>%
    filter(SYMBOL %in% protein_coding_genes[[symbol_col]]) # Use dynamic symbol column name

  if (nrow(peak_anno_df_coding) == 0) {
    stop("No peaks were annotated to protein-coding genes.")
  }

  # --- Select Target Genes based on 'link' method ---
  message("Selecting target genes based on method: ", link)
  genes_selected <- character(0) # Initialize empty vector

  if (link == "pathway") {
    # Identify genes near promoters
    promoter_genes_entrez <- peak_anno_df_coding %>%
      filter(grepl("Promoter", annotation)) %>% # annotation is standard ChIPseeker output column
      pull(geneId) %>% # geneId is Entrez ID
      unique()

    if (length(promoter_genes_entrez) > 0) {
      # Perform GO Enrichment Analysis
      go_results <- enrichGO(gene = promoter_genes_entrez,
                             OrgDb = anno_db,
                             ont = "BP", # Biological Process
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2, # Relaxed q-value often needed
                             readable = TRUE)     # Converts Entrez IDs to Gene Symbols

      if (!is.null(go_results) && nrow(go_results@result) > 0) {
        # Get genes from the top significant pathway
        top_pathway_genes <- strsplit(go_results@result$geneID[1], "/")[[1]]
        genes_selected <- unique(top_pathway_genes)
        message("Selected ", length(genes_selected), " genes from top GO pathway: ", go_results@result$Description[1])
      } else {
        warning("GO enrichment analysis yielded no significant results. No pathway genes selected.")
      }
    } else {
      warning("No peaks found in promoter regions for GO analysis.")
    }

  } else if (link == "number") {
    # Calculate frequency of genes associated with peaks
    gene_freq <- table(peak_anno_df_coding$SYMBOL)

    # Filter for valid protein-coding gene symbols
    gene_freq <- gene_freq[names(gene_freq) %in% protein_coding_genes[[symbol_col]]]

    if (length(gene_freq) > 0) {
      # Determine cutoff based on quantile
      cutoff <- quantile(gene_freq, probs = (1 - number_proportion), na.rm = TRUE)
      # Select genes above or equal to the cutoff
      genes_selected <- names(gene_freq[gene_freq >= cutoff])
      message("Selected ", length(genes_selected), " genes based on peak frequency (top ", number_proportion * 100, "%)")
    } else {
      warning("No genes found after filtering for protein-coding symbols in frequency analysis.")
    }
  }

  # Combine TRs and selected target genes for network nodes
  all_selected_genes <- unique(c(trs_selected, genes_selected))

  if (length(all_selected_genes) < 2) {
    stop("Less than 2 genes/TRs selected. Cannot build a network.")
  }

  # --- Fetch Interactions from STRINGdb ---
  message("Fetching interactions from STRINGdb (threshold: ", string_score_threshold, ")..")
  string_db <- STRINGdb$new(version = "11.5", # Or latest version
                            species = species_id,
                            score_threshold = string_score_threshold,
                            input_directory = "") # Optional: path to local STRING data

  # Map gene symbols to STRING IDs
  mapped_genes <- string_db$map(data.frame(gene = all_selected_genes), "gene", removeUnmappedRows = TRUE)

  if (nrow(mapped_genes) > 0) {
    # Get interactions between the mapped genes
    interactions_string <- string_db$get_interactions(mapped_genes$STRING_id)

    # Map STRING IDs back to gene symbols
    interactions_string <- interactions_string %>%
      left_join(select(mapped_genes, STRING_id, gene_from = gene), by = c("from" = "STRING_id")) %>%
      left_join(select(mapped_genes, STRING_id, gene_to = gene), by = c("to" = "STRING_id")) %>%
      filter(!is.na(gene_from) & !is.na(gene_to)) %>% # Ensure both partners are in our list
      select(from = gene_from, to = gene_to, score = combined_score) %>%
      mutate(source = "STRING") %>%
      distinct() # Avoid duplicates if mapping wasn't 1:1
  } else {
    warning("None of the selected genes could be mapped by STRINGdb.")
    interactions_string <- data.frame(from=character(), to=character(), score=numeric(), source=character())
  }

  # --- Fetch Interactions from TFLink ---
  message("Fetching interactions from TFLink...")

  # Helper function to safely read TFLink CSV
  read_tflink <- function(file_path) {
    if (!file.exists(file_path)) {
      warning("TFLink file not found: ", file_path)
      return(data.frame(Name.TF = character(), Name.Target = character())) # Return empty frame
    }
    tryCatch({
      data.table::fread(file_path)
    }, error = function(e) {
      warning("Error reading TFLink file ", file_path, ": ", e$message)
      return(data.frame(Name.TF = character(), Name.Target = character())) # Return empty frame on error
    })
  }

  tflink_small <- read_tflink(tflink_ss_file)
  tflink_large <- read_tflink(tflink_ls_file)

  # Filter TFLink interactions for selected genes
  interactions_tflink_small <- tflink_small %>%
    filter(Name.TF %in% all_selected_genes & Name.Target %in% all_selected_genes) %>%
    select(from = Name.TF, to = Name.Target) %>%
    mutate(source = "TFLink_SS", score = NA) # Add source, score placeholder

  interactions_tflink_large <- tflink_large %>%
    filter(Name.TF %in% all_selected_genes & Name.Target %in% all_selected_genes) %>%
    select(from = Name.TF, to = Name.Target) %>%
    mutate(source = "TFLink_LS", score = NA) # Add source, score placeholder

  # --- Combine Interactions ---
  message("Combining interactions...")
  all_interactions <- bind_rows(interactions_string, interactions_tflink_small, interactions_tflink_large) %>%
    # Ensure interaction is only between nodes present in the final node list
    filter(from %in% all_selected_genes & to %in% all_selected_genes) %>%
    # Remove self-loops for clarity in visualization
    filter(from != to) %>%
    # Keep only one edge between two nodes regardless of direction for undirected graph
    # Sort node pairs alphabetically to treat A-B and B-A as the same
    rowwise() %>%
    mutate(node1 = min(from, to), node2 = max(from, to)) %>%
    ungroup() %>%
    distinct(node1, node2, .keep_all = TRUE) %>% # Keep strongest score or first source if duplicates exist
    select(from, to, source, score)


  if (nrow(all_interactions) == 0) {
    warning("No interactions found between the selected genes/TRs after combining sources.")
    # Optionally, create a plot with disconnected nodes or stop
  }

  # --- Prepare Nodes Data Frame ---
  nodes_df <- data.frame(
    id = all_selected_genes,
    type = ifelse(all_selected_genes %in% trs_selected, "TR", "Target Gene")
  ) %>%
    # Ensure nodes included are actually part of an interaction or are TRs
    filter(id %in% unique(c(all_interactions$from, all_interactions$to)) | type == "TR")

  # Refilter interactions to only include edges where both nodes are in the final nodes_df
  final_nodes_ids <- nodes_df$id
  edges_df <- all_interactions %>%
    filter(from %in% final_nodes_ids & to %in% final_nodes_ids)

  if (nrow(edges_df) == 0 && nrow(nodes_df) > 0) {
    warning("No edges connect the final set of nodes. The plot will show disconnected nodes.")
  } else if (nrow(nodes_df) == 0) {
    stop("Cannot create graph: No nodes remaining after filtering.")
  }

  # --- Create igraph object ---
  message("Building graph object...")
  graph <- graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = FALSE)

  # --- Create Network Plot using ggraph ---
  message("Generating network plot using ggraph with layout: ", layout_algorithm)

  # Define colors
  edge_colors <- c("STRING" = "red", "TFLink_SS" = "orange", "TFLink_LS" = "grey")
  node_fill_colors <- c("TR" = "pink", "Target Gene" = "lightblue")
  # Create layout *before* plotting - helps reproducibility and access to coordinates if needed
  # Using create_layout is preferred over passing layout string directly to ggraph sometimes
  set.seed(123) # Set seed for reproducible layout
  graph_layout <- create_layout(graph, layout = "nicely")

  # Build the plot
  gg_plot <- ggraph(graph_layout) +
    # Edges: Draw links first (UNDIRECTED)
    geom_edge_link(aes(color = source),          # Map edge color to the source
                   alpha = 0.7,
                   # arrow = arrow(length = unit(1.5, 'mm')), # REMOVED this line
                   end_cap = circle(1.5, 'mm'),     # Keep edges from overlapping node center
                   show.legend = TRUE) +
    scale_edge_color_manual(values = edge_colors, name = "Interaction Source") +

    # Nodes: Use geom_node_label
    geom_node_label(aes(label = name, fill = type),
                    color = "black",
                    size = label_size,
                    label.padding = unit(0.2, "lines"),
                    label.r = unit(0.15, "lines"),
                    label.size = 0.1,
                    show.legend = FALSE
    ) +
    scale_fill_manual(values = node_fill_colors, name = "Node Type") +

    # Theme
    theme_graph(base_family = 'sans') +
    labs(
      title = "Gene Regulatory Network",
      subtitle = paste("Top", thres_TRs, "TRs and associated genes (Link:", link, ", Genome:", genome, ")")
    ) +
    theme(
      legend.position = "right"
    )

  message("Plot generated successfully.")
  return(gg_plot)
}

