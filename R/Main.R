vBIT <- function(file, output_path, show=TRUE, plot.bar=TRUE, format=NULL, N = 5000 ,bin_width = 1000, genome=c("hg38","mm10")){
  print("Load and map peaks to bins...")

  output_path = R.utils::getAbsolutePath(output_path)

  peak_inds <- import_input_regions(file = file, format = format, bin_width = bin_width, genome = genome)
  filtered_peak_inds <- filter_by_distal(peak_inds, genome=genome)

  print("Done.")
  print(paste0("compare the input regions with the pre-compiled reference ChIP-seq data, bin width used: ",bin_width," bps"))

  alignment_results <- alignment_wrapper(filtered_peak_inds, bin_width = bin_width, genome = genome)

  print("Done.")

  xct <- alignment_results$GOOD
  nct <- alignment_results$TOTAL

  tr_labels <- as.numeric(factor(alignment_results$TR))

  print(paste0("Start vBIT inference, iterations: ",N))

  VI_chain_results <- Main_Function(N, xct, nct, tr_labels)

  VI_chain_results[["TR_names"]] <- alignment_results$TR

  print("Done.")

  file_name<-paste0(output_path,"/",tools::file_path_sans_ext(basename(file)),".rds")
  saveRDS(VI_chain_results,file_name)

  print(paste0("Output data saved as ",file_name))

  if(show==TRUE){
    display_tables(file_path = file_name, output_path = output_path)
  }

#  if(plot.bar==TRUE){
#    rank_plot(file_path = file_name, output_path = output_path, burnin = burnin)
#  }

  return()
}













