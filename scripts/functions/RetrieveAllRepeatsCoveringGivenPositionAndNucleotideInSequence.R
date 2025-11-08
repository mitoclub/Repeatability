## KP function to Get All Repeats for All Motivs Overlapping The Position with a Given Nucleotide
# Required libraries
library(Biostrings)  # For DNAString and sequence manipulation
library(stringdist)  # For hamming distance calculation

# Complete function including sequence editing, motif extraction, repeat search, and output writing
RetrieveAllRepeatsCoveringGivenPositionAndNucleotideInSequence <- function(seq, pos, fixed_nuc, output_folder,
                                        max_flank_left = 20, max_flank_right = 20,
                                        min_length = 5, max_length = 41,
                                        mismatch_fraction = 0.2) 
  
  {  # Default 20% mismatches
  # ---- Step 1: Ensure nucleotide at 'pos' is 'fixed_nuc', modify if necessary ----
  current_nuc <- as.character(subseq(seq, start=pos, width=1))
  if (current_nuc != fixed_nuc) {
    seq_char <- unlist(strsplit(as.character(seq), ""))
    seq_char[pos] <- fixed_nuc
    seq <- DNAString(paste(seq_char, collapse=""))
  }
  
  # ---- Helper function: Extract motif around position with asymmetric flanks ----
  get_motif_around_position <- function(seq, pos, left_flank, right_flank) {
    start_pos <- max(pos - left_flank, 1)
    end_pos <- min(pos + right_flank, length(seq))
    subseq(seq, start=start_pos, end=end_pos)
  }
  
  # ---- Helper function: Find approximate repeats allowing mismatches (Hamming distance) ----
  find_approximate_repeats <- function(seq, motif, max_mismatch) {
    motif_str <- as.character(motif)
    motif_length <- nchar(motif_str)
    seq_length <- length(seq)
    # Extract all k-mers of motif length from sequence
    all_kmers <- substring(as.character(seq), 1:(seq_length - motif_length + 1), motif_length:seq_length)
    # Calculate Hamming distances of motif vs all k-mers
    distances <- stringdist(motif_str, all_kmers, method = "hamming")
    # Positions where mismatch count is within allowed threshold
    matches_pos <- which(distances <= max_mismatch)
    matched_kmers <- all_kmers[matches_pos]
    # Return data frame with repeat info
    data.frame(
      repeat.seq = matched_kmers,
      repeat.start = matches_pos,
      repeat.end = matches_pos + motif_length - 1,
      repeat.hamming.distance = distances[matches_pos],
      stringsAsFactors = FALSE
    )
  }
  
  # ---- Main search: Find repeats for motifs of varying lengths and flanks ----
  all_results <- data.frame()
  for (motif_length in min_length:max_length) {
    for (left_flank in 0:(motif_length - 1)) {
      right_flank <- motif_length - left_flank - 1
      # Only consider flanks within max allowed size
      if (left_flank <= max_flank_left && right_flank <= max_flank_right) {
        motif_seq <- get_motif_around_position(seq, pos, left_flank, right_flank)
        motif_string <- as.character(motif_seq)
        this_length <- nchar(motif_string)
        max_mismatch_allowed <- floor(mismatch_fraction * this_length)  # User-defined fraction
        repeats_df <- find_approximate_repeats(seq, motif_seq, max_mismatch_allowed)
        if (nrow(repeats_df) > 0) {
          # Add motif and position metadata to repeats dataframe
          repeats_df$motif.seq <- motif_string
          repeats_df$motif.length <- this_length
          repeats_df$motif.start <- max(pos - left_flank, 1)
          repeats_df$motif.end <- min(pos + right_flank, length(seq))
          repeats_df$nuc <- fixed_nuc
          repeats_df$effective.length <- repeats_df$motif.length - repeats_df$repeat.hamming.distance
          # Select and order columns
          repeats_df <- repeats_df[, c(
            "nuc", "motif.seq", "motif.length",
            "motif.start", "motif.end",
            "repeat.seq", "repeat.start", "repeat.end",
            "repeat.hamming.distance", "effective.length"
          )]
          # Aggregate results
          all_results <- rbind(all_results, repeats_df)
        }
      }
    }
  }
  
  # ---- Sort results and write to output file if any repeats found ----
  if (nrow(all_results) > 0) {
    all_results <- cbind(pos = pos, all_results)  # Add 'pos' as first column
    all_results <- all_results[order(-all_results$effective.length), ]  # Sort descending by effective length
    
    # Construct output filename automatically using pos and fixed_nuc
    output_filename <- paste0("01KP.", pos, ".", fixed_nuc, ".txt")
    output_path <- file.path(output_folder, output_filename)
    
    # Write results as tab-delimited text with header, no row names or quotes
    write.table(all_results, file = output_path, sep = "\t",
                row.names = FALSE, quote = FALSE)
  } else {
    message("No repeats found; no file created.")
  }
  
  # Return results as data frame
  return(all_results)
}

### how to source it and use from R code
# source("/home/popadin/Documents/REPEATABILITY/scripts/functions/RetrieveAllRepeatsCoveringGivenPositionAndNucleotideInSequence.R")
# RetrieveAllRepeatsCoveringGivenPositionAndNucleotideInSequence(seq, pos, fixed_nuc, output_folder) # other parameters keep by default
