# ==== Load required package ====
if (!require("data.table")) install.packages("data.table")
library(data.table)

# ==== Define folders ====
readcount_dir <- "./readcount"
nexus_dir <- "./nexus"
ref_path <- "/home/yil218/data/10K.fna"
output_base <- "./pileup"
# Remove and recreate top-level output folder if exists
if (dir.exists(output_base)) {
  cat("Removing existing pileup folder:", output_base, "\n")
  unlink(output_base, recursive = TRUE)
}
dir.create(output_base)

# ==== Load reference sequence ====
ref_file <- readLines(ref_path)
ref_seq <- toupper(paste(ref_file[!grepl("^>", ref_file)], collapse = ""))

# ==== Get list of readcount files ====
readcount_files <- list.files(readcount_dir, full.names = TRUE, pattern = "_readcount\\.log$")

# ==== Process each readcount file ====
for (file_path in readcount_files) {
  start_time_all <- Sys.time()
  
  # Extract r number, e.g., r0
  file_name <- basename(file_path)
  r_match <- regmatches(file_name, regexpr("r[0-9]+", file_name))
  if (length(r_match) == 0) {
    cat("Skipping file (no r number):", file_name, "\n")
    next
  }
  r_number <- r_match
  
  # Load corresponding NEXUS file
  nexus_file <- file.path(nexus_dir, paste0("sequenceGeneration_", r_number, "_scA.nexus"))
  if (!file.exists(nexus_file)) {
    cat("NEXUS file not found for", r_number, "\n")
    next
  }
  nexus_lines <- readLines(nexus_file)
  taxlabel_line_index <- grep("taxlabels", nexus_lines)
  if (length(taxlabel_line_index) == 0) {
    cat("No 'taxlabels' line in NEXUS file for", r_number, "\n")
    next
  }
  
  # Extract cell names
  taxlabel_line <- nexus_lines[taxlabel_line_index]
  taxlabel_line <- gsub("^\\s*taxlabels\\s+", "", taxlabel_line)
  taxlabel_line <- gsub(";", "", taxlabel_line)
  cell_names <- unlist(strsplit(taxlabel_line, "\\s+"))
  
  # Load readcount lines
  lines <- readLines(file_path)
  if (length(cell_names) != length(lines)) {
    cat("Cell name count", length(cell_names), "does not match line count",length(lines),"in", r_number, "\n")
    next
  }
  
  # Output folder for this r
  output_folder <- file.path(output_base, r_number)
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Bases
  bases <- c("A", "C", "G", "T")
  
  # Process each cell
  for (i in seq_along(lines)) {
    line <- lines[i]
    entries <- unlist(strsplit(line, ","))
    
    rows_list <- vector("list", length(entries))
    row_idx <- 1
    
    for (j in seq_along(entries)) {
      counts <- as.integer(unlist(strsplit(entries[j], ":")))
      if (length(counts) != 4) next
      count_sum <- sum(counts)
      if (count_sum == 0) next
      
      pos <- j
      ref_base <- substr(ref_seq, pos, pos)
      read_result <- ""
      
      for (k in seq_along(bases)) {
        base_count <- counts[k]
        if (base_count > 0) {
          if (bases[k] == ref_base) {
            read_result <- paste0(read_result, strrep(".", base_count))
          } else {
            read_result <- paste0(read_result, strrep(bases[k], base_count))
          }
        }
      }
      
      rows_list[[row_idx]] <- list(
        Name = "chr1",
        Position = pos,
        Ref = ref_base,
        Count = count_sum,
        ReadResults = read_result
      )
      row_idx <- row_idx + 1
    }
    
    rows_list <- rows_list[1:(row_idx - 1)]
    if (length(rows_list) == 0) {
      cat("No valid sites for cell", i, "(", cell_names[i], "), skipping.\n")
      next
    }
    
    cell_pileup <- rbindlist(rows_list)
    
    # Output file
    cell_name <- cell_names[i]
    output_file <- file.path(output_folder, paste0(cell_name, ".pileup"))
    fwrite(cell_pileup, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    cat(r_number, "cell", i, "(", cell_name, ") written to:", output_file, "\n")
  }
  
  cat("All pileups for", r_number, "written in",
      round(difftime(Sys.time(), start_time_all, units = "secs"), 2), "seconds.\n\n")
}

cat("All readcount files processed.\n")