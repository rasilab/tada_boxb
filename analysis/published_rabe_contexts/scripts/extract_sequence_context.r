suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(Rsamtools)
  library(GenomicRanges)
  library(Biostrings)
})

# Snakemake supplies all analysis choices and repository-relative paths.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 9)
genome <- args[1]
counts_file <- args[2]
sites_file <- args[3]
summary_file <- args[4]
qc_file <- args[5]
sample_name <- args[6]
read_cutoff <- as.numeric(args[7])
chunk_size <- as.integer(args[8])
converted_base <- args[9]
bases <- c("A", "C", "G", "T")
count_columns <- c("ref", "pos", "strand", "convertedBaseCount", "unconvertedBaseCount")

stopifnot(file.exists(paste0(genome, ".fai")), read_cutoff >= 0, chunk_size > 0)
fa <- FaFile(genome)
open(fa)
reference_index <- scanFaIndex(fa)
reference_lengths <- setNames(width(reference_index), as.character(seqnames(reference_index)))
site_output <- gzfile(sites_file, "wt")
writeLines(paste(c(count_columns, "total_count", "percent_edited",
                  "upstream_base", "downstream_base"), collapse = ","), site_output)
summary <- expand.grid(upstream_base = bases, downstream_base = bases,
                       stringsAsFactors = FALSE) %>%
  mutate(n_sites = 0, sum_editing_fraction = 0)
qc <- c(input_sites = 0, covered_sites = 0, edited_covered_sites = 0,
        ambiguous_flank_sites = 0, summarized_sites = 0)

# One callback bounds memory while preserving the mean over site observations.
process_chunk <- function(data, position) {
  if (nrow(problems(data))) stop("Malformed count CSV at row ", position)
  stopifnot(identical(names(data), count_columns), !anyNA(data),
            all(data$strand %in% c("+", "-")),
            all(data$pos > 0 & data$pos == floor(data$pos)),
            all(data$convertedBaseCount >= 0), all(data$unconvertedBaseCount >= 0),
            all(data$convertedBaseCount == floor(data$convertedBaseCount)),
            all(data$unconvertedBaseCount == floor(data$unconvertedBaseCount)))
  data <- mutate(data, total_count = convertedBaseCount + unconvertedBaseCount)
  qc["input_sites"] <<- qc["input_sites"] + nrow(data)
  qc["covered_sites"] <<- qc["covered_sites"] + sum(data$total_count > read_cutoff)
  data <- filter(data, convertedBaseCount > 0, total_count > read_cutoff)
  qc["edited_covered_sites"] <<- qc["edited_covered_sites"] + nrow(data)
  if (!nrow(data)) return(invisible(NULL))

  # Reference coordinates are one-based. A strand-oriented triplet supplies
  # both flanks in one lookup, including the reverse complement on minus strands.
  stopifnot(all(data$ref %in% names(reference_lengths)), all(data$pos > 1),
            all(data$pos < reference_lengths[data$ref]))
  regions <- GRanges(data$ref, IRanges(data$pos - 1, data$pos + 1), data$strand)
  context <- as.character(getSeq(fa, regions))
  stopifnot(all(nchar(context) == 3), all(substr(context, 2, 2) == converted_base))
  data <- mutate(data, percent_edited = convertedBaseCount / total_count,
                 upstream_base = substr(context, 1, 1),
                 downstream_base = substr(context, 3, 3))
  write.table(data, site_output, sep = ",", row.names = FALSE,
              col.names = FALSE, quote = TRUE)
  ambiguous <- data$upstream_base == "N" | data$downstream_base == "N"
  qc["ambiguous_flank_sites"] <<- qc["ambiguous_flank_sites"] + sum(ambiguous)
  data <- data[!ambiguous, ]
  stopifnot(all(data$upstream_base %in% bases), all(data$downstream_base %in% bases))
  qc["summarized_sites"] <<- qc["summarized_sites"] + nrow(data)
  if (!nrow(data)) return(invisible(NULL))

  chunk_summary <- data %>%
    group_by(upstream_base, downstream_base) %>%
    summarise(n = n(), total = sum(percent_edited), .groups = "drop")
  index <- match(paste(chunk_summary$upstream_base, chunk_summary$downstream_base),
                 paste(summary$upstream_base, summary$downstream_base))
  summary$n_sites[index] <<- summary$n_sites[index] + chunk_summary$n
  summary$sum_editing_fraction[index] <<-
    summary$sum_editing_fraction[index] + chunk_summary$total
  invisible(NULL)
}

read_csv_chunked(counts_file, SideEffectChunkCallback$new(process_chunk),
                 chunk_size = chunk_size, progress = FALSE,
                 col_types = cols(ref = col_character(), pos = col_double(),
                                  strand = col_character(),
                                  convertedBaseCount = col_double(),
                                  unconvertedBaseCount = col_double()))
close(site_output)
close(fa)
stopifnot(qc["input_sites"] > 0, sum(summary$n_sites) == qc["summarized_sites"])
summary %>%
  mutate(sample_name = sample_name,
         mean_percent = if_else(n_sites > 0, 100 * sum_editing_fraction / n_sites, NA_real_)) %>%
  select(sample_name, everything()) %>%
  write_csv(summary_file, na = "")
tibble::as_tibble_row(c(list(sample_name = sample_name), as.list(qc))) %>% write_csv(qc_file)
