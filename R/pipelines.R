# FUNCTION: find_startstop_ORFs ---------------------------------------------------------------

#' Find Start-Stop ORFs
#'
#' @param file_name Path to the input FASTA file.
#' @importFrom magrittr %>%
#' @export
find_startstop_ORFs <- function(file_name) {

  # Sanity check: Ensure file exists before starting
  if (!file.exists(file_name)) {
    stop(paste("Error: The file '", file_name, "' could not be found in the current working directory.", sep=""))
  }

  # 1. INTERNAL LOGIC TABLES ------------------------------------------------------------------
  kozak_reference <- tibble::tibble(
    dinuc = c("AG", "GG", "AA", "GA", "AC", "AT", "GC", "GT", "CG", "TG", "CA", "CC", "CT", "TA", "TC", "TT"),
    kozak_strength = c("strong", "strong", "moderate", "moderate", "moderate", "moderate",
                       "moderate", "moderate", "moderate", "moderate", "weak", "weak",
                       "weak", "weak", "weak", "weak")
  )

  frame_correction_table <- tibble::tibble(
    condition = c(
      rep("plus_0", 3), rep("plus_1", 3), rep("plus_2", 3),
      rep("revcom_plus_0", 3), rep("revcom_plus_1", 3), rep("revcom_plus_2", 3)
    ),
    frame_to_change = c(
      "plus_0", "plus_1", "plus_2",
      "plus_0", "plus_1", "plus_2",
      "plus_0", "plus_1", "plus_2",
      "revcom_plus_0", "revcom_plus_1", "revcom_plus_2",
      "revcom_plus_0", "revcom_plus_1", "revcom_plus_2",
      "revcom_plus_0", "revcom_plus_1", "revcom_plus_2"
    ),
    correct_frame = c(
      "plus_0", "plus_1", "plus_2",
      "plus_2", "plus_0", "plus_1",
      "plus_1", "plus_2", "plus_0",
      "revcom_plus_0", "revcom_plus_1", "revcom_plus_2",
      "revcom_plus_2", "revcom_plus_0", "revcom_plus_1",
      "revcom_plus_1", "revcom_plus_2", "revcom_plus_0"
    )
  )

  # 2. INTERNAL HELPER FUNCTIONS --------------------------------------------------------------
  strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

  complem <- function(x) {
    x <- tolower(x)
    x <- gsub("a","T",x)
    x <- gsub("c","G",x)
    x <- gsub("t","A",x)
    x <- gsub("g","C",x)
    return(tolower(x))
  }

  has_ambiguous_bases <- function(sequence_str) {
    grepl("[^ATCG]", toupper(sequence_str))
  }

  get_kozak_dinuc <- function(full_sequence, orf_sequence, orf_start_in_full_seq) {
    kozak_dinucs <- character(length(orf_sequence))

    for (k in seq_along(orf_sequence)) {
      current_full_seq <- full_sequence[k]
      current_orf_seq <- orf_sequence[k]
      current_orf_start <- orf_start_in_full_seq[k]

      full_seq_len <- nchar(current_full_seq)
      orf_seq_len <- nchar(current_orf_seq)

      pos_minus3_in_full <- current_orf_start - 3
      if (pos_minus3_in_full >= 1 & pos_minus3_in_full <= full_seq_len) {
        dinuc1 <- stringr::str_sub(current_full_seq, pos_minus3_in_full, pos_minus3_in_full)
      } else {
        dinuc1 <- "_"
      }

      pos_plus4_in_orf <- 4
      if (pos_plus4_in_orf >= 1 & pos_plus4_in_orf <= orf_seq_len) {
        dinuc2 <- stringr::str_sub(current_orf_seq, pos_plus4_in_orf, pos_plus4_in_orf)
      } else {
        dinuc2 <- NA_character_
      }

      kozak_dinucs[k] <- paste0(as.character(dinuc1)[1], as.character(dinuc2)[1])
    }
    return(kozak_dinucs)
  }

  process_strand_ORFs <- function(sequences_df, strand_type, frame_prefix) {
    sequences_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        processed_sequence = if (strand_type == "negative") {
          toupper(complem(strReverse(sequence)))
        } else {
          toupper(sequence)
        },
        ORFs = list(as.data.frame(ORFik::findORFs(processed_sequence, longestORF = FALSE, "ATG")))
      ) %>%
      tidyr::unnest(ORFs) %>%
      dplyr::ungroup() %>%
      tibble::as_tibble() %>%
      dplyr::rename(start = start, end = end, width = width) %>%
      dplyr::mutate(
        ORF_ID = paste(seq_name, dplyr::row_number(), sep = "-"),
        length = width,
        nuc_sequence = stringr::str_sub(processed_sequence, start, end),
        kozak_dinuc = get_kozak_dinuc(processed_sequence, nuc_sequence, start)
      ) %>%
      dplyr::select(seq = seq_name, ORF_ID, start, end, length, nuc_sequence, kozak_dinuc)
  }

  # 3. CORE PIPELINE PROCESSING ----------------------------------------------------------------

  # Generate base segment name for files
  segment_name <- sub("\\.[^.]+$", "", file_name)
  segment_name <- sub("_complete$", "", segment_name)

  # PASS 1: Check for sequences containing ambiguous bases
  all_seqs_raw_list <- seqinr::read.fasta(file = file_name, set.attributes = FALSE, as.string = TRUE)
  has_ambiguous_lgl_vec <- sapply(all_seqs_raw_list, has_ambiguous_bases)
  ambiguous_names <- names(all_seqs_raw_list[has_ambiguous_lgl_vec])

  # Export excluded sequences to a text file
  writeLines(ambiguous_names, paste0(segment_name, "_excluded_sequences.txt"))

  if(length(ambiguous_names) > 0) {
    warning(paste("Filtered out", length(ambiguous_names), "sequences due to ambiguous bases. Names written to 'excluded_sequences.txt'"))
  }

  # PASS 2: Import clean sequences using Biostrings
  fasta_file <- Biostrings::readDNAStringSet(file_name)
  all_biostring_names <- names(fasta_file)
  clean_names <- all_biostring_names[!(all_biostring_names %in% ambiguous_names)]
  sequences_without_ambiguities <- fasta_file[clean_names]

  if(length(sequences_without_ambiguities) == 0) {
    stop("All sequences inside the file were filtered out due to ambiguities. Execution halted.")
  }

  temp_filtered_fasta <- tempfile(pattern = "filtered_sequences_", fileext = ".fasta")
  Biostrings::writeXStringSet(sequences_without_ambiguities, temp_filtered_fasta)
  all_sequences_list <- seqinr::read.fasta(file = temp_filtered_fasta, set.attributes = FALSE, as.string = TRUE)

  all_sequences_df <- tibble::tibble(
    seq_name = names(all_sequences_list),
    sequence = as.character(all_sequences_list)
  )

  # Process positive strand
  positive_ORF_1 <- process_strand_ORFs(all_sequences_df, "positive", "plus_")

  # Frame correction (positive)
  max_pos_start <- max(as.numeric(positive_ORF_1$start), na.rm = TRUE)
  placeholder_frame_table <- tibble::tibble(
    nuc_position = as.character(1:(max_pos_start + 1000)),
    placeholder_frame = rep(c("plus_0", "plus_1", "plus_2"), length.out = max_pos_start + 1000)
  )

  positive_ORF_2 <- positive_ORF_1 %>%
    dplyr::mutate(start = as.character(start)) %>%
    dplyr::left_join(placeholder_frame_table, by = c("start" = "nuc_position"))

  positive_ORF_3 <- positive_ORF_2 %>%
    dplyr::group_by(seq) %>%
    dplyr::mutate(
      length = as.numeric(length),
      longest_orf_placeholder = placeholder_frame[which.max(length)],
      correct_frame_lookup = purrr::map(longest_orf_placeholder, ~frame_correction_table %>%
                                   dplyr::filter(condition == .x)),
      frame = purrr::map2_chr(placeholder_frame, correct_frame_lookup,
                       ~dplyr::filter(.y, frame_to_change == .x)$correct_frame)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-placeholder_frame, -longest_orf_placeholder, -correct_frame_lookup) %>%
    dplyr::mutate(sense = "positive")

  # Process negative strand
  negative_ORF_1 <- process_strand_ORFs(all_sequences_df, "negative", "revcom_plus_")

  # Frame correction (negative)
  max_neg_start <- max(as.numeric(negative_ORF_1$start), na.rm = TRUE)
  placeholder_frame_table_neg <- tibble::tibble(
    nuc_position = as.character(1:(max_neg_start + 1000)),
    placeholder_frame = rep(c("revcom_plus_0", "revcom_plus_1", "revcom_plus_2"), length.out = max_neg_start + 1000)
  )

  negative_ORF_2 <- negative_ORF_1 %>%
    dplyr::mutate(start = as.character(start)) %>%
    dplyr::left_join(placeholder_frame_table_neg, by = c("start" = "nuc_position"))

  negative_ORF_3 <- negative_ORF_2 %>%
    dplyr::group_by(seq) %>%
    dplyr::mutate(
      length = as.numeric(length),
      longest_orf_placeholder = placeholder_frame[which.max(length)],
      correct_frame_lookup = purrr::map(longest_orf_placeholder, ~frame_correction_table %>%
                                   dplyr::filter(condition == .x)),
      frame = purrr::map2_chr(placeholder_frame, correct_frame_lookup,
                       ~dplyr::filter(.y, frame_to_change == .x)$correct_frame)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-placeholder_frame, -longest_orf_placeholder, -correct_frame_lookup) %>%
    dplyr::mutate(sense = "negative")

  # Compiling data
  all_ORFs <- dplyr::bind_rows(positive_ORF_3, negative_ORF_3)

  # Adding Kozak strength calculations
  all_ORFs <- all_ORFs %>%
    dplyr::left_join(kozak_reference, by = c("kozak_dinuc" = "dinuc"))

  output_file_name <- paste0(segment_name, "_start-stop_ORFs.csv")
  write.csv(all_ORFs, file = output_file_name, row.names = FALSE)

  # Translating longest ORFs and exporting as FASTA file
  sequences_for_translation <- all_ORFs %>%
    dplyr::filter(frame == "plus_0") %>%
    dplyr::group_by(seq) %>%
    dplyr::slice_max(order_by = length, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  all_translated_sequences <- sequences_for_translation %>%
    dplyr::mutate(
      nuc_sequence_lower = tolower(nuc_sequence),
      translated_aa = map_chr(nuc_sequence_lower, ~{
        seq_chars <- strsplit(.x, "")[[1]]
        paste(seqinr::translate(seq_chars, frame = 0), collapse = "")
      })
    ) %>%
    pull(translated_aa, name = seq)

  translation_file_name <- paste(segment_name, "translated.fasta", sep = "-")
  seqinr::write.fasta(as.list(all_translated_sequences), names = names(all_translated_sequences), file.out = translation_file_name)

  # Diagnostic Messages
  message("--- Run Complete ---")
  message(paste("1. Excluded headers list exported to: 'excluded_sequences.txt'"))
  message(paste("2. Detailed ORF mapping spreadsheet created: '", output_file_name, "'", sep=""))
  message(paste("3. Canonical protein translation file created: '", translation_file_name, "'", sep=""))

  # Return dataframe to environment
  return(all_ORFs)
}

# FUNCTION: find_stopstop_ORFs ----------------------------------------------------------------

#' Find Stop-Stop ORFs
#'
#' @param file_name Path to the input FASTA file.
#' @importFrom magrittr %>%
#' @export
find_stopstop_ORFs <- function(file_name) {

  # Sanity check: Ensure file exists before running pipeline
  if (!file.exists(file_name)) {
    stop(paste("Error: The file '", file_name, "' could not be found in the current working directory.", sep=""))
  }

  # 1. INTERNAL LOGIC TABLES ------------------------------------------------------------------
  frame_correction_table <- tibble::tibble(
    condition = c(
      rep("plus_0", 3), rep("plus_1", 3), rep("plus_2", 3),
      rep("revcom_plus_0", 3), rep("revcom_plus_1", 3), rep("revcom_plus_2", 3)
    ),
    frame_to_change = c(
      "plus_0", "plus_1", "plus_2",
      "plus_0", "plus_1", "plus_2",
      "plus_0", "plus_1", "plus_2",
      "revcom_plus_0", "revcom_plus_1", "revcom_plus_2",
      "revcom_plus_0", "revcom_plus_1", "revcom_plus_2",
      "revcom_plus_0", "revcom_plus_1", "revcom_plus_2"
    ),
    correct_frame = c(
      "plus_0", "plus_1", "plus_2",
      "plus_2", "plus_0", "plus_1",
      "plus_1", "plus_2", "plus_0",
      "revcom_plus_0", "revcom_plus_1", "revcom_plus_2",
      "revcom_plus_2", "revcom_plus_0", "revcom_plus_1",
      "revcom_plus_1", "revcom_plus_2", "revcom_plus_0"
    )
  )

  # 2. INTERNAL HELPER FUNCTIONS --------------------------------------------------------------
  strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

  complem <- function(x) {
    x <- tolower(x)
    x <- gsub("a","T",x)
    x <- gsub("c","G",x)
    x <- gsub("t","A",x)
    x <- gsub("g","C",x)
    return(tolower(x))
  }

  has_ambiguous_bases <- function(sequence_str) {
    grepl("[^ATCG]", toupper(sequence_str))
  }

  # Captures absolute continuous stop-free blocks, including terminal ends
  find_stop_free_blocks <- function(sequence_str) {
    seq_len <- nchar(sequence_str)
    all_orfs_list <- list()

    for (f in 0:2) {
      sub_seq <- stringr::str_sub(sequence_str, f + 1, seq_len)
      codons <- strsplit(sub_seq, "(?<=...)", perl = TRUE)[[1]]

      # Fixed the assignment operator typo here:
      is_stop <- codons %in% c("TAA", "TAG", "TGA")
      stop_indices <- which(is_stop)

      start_chunk_idx <- c(1, stop_indices + 1)
      end_chunk_idx <- c(stop_indices - 1, length(codons))

      for (i in seq_along(start_chunk_idx)) {
        s_idx <- start_chunk_idx[i]
        e_idx <- end_chunk_idx[i]

        if (s_idx <= e_idx) {
          orf_codons <- codons[s_idx:e_idx]
          nuc_seq <- paste(orf_codons, collapse = "")

          start_pos <- f + ((s_idx - 1) * 3) + 1
          end_pos <- f + (e_idx * 3)

          if (nchar(nuc_seq) > 0) {
            all_orfs_list[[length(all_orfs_list) + 1]] <- tibble::tibble(
              start = start_pos,
              end = end_pos,
              width = nchar(nuc_seq)
            )
          }
        }
      }
    }
    return(dplyr::bind_rows(all_orfs_list))
  }

  process_strand_ORFs <- function(sequences_df, strand_type, frame_prefix) {
    sequences_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        processed_sequence = if (strand_type == "negative") {
          toupper(complem(strReverse(sequence)))
        } else {
          toupper(sequence)
        },
        ORFs = list(find_stop_free_blocks(processed_sequence))
      ) %>%
      tidyr::unnest(ORFs) %>%
      dplyr::ungroup() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        ORF_ID = paste(seq_name, dplyr::row_number(), sep = "-"),
        length = width,
        nuc_sequence = stringr::str_sub(processed_sequence, start, end)
      ) %>%
      dplyr::select(seq = seq_name, ORF_ID, start, end, length, nuc_sequence)
  }

  # 3. CORE PIPELINE PROCESSING ----------------------------------------------------------------
  segment_name <- sub("\\.[^.]+$", "", file_name)
  segment_name <- sub("_complete$", "", segment_name)

  # PASS 1: Check for sequences containing ambiguous bases
  all_seqs_raw_list <- seqinr::read.fasta(file = file_name, set.attributes = FALSE, as.string = TRUE)
  has_ambiguous_lgl_vec <- sapply(all_seqs_raw_list, has_ambiguous_bases)
  ambiguous_names <- names(all_seqs_raw_list[has_ambiguous_lgl_vec])

  writeLines(ambiguous_names, "excluded_sequences.txt")

  if(length(ambiguous_names) > 0) {
    warning(paste("Filtered out", length(ambiguous_names), "sequences due to ambiguous bases. Names written to 'excluded_sequences.txt'"))
  }

  # PASS 2: Import clean sequences using Biostrings
  fasta_file <- Biostrings::readDNAStringSet(file_name)
  all_biostring_names <- names(fasta_file)
  clean_names <- all_biostring_names[!(all_biostring_names %in% ambiguous_names)]
  sequences_without_ambiguities <- fasta_file[clean_names]

  if(length(sequences_without_ambiguities) == 0) {
    stop("All sequences inside the file were filtered out due to ambiguities. Execution halted.")
  }

  temp_filtered_fasta <- tempfile(pattern = "filtered_sequences_", fileext = ".fasta")
  Biostrings::writeXStringSet(sequences_without_ambiguities, temp_filtered_fasta)
  all_sequences_list <- seqinr::read.fasta(file = temp_filtered_fasta, set.attributes = FALSE, as.string = TRUE)

  all_sequences_df <- tibble::tibble(
    seq_name = names(all_sequences_list),
    sequence = as.character(all_sequences_list)
  )

  # Process positive strand
  positive_ORF_1 <- process_strand_ORFs(all_sequences_df, "positive", "plus_")

  # Frame correction (positive)
  max_pos_start <- max(as.numeric(positive_ORF_1$start), na.rm = TRUE)
  placeholder_frame_table <- tibble::tibble(
    nuc_position = as.character(1:(max_pos_start + 1000)),
    placeholder_frame = rep(c("plus_0", "plus_1", "plus_2"), length.out = max_pos_start + 1000)
  )

  positive_ORF_2 <- positive_ORF_1 %>%
    dplyr::mutate(start = as.character(start)) %>%
    dplyr::left_join(placeholder_frame_table, by = c("start" = "nuc_position"))

  positive_ORF_3 <- positive_ORF_2 %>%
    dplyr::group_by(seq) %>%
    dplyr::mutate(
      length = as.numeric(length),
      longest_orf_placeholder = placeholder_frame[which.max(length)],
      correct_frame_lookup = purrr::map(longest_orf_placeholder, ~frame_correction_table %>%
                                   dplyr::filter(condition == .x)),
      frame = purrr::map2_chr(placeholder_frame, correct_frame_lookup,
                       ~dplyr::filter(.y, frame_to_change == .x)$correct_frame)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-placeholder_frame, -longest_orf_placeholder, -correct_frame_lookup) %>%
    dplyr::mutate(sense = "positive")

  # Process negative strand
  negative_ORF_1 <- process_strand_ORFs(all_sequences_df, "negative", "revcom_plus_")

  # Frame correction (negative)
  max_neg_start <- max(as.numeric(negative_ORF_1$start), na.rm = TRUE)
  placeholder_frame_table_neg <- tibble::tibble(
    nuc_position = as.character(1:(max_neg_start + 1000)),
    placeholder_frame = rep(c("revcom_plus_0", "revcom_plus_1", "revcom_plus_2"), length.out = max_neg_start + 1000)
  )

  negative_ORF_2 <- negative_ORF_1 %>%
    dplyr::mutate(start = as.character(start)) %>%
    dplyr::left_join(placeholder_frame_table_neg, by = c("start" = "nuc_position"))

  negative_ORF_3 <- negative_ORF_2 %>%
    dplyr::group_by(seq) %>%
    dplyr::mutate(
      length = as.numeric(length),
      longest_orf_placeholder = placeholder_frame[which.max(length)],
      correct_frame_lookup = purrr::map(longest_orf_placeholder, ~frame_correction_table %>%
                                   dplyr::filter(condition == .x)),
      frame = purrr::map2_chr(placeholder_frame, correct_frame_lookup,
                       ~dplyr::filter(.y, frame_to_change == .x)$correct_frame)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-placeholder_frame, -longest_orf_placeholder, -correct_frame_lookup) %>%
    dplyr::mutate(sense = "negative")

  # Compiling data
  all_ORFs <- dplyr::bind_rows(positive_ORF_3, negative_ORF_3)

  output_file_name <- paste0(segment_name, "_stop-stop_ORFs.csv")
  write.csv(all_ORFs, file = output_file_name, row.names = FALSE)

  # Translating longest ORFs and exporting as FASTA file
  sequences_for_translation <- all_ORFs %>%
    dplyr::filter(frame == "plus_0") %>%
    dplyr::group_by(seq) %>%
    dplyr::slice_max(order_by = length, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  all_translated_sequences <- sequences_for_translation %>%
    dplyr::mutate(
      nuc_sequence_lower = tolower(nuc_sequence),
      translated_aa = map_chr(nuc_sequence_lower, ~{
        seq_chars <- strsplit(.x, "")[[1]]
        paste(seqinr::translate(seq_chars, frame = 0), collapse = "")
      })
    ) %>%
    pull(translated_aa, name = seq)

  translation_file_name <- paste(segment_name, "translated.fasta", sep = "-")
  seqinr::write.fasta(as.list(all_translated_sequences), names = names(all_translated_sequences), file.out = translation_file_name)

  # Diagnostic Messages
  message("--- Run Complete ---")
  message(paste("1. Excluded headers list exported to: 'excluded_sequences.txt'"))
  message(paste("2. Detailed Stop-Stop ORF mapping spreadsheet created: '", output_file_name, "'", sep=""))
  message(paste("3. Canonical protein translation file created: '", translation_file_name, "'", sep=""))

  return(all_ORFs)
}

# FUNCTION: graph_startstop_ORFs (WITH A4 PANEL EXPORT) --------------------------------------------

#' Plot Start-Stop ORFs
#'
#' @param orf_data Path to the input CSV file generated by find_startstop_ORFs.
#' @importFrom magrittr %>%
#' @export
graph_startstop_ORFs <- function(orf_data, output_prefix = "ORF_plot") {

  # 1. INPUT HANDLING -------------------------------------------------------------------------
  if (is.character(orf_data)) {
    if (!file.exists(orf_data)) {
      stop(paste("Error: The file '", orf_data, "' could not be found.", sep=""))
    }
    orf_data <- read.csv(orf_data, header = TRUE, stringsAsFactors = FALSE)
  }

  if (!"kozak_strength" %in% colnames(orf_data)) {
    stop("Error: Input data must be the output of the start-stop pipeline (requires a 'kozak_strength' column).")
  }

  # 2. DATA CLEANING & TYPE FORCING -----------------------------------------------------------
  orf_data$start <- suppressWarnings(as.numeric(orf_data$start))
  orf_data$end   <- suppressWarnings(as.numeric(orf_data$end))
  orf_data <- orf_data[!is.na(orf_data$start) & !is.na(orf_data$end), ]

  if (nrow(orf_data) == 0) {
    stop("Error: Zero rows remained after forcing numeric conversion.")
  }

  # 3. INTERNAL HELPER ------------------------------------------------------------------------
  ceiling_to_nearest_100 <- function(x) { ceiling(x / 100) * 100 }

  # 4. GRAPHING DATA PREPARATION --------------------------------------------------------------
  max_seq_len <- max(orf_data$end, na.rm = TRUE)
  plot_x_limit <- ceiling_to_nearest_100(max_seq_len)

  visualisation_dataset <- orf_data %>%
    dplyr::mutate(block_seq_number = as.numeric(factor(seq, levels = unique(seq))))

  all_frames <- c("plus_0", "plus_1", "plus_2", "revcom_plus_0", "revcom_plus_1", "revcom_plus_2")

  # Environment container to hold plots for the composite A4 panel
  plot_list <- list()

  # 5. ITERATIVE PLOT GENERATION --------------------------------------------------------------
  for (current_frame in all_frames) {

    block_dataset_for_graph <- visualisation_dataset %>%
      dplyr::filter(frame == current_frame)

    # If a frame has no data, create a blank placeholder plot to keep the layout grid aligned
    if (nrow(block_dataset_for_graph) == 0) {
      p_blank <- ggplot() +
        theme_void() +
        ggtitle(paste(output_prefix, current_frame, "(No ORFs found)", sep = " - ")) +
        theme(plot.title = element_text(size = 10, face = "italic", hjust = 0.5))

      plot_list[[current_frame]] <- p_blank
      next
    }

    block_graph_dataset <- block_dataset_for_graph %>%
      dplyr::mutate(start = as.character(start), end = as.character(end)) %>%
      tidyr::pivot_longer(
        cols = c(start, end),
        names_to = "position_type",
        values_to = "position"
      ) %>%
      dplyr::mutate(
        category = case_when(
          position_type == "end" ~ "stop",
          is.na(kozak_strength) | kozak_strength == "" ~ "unknown",
          TRUE ~ kozak_strength
        )
      ) %>%
      dplyr::select(V1 = block_seq_number, position, category)

    block_graph_dataset$position <- as.numeric(block_graph_dataset$position)

    all_categories_values <- c("strong", "moderate", "weak", "stop", "unknown")
    all_categories_labels <- c("strong Kozak", "moderate Kozak", "weak Kozak", "stop codon", "unknown Kozak")

    block_graph_dataset$category <- factor(
      block_graph_dataset$category,
      levels = all_categories_values,
      labels = all_categories_labels
    )

    dummy_data <- tibble::tibble(
      V1 = NA_real_,
      position = NA_real_,
      category = factor(
        all_categories_values,
        levels = all_categories_values,
        labels = all_categories_labels
      )
    )

    plot_data_with_dummy <- dplyr::bind_rows(block_graph_dataset, dummy_data)
    block_graph_title <- paste(output_prefix, current_frame, sep = " - ")

    block_graph <- ggplot(plot_data_with_dummy,
                          aes(x = position, y = V1, col = category, shape = category)) +
      geom_point(size = 2.0, na.rm = TRUE) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(size = 10, face = "bold")) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      ggtitle(block_graph_title) +

      scale_color_manual(
        name = "Kozak strength/stop codon",
        values = c(
          "strong Kozak" = "chartreuse3",
          "moderate Kozak" = "darkgoldenrod2",
          "weak Kozak" = "deeppink",
          "unknown Kozak" = "blue",
          "stop codon" = "black"
        ),
        drop = FALSE
      ) +
      scale_shape_manual(
        name = "Kozak strength/stop codon",
        values = c(
          "strong Kozak" = 16,
          "moderate Kozak" = 16,
          "weak Kozak" = 16,
          "unknown Kozak" = 16,
          "stop codon" = 16
        ),
        drop = FALSE
      ) +
      coord_cartesian(xlim = c(0, plot_x_limit)) +
      ylab("Unique Sequences") + xlab("nucleotide position")

    # Save the individual standalone image
    block_graph_name <- paste0(output_prefix, "_", current_frame, "_start-stop_ORF.png")
    ggsave(filename = block_graph_name, plot = block_graph, width = 10, height = 6, dpi = 300)
    message(paste("Saved individual plot: '", block_graph_name, "'", sep=""))

    # Store plot for composite mapping
    plot_list[[current_frame]] <- block_graph
  }

  # 6. ASSEMBLE A4 PORTRAIT PANEL SHEET (USING PATCHWORK) -------------------------------------
  message("Assembling consolidated A4 portrait panel layout...")

  # Extract individual plots for programmatic readability
  p1 <- plot_list[["plus_0"]]
  p2 <- plot_list[["plus_1"]]
  p3 <- plot_list[["plus_2"]]
  p4 <- plot_list[["revcom_plus_0"]]
  p5 <- plot_list[["revcom_plus_1"]]
  p6 <- plot_list[["revcom_plus_2"]]

  # Clean up internal axes text for tidy paneling (keeps outer labels intact)
  p1 <- p1 + xlab(NULL); p4 <- p4 + xlab(NULL); p4 <- p4 + ylab(NULL)
  p2 <- p2 + xlab(NULL); p5 <- p5 + xlab(NULL); p5 <- p5 + ylab(NULL)
  p6 <- p6 + ylab(NULL)

  # Design the 3x2 matrix math grid layout:
  # Left column (top to bottom): plus_0, plus_1, plus_2
  # Right column (top to bottom): revcom_plus_0, revcom_plus_1, revcom_plus_2
  a4_composite_panel <- (p1 + p4) / (p2 + p5) / (p3 + p6)

  # Format consolidated multi-plot legend attributes
  a4_composite_panel <- a4_composite_panel +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.box = "horizontal",
          legend.title = element_text(size = 9, face = "bold"),
          legend.text = element_text(size = 8))

  # Calculate dimensions matching to A4 paper ratio (8.27 x 11.69 inches)
  panel_file_name <- paste0(output_prefix, "_A4_composite_panel.png")
  ggsave(filename = panel_file_name, plot = a4_composite_panel,
         width = 8.27, height = 11.69, dpi = 300)

  message(paste("Successfully exported panel sheet: '", panel_file_name, "'", sep=""))
}

# FUNCTION: graph_stopstop_ORFs --------------------------------------------------------------------

#' Plot Stop-Stop ORFs
#'
#' @param orf_data Path to the input CSV file generated by find_startstop_ORFs.
#' @importFrom magrittr %>%
#' @export
graph_stopstop_ORFs <- function(orf_data, output_prefix = "ORF_plot") {

  # 1. INPUT HANDLING -------------------------------------------------------------------------
  if (is.character(orf_data)) {
    if (!file.exists(orf_data)) {
      stop(paste("Error: The file '", orf_data, "' could not be found.", sep=""))
    }
    message("Loading ORF dataset from CSV file...")
    orf_data <- read.csv(orf_data, header = TRUE, stringsAsFactors = FALSE)
  }

  if ("kozak_strength" %in% colnames(orf_data)) {
    stop("Error: Input data appears to be from the start-stop pipeline. Use startstop_graph() instead.")
  }

  # 2. DATA CLEANING & TYPE FORCING -----------------------------------------------------------
  orf_data$start <- suppressWarnings(as.numeric(orf_data$start))
  orf_data$end   <- suppressWarnings(as.numeric(orf_data$end))
  orf_data <- orf_data[!is.na(orf_data$start) & !is.na(orf_data$end), ]

  if (nrow(orf_data) == 0) {
    stop("Error: Zero rows remained after forcing numeric conversion.")
  }

  # 3. INTERNAL HELPER ------------------------------------------------------------------------
  ceiling_to_nearest_100 <- function(x) { ceiling(x / 100) * 100 }

  # 4. GRAPHING DATA PREPARATION --------------------------------------------------------------
  max_seq_len <- max(orf_data$end, na.rm = TRUE)
  plot_x_limit <- ceiling_to_nearest_100(max_seq_len)

  visualisation_dataset <- orf_data %>%
    dplyr::mutate(block_seq_number = as.numeric(factor(seq, levels = unique(seq))))

  all_frames <- c("plus_0", "plus_1", "plus_2", "revcom_plus_0", "revcom_plus_1", "revcom_plus_2")

  # Container to hold plots for the final composite A4 panel
  plot_list <- list()

  # 5. ITERATIVE PLOT GENERATION --------------------------------------------------------------
  for (current_frame in all_frames) {

    block_dataset_for_graph <- visualisation_dataset %>%
      dplyr::filter(frame == current_frame)

    # Handle frames with missing data to keep the 3x2 matrix layout structurally aligned
    if (nrow(block_dataset_for_graph) == 0) {
      p_blank <- ggplot() +
        theme_void() +
        ggtitle(paste(output_prefix, current_frame, "(No stop codons found)", sep = " - ")) +
        theme(plot.title = element_text(size = 10, face = "italic", hjust = 0.5))

      plot_list[[current_frame]] <- p_blank
      next
    }

    # Modifying data map to plot ONLY the 'end' column (the position of the stop codon)
    block_graph_dataset <- block_dataset_for_graph %>%
      dplyr::select(V1 = block_seq_number, position = end) %>%
      dplyr::mutate(category = factor("stop codon")) # Forces a unified categorical value for the legend

    block_graph_title <- paste(output_prefix, current_frame, sep = " - ")

    # Build the clean dot map
    block_graph <- ggplot(block_graph_dataset,
                          aes(x = position, y = V1, col = category, shape = category)) +
      geom_point(size = 1.5, na.rm = TRUE) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(size = 10, face = "bold")) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      ggtitle(block_graph_title) +

      # Enforce a single black dot classification style matching formatting rules
      scale_color_manual(
        name = "Classification",
        values = c("stop codon" = "black")
      ) +
      scale_shape_manual(
        name = "Classification",
        values = c("stop codon" = 16)
      ) +
      coord_cartesian(xlim = c(0, plot_x_limit)) +
      ylab("Unique Sequences") + xlab("nucleotide position")

    # Save standalone plot
    block_graph_name <- paste0(output_prefix, "_", current_frame, "_stop-stop_ORF.png")
    ggsave(filename = block_graph_name, plot = block_graph, width = 10, height = 6, dpi = 300)
    message(paste("Saved individual plot: '", block_graph_name, "'", sep=""))

    # Store plot object
    plot_list[[current_frame]] <- block_graph
  }

  # 6. ASSEMBLE A4 PORTRAIT PANEL SHEET (PATCHWORK) -------------------------------------------
  message("Assembling consolidated A4 portrait panel layout...")

  p1 <- plot_list[["plus_0"]]
  p2 <- plot_list[["plus_1"]]
  p3 <- plot_list[["plus_2"]]
  p4 <- plot_list[["revcom_plus_0"]]
  p5 <- plot_list[["revcom_plus_1"]]
  p6 <- plot_list[["revcom_plus_2"]]

  # Eliminate repetitive interior labels
  p1 <- p1 + xlab(NULL); p4 <- p4 + xlab(NULL); p4 <- p4 + ylab(NULL)
  p2 <- p2 + xlab(NULL); p5 <- p5 + xlab(NULL); p5 <- p5 + ylab(NULL)
  p6 <- p6 + ylab(NULL)

  # Construct the 3x2 grid
  a4_composite_panel <- (p1 + p4) / (p2 + p5) / (p3 + p6)

  # Consolidate labels and place a single legend along the bottom-right frame margins
  a4_composite_panel <- a4_composite_panel +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.title = element_text(size = 9, face = "bold"),
          legend.text = element_text(size = 8))

  # Export at the exact dimensions of an A4 sheet
  panel_file_name <- paste0(output_prefix, "_A4_composite_panel.png")
  ggsave(filename = panel_file_name, plot = a4_composite_panel,
         width = 8.27, height = 11.69, dpi = 300)

  message(paste("Successfully exported panel sheet: '", panel_file_name, "'", sep=""))
}
