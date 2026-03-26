# ============================================================
# Script: build_phyloseq_zOTU_genus_prev2_multiNMDS.R
# Purpose:
#   1) Build a phyloseq object from zOTUs
#   2) Apply the most suitable prevalence filter: prev2
#   3) Collapse taxa at the genus level
#   4) Run NMDS, PERMANOVA, and betadisper
#   5) Allow multiple NMDS runs colored by different metadata variables
#   6) Automatically display figures in RStudio
#   7) Export all outputs to the output/ folder
#
# How to use:
#   - Edit only the "color_vars" vector near the end of the script
#     to choose which metadata variables will be used to color the NMDS.
#
# Notes:
#   - NMDS figures are automatically exported as PNG and PDF
#   - NMDS figures are also automatically printed to the RStudio Plots panel
#   - If a metadata variable does not exist, the script skips it with a warning
#
# Methodological note incorporated:
#   Previous analyses showed that OTU98 clustering did not sufficiently
#   remove year/platform-associated bias. The best compromise was obtained
#   using prevalence filtering >= 2 samples and genus-level taxonomic collapse,
#   rather than relying on OTU98 clustering.
# ============================================================

# ---------------------------
# 1. Working directories
# ---------------------------

args <- commandArgs(trailingOnly = TRUE)


input_dir  <- file.path(args[1])
output_dir <- file.path(args[2])
color_vars <- args[3]

if (is.na(color_vars)) {
    text_lines <- paste("",
                        "ATLASMX ITS — DB BUILDER V9.1",
                        "",
                        "build_phyloseq_zOTU_genus_prev2_multiNMDS.R",
                        "",
                        "Usage: Rscript build_phyloseq_zOTU_genus_prev2_multiNMDS.R <input_directory> <output_directory> <color_vars>",
                        "",
                        "Input directory",
                        '    String indicating the path of the input directory, for example "input"',
                        "Output directory",
                        '    String indicating the path of the output directory, for example "output"',
                        "Color variables",
                        '    String or list separated by commas to color NMDS, for example "ecoregion_WWF" or "ecoregion_WWF,vegetation_CONABIO"',
                        "",
                        "Example: Rscript build_phyloseq_zOTU_genus_prev2_multiNMDS.R input output ecoregion_WWF,vegetation_CONABIO",
                        "",
                        
                        sep = "\n")
    cat(text_lines)
    quit()
}
         
# process color vars
color_vars <- strsplit(color_vars, split = ",")
# color_vars <- c("año_atlas","ecoregion_WWF","vegetacion_CONABIO")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------
# 2. Required packages
# ---------------------------

library(phyloseq)
library(tidyverse)
library(data.table)
library(vegan)
library(ggplot2)


# ---------------------------
# 3. Input files
# ---------------------------
abundance_file <- file.path(input_dir, "abundance_table.csv")
traits_file    <- file.path(input_dir, "fungal_traits_table.csv")
taxonomy_file  <- file.path(input_dir, "taxonomy.csv")
metadata_file  <- file.path(input_dir, "sample_metadata.csv")

files_to_check <- c(abundance_file, traits_file, taxonomy_file, metadata_file)
missing_files <- files_to_check[!file.exists(files_to_check)]

if (length(missing_files) > 0) {
  stop("The following files are missing:\n", paste(missing_files, collapse = "\n"))
}

# ---------------------------
# 4. Read input tables
# ---------------------------
abundance_df <- fread(abundance_file, encoding = "UTF-8", data.table = FALSE)
traits_df    <- fread(traits_file, encoding = "UTF-8", data.table = FALSE)
taxonomy_df  <- fread(taxonomy_file, encoding = "UTF-8", data.table = FALSE)
metadata_df  <- fread(metadata_file, encoding = "UTF-8", data.table = FALSE)

# ---------------------------
# 5. Helper function to detect ID columns
# ---------------------------
detect_id_column <- function(
    df,
    candidate_names = c(
      "OTU", "otu", "OTU_ID", "otu_id", "FeatureID", "featureid",
      "feature", "ASV", "asv", "zOTU", "zotu", "SequenceID",
      "sequence_id", "sample_ID", "SampleID", "sampleid"
    )
) {
  cn <- colnames(df)
  hit <- cn[cn %in% candidate_names]
  if (length(hit) > 0) return(hit[1])
  return(cn[1])
}

otu_id_col_tax    <- detect_id_column(taxonomy_df)
otu_id_col_traits <- detect_id_column(traits_df)
otu_id_col_abund  <- detect_id_column(abundance_df)
sample_id_col_md  <- detect_id_column(
  metadata_df,
  candidate_names = c("sample_ID", "SampleID", "sampleid", "sample", "Sample")
)

message("Detected ID column in taxonomy: ", otu_id_col_tax)
message("Detected ID column in traits: ", otu_id_col_traits)
message("Detected ID column in abundance: ", otu_id_col_abund)
message("Detected ID column in metadata: ", sample_id_col_md)

# ---------------------------
# 6. Prepare abundance table
# ---------------------------
abundance_df[[otu_id_col_abund]] <- as.character(abundance_df[[otu_id_col_abund]])
rownames(abundance_df) <- abundance_df[[otu_id_col_abund]]
abundance_df[[otu_id_col_abund]] <- NULL

abundance_mat <- as.matrix(abundance_df)
mode(abundance_mat) <- "numeric"
abundance_mat[is.na(abundance_mat)] <- 0

# Remove completely empty zOTUs
abundance_mat <- abundance_mat[rowSums(abundance_mat) > 0, , drop = FALSE]

# ---------------------------
# 7. Prepare metadata
# ---------------------------
metadata_df[[sample_id_col_md]] <- as.character(metadata_df[[sample_id_col_md]])
rownames(metadata_df) <- metadata_df[[sample_id_col_md]]

# ---------------------------
# 8. Prepare taxonomy + traits
# ---------------------------
taxonomy_df[[otu_id_col_tax]] <- as.character(taxonomy_df[[otu_id_col_tax]])
traits_df[[otu_id_col_traits]] <- as.character(traits_df[[otu_id_col_traits]])

colnames(taxonomy_df)[colnames(taxonomy_df) == otu_id_col_tax] <- "OTU_ID"
colnames(traits_df)[colnames(traits_df) == otu_id_col_traits]  <- "OTU_ID"

shared_cols <- intersect(colnames(taxonomy_df), colnames(traits_df))
shared_cols <- setdiff(shared_cols, "OTU_ID")

if (length(shared_cols) > 0) {
  colnames(traits_df)[colnames(traits_df) %in% shared_cols] <- paste0(
    colnames(traits_df)[colnames(traits_df) %in% shared_cols],
    "_traits"
  )
}

otu_metadata_df <- taxonomy_df %>%
  left_join(traits_df, by = "OTU_ID")

rownames(otu_metadata_df) <- otu_metadata_df$OTU_ID
otu_metadata_df$OTU_ID <- NULL

otu_metadata_mat <- as.matrix(
  otu_metadata_df %>%
    mutate(across(everything(), as.character))
)

otu_metadata_mat[is.na(otu_metadata_mat)] <- "NA"

# ---------------------------
# 9. Intersect IDs across tables
# ---------------------------
common_otus <- Reduce(intersect, list(
  rownames(abundance_mat),
  rownames(otu_metadata_mat)
))

common_samples <- intersect(
  colnames(abundance_mat),
  rownames(metadata_df)
)

if (length(common_otus) == 0) {
  stop("No shared zOTUs were found between abundance_table and taxonomy/traits.")
}

if (length(common_samples) == 0) {
  stop("No shared samples were found between abundance_table and sample_metadata.")
}

message("Shared zOTUs: ", length(common_otus))
message("Shared samples: ", length(common_samples))

abundance_mat    <- abundance_mat[common_otus, common_samples, drop = FALSE]
otu_metadata_mat <- otu_metadata_mat[common_otus, , drop = FALSE]
metadata_df      <- metadata_df[common_samples, , drop = FALSE]

# ---------------------------
# 10. Build the base phyloseq object
# ---------------------------
OTU  <- otu_table(abundance_mat, taxa_are_rows = TRUE)
TAX  <- tax_table(otu_metadata_mat)
SAMP <- sample_data(metadata_df)

ps <- phyloseq(OTU, TAX, SAMP)
ps

# ---------------------------
# 11. Save the base phyloseq object
# ---------------------------
saveRDS(ps, file = file.path(output_dir, "ps_zotu.rds"))
save(ps, file = file.path(output_dir, "ps_zotu.RData"))

# ---------------------------
# 12. General summary of the base object
# ---------------------------
summary_lines <- c(
  "=======================================================",
  "Summary of the base phyloseq object (zOTUs)",
  "=======================================================",
  paste("Number of zOTUs:", ntaxa(ps)),
  paste("Number of samples:", nsamples(ps)),
  paste("Total read count:", sum(otu_table(ps))),
  "",
  "Reads per sample (summary):",
  capture.output(summary(sample_sums(ps))),
  "",
  "zOTU prevalence (number of samples where each zOTU is present):",
  capture.output(summary(apply(otu_table(ps), 1, function(x) sum(x > 0)))),
  "======================================================="
)

writeLines(summary_lines, con = file.path(output_dir, "00_phyloseq_summary_zotu.txt"))

# ---------------------------
# 13. Export base exploratory tables
# ---------------------------
reads_per_sample <- data.frame(
  sample_ID = sample_names(ps),
  total_reads = sample_sums(ps)
) %>%
  left_join(
    metadata_df %>% rownames_to_column("sample_ID"),
    by = "sample_ID"
  )

write.csv(
  reads_per_sample,
  file = file.path(output_dir, "01_reads_per_sample_zotu.csv"),
  row.names = FALSE
)

otu_prevalence <- data.frame(
  OTU_ID = taxa_names(ps),
  total_reads = taxa_sums(ps),
  prevalence = apply(otu_table(ps), 1, function(x) sum(x > 0))
)

otu_prevalence <- otu_prevalence %>%
  left_join(
    otu_metadata_df %>% rownames_to_column("OTU_ID"),
    by = "OTU_ID"
  )

write.csv(
  otu_prevalence,
  file = file.path(output_dir, "02_zotu_prevalence.csv"),
  row.names = FALSE
)

alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson")) %>%
  rownames_to_column("sample_ID") %>%
  left_join(
    metadata_df %>% rownames_to_column("sample_ID"),
    by = "sample_ID"
  )

write.csv(
  alpha_div,
  file = file.path(output_dir, "03_alpha_diversity_zotu.csv"),
  row.names = FALSE
)

tax_df_export <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE) %>%
  rownames_to_column("OTU_ID")

if ("phylum" %in% colnames(tax_df_export)) {
  phylum_abund <- data.frame(
    OTU_ID = taxa_names(ps),
    reads = taxa_sums(ps)
  ) %>%
    left_join(tax_df_export, by = "OTU_ID") %>%
    group_by(phylum) %>%
    summarise(total_reads = sum(reads), .groups = "drop") %>%
    arrange(desc(total_reads)) %>%
    mutate(relative_abundance = total_reads / sum(total_reads))
  
  write.csv(
    phylum_abund,
    file = file.path(output_dir, "04_abundance_by_phylum_zotu.csv"),
    row.names = FALSE
  )
}

possible_primary_cols <- c(
  "primary_lifestyle", "primary_lifestyle_traits", "primary",
  "Primary", "primary_lifestyle_FT"
)
found_primary_col <- possible_primary_cols[possible_primary_cols %in% colnames(tax_df_export)]

if (length(found_primary_col) > 0) {
  primary_col <- found_primary_col[1]
  
  lifestyle_abund <- data.frame(
    OTU_ID = taxa_names(ps),
    reads = taxa_sums(ps)
  ) %>%
    left_join(tax_df_export, by = "OTU_ID") %>%
    group_by(.data[[primary_col]]) %>%
    summarise(total_reads = sum(reads), .groups = "drop") %>%
    arrange(desc(total_reads)) %>%
    mutate(relative_abundance = total_reads / sum(total_reads))
  
  colnames(lifestyle_abund)[1] <- "primary_lifestyle"
  
  write.csv(
    lifestyle_abund,
    file = file.path(output_dir, "05_abundance_by_primary_lifestyle_zotu.csv"),
    row.names = FALSE
  )
}

# ---------------------------
# 14. Base exploratory plots
# ---------------------------
p_reads <- ggplot(reads_per_sample, aes(x = total_reads)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  labs(
    title = "Distribution of reads per sample (zOTUs)",
    x = "Reads per sample",
    y = "Frequency"
  )

ggsave(
  filename = file.path(output_dir, "06_hist_reads_per_sample_zotu.png"),
  plot = p_reads,
  width = 8,
  height = 5,
  dpi = 300
)

if ("sistema" %in% colnames(alpha_div)) {
  p_shannon <- ggplot(alpha_div, aes(x = sistema, y = Shannon)) +
    geom_boxplot() +
    theme_bw() +
    labs(
      title = "Shannon diversity by system (zOTUs)",
      x = "System",
      y = "Shannon"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = file.path(output_dir, "07_shannon_by_sistema_zotu.png"),
    plot = p_shannon,
    width = 8,
    height = 5,
    dpi = 300
  )
}

top20_otus <- otu_prevalence %>%
  arrange(desc(total_reads)) %>%
  slice(1:20)

p_top20 <- ggplot(top20_otus, aes(x = reorder(OTU_ID, total_reads), y = total_reads)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top 20 most abundant zOTUs",
    x = "zOTU",
    y = "Total reads"
  )

ggsave(
  filename = file.path(output_dir, "08_top20_zotus.png"),
  plot = p_top20,
  width = 8,
  height = 7,
  dpi = 300
)

if (exists("phylum_abund")) {
  p_phylum <- phylum_abund %>%
    slice(1:15) %>%
    ggplot(aes(x = reorder(phylum, total_reads), y = total_reads)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = "Top phyla by abundance (zOTUs)",
      x = "Phylum",
      y = "Total reads"
    )
  
  ggsave(
    filename = file.path(output_dir, "09_top_phyla_zotu.png"),
    plot = p_phylum,
    width = 8,
    height = 6,
    dpi = 300
  )
}

if (exists("lifestyle_abund")) {
  p_lifestyle <- lifestyle_abund %>%
    slice(1:15) %>%
    ggplot(aes(x = reorder(primary_lifestyle, total_reads), y = total_reads)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = "Top primary lifestyles by abundance (zOTUs)",
      x = "Primary lifestyle",
      y = "Total reads"
    )
  
  ggsave(
    filename = file.path(output_dir, "10_top_primary_lifestyles_zotu.png"),
    plot = p_lifestyle,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ---------------------------
# 15. Helper functions for filtering and cleaning
# ---------------------------
filter_phyloseq_for_beta <- function(ps_obj) {
  ps_obj <- prune_samples(sample_sums(ps_obj) > 0, ps_obj)
  ps_obj <- prune_taxa(taxa_sums(ps_obj) > 0, ps_obj)
  return(ps_obj)
}

filter_by_prevalence_count <- function(ps_obj, min_samples = 2) {
  prev <- apply(otu_table(ps_obj), 1, function(x) sum(x > 0))
  ps_f <- prune_taxa(prev >= min_samples, ps_obj)
  ps_f <- filter_phyloseq_for_beta(ps_f)
  return(ps_f)
}

clean_tax_rank <- function(x, unknown_label) {
  x <- as.character(x)
  x[is.na(x)] <- unknown_label
  x[x == ""] <- unknown_label
  x[x == "NA"] <- unknown_label
  x[x == "NaN"] <- unknown_label
  x
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  return(x)
}

# ---------------------------
# 16. Methodological decision
# ---------------------------
# Previous analyses indicated that the most useful prevalence filter was:
#   presence in at least 2 samples (prev2)
#
# Additionally, it was observed that:
#   - OTU98 clustering did not sufficiently remove year/platform bias
#   - genus-level collapse reduced technical heterogeneity better than
#     family-level collapse, without losing as much biological structure
#
# Therefore, the main workflow adopted here is:
#   zOTUs -> prev2 filter -> genus-level collapse

# ---------------------------
# 17. Create ps_prev2 (zOTUs with prevalence >= 2)
# ---------------------------
ps_prev2 <- filter_by_prevalence_count(ps, min_samples = 2)

saveRDS(ps_prev2, file = file.path(output_dir, "ps_prev2_zotu.rds"))
save(ps_prev2, file = file.path(output_dir, "ps_prev2_zotu.RData"))

# ---------------------------
# 18. Collapse ps_prev2 at the genus level
# ---------------------------
tax_table(ps_prev2)[, "genus"] <- clean_tax_rank(
  tax_table(ps_prev2)[, "genus"],
  "Unassigned_Genus"
)

ps_genus_prev2 <- tax_glom(ps_prev2, taxrank = "genus", NArm = FALSE)
ps_genus_prev2 <- filter_phyloseq_for_beta(ps_genus_prev2)

saveRDS(ps_genus_prev2, file = file.path(output_dir, "ps_genus_prev2_from_zotu.rds"))
save(ps_genus_prev2, file = file.path(output_dir, "ps_genus_prev2_from_zotu.RData"))

# ---------------------------
# 19. Summary of the main working objects
# ---------------------------
workflow_summary <- tibble(
  dataset = c("ps_zotu", "ps_prev2_zotu", "ps_genus_prev2_from_zotu"),
  n_samples = c(nsamples(ps), nsamples(ps_prev2), nsamples(ps_genus_prev2)),
  n_taxa = c(ntaxa(ps), ntaxa(ps_prev2), ntaxa(ps_genus_prev2)),
  total_reads = c(sum(otu_table(ps)), sum(otu_table(ps_prev2)), sum(otu_table(ps_genus_prev2)))
)

write.csv(
  workflow_summary,
  file = file.path(output_dir, "11_workflow_summary_zotu_to_genus.csv"),
  row.names = FALSE
)

# ---------------------------
# 20. NMDS function
# ---------------------------
run_nmds_phyloseq <- function(ps_obj, dataset_name, output_dir, color_by = "año_atlas") {
  
  ps_obj <- filter_phyloseq_for_beta(ps_obj)
  
  cat("\n====================================\n")
  cat("Running NMDS for:", dataset_name, "\n")
  cat("Samples:", nsamples(ps_obj), "\n")
  cat("Taxa:", ntaxa(ps_obj), "\n")
  cat("Color variable:", color_by, "\n")
  cat("====================================\n")
  
  if (nsamples(ps_obj) < 3) {
    warning(paste0("Dataset ", dataset_name, ": fewer than 3 samples after filtering. NMDS skipped."))
    return(NULL)
  }
  
  ps_relabund <- transform_sample_counts(ps_obj, function(x) {
    if (sum(x) == 0) return(x)
    x / sum(x)
  })
  
  zero_samples_after <- sample_sums(ps_relabund) == 0
  if (any(zero_samples_after)) {
    warning(paste0("Dataset ", dataset_name, ": empty samples still present after transformation. NMDS skipped."))
    return(NULL)
  }
  
  ord_nmds <- tryCatch({
    ordinate(
      ps_relabund,
      method = "NMDS",
      distance = "bray",
      trymax = 300
    )
  }, error = function(e) {
    warning(paste0("NMDS failed for ", dataset_name, ": ", e$message))
    return(NULL)
  })
  
  if (is.null(ord_nmds)) return(NULL)
  
  ord_scores <- as.data.frame(vegan::scores(ord_nmds, display = "sites"))
  ord_scores$sample_ID <- rownames(ord_scores)
  
  if (ncol(ord_scores) < 2) {
    warning(paste0("Could not extract two ordination axes for ", dataset_name))
    return(NULL)
  }
  
  colnames(ord_scores)[1:2] <- c("NMDS1", "NMDS2")
  
  metadata_df_nmds <- as(sample_data(ps_obj), "data.frame")
  metadata_df_nmds$sample_ID <- rownames(metadata_df_nmds)
  
  ord_df <- ord_scores %>%
    left_join(metadata_df_nmds, by = "sample_ID")
  
  safe_color_by <- sanitize_filename(color_by)
  
  write.csv(
    ord_df,
    file = file.path(output_dir, paste0("NMDS_coordinates_", dataset_name, "_colored_by_", safe_color_by, ".csv")),
    row.names = FALSE
  )
  
  if (color_by %in% colnames(ord_df)) {
    
    ord_df[[color_by]] <- as.factor(ord_df[[color_by]])
    
    # Hide the legend automatically when too many categories are present
    n_levels <- length(unique(ord_df[[color_by]]))
    legend_pos <- ifelse(n_levels > 20, "none", "right")
    
    p_nmds <- ggplot(
      ord_df,
      aes(x = NMDS1, y = NMDS2, color = .data[[color_by]])
    ) +
      geom_point(size = 3, alpha = 0.85) +
      theme_bw() +
      labs(
        title = paste0("NMDS (Bray-Curtis) - ", dataset_name),
        subtitle = paste0(
          "Stress = ", round(ord_nmds$stress, 4),
          " | Colored by: ", color_by
        ),
        x = "NMDS1",
        y = "NMDS2",
        color = color_by
      ) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.title = element_text(face = "bold"),
        legend.position = legend_pos
      )
    
  } else {
    
    warning(paste0("Variable '", color_by, "' was not found in metadata for ", dataset_name, ". Figure will be exported without valid color mapping."))
    
    p_nmds <- ggplot(ord_df, aes(x = NMDS1, y = NMDS2)) +
      geom_point(size = 3, alpha = 0.85) +
      theme_bw() +
      labs(
        title = paste0("NMDS (Bray-Curtis) - ", dataset_name),
        subtitle = paste0(
          "Stress = ", round(ord_nmds$stress, 4),
          " | No valid color variable found"
        ),
        x = "NMDS1",
        y = "NMDS2"
      ) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = "black")
      )
  }
  
  # Automatically display the figure in RStudio
  print(p_nmds)
  
  # Automatically export the figure
  ggsave(
    filename = file.path(output_dir, paste0("NMDS_bray_", dataset_name, "_colored_by_", safe_color_by, ".png")),
    plot = p_nmds,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, paste0("NMDS_bray_", dataset_name, "_colored_by_", safe_color_by, ".pdf")),
    plot = p_nmds,
    width = 8,
    height = 6
  )
  
  return(data.frame(
    dataset = dataset_name,
    color_by = color_by,
    stress = ord_nmds$stress,
    n_samples = nsamples(ps_obj),
    n_taxa = ntaxa(ps_obj),
    stringsAsFactors = FALSE
  ))
}

# ---------------------------
# 21. PERMANOVA + betadisper function
# ---------------------------
run_permanova_betadisper <- function(ps_obj, dataset_name, output_dir, group_var = "año_atlas", permutations = 999) {
  
  ps_obj <- filter_phyloseq_for_beta(ps_obj)
  
  if (nsamples(ps_obj) < 3) {
    warning(paste0("Dataset ", dataset_name, ": fewer than 3 samples. PERMANOVA/betadisper skipped."))
    return(NULL)
  }
  
  metadata_df_stat <- as(sample_data(ps_obj), "data.frame")
  metadata_df_stat$sample_ID <- rownames(metadata_df_stat)
  
  if (!group_var %in% colnames(metadata_df_stat)) {
    warning(paste0("Variable ", group_var, " does not exist in ", dataset_name, "."))
    return(NULL)
  }
  
  metadata_df_stat[[group_var]] <- as.factor(metadata_df_stat[[group_var]])
  
  if (length(unique(metadata_df_stat[[group_var]])) < 2) {
    warning(paste0("Variable ", group_var, " in ", dataset_name, " has fewer than 2 levels."))
    return(NULL)
  }
  
  ps_relabund <- transform_sample_counts(ps_obj, function(x) {
    if (sum(x) == 0) return(x)
    x / sum(x)
  })
  
  dist_mat <- tryCatch({
    phyloseq::distance(ps_relabund, method = "bray")
  }, error = function(e) {
    warning(paste0("Could not calculate Bray distance for ", dataset_name, ": ", e$message))
    return(NULL)
  })
  
  if (is.null(dist_mat)) return(NULL)
  
  permanova_res <- tryCatch({
    vegan::adonis2(
      dist_mat ~ metadata_df_stat[[group_var]],
      data = metadata_df_stat,
      permutations = permutations
    )
  }, error = function(e) {
    warning(paste0("PERMANOVA failed for ", dataset_name, ": ", e$message))
    return(NULL)
  })
  
  betadisper_res <- tryCatch({
    bd <- vegan::betadisper(dist_mat, metadata_df_stat[[group_var]])
    perm <- vegan::permutest(bd, permutations = permutations)
    list(betadisper = bd, permutest = perm)
  }, error = function(e) {
    warning(paste0("betadisper failed for ", dataset_name, ": ", e$message))
    return(NULL)
  })
  
  safe_group_var <- sanitize_filename(group_var)
  
  if (!is.null(permanova_res)) {
    permanova_df <- as.data.frame(permanova_res)
    permanova_df$term <- rownames(permanova_df)
    
    write.csv(
      permanova_df,
      file = file.path(output_dir, paste0("PERMANOVA_", dataset_name, "_by_", safe_group_var, ".csv")),
      row.names = FALSE
    )
  }
  
  if (!is.null(betadisper_res)) {
    bd_scores <- data.frame(
      sample_ID = names(betadisper_res$betadisper$distances),
      distance_to_centroid = betadisper_res$betadisper$distances,
      group = metadata_df_stat[[group_var]]
    )
    
    write.csv(
      bd_scores,
      file = file.path(output_dir, paste0("betadisper_distances_", dataset_name, "_by_", safe_group_var, ".csv")),
      row.names = FALSE
    )
    
    p_bd <- ggplot(bd_scores, aes(x = group, y = distance_to_centroid)) +
      geom_boxplot() +
      theme_bw() +
      labs(
        title = paste0("betadisper - ", dataset_name),
        subtitle = paste0("Grouped by: ", group_var),
        x = group_var,
        y = "Distance to centroid"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      filename = file.path(output_dir, paste0("betadisper_boxplot_", dataset_name, "_by_", safe_group_var, ".png")),
      plot = p_bd,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  permanova_summary <- NULL
  if (!is.null(permanova_res)) {
    permanova_tab <- as.data.frame(permanova_res)
    permanova_summary <- data.frame(
      dataset = dataset_name,
      method = "PERMANOVA",
      variable = group_var,
      R2 = permanova_tab$R2[1],
      F = permanova_tab$F[1],
      p_value = permanova_tab$`Pr(>F)`[1]
    )
  }
  
  betadisper_summary <- NULL
  if (!is.null(betadisper_res)) {
    perm_tab <- as.data.frame(betadisper_res$permutest$tab)
    betadisper_summary <- data.frame(
      dataset = dataset_name,
      method = "betadisper",
      variable = group_var,
      R2 = NA,
      F = perm_tab$F[1],
      p_value = perm_tab$`Pr(>F)`[1]
    )
  }
  
  return(bind_rows(permanova_summary, betadisper_summary))
}

# ---------------------------
# 22. Run the base NMDS colored by año_atlas
# ---------------------------
color_by <- "año_atlas"

nmds_results <- bind_rows(
  run_nmds_phyloseq(ps_prev2, "ps_prev2_zotu", output_dir, color_by = color_by),
  run_nmds_phyloseq(ps_genus_prev2, "ps_genus_prev2_from_zotu", output_dir, color_by = color_by)
)

write.csv(
  nmds_results,
  file = file.path(output_dir, "12_NMDS_summary_zotu_vs_genus.csv"),
  row.names = FALSE
)

# ---------------------------
# 23. Run base PERMANOVA + betadisper by año_atlas
# ---------------------------
group_var <- "año_atlas"

stat_results <- bind_rows(
  run_permanova_betadisper(ps_prev2, "ps_prev2_zotu", output_dir, group_var = group_var),
  run_permanova_betadisper(ps_genus_prev2, "ps_genus_prev2_from_zotu", output_dir, group_var = group_var)
)

write.csv(
  stat_results,
  file = file.path(output_dir, "13_PERMANOVA_betadisper_summary_zotu_vs_genus.csv"),
  row.names = FALSE
)

# ---------------------------
# 24. Export collapsed tables for inspection
# ---------------------------
genus_abund_table <- as.data.frame(otu_table(ps_genus_prev2)) %>%
  rownames_to_column("genus")

write.csv(
  genus_abund_table,
  file = file.path(output_dir, "14_genus_abundance_table_from_zotu.csv"),
  row.names = FALSE
)

genus_tax_table <- as.data.frame(tax_table(ps_genus_prev2), stringsAsFactors = FALSE)
genus_tax_table$Collapsed_Label <- rownames(genus_tax_table)

write.csv(
  genus_tax_table,
  file = file.path(output_dir, "15_genus_tax_table_from_zotu.csv"),
  row.names = TRUE
)

write.csv(
  as(sample_data(ps_prev2), "data.frame") %>% rownames_to_column("sample_ID"),
  file = file.path(output_dir, "16_sample_data_ps_prev2_zotu.csv"),
  row.names = FALSE
)

write.csv(
  as(sample_data(ps_genus_prev2), "data.frame") %>% rownames_to_column("sample_ID"),
  file = file.path(output_dir, "17_sample_data_ps_genus_prev2_from_zotu.csv"),
  row.names = FALSE
)

# ---------------------------
# 25. Save the tax_table and sample_data used
# ---------------------------
write.csv(
  as.data.frame(otu_metadata_mat) %>% rownames_to_column("OTU_ID"),
  file = file.path(output_dir, "18_taxonomy_traits_used_in_phyloseq_zotu.csv"),
  row.names = FALSE
)

write.csv(
  metadata_df %>% rownames_to_column("sample_ID"),
  file = file.path(output_dir, "19_sample_metadata_used_in_phyloseq_zotu.csv"),
  row.names = FALSE
)

# ---------------------------
# 26. Run multiple NMDS analyses by different metadata variables
# ---------------------------
# ============================================================
# Add or remove variables depending on how you want to color NMDS
# ============================================================


multi_nmds_results <- list()

for (cv in color_vars) {
  cat("\n####################################################\n")
  cat("Running multi-NMDS block for variable:", cv, "\n")
  cat("####################################################\n")
  
  res_prev2 <- run_nmds_phyloseq(
    ps_prev2,
    "ps_prev2_zotu",
    output_dir,
    color_by = cv
  )
  
  res_genus <- run_nmds_phyloseq(
    ps_genus_prev2,
    "ps_genus_prev2_from_zotu",
    output_dir,
    color_by = cv
  )
  
  multi_nmds_results[[cv]] <- bind_rows(res_prev2, res_genus)
}

multi_nmds_summary <- bind_rows(multi_nmds_results)

write.csv(
  multi_nmds_summary,
  file = file.path(output_dir, "20_multiNMDS_summary_by_metadata.csv"),
  row.names = FALSE
)

# ---------------------------
# 27. Final message
# ---------------------------
cat("\nProcess completed successfully.\n")
cat("Main objects saved in:\n")
cat(" - ", file.path(output_dir, "ps_zotu.rds"), "\n")
cat(" - ", file.path(output_dir, "ps_prev2_zotu.rds"), "\n")
cat(" - ", file.path(output_dir, "ps_genus_prev2_from_zotu.rds"), "\n")
cat("The following outputs were also exported:\n")
cat(" - base exploratory summaries for zOTUs\n")
cat(" - base NMDS colored by año_atlas\n")
cat(" - PERMANOVA and betadisper by año_atlas\n")
cat(" - multiple NMDS plots for all variables listed in color_vars\n")
cat(" - genus-level collapsed abundance and taxonomy tables\n")
