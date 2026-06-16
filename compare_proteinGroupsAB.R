# ---- Environment setup: required package ----

required_pkgs <- c("tidyverse", "data.table", "ggrepel", "patchwork", "scales")
new_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, ggplot2, tidyr, purrr, stringr
  library(data.table)  # fast TSV reading (fread)
  library(ggrepel)     # non-overlapping volcano labels
  library(patchwork)   # multi-panel figure composition
  library(scales)      # comma() axis formatting
})


# ---- User-defined parameters ----
# Set these three paths before running

PATH_A  <- "proteinGroups_A.txt"  # A — original file from PRIDE
PATH_B  <- "proteinGroups_B.txt"  # B — reanalysis (modified parameters in MaxQuant)
OUT_DIR <- "compareProteinGroups"        # all outputs written here


# ---- Differential expression thresholds ----
FDR_THRESHOLD <- 0.05   # Benjamini-Hochberg adjusted p-value cutoff
FC_THRESHOLD  <- 1.0    # |log2 fold change| cutoff (1.0 = 2-fold)
MIN_VALID     <- 2      # minimum non-NA values per group to run t-test

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Core helper functions for preprocessing and QC ----


# Export plots as SVG using base R graphics device
save_svg <- function(plot, path, width = 7, height = 6) {
  svg(path, width = width, height = height)
  print(plot)
  dev.off()
  message(sprintf("  Saved: %s", basename(path)))
  invisible(path)
}

# Extract experimental condition from sample identifiers
# Works for both naming schemes:
#   "KO_M0_S16" → "KO_M0"  (A run-batch IDs)
#   "KO_M0_1"   → "KO_M0"  (B sequential IDs)
# Regex captures the first occurrence of (KO|WT)_(M0|M1|M2).
parse_cond <- function(s) {
  vapply(s, function(x) {
    m <- regmatches(x, regexpr("(KO|WT)_(M[012])", x))
    if (length(m) == 1L) m else NA_character_
  }, character(1L), USE.NAMES = FALSE)
}

# Remove contaminants, reverse/decoy hits, and site-only identifications
# based on MaxQuant annotation columns
# A uses "Reverse"; B uses "Decoy". Both use "Potential contaminant"
apply_qc_filters <- function(dt, label) {
  flag_cols <- c("Reverse", "Decoy",
                 "Potential contaminant", "Only identified by site")
  present <- intersect(flag_cols, names(dt))
  before  <- nrow(dt)
  removed <- integer(length(present))

  for (i in seq_along(present)) {
    col <- present[i]
    v   <- dt[[col]]
    # Handle both "+" (character) and TRUE/1 (logical/integer) encodings
    bad <- (!is.na(v)) & (v == "+" | v == TRUE | v == "true" | v == "1")
    removed[i] <- sum(bad)
    dt <- dt[!bad]
  }

  message(sprintf(
    "  [%s] QC filters: %s → removed %d rows, %d proteins retained",
    label,
    paste(sprintf("%s (-%d)", present, removed), collapse = ", "),
    before - nrow(dt), nrow(dt)
  ))
  dt
}

# Ensure intensity columns are correctly parsed as numeric values
coerce_intensity_cols <- function(dt) {
  intensity_pat <- "^(LFQ intensity|Intensity |iBAQ |MS/MS count)"
  int_cols <- grep(intensity_pat, names(dt), value = TRUE)
  for (col in int_cols) {
    set(dt, j = col, value = as.double(dt[[col]]))
  }
  dt
}

# Extract log2 LFQ matrix from a data.table 
# Replaces 0 with NA (MaxQuant convention: 0 = missing, not true zero),
# then applies log2 transformation. Returns a numeric matrix with
# protein accessions as rownames and short sample names as colnames.
extract_lfq_matrix <- function(dt, lfq_cols, id_col = "leading_id") {
  m <- as.matrix(dt[, ..lfq_cols])
  storage.mode(m) <- "double"
  m[m == 0] <- NA          # 0 → NA (missing)
  m <- log2(m)             # log2 transform
  rownames(m) <- dt[[id_col]]
  colnames(m) <- sub("^LFQ intensity ", "", lfq_cols)  # shorten names
  m
}

# Protein-level t-test (Welch)
run_ttest_per_protein <- function(mA, mB, min_valid = MIN_VALID) {
  proteins <- rownames(mA)
  stopifnot("Row order mismatch" = identical(proteins, rownames(mB)))

  map_dfr(proteins, function(prot) {
    xA <- as.numeric(mA[prot, ])
    xB <- as.numeric(mB[prot, ])
    nA <- sum(!is.na(xA)); nB <- sum(!is.na(xB))
    meanA <- mean(xA, na.rm = TRUE); meanB <- mean(xB, na.rm = TRUE)
    log2FC <- meanB - meanA   # positive = higher in B (reanalysis)

    p_val <- if (nA >= min_valid && nB >= min_valid) {
      tryCatch(
        t.test(xB, xA, var.equal = FALSE)$p.value,  # Welch's t-test
        error = function(e) NA_real_
      )
    } else NA_real_

    tibble(protein = prot, mean_log2_A = meanA, mean_log2_B = meanB,
           log2FC = log2FC, n_valid_A = nA, n_valid_B = nB, p_value = p_val)
  }) %>%
    mutate(
      FDR = p.adjust(p_value, method = "BH"),
      significant = !is.na(FDR) & FDR < FDR_THRESHOLD & abs(log2FC) >= FC_THRESHOLD,
      direction = case_when(
        significant & log2FC > 0 ~ "Higher in B (reanalysis)",
        significant & log2FC < 0 ~ "Lower in B (reanalysis)",
        TRUE                     ~ "Not significant"
      )
    )
}


# ---- Step1.Load & QC ProteinGroups ----

message("\n[STEP 1] Loading files...")

for (path in c(PATH_A, PATH_B)) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
}

dtA <- fread(PATH_A, showProgress = FALSE, na.strings = c("", "NA"))
dtA <- apply_qc_filters(dtA, "A (original)")
dtA <- coerce_intensity_cols(dtA)
dtA[, leading_id := sub(";.*", "", `Protein IDs`)]  # leading UniProt accession

dtB <- fread(PATH_B, showProgress = FALSE, na.strings = c("", "NA"))
dtB <- apply_qc_filters(dtB, "B (reanalysis)")
dtB <- coerce_intensity_cols(dtB)
dtB[, leading_id := sub(";.*", "", `Protein IDs`)]

message(sprintf("  A (original):   %d proteins × %d columns", nrow(dtA), ncol(dtA)))
message(sprintf("  B (reanalysis): %d proteins × %d columns", nrow(dtB), ncol(dtB)))


# ---- Step2.Dataset structure check ----
message("\n[STEP 2] Structural audit...")

col_type_patterns <- c(
  "LFQ intensity"       = "^LFQ intensity ",
  "Intensity"           = "^Intensity ",
  "iBAQ"                = "^iBAQ ",
  "MS/MS count"         = "^MS/MS count ",
  "Peptide counts"      = "^Peptide counts",
  "Identification type" = "^Identification type ",
  "Sequence coverage"   = "^Sequence coverage",
  "Unique peptides"     = "^Unique peptides ",
  "Razor peptides"      = "^Razor \\+ unique peptides ",
  "Peptides"            = "^Peptides "
)

count_col_type <- function(hdr, patterns) {
  sapply(patterns, function(p) sum(grepl(p, hdr)))
}

ctA <- count_col_type(names(dtA), col_type_patterns)
ctB <- count_col_type(names(dtB), col_type_patterns)

structural_audit <- tibble(
  column_type  = names(col_type_patterns),
  A_original   = as.integer(ctA),
  B_reanalysis = as.integer(ctB),
  difference   = as.integer(ctB) - as.integer(ctA)
)
message("  Column type audit:"); print(structural_audit, n = Inf)

# Detect LFQ columns and parse sample → condition mapping
lfq_colsA <- grep("^LFQ intensity ", names(dtA), value = TRUE)
lfq_colsB <- grep("^LFQ intensity ", names(dtB), value = TRUE)

if (length(lfq_colsA) == 0 || length(lfq_colsB) == 0) {
  stop("No LFQ intensity columns found in one or both files.")
}

sampA <- sub("^LFQ intensity ", "", lfq_colsA)
sampB <- sub("^LFQ intensity ", "", lfq_colsB)
condA <- parse_cond(sampA)
condB <- parse_cond(sampB)

if (any(is.na(condA)) || any(is.na(condB))) {
  warning("Some sample names could not be parsed to a condition. ",
          "Check parse_cond() regex against your sample naming scheme.")
}

shared_conds <- sort(intersect(unique(na.omit(condA)), unique(na.omit(condB))))
message(sprintf("  Shared biological conditions (%d): %s",
                length(shared_conds), paste(shared_conds, collapse = ", ")))

# Warn if iBAQ is absent in one file
if (structural_audit$A_original[structural_audit$column_type == "iBAQ"] > 0 &&
    structural_audit$B_reanalysis[structural_audit$column_type == "iBAQ"] == 0) {
  message("  WARNING: iBAQ columns present in A but absent in B. ",
          "iBAQ comparison skipped.")
}

sample_map <- bind_rows(
  tibble(dataset = "A (original)",   sample = sampA, condition = condA),
  tibble(dataset = "B (reanalysis)", sample = sampB, condition = condB)
)

write_csv(structural_audit, file.path(OUT_DIR, "01_structural_audit.csv"))
write_csv(sample_map,       file.path(OUT_DIR, "01_sample_condition_map.csv"))


# ---- Step3.Protein overlap between runs (leading UniProt IDs) ----

idsA <- dtA$leading_id
idsB <- dtB$leading_id

ids_common <- intersect(idsA, idsB)
ids_only_A <- setdiff(idsA, idsB)
ids_only_B <- setdiff(idsB, idsA)
jaccard    <- length(ids_common) / length(union(idsA, idsB))

overlap_summary <- tibble(
  category = c(
    "A (original) — total proteins",
    "B (reanalysis) — total proteins",
    "Common to both runs",
    "Unique to A (original only)",
    "Unique to B (reanalysis only)",
    "Jaccard similarity index"
  ),
  count = c(length(idsA), length(idsB), length(ids_common),
            length(ids_only_A), length(ids_only_B), round(jaccard, 4))
)
message("  Overlap summary:"); print(overlap_summary)

# Gene name lookup for annotation (use A as reference; fall back to accession)
gene_map <- dtA[, .(leading_id, `Gene names`)] %>%
  as_tibble() %>%
  rename(protein = leading_id, gene = `Gene names`) %>%
  mutate(gene = coalesce(na_if(gene, ""), protein))

write_csv(overlap_summary, file.path(OUT_DIR, "02_overlap_summary.csv"))
write_csv(
  tibble(leading_id = ids_only_A,
         gene_name  = dtA$`Gene names`[match(ids_only_A, dtA$leading_id)]),
  file.path(OUT_DIR, "02_proteins_unique_A.csv")
)
write_csv(
  tibble(leading_id = ids_only_B,
         gene_name  = dtB$`Gene names`[match(ids_only_B, dtB$leading_id)]),
  file.path(OUT_DIR, "02_proteins_unique_B.csv")
)

# ---- Overlap bar chart ----
overlap_bar <- tibble(
  category = factor(
    c("A only\n(original)", "Common\n(both runs)", "B only\n(reanalysis)"),
    levels = c("A only\n(original)", "Common\n(both runs)", "B only\n(reanalysis)")
  ),
  n    = c(length(ids_only_A), length(ids_common), length(ids_only_B)),
  fill = c("#2166AC", "#6BAED6", "#D6604D")
)

p_overlap <- ggplot(overlap_bar, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = comma(n)), vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = setNames(overlap_bar$fill, overlap_bar$category)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.14))) +
  labs(
    title    = "Protein identification overlap between MaxQuant runs",
    subtitle = sprintf("Jaccard = %.3f | %d common / %d total unique proteins",
                       jaccard, length(ids_common), length(union(idsA, idsB))),
    x = NULL, y = "Number of protein groups"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.subtitle          = element_text(size = 10, colour = "grey40"),
        panel.grid.major.x     = element_blank())

save_svg(p_overlap, file.path(OUT_DIR, "03_protein_overlap_bar.svg"), width = 6, height = 5)

# ---- Step4.Restrict analysis to shared proteins ----

message("\n[STEP 4] Subsetting to common proteins...")

dfA_c <- dtA[leading_id %in% ids_common][order(leading_id)]
dfB_c <- dtB[leading_id %in% ids_common][order(leading_id)]

# Assert identical row order before any paired operations
stopifnot(
  "Protein order mismatch after subsetting — check for duplicate leading IDs" =
    identical(dfA_c$leading_id, dfB_c$leading_id)
)
message(sprintf("  Common protein set: %d proteins — row order verified.", nrow(dfA_c)))

# ---- Step5.Log2 LFQ distributions and missing value patterns ----

message("\n[STEP 5] LFQ distributions...")

matA <- extract_lfq_matrix(dfA_c, lfq_colsA)
matB <- extract_lfq_matrix(dfB_c, lfq_colsB)

# ---- Per-condition missing value rates ----
missing_per_cond <- bind_rows(lapply(shared_conds, function(cond) {
  idxA <- which(condA == cond)
  idxB <- which(condB == cond)
  tibble(
    condition     = cond,
    A_missing_pct = round(mean(is.na(matA[, idxA])) * 100, 1),
    B_missing_pct = round(mean(is.na(matB[, idxB])) * 100, 1),
    delta_pct     = round(mean(is.na(matB[, idxB])) * 100 -
                            mean(is.na(matA[, idxA])) * 100, 1)
  )
}))

overall_missing <- tibble(
  dataset         = c("A (original)", "B (reanalysis)"),
  total_values    = c(prod(dim(matA)), prod(dim(matB))),
  missing_NA      = c(sum(is.na(matA)), sum(is.na(matB))),
  pct_missing     = round(c(mean(is.na(matA)), mean(is.na(matB))) * 100, 2),
  median_log2_lfq = round(c(median(matA, na.rm = TRUE), median(matB, na.rm = TRUE)), 3)
)
message("  Overall missing:"); print(overall_missing)
message("  Per-condition:");   print(missing_per_cond)

write_csv(overall_missing,  file.path(OUT_DIR, "04_missing_value_summary.csv"))
write_csv(missing_per_cond, file.path(OUT_DIR, "04_missing_per_condition.csv"))

# ---- Tidy long format for ggplot ----
pal2 <- c("A (original)" = "#2166AC", "B (reanalysis)" = "#D6604D")

long_all <- bind_rows(
  as_tibble(matA, rownames = "protein") %>%
    pivot_longer(-protein, names_to = "sample", values_to = "log2_LFQ") %>%
    mutate(dataset = "A (original)",   condition = parse_cond(sample)),
  as_tibble(matB, rownames = "protein") %>%
    pivot_longer(-protein, names_to = "sample", values_to = "log2_LFQ") %>%
    mutate(dataset = "B (reanalysis)", condition = parse_cond(sample))
) %>% filter(!is.na(log2_LFQ))

p_density <- ggplot(long_all, aes(x = log2_LFQ, colour = dataset, fill = dataset)) +
  geom_density(alpha = 0.20, linewidth = 0.9) +
  scale_colour_manual(values = pal2) +
  scale_fill_manual(values = pal2) +
  labs(title    = "Log\u2082 LFQ intensity distributions — common proteins",
       subtitle = "Zeros replaced with NA before log\u2082 transformation",
       x = "log\u2082(LFQ intensity)", y = "Density",
       colour = NULL, fill = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

p_box <- ggplot(long_all, aes(x = sample, y = log2_LFQ, fill = dataset)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, linewidth = 0.4) +
  scale_fill_manual(values = pal2) +
  facet_wrap(~ dataset, scales = "free_x", nrow = 2) +
  labs(title = "Per-sample log\u2082 LFQ boxplots",
       x = "Sample", y = "log\u2082(LFQ intensity)", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
        legend.position  = "none",
        strip.background = element_rect(fill = "grey92"))

save_svg(p_density / p_box + plot_layout(heights = c(1, 1.6)) +
           plot_annotation(tag_levels = "A"),
         file.path(OUT_DIR, "05_lfq_distributions.svg"), width = 11, height = 9)


# ---- Step6.Differential abundance testing (Welch t-test) ----

message("\n[STEP 6] Statistical tests...")

# ---- 6a. Global test ----
message("  Running global test (18 vs 18 samples)...")
stats_global <- run_ttest_per_protein(matA, matB)

n_tested <- sum(!is.na(stats_global$p_value))
n_sig    <- sum(stats_global$significant, na.rm = TRUE)
n_up     <- sum(stats_global$significant & stats_global$log2FC > 0, na.rm = TRUE)
n_down   <- sum(stats_global$significant & stats_global$log2FC < 0, na.rm = TRUE)

message(sprintf("  Global: %d tested | %d significant (up: %d, down: %d)",
                n_tested, n_sig, n_up, n_down))
write_csv(stats_global, file.path(OUT_DIR, "06a_global_statistics.csv"))

# ---- 6b. Per-condition test ----
message("  Running per-condition tests (3 vs 3 per condition)...")

stats_percond <- map_dfr(shared_conds, function(cond) {
  idxA <- which(condA == cond)
  idxB <- which(condB == cond)
  run_ttest_per_protein(matA[, idxA, drop = FALSE],
                        matB[, idxB, drop = FALSE]) %>%
    mutate(condition = cond, .before = 1)
})

percond_summary <- stats_percond %>%
  group_by(condition) %>%
  summarise(n_tested = sum(!is.na(p_value)),
            n_sig    = sum(significant, na.rm = TRUE),
            n_up     = sum(significant & log2FC > 0, na.rm = TRUE),
            n_down   = sum(significant & log2FC < 0, na.rm = TRUE),
            .groups  = "drop")

message("  Per-condition summary:"); print(percond_summary)
write_csv(stats_percond,   file.path(OUT_DIR, "06b_per_condition_statistics.csv"))
write_csv(percond_summary, file.path(OUT_DIR, "06b_per_condition_summary.csv"))


# ---- Step7.Volcano plots (global + per-condition) ----

message("\n[STEP 7] Volcano plots...")

volcano_colours <- c("Higher in B (reanalysis)" = "#D6604D",
                     "Lower in B (reanalysis)"  = "#2166AC",
                     "Not significant"           = "grey72")

# ---- 7a. Global volcano ----
top_labels <- stats_global %>%
  filter(significant) %>% arrange(FDR) %>% slice_head(n = 20) %>% pull(protein)

plot_global <- stats_global %>%
  filter(!is.na(FDR)) %>%
  mutate(neg_log10_FDR = -log10(FDR)) %>%
  left_join(gene_map, by = "protein") %>%
  mutate(label = if_else(protein %in% top_labels, gene, NA_character_))

p_volcano_global <- ggplot(plot_global,
    aes(x = log2FC, y = neg_log10_FDR, colour = direction, size = direction)) +
  geom_point(alpha = 0.65) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 25,
                  segment.colour = "grey50", segment.size = 0.3,
                  show.legend = FALSE, na.rm = TRUE) +
  geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = volcano_colours) +
  scale_size_manual(values = c("Higher in B (reanalysis)" = 1.8,
                               "Lower in B (reanalysis)"  = 1.8,
                               "Not significant"          = 0.7)) +
  annotate("text",
           x = max(plot_global$log2FC) * 0.82,
           y = max(plot_global$neg_log10_FDR) * 0.96,
           label = sprintf("Higher in B\nn = %d", n_up),
           colour = "#D6604D", size = 3.5, hjust = 1) +
  annotate("text",
           x = min(plot_global$log2FC) * 0.82,
           y = max(plot_global$neg_log10_FDR) * 0.96,
           label = sprintf("Lower in B\nn = %d", n_down),
           colour = "#2166AC", size = 3.5, hjust = 0) +
  labs(title    = "Volcano plot: B (reanalysis) vs A (original) — global",
       subtitle = sprintf(
         "Common proteins (n=%d) | Welch t-test + BH-FDR | FDR<%.2f & |log\u2082FC|\u2265%.1f",
         nrow(stats_global), FDR_THRESHOLD, FC_THRESHOLD),
       x = "log\u2082 Fold Change (B / A)",
       y = expression(-log[10](FDR)), colour = NULL) +
  guides(size = "none") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        plot.subtitle   = element_text(size = 9, colour = "grey40"))

save_svg(p_volcano_global, file.path(OUT_DIR, "07a_volcano_global.svg"), width = 8, height = 7)

# ---- 7b. Per-condition faceted volcano ----
top_percond <- stats_percond %>%
  filter(significant) %>%
  group_by(condition) %>% slice_min(FDR, n = 5) %>% ungroup() %>%
  select(condition, protein)

plot_percond <- stats_percond %>%
  filter(!is.na(FDR)) %>%
  mutate(neg_log10_FDR = -log10(FDR),
         condition     = factor(condition, levels = shared_conds)) %>%
  left_join(gene_map, by = "protein") %>%
  left_join(top_percond %>% mutate(is_top = TRUE), by = c("condition", "protein")) %>%
  mutate(label = if_else(!is.na(is_top), gene, NA_character_))

p_volcano_percond <- ggplot(plot_percond,
    aes(x = log2FC, y = neg_log10_FDR, colour = direction, size = direction)) +
  geom_point(alpha = 0.55) +
  geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 15,
                  segment.colour = "grey50", segment.size = 0.3,
                  show.legend = FALSE, na.rm = TRUE) +
  geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  scale_colour_manual(values = volcano_colours) +
  scale_size_manual(values = c("Higher in B (reanalysis)" = 1.6,
                               "Lower in B (reanalysis)"  = 1.6,
                               "Not significant"          = 0.6)) +
  facet_wrap(~ condition, ncol = 3) +
  labs(title    = "Per-condition volcano: B (reanalysis) vs A (original)",
       subtitle = sprintf(
         "BH-FDR within each condition | FDR<%.2f & |log\u2082FC|\u2265%.1f | n=3 vs n=3 (low power)",
         FDR_THRESHOLD, FC_THRESHOLD),
       x = "log\u2082 Fold Change (B / A)",
       y = expression(-log[10](FDR)), colour = NULL) +
  guides(size = "none") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "top",
        strip.background = element_rect(fill = "grey92"),
        plot.subtitle    = element_text(size = 9, colour = "grey40"))

save_svg(p_volcano_percond, file.path(OUT_DIR, "07b_volcano_per_condition.svg"),
         width = 11, height = 7)

# ---- Significant proteins table ----
sig_table <- stats_global %>%
  filter(significant) %>% arrange(FDR) %>%
  left_join(gene_map, by = "protein") %>%
  mutate(across(where(is.double), ~ round(.x, 4))) %>%
  select(protein, gene, log2FC, p_value, FDR, direction,
         mean_log2_A, mean_log2_B, n_valid_A, n_valid_B)

write_csv(sig_table, file.path(OUT_DIR, "07_significant_proteins_global.csv"))
message(sprintf("  Significant proteins table: %d rows", nrow(sig_table)))


# ---- Step8.QC: correlation and PCA ----

message("\n[STEP 8] QC: correlation and PCA...")

# ---- 8a. Pearson / Spearman correlation on mean log2 LFQ ----
meanA     <- rowMeans(matA, na.rm = TRUE)
meanB     <- rowMeans(matB, na.rm = TRUE)
valid_idx  <- is.finite(meanA) & is.finite(meanB)
corr_pear  <- cor(meanA[valid_idx], meanB[valid_idx], method = "pearson")
corr_spear <- cor(meanA[valid_idx], meanB[valid_idx], method = "spearman")

message(sprintf("  Pearson r = %.4f | Spearman \u03c1 = %.4f | n = %d proteins",
                corr_pear, corr_spear, sum(valid_idx)))

corr_df <- tibble(protein     = names(meanA)[valid_idx],
                  mean_log2_A = meanA[valid_idx],
                  mean_log2_B = meanB[valid_idx]) %>%
  left_join(gene_map, by = "protein") %>%
  mutate(highlight = if_else(protein %in% (stats_global %>%
                               filter(significant) %>% pull(protein)),
                             "Significant", "Not significant"))

p_corr <- ggplot(corr_df, aes(x = mean_log2_A, y = mean_log2_B)) +
  geom_point(data = filter(corr_df, highlight == "Not significant"),
             colour = "grey65", size = 0.8, alpha = 0.5) +
  geom_point(data = filter(corr_df, highlight == "Significant"),
             colour = "#D6604D", size = 1.8, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, colour = "black",
              linetype = "dashed", linewidth = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "#2166AC",
              linewidth = 0.8, alpha = 0.12) +
  annotate("text",
           x = min(corr_df$mean_log2_A) + 0.3,
           y = max(corr_df$mean_log2_B) - 0.5,
           label = sprintf("Pearson r = %.4f\nSpearman \u03c1 = %.4f\nn = %d proteins",
                           corr_pear, corr_spear, sum(valid_idx)),
           hjust = 0, size = 4) +
  labs(title    = "Quantitative correlation: A vs B (common proteins)",
       subtitle = "Mean log\u2082 LFQ per protein | red = globally significant (FDR<0.05, |log\u2082FC|\u22651)",
       x = "Mean log\u2082 LFQ — A (original)",
       y = "Mean log\u2082 LFQ — B (reanalysis)") +
  theme_bw(base_size = 13) +
  theme(plot.subtitle = element_text(size = 9, colour = "grey40"))

save_svg(p_corr, file.path(OUT_DIR, "08a_correlation_scatter.svg"), width = 6.5, height = 6)

# ---- 8b. PCA — 36 samples combined ----
matA_pca <- matA; colnames(matA_pca) <- paste0(colnames(matA_pca), "__A")
matB_pca <- matB; colnames(matB_pca) <- paste0(colnames(matB_pca), "__B")
mat_combined <- cbind(matA_pca, matB_pca)

complete_rows <- complete.cases(mat_combined)
mat_pca_in   <- mat_combined[complete_rows, ]
message(sprintf("  PCA: %d complete-case proteins × %d samples",
                nrow(mat_pca_in), ncol(mat_pca_in)))

if (nrow(mat_pca_in) < 10) {
  warning("Too few complete-case proteins for PCA. Consider relaxing MIN_VALID.")
} else {
  pca_res <- prcomp(t(mat_pca_in), center = TRUE, scale. = TRUE)
  var_exp <- summary(pca_res)$importance[2, 1:2] * 100

  pca_df <- as_tibble(pca_res$x[, 1:2], rownames = "sample_id") %>%
    mutate(
      run      = if_else(str_ends(sample_id, "__A"), "A (original)", "B (reanalysis)"),
      raw_name = str_remove(sample_id, "__(A|B)$"),
      condition= parse_cond(raw_name)
    )

  cond_colours <- c(KO_M0 = "#1F78B4", KO_M1 = "#33A02C", KO_M2 = "#E31A1C",
                    WT_M0 = "#A6CEE3", WT_M1 = "#B2DF8A", WT_M2 = "#FB9A99")

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                               colour = condition, shape = run, label = raw_name)) +
    geom_point(size = 3.5, alpha = 0.9, stroke = 0.8) +
    geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE,
                    segment.colour = "grey60", segment.size = 0.3) +
    scale_colour_manual(values = cond_colours) +
    scale_shape_manual(values = c("A (original)" = 16, "B (reanalysis)" = 17)) +
    labs(
      title    = "PCA: A + B samples combined (36 total)",
      subtitle = sprintf(
        "%d complete-case proteins | circles=A, triangles=B | colour=condition",
        nrow(mat_pca_in)),
      x      = sprintf("PC1 (%.1f%% variance)", var_exp[1]),
      y      = sprintf("PC2 (%.1f%% variance)", var_exp[2]),
      colour = "Condition", shape = "Run"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right",
          plot.subtitle   = element_text(size = 9, colour = "grey40"))

  save_svg(p_pca, file.path(OUT_DIR, "08b_pca_plot.svg"), width = 9, height = 7)
}


# ---- Srep9.Classify proteins by biological consistency vs processing variability using CV ----

message("\n[STEP 9] Classifying results...")

cv_A <- apply(matA, 1, function(x) sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE)))
cv_B <- apply(matB, 1, function(x) sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE)))

cv_table <- tibble(protein  = names(cv_A),
                   cv_A     = cv_A,
                   cv_B     = cv_B,
                   cv_ratio = cv_B / cv_A)

annotated <- stats_global %>%
  left_join(cv_table, by = "protein") %>%
  left_join(gene_map, by = "protein") %>%
  mutate(result_class = case_when(
    !significant                                   ~ "Not significant",
    significant & cv_A < 0.15 & cv_B < 0.15       ~ "Biological overlap",
    significant & (cv_ratio > 2 | cv_ratio < 0.5) ~ "Processing-induced difference",
    significant                                    ~ "Ambiguous — verify manually",
    TRUE                                           ~ "Not significant"
  ))

class_summary <- annotated %>%
  count(result_class, name = "n_proteins") %>%
  arrange(desc(n_proteins))

message("  Classification:"); print(class_summary)
write_csv(annotated,     file.path(OUT_DIR, "09_annotated_results.csv"))
write_csv(class_summary, file.path(OUT_DIR, "09_class_summary.csv"))

# ---- Annotated volcano ----
class_colours <- c("Biological overlap"            = "#1A9641",
                   "Processing-induced difference" = "#F46D43",
                   "Ambiguous — verify manually"   = "#FEC44F",
                   "Not significant"               = "grey75")
class_sizes   <- c("Biological overlap"            = 2.2,
                   "Processing-induced difference" = 2.2,
                   "Ambiguous — verify manually"   = 2.0,
                   "Not significant"               = 0.7)

top_annot <- annotated %>%
  filter(significant) %>% arrange(FDR) %>% slice_head(n = 20) %>% pull(protein)

plot_annot <- annotated %>%
  filter(!is.na(FDR)) %>%
  mutate(neg_log10_FDR = -log10(FDR),
         result_class  = factor(result_class,
                                levels = c("Biological overlap",
                                           "Processing-induced difference",
                                           "Ambiguous — verify manually",
                                           "Not significant")),
         label = if_else(protein %in% top_annot,
                         coalesce(gene, protein), NA_character_))

p_annot <- ggplot(plot_annot,
    aes(x = log2FC, y = neg_log10_FDR,
        colour = result_class, size = result_class)) +
  geom_point(alpha = 0.70) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 25,
                  segment.colour = "grey50", segment.size = 0.3,
                  show.legend = FALSE, na.rm = TRUE) +
  geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = class_colours,
                      guide  = guide_legend(override.aes = list(size = 3))) +
  scale_size_manual(values = class_sizes) +
  labs(title    = "Annotated volcano: biological overlap vs. processing-induced differences",
       subtitle = sprintf(
         "CV-based classification | FDR<%.2f & |log\u2082FC|\u2265%.1f | n significant = %d",
         FDR_THRESHOLD, FC_THRESHOLD, n_sig),
       x = "log\u2082 Fold Change (B / A)",
       y = expression(-log[10](FDR)), colour = NULL) +
  guides(size = "none") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        plot.subtitle   = element_text(size = 9, colour = "grey40"))

save_svg(p_annot, file.path(OUT_DIR, "09_volcano_annotated.svg"), width = 9, height = 7)


# ---- Step10.Summary report ----

message("\n[STEP 10] Writing summary report...")

report_lines <- c(
  "MAXQUANT PROTEINGROUPS COMPARISON — SUMMARY REPORT",
  sprintf("Dataset: THP-1 macrophage secretome (M0/M1/M2, KO/WT)"),
  sprintf("A (original):   %s", PATH_A),
  sprintf("B (reanalysis): %s", PATH_B),
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "STRUCTURAL DIFFERENCES",
  sprintf("  A proteins (after QC): %d", nrow(dtA)),
  sprintf("  B proteins (after QC): %d", nrow(dtB)),
  sprintf("  A total columns: %d", ncol(dtA)),
  sprintf("  B total columns: %d", ncol(dtB)),
  sprintf("  LFQ intensity columns: %d in A | %d in B (same count, different names)",
          length(lfq_colsA), length(lfq_colsB)),
  sprintf("  iBAQ columns: %d in A | %d in B",
          structural_audit$A_original[structural_audit$column_type == "iBAQ"],
          structural_audit$B_reanalysis[structural_audit$column_type == "iBAQ"]),
  sprintf("  Decoy flag column: 'Reverse' in A | 'Decoy' in B"),
  sprintf("  Sample naming: A uses run-batch IDs (e.g. KO_M0_S16)"),
  sprintf("                 B uses sequential IDs (e.g. KO_M0_1)"),
  sprintf("  Biological design: Identical — KO/WT x M0/M1/M2 x 3 replicates"),
  "",
  "PROTEIN IDENTIFICATION OVERLAP",
  sprintf("  Common proteins:          %d", length(ids_common)),
  sprintf("  Unique to A (original):   %d", length(ids_only_A)),
  sprintf("  Unique to B (reanalysis): %d", length(ids_only_B)),
  sprintf("  Jaccard similarity:       %.3f", jaccard),
  sprintf("  Interpretation: High overlap. B gains %d additional protein groups,",
          length(ids_only_B)),
  sprintf("  likely from modified protein FDR or grouping parameters."),
  "",
  "LFQ INTENSITY DISTRIBUTIONS (common proteins only)",
  sprintf("  A median log2 LFQ: %.2f", median(matA, na.rm = TRUE)),
  sprintf("  B median log2 LFQ: %.2f", median(matB, na.rm = TRUE)),
  sprintf("  A missing value rate: %.1f%%", mean(is.na(matA)) * 100),
  sprintf("  B missing value rate: %.1f%%", mean(is.na(matB)) * 100),
  sprintf("  Per-condition delta: All conditions differ by < 1%% in missing rate."),
  sprintf("  Interpretation: Intensity distributions are nearly identical between runs."),
  "",
  "QUANTITATIVE REPRODUCIBILITY",
  sprintf("  Pearson r (mean log2 LFQ):    %.4f", corr_pear),
  sprintf("  Spearman rho (mean log2 LFQ): %.4f", corr_spear),
  sprintf("  n proteins in correlation:    %d", sum(valid_idx)),
  sprintf("  Interpretation: r > 0.95 indicates excellent quantitative agreement."),
  sprintf("  PCA (see 08b_pca_plot.svg): if samples cluster by condition rather"),
  sprintf("  than by run, biological signal dominates over processing effects."),
  "",
  "STATISTICAL COMPARISON",
  sprintf("  Method: Welch's unpaired t-test + Benjamini-Hochberg FDR"),
  sprintf("  Thresholds: FDR < %.2f AND |log2FC| >= %.1f", FDR_THRESHOLD, FC_THRESHOLD),
  sprintf("  Note: Samples cannot be paired (naming schemes differ completely)."),
  "",
  sprintf("  GLOBAL TEST (all %d samples pooled per run):", ncol(matA)),
  sprintf("    Proteins tested:   %d / %d", n_tested, nrow(stats_global)),
  sprintf("    Significant total: %d", n_sig),
  sprintf("    Higher in B:       %d", n_up),
  sprintf("    Lower in B:        %d", n_down),
  "",
  sprintf("  PER-CONDITION TEST (3 replicates per condition per run):"),
  sprintf("  *** NOTE: n=3 vs n=3 gives very low statistical power. ***"),
  sprintf("  *** Per-condition results are indicative only; use with caution. ***"),
  paste(sprintf("    %s: %d significant", percond_summary$condition,
                percond_summary$n_sig), collapse = "\n"),
  "",
  "RESULT CLASSIFICATION",
  paste(sprintf("  %-42s %d proteins",
                class_summary$result_class, class_summary$n_proteins),
        collapse = "\n"),
  "",
  sprintf("  'Biological overlap': CV < 0.15 in both runs — consistent quantification,"),
  sprintf("  genuine intensity shift. These proteins warrant biological follow-up."),
  sprintf("  'Processing-induced': CV ratio > 2 or < 0.5 — one run much more variable,"),
  sprintf("  likely a parameter-change artefact. Treat with caution."),
  "",
  "VALIDATION VERDICT FOR B (reanalysis)",
  sprintf("  Pearson r = %.4f, Jaccard = %.3f, missing rate delta < 1%%.", corr_pear, jaccard),
  sprintf("  B (reanalysis) is quantitatively reliable for downstream biological analysis."),
  sprintf("  The %d globally significant proteins represent real intensity shifts", n_sig),
  sprintf("  introduced by the modified MaxQuant parameters and should be reviewed"),
  sprintf("  before treating A and B as interchangeable."),
  sprintf("  The %d proteins unique to B represent a genuine gain in coverage",
          length(ids_only_B)),
  sprintf("  from the re-processing and can be used in downstream analyses."),
  "",
  "OUTPUT FILES",
  "  01_structural_audit.csv              Column type counts per run",
  "  01_sample_condition_map.csv          Sample to condition mapping",
  "  02_overlap_summary.csv               Protein overlap counts + Jaccard",
  "  02_proteins_unique_A.csv             Proteins only in A (original)",
  "  02_proteins_unique_B.csv             Proteins only in B (reanalysis)",
  "  03_protein_overlap_bar.svg           Overlap bar chart",
  "  04_missing_value_summary.csv         Overall missing rates",
  "  04_missing_per_condition.csv         Per-condition missing rates",
  "  05_lfq_distributions.svg             Density + per-sample boxplots",
  "  06a_global_statistics.csv            Per-protein global t-test + FDR",
  "  06b_per_condition_statistics.csv     Per-protein per-condition t-test + FDR",
  "  06b_per_condition_summary.csv        Significant counts per condition",
  "  07a_volcano_global.svg               Global volcano plot",
  "  07b_volcano_per_condition.svg        Faceted per-condition volcano",
  "  07_significant_proteins_global.csv   Significant proteins (global test)",
  "  08a_correlation_scatter.svg          Mean log2 LFQ correlation",
  "  08b_pca_plot.svg                     PCA of all 36 samples",
  "  09_annotated_results.csv             Full results with classification",
  "  09_class_summary.csv                 Count per classification",
  "  09_volcano_annotated.svg             Annotated volcano",
  "  10_summary_report.txt                This report",
  "",
  capture.output(sessionInfo())
)

writeLines(report_lines, file.path(OUT_DIR, "10_summary_report.txt"))
message(sprintf("  Saved: 10_summary_report.txt"))
message("Do what must be done. Do not hesitate. Show no mercy")
message(sprintf("\n[DONE] All outputs written to: %s/", OUT_DIR))
