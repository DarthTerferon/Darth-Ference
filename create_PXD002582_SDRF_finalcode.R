# SDRF-Proteomics file generator for PXD002582 — PROTEOME ONLY

# ---- 0. Setup ----------------------------------------------------------------
output_path <- "SDRF_PXD002582_proteome.tsv"
ftp_base    <- "https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/08/PXD002582/"

# ---- 1. Define 48 proteome RAW filenames --------------------------------

# Set1
proteome_set1 <- c(
  "20131128_Macrophage_set1_FN_1_01.raw",
  "20131128_Macrophage_set1_FN_1_02.raw",
  "20131128_Macrophage_set1_FN_2_01.raw",
  "20131128_Macrophage_set1_FN_2_02.raw",
  "20131128_Macrophage_set1_FN_3_01.raw",
  "20131128_Macrophage_set1_FN_3_02.raw",
  "20131128_Macrophage_set1_FN_4_01.raw",
  "20131128_Macrophage_set1_FN_4_02.raw",
  "20131128_Macrophage_set1_FN_5_01.raw",
  "20131128_Macrophage_set1_FN_5_02.raw",
  "20131128_Macrophage_set1_FN_6_01.raw",
  "20131128_Macrophage_set1_FN_6_02.raw",
  "20131128_Macrophage_set1_FN_7_01.raw",
  "20131128_Macrophage_set1_FN_7_02.raw",
  "20131128_Macrophage_set1_FN_8_01.raw",
  "20131128_Macrophage_set1_FN_8_02.raw",
  "20131128_Macrophage_set1_FN_9_01.raw",
  "20131128_Macrophage_set1_FN_9_02.raw",
  "20131128_Macrophage_set1_FN_10_01.raw",
  "20131128_Macrophage_set1_FN_10_02.raw",
  "20131128_Macrophage_set1_FN_11_01.raw",
  "20131128_Macrophage_set1_FN_11_02.raw",
  "20131128_Macrophage_set1_FN_12_01.raw",
  "20131128_Macrophage_set1_FN_12_02.raw"
)

# Set2
proteome_set2 <- paste0("20140203_Macrophage_", 1:12, ".raw")

# Set3
proteome_set3 <- paste0("20140205_Macrophage_", 1:12, ".raw")

all_proteome_files <- c(proteome_set1, proteome_set2, proteome_set3)

stopifnot(length(proteome_set1) == 24)
stopifnot(length(proteome_set2) == 12)
stopifnot(length(proteome_set3) == 12)
stopifnot(length(all_proteome_files) == 48)

# ---- 2. Function to parse metadata from filename --------------------------------

parse_file_metadata <- function(filename) {
  base <- sub("\\.raw$", "", filename, ignore.case = TRUE)

  if (grepl("set1_FN", base)) {
    m     <- regmatches(base, regexpr("FN_(\\d+)_(\\d+)$", base))
    parts <- strsplit(m, "_")[[1]]
    list(biorep   = 1L,
         fraction = as.integer(parts[2]),
         tech_rep = as.integer(parts[3]))

  } else if (grepl("20140203_Macrophage", base)) {
    list(biorep   = 2L,
         fraction = as.integer(sub(".*Macrophage_(\\d+)$", "\\1", base)),
         tech_rep = 1L)

  } else if (grepl("20140205_Macrophage", base)) {
    list(biorep   = 3L,
         fraction = as.integer(sub(".*Macrophage_(\\d+)$", "\\1", base)),
         tech_rep = 1L)

  } else {
    stop(paste("Unrecognised filename pattern:", filename))
  }
}

# ---- 3. Condition mapping ------------------------------------

tmt_channels <- data.frame(
  label           = c("TMT126",             "TMT127",             "TMT128",           "TMT129"),
  condition       = c("GM-BMM",             "GM-BMM+LPS",         "M-BMM",            "M-BMM+LPS"),
  macrophage_type = c("GM-CSF macrophage",  "GM-CSF macrophage",  "M-CSF macrophage", "M-CSF macrophage"),
  lps_treatment   = c("untreated",          "LPS treated",        "untreated",        "LPS treated"),
  stringsAsFactors = FALSE
)

# ---- 4. Build SDRF rows for one RAW file ------------------------------------

build_rows <- function(filename, meta) {
  assay_name <- sub("\\.raw$", "", filename, ignore.case = TRUE)
  file_uri   <- paste0(ftp_base, filename)

  rows <- lapply(seq_len(nrow(tmt_channels)), function(i) {
    ch <- tmt_channels[i, ]
    source_name <- paste0(
      gsub("[^A-Za-z0-9]", "_", ch$condition),
      "_biorep", meta$biorep
    )
    data.frame(
      # --- Sample metadata (characteristics) ---
      `source name`                              = source_name,
      `characteristics[organism]`               = "Mus musculus",
      `characteristics[organism part]`          = "bone marrow",
      `characteristics[sex]`                    = "male",
      `characteristics[age]`                    = "5-10 weeks",
      `characteristics[developmental stage]`   = "adult",
      `characteristics[biological replicate]`   = meta$biorep,
      `characteristics[material type]`          = "cell",
      `characteristics[cell type]`              = "macrophage",
      `characteristics[disease]`                = "normal",
      `characteristics[macrophage type]`        = ch$macrophage_type,
      `characteristics[LPS stimulation]`        = ch$lps_treatment,
      # --- Assay metadata ---
      `assay name`                              = assay_name,
      `technology type`                         = "proteomic profiling by mass spectrometry",
      # --- Comment fields ---
      `comment[label]`                          = ch$label,
      `comment[file uri]`                       = file_uri,
      `comment[fraction identifier]`            = meta$fraction,
      `comment[technical replicate]`            = meta$tech_rep,
      `comment[instrument]`                     = "NT=Q Exactive;AC=MS:1001911",
      `comment[modification parameters]`        = "NT=TMT4plex;PP=Any N-term;AC=UNIMOD:214;MT=fixed",
      `comment[modification parameters].1`      = "NT=TMT4plex;TA=K;AC=UNIMOD:214;MT=fixed",
      `comment[modification parameters].2`      = "NT=Carbamidomethyl;TA=C;AC=UNIMOD:4;MT=fixed",
      `comment[modification parameters].3`      = "NT=Oxidation;TA=M;AC=UNIMOD:35;MT=variable",
      `comment[cleavage agent details]`         = "NT=Trypsin",
      `comment[fragment mass tolerance]`        = "0.8 Da",
      `comment[precursor mass tolerance]`       = "10 ppm",
      `comment[proteomics data acquisition method]` = "NT=Data-dependent acquisition;AC=PRIDE:0000627",
      `comment[dissociation method]`            = "NT=HCD;AC=PRIDE:0000590",
      `comment[fractionation method]`           = "peptide isoelectric focusing",
      `comment[data file]`                      = filename,
      `comment[sdrf version]`                   = "v1.1.0",
      `comment[sdrf template]`                  = "NT=ms-proteomics;VV=v1.1.0",
      # --- Factor values (study variables) ---
      `factor value[macrophage type]`           = ch$macrophage_type,
      `factor value[LPS treatment]`             = ch$lps_treatment,

      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
}

# ---- 5. Generate all rows ---------------------------------------------------

sdrf_rows <- lapply(all_proteome_files, function(f) {
  meta <- parse_file_metadata(f)
  build_rows(f, meta)
})

sdrf <- do.call(rbind, sdrf_rows)
rownames(sdrf) <- NULL

cat("SDRF dimensions:", nrow(sdrf), "rows x", ncol(sdrf), "columns\n")
cat("Expected: 48 files x 4 TMT channels =", 48 * 4, "rows\n")

# ---- 6. Rename duplicate comment[modification parameters] columns -----------

col_names <- colnames(sdrf)
col_names[col_names == "comment[modification parameters].1"] <- "comment[modification parameters]"
col_names[col_names == "comment[modification parameters].2"] <- "comment[modification parameters]"
col_names[col_names == "comment[modification parameters].3"] <- "comment[modification parameters]"
colnames(sdrf) <- col_names

# ---- 7. Validate -----------------------------------------------------------

cat("\n--- Validation ---\n")
cat("Unique source names:", length(unique(sdrf[["source name"]])), "\n")
print(sort(unique(sdrf[["source name"]])))
cat("Unique assay names:", length(unique(sdrf[["assay name"]])), "\n")
cat("TMT labels:", paste(sort(unique(sdrf[["comment[label]"]])), collapse = ", "), "\n")
cat("Biological replicates:", paste(sort(unique(sdrf[["characteristics[biological replicate]"]])), collapse = ", "), "\n")
cat("Fraction range:", range(sdrf[["comment[fraction identifier]"]]), "\n")
cat("Technical replicates:", paste(sort(unique(sdrf[["comment[technical replicate]"]])), collapse = ", "), "\n")
cat("Sex:", unique(sdrf[["characteristics[sex]"]]), "\n")
cat("Age:", unique(sdrf[["characteristics[age]"]]), "\n")
cat("Material type:", unique(sdrf[["characteristics[material type]"]]), "\n")
cat("Instrument:", unique(sdrf[["comment[instrument]"]]), "\n")
cat("Cleavage agent:", unique(sdrf[["comment[cleavage agent details]"]]), "\n")
cat("DDA method:", unique(sdrf[["comment[proteomics data acquisition method]"]]), "\n")
cat("Dissociation:", unique(sdrf[["comment[dissociation method]"]]), "\n")
cat("Precursor tolerance:", unique(sdrf[["comment[precursor mass tolerance]"]]), "\n")
cat("Fragment tolerance:", unique(sdrf[["comment[fragment mass tolerance]"]]), "\n")
cat("SDRF version:", unique(sdrf[["comment[sdrf version]"]]), "\n")
cat("SDRF template:", unique(sdrf[["comment[sdrf template]"]]), "\n")

mod_cols <- which(colnames(sdrf) == "comment[modification parameters]")
cat("Modification parameter columns:", length(mod_cols), "\n")
for (i in seq_along(mod_cols)) {
  cat("  [", i, "]", unique(sdrf[[mod_cols[i]]]), "\n")
}

cat("File URI example:", sdrf[["comment[file uri]"]][1], "\n")

required_cols <- c("source name", "characteristics[organism]", "assay name",
                   "comment[label]", "comment[instrument]", "comment[data file]",
                   "comment[fraction identifier]", "comment[technical replicate]",
                   "characteristics[biological replicate]", "comment[file uri]",
                   "comment[sdrf version]", "comment[sdrf template]",
                   "characteristics[material type]")
cat("\nNA/empty check in required columns:\n")
all_ok <- TRUE
for (col in required_cols) {
  n_na <- sum(is.na(sdrf[[col]]) | sdrf[[col]] == "")
  status <- if (n_na == 0) "OK" else paste("PROBLEM:", n_na, "empty")
  cat(sprintf("  %-55s %s\n", col, status))
  if (n_na > 0) all_ok <- FALSE
}
cat("All required fields complete:", all_ok, "\n")

# ---- 8. Write SDRF -----------------------------------------------------------

out_dir <- file.path(getwd(), "results")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

file_path <- file.path(out_dir, output_path)

write.table(
  sdrf,
  file = file_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  fileEncoding = "UTF-8"
)

cat("\nSDRF written to:", file_path, "\n")
cat("File size:", file.size(file_path), "bytes\n")
cat("Columns:", ncol(sdrf), "\n")
cat("Rows:", nrow(sdrf), "\n")

if (nrow(sdrf) == 192) {message ("GOOD!.. GOOOOOOD!..")}
