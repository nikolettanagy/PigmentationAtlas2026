library(data.table)

input_file <- "results/GCST90691600_harmonized.tsv.gz"
output_file <- "results/GCST90691600_clumping_input_uniqueID.tsv"

cat("Reading harmonized GWAS...\n")

gwas <- fread(
  input_file,
  select = c("variant_id", "p_value")
)

cat("Rows read:", nrow(gwas), "\n")

clump_input <- gwas[
  !is.na(variant_id) &
  variant_id != "" &
  !is.na(p_value) &
  p_value > 0 &
  p_value <= 0.05
]

# GWAS format: 1_11063_T_G
# Reference format: 1:11063:T:G
clump_input[, ID := gsub("_", ":", variant_id, fixed = TRUE)]

clump_input <- clump_input[, .(
  ID,
  P = p_value
)]

setorder(clump_input, ID, P)
clump_input <- unique(clump_input, by = "ID")

fwrite(
  clump_input,
  output_file,
  sep = "\t"
)

cat("Variants retained for clumping:", nrow(clump_input), "\n")
cat("Unique allele-specific IDs:", uniqueN(clump_input$ID), "\n")
cat("Saved to:", output_file, "\n")