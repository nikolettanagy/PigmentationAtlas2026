#!/usr/bin/env Rscript
# Run from the extracted candidate root, not from scripts/.
base <- "examples/CLUMP_0001"
g <- read.delim(gzfile(file.path(base, "gwas/gwas_region.tsv.gz")))
q <- read.delim(gzfile(file.path(base, "apaqtl/Skin_Not_Sun_Exposed_Suprapubic",
                                 "CLUMP_0001_chr1_1122649_1461195.apaqtl.tsv.gz")))
h <- read.delim(gzfile(file.path(base, "coloc/COLOC_0000001_harmonized_input.tsv.gz")))
q <- q[q$phenotype_id == "ENSG00000008130.15_chr1:1751232-1753060" &
         q$gene_id == "ENSG00000008130.15", ]
g <- g[!is.na(g$base_pair_location) & is.finite(g$beta) &
         is.finite(g$standard_error) & g$standard_error > 0 &
         is.finite(g$p_value) & g$p_value > 0 & g$p_value <= 1, ]
q <- q[!is.na(q$variant_position) & is.finite(q$slope) &
         is.finite(q$slope_se) & q$slope_se > 0 &
         is.finite(q$pval_nominal) & q$pval_nominal > 0 & q$pval_nominal <= 1, ]
# Reproduce Module 12.2's lowest-p row rule for duplicate genomic positions.
g <- g[order(g$base_pair_location, g$p_value), ]
q <- q[order(q$variant_position, q$pval_nominal), ]
g <- g[!duplicated(g$base_pair_location), ]
q <- q[!duplicated(q$variant_position), ]
stopifnot(nrow(g) == 3030L, nrow(q) == 1745L)
m <- merge(g, q, by.x = "base_pair_location", by.y = "variant_position")
stopifnot(nrow(m) == 1532L)
complement <- function(z) chartr("ACGT", "TGCA", z)
ge <- as.character(m$effect_allele.x); go <- as.character(m$other_allele.x)
ref <- as.character(m$ref); alt <- as.character(m$alt)
direct <- alt == ge & ref == go
swapped <- alt == go & ref == ge
c_direct <- complement(alt) == ge & complement(ref) == go
c_swapped <- complement(alt) == go & complement(ref) == ge
relationship <- ifelse(direct, "DIRECT", ifelse(swapped, "SWAPPED",
                   ifelse(c_direct, "COMPLEMENT_DIRECT",
                     ifelse(c_swapped, "COMPLEMENT_SWAPPED", "INCOMPATIBLE"))))
palindromic <- paste(pmin(ge, go), pmax(ge, go)) %in% c("A T", "C G")
af_g <- m$effect_allele_frequency.x
af_q <- m$af
resolvable <- !palindromic | (is.finite(af_g) & is.finite(af_q) &
                  af_g > 0 & af_g < 1 & af_q > 0 & af_q < 1 &
                  (abs(af_g - af_q) <= .10 | abs(af_g - (1 - af_q)) <= .10))
keep <- relationship != "INCOMPATIBLE" & resolvable
m <- m[keep, ]; relationship <- relationship[keep]
stopifnot(nrow(m) == 1524L)
# Check the complete selected variant set, allele orientation and effects.
ix <- match(m$variant_id.x, h$gwas_variant_id)
stopifnot(!anyNA(ix), !anyDuplicated(ix),
          setequal(m$variant_id.y, h$apaqtl_variant_id),
          identical(as.character(m$variant_id.y), as.character(h$apaqtl_variant_id[ix])),
          identical(as.character(relationship), as.character(h$allele_relationship[ix])),
          all(abs(m$beta - h$gwas_beta[ix]) < 1e-9),
          all(abs(m$slope * ifelse(relationship %in% c("SWAPPED", "COMPLEMENT_SWAPPED"),
                                  -1, 1) - h$apaqtl_beta_harmonized[ix]) < 1e-9))
cat("PASS: reconstructed 1524 selected and allele-harmonized variants from raw slices.\n")
source("scripts/check_coloc_example.R", local = TRUE)
