#Script for converting HIBAG to a format that can be read by HATK/HLA2HPED

debug <- F

if (debug) {
    input <- "/home/oystein/hla_imputation_pipeout/2026.02.24/HIBAG/hibag.A"
    hla <- "A"
    output <- "/home/oystein/test/hibag.hla2hped_input.A"
} else {
    args <- commandArgs(TRUE)
    input <- args[1]
    hla <- args[2]
    output <- args[3]
}

library(dplyr)
library(tidyr)

hibag <- read.table(input, header=T, sep = "\t")

hibag$fid <- hibag$sample.id
hibag$hla <- hla

hibag$alleles_1field <- paste0(sapply(hibag$allele1, function(x) strsplit(x, split = ":")[[1]][1]), ",", sapply(hibag$allele2, function(x) strsplit(x, split = ":")[[1]][1]))
hibag$alleles_2field <- paste0(sapply(hibag$allele1, function(x) gsub(":", "", x)), ",", sapply(hibag$allele2, function(x) gsub(":", "", x)))
hla2hped_input <- hibag %>% select(fid, sample.id, hla, alleles_1field, alleles_2field, prob, matching)
write.table(hla2hped_input, file=output, col.names=F, row.names=F, sep = "\t", quote=F)