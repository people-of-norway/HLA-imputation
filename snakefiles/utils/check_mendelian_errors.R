library(dplyr)
library(tidyr)

debug <- F

if (debug) {
    args <- c(
        "/home/oystein/hla_imputation_pipeout/2025.08.23/HIBAG/hibag",
        "/home/oystein/hla_imputation_pipeout/2025.08.23/CookHLA/cookhla_output.MHC.HLA_IMPUTATION_OUT.alleles",
        "/home/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fam",
        "/home/oystein/hla_imputation_pipeout/2025.08.23/docs/error_rates"
    )
} else {
    args <- commandArgs(TRUE)
}

hibag_trunk_file <- args[1]
cookhla_file <- args[2]
fam_file <- args[3]
out_file <- args[4]


fam <- read.table(
    file = fam_file,
    header = F,
    col.names = c("fid", "iid", "pat", "mat", "sex", "phen")
)

trios <- subset(fam, !is.na(pat) & pat != "0" & !is.na(mat) & mat != "0")

hibag_a <- read.table(paste0(hibag_trunk_file, ".A"), sep = "\t", header = T) %>% rename(iid = sample.id)
hibag_b <- read.table(paste0(hibag_trunk_file, ".B"), sep = "\t", header = T) %>% rename(iid = sample.id)
hibag_c <- read.table(paste0(hibag_trunk_file, ".C"), sep = "\t", header = T) %>% rename(iid = sample.id)
hibag_dpb1 <- read.table(paste0(hibag_trunk_file, ".DPB1"), sep = "\t", header = T) %>% rename(iid = sample.id)
hibag_dqb1 <- read.table(paste0(hibag_trunk_file, ".DQB1"), sep = "\t", header = T) %>% rename(iid = sample.id)
hibag_drb1 <- read.table(paste0(hibag_trunk_file, ".DRB1"), sep = "\t", header = T) %>% rename(iid = sample.id)


cookhla <- read.table(cookhla_file, sep = " ", header = F, col.names = c("fid", "iid", "hla", "dig2", "dig4", "post1", "post2", "confidence"))
cookhla <- cookhla %>%
    separate(dig4, into = c("allele1", "allele2"), sep = ",") %>%
    mutate(
        allele1 = sub("(\\d{2})(\\d{2})", "\\1:\\2", allele1),
        allele2 = sub("(\\d{2})(\\d{2})", "\\1:\\2", allele2)
    )

cookhla_a <- subset(cookhla, hla == "A")
cookhla_b <- subset(cookhla, hla == "B")
cookhla_c <- subset(cookhla, hla == "C")
cookhla_dqa1 <- subset(cookhla, hla == "DQA1")
cookhla_dqb1 <- subset(cookhla, hla == "DQB1")
cookhla_drb1 <- subset(cookhla, hla == "DRB1")


get_trios_alleles <- function(imp, trios) {
    trios_alleles <- trios %>%
        left_join(imp %>% select(iid, child_allele1 = allele1, child_allele2 = allele2), by = "iid") %>%
        left_join(imp %>% select(pat = iid, pat_allele1 = allele1, pat_allele2 = allele2), by = "pat") %>%
        left_join(imp %>% select(mat = iid, mat_allele1 = allele1, mat_allele2 = allele2), by = "mat")
    return(na.omit(trios_alleles))
}


trios_hibag_a <- get_trios_alleles(hibag_a, trios)
trios_hibag_b <- get_trios_alleles(hibag_b, trios)
trios_hibag_c <- get_trios_alleles(hibag_c, trios)
trios_hibag_dpb1 <- get_trios_alleles(hibag_dpb1, trios)
trios_hibag_dqb1 <- get_trios_alleles(hibag_dqb1, trios)
trios_hibag_drb1 <- get_trios_alleles(hibag_drb1, trios)

trios_cookhla_a <- get_trios_alleles(cookhla_a, trios)
trios_cookhla_b <- get_trios_alleles(cookhla_b, trios)
trios_cookhla_c <- get_trios_alleles(cookhla_c, trios)
trios_cookhla_dqa1 <- get_trios_alleles(cookhla_dqa1, trios)
trios_cookhla_dqb1 <- get_trios_alleles(cookhla_dqb1, trios)
trios_cookhla_drb1 <- get_trios_alleles(cookhla_drb1, trios)

mendelian_check <- function(row) {
    child_allele1 <- row["child_allele1"]
    child_allele2 <- row["child_allele2"]
    paternal_alleles <- c(row["pat_allele1"], row["pat_allele2"])
    maternal_alleles <- c(row["mat_allele1"], row["mat_allele2"])

    possible <- (
        child_allele1 %in% paternal_alleles && child_allele2 %in% maternal_alleles ||
            child_allele2 %in% paternal_alleles && child_allele1 %in% maternal_alleles
    )

    return(possible)
}

trios_hibag_a$mendel <- apply(trios_hibag_a, 1, mendelian_check)
trios_hibag_b$mendel <- apply(trios_hibag_b, 1, mendelian_check)
trios_hibag_c$mendel <- apply(trios_hibag_c, 1, mendelian_check)
trios_hibag_dpb1$mendel <- apply(trios_hibag_dpb1, 1, mendelian_check)
trios_hibag_dqb1$mendel <- apply(trios_hibag_dqb1, 1, mendelian_check)
trios_hibag_drb1$mendel <- apply(trios_hibag_drb1, 1, mendelian_check)

trios_cookhla_a$mendel <- apply(trios_cookhla_a, 1, mendelian_check)
trios_cookhla_b$mendel <- apply(trios_cookhla_b, 1, mendelian_check)
trios_cookhla_c$mendel <- apply(trios_cookhla_c, 1, mendelian_check)
trios_cookhla_dqa1$mendel <- apply(trios_cookhla_dqa1, 1, mendelian_check)
trios_cookhla_dqb1$mendel <- apply(trios_cookhla_dqb1, 1, mendelian_check)
trios_cookhla_drb1$mendel <- apply(trios_cookhla_drb1, 1, mendelian_check)


write(
    x = "hla\tsoftware\tn_trios\tn_erros\terror_rate",
    file = out_file,
    append = F
)


add_table_row <- function(trios_alleles, software, hla) {
    n_errors <- sum(!trios_alleles$mendel)
    n_rows <- nrow(trios_alleles)
    error_rate <- n_errors / n_rows
    write(
        x = paste0(hla, "\t", software, "\t", n_rows, "\t", n_errors, "\t", error_rate),
        file = out_file,
        append = T
    )
}

add_table_row(trios_hibag_a, "HIBAG", "A")
add_table_row(trios_cookhla_a, "CookHLA", "A")
add_table_row(trios_hibag_b, "HIBAG", "B")
add_table_row(trios_cookhla_b, "CookHLA", "B")
add_table_row(trios_hibag_c, "HIBAG", "C")
add_table_row(trios_cookhla_c, "CookHLA", "C")
add_table_row(trios_hibag_dpb1, "HIBAG", "DPB1")
add_table_row(trios_cookhla_dqa1, "CookHLA", "DQA1")
add_table_row(trios_hibag_dqb1, "HIBAG", "DQB1")
add_table_row(trios_cookhla_dqb1, "CookHLA", "DQB1")
add_table_row(trios_hibag_drb1, "HIBAG", "DRB1")
add_table_row(trios_cookhla_drb1, "CookHLA", "DRB1")
