library(dplyr)
library(tidyr)

debug <- T

if (debug) {
    args <- c(
        "/home/oystein/hla_imputation_pipeout/2025.08.23/HIBAG/hibag",
        "/home/oystein/hla_imputation_pipeout/2025.08.23/CookHLA/cookhla_output.MHC.HLA_IMPUTATION_OUT.alleles",
        "/home/oystein/github/HLA-imputation/snakefiles/resources/norwegian_allele_frequencies/common/HLA",
        "/home/oystein/hla_imputation_pipeout/2025.08.23/docs/"
    )
} else {
    args <- commandArgs(TRUE)
}

hibag_trunk_file <- args[1]
cookhla_file <- args[2]
ref_trunk_file <- args[3]
out_folder <- args[4]

hibag_a <- read.table(paste0(hibag_trunk_file, ".A"), sep = "\t", header = T)
hibag_b <- read.table(paste0(hibag_trunk_file, ".B"), sep = "\t", header = T)
hibag_c <- read.table(paste0(hibag_trunk_file, ".C"), sep = "\t", header = T)
hibag_dpb1 <- read.table(paste0(hibag_trunk_file, ".DPB1"), sep = "\t", header = T)
hibag_dqb1 <- read.table(paste0(hibag_trunk_file, ".DQB1"), sep = "\t", header = T)
hibag_drb1 <- read.table(paste0(hibag_trunk_file, ".DRB1"), sep = "\t", header = T)


cookhla <- read.table(cookhla_file, sep = " ", header = F, col.names = c("fid", "iid", "hla", "dig2", "dig4", "post1", "post2", "confidence"))
cookhla <- cookhla %>%
    separate(dig4, into = c("allele1", "allele2"), sep = ",") %>%
    mutate(
        allele1 = sub("(\\d{2})(\\d{2})", "\\1:\\2", allele1),
        allele2 = sub("(\\d{2})(\\d{2})", "\\1:\\2", allele2)
    )


extract_common_allele_freqs <- function(allele_table, ref) {
    allele_freqs <- table(c(allele_table$allele1, allele_table$allele2))
    norm_freqs <- allele_freqs / sum(allele_freqs)
    filtered_freqs <- norm_freqs[names(norm_freqs) %in% ref$allele]
    filtered_freqs <- c(filtered_freqs, "rare" = 1 - sum(filtered_freqs, na.rm = TRUE))
    return(filtered_freqs)
}


plot_frequencies <- function(cookhla, hibag, ref, hla) {
    freq_matrix <- rbind(cookhla[ref$allele], hibag[ref$allele], ref$freq)
    png(paste0(out_folder, "frequencies_HLA_", hla, ".png"), width = 800, height = 600)
    barplot(freq_matrix,
        beside = TRUE,
        main = paste0("Allele frequencies for HLA-", hla),
        xlab = "Allele", ylab = "Frequency", names.arg = ref$allele,
        col = c("lightblue", "pink", "darkgreen"), las = 2
    )
    legend("topright", legend = c("CookHLA", "HIBAG", "Reference"), fill = c("lightblue", "pink", "darkgreen"))
    dev.off()
}

plot_hibag_ref <- function(hibag, ref, hla) {
    freq_matrix <- rbind(hibag[ref$allele], ref$freq)
    png(paste0(out_folder, "frequencies_HLA_", hla, ".png"), width = 800, height = 600)
    barplot(freq_matrix,
        beside = TRUE,
        main = paste0("Allele frequencies for HLA-", hla),
        xlab = "Allele", ylab = "Frequency", names.arg = ref$allele,
        col = c("pink", "darkgreen"), las = 2
    )
    legend("topright", legend = c("HIBAG", "Reference"), fill = c("pink", "darkgreen"))
    dev.off()
}

add_rare_row <- function(ref) {
    rare_freq <- 1 - sum(ref$freq)
    rare_row <- data.frame(allele = "rare", freq = rare_freq)
    ref_rare <- rbind(ref, rare_row)
    return(ref_rare)
}

ref_a <- add_rare_row(read.table(paste0(ref_trunk_file, ".A"), header = T))
ref_b <- add_rare_row(read.table(paste0(ref_trunk_file, ".B"), header = T))
ref_c <- add_rare_row(read.table(paste0(ref_trunk_file, ".C"), header = T))
ref_dpb1 <- add_rare_row(read.table(paste0(ref_trunk_file, ".DPB1"), header = T))
ref_dqb1 <- add_rare_row(read.table(paste0(ref_trunk_file, ".DQB1"), header = T))
ref_drb1 <- add_rare_row(read.table(paste0(ref_trunk_file, ".DRB1"), header = T))


common_hibag_a <- extract_common_allele_freqs(hibag_a, ref_a)
common_hibag_b <- extract_common_allele_freqs(hibag_b, ref_b)
common_hibag_c <- extract_common_allele_freqs(hibag_c, ref_c)
common_hibag_dpb1 <- extract_common_allele_freqs(hibag_dpb1, ref_dpb1)
common_hibag_dqb1 <- extract_common_allele_freqs(hibag_dqb1, ref_dqb1)
common_hibag_drb1 <- extract_common_allele_freqs(hibag_drb1, ref_drb1)


common_cookhla_a <- extract_common_allele_freqs(subset(cookhla, hla == "A"), ref_a)
common_cookhla_b <- extract_common_allele_freqs(subset(cookhla, hla == "B"), ref_b)
common_cookhla_c <- extract_common_allele_freqs(subset(cookhla, hla == "C"), ref_c)
common_cookhla_dqb1 <- extract_common_allele_freqs(subset(cookhla, hla == "DQB1"), ref_dqb1)
common_cookhla_drb1 <- extract_common_allele_freqs(subset(cookhla, hla == "DRB1"), ref_drb1)


plot_frequencies(common_cookhla_a, common_hibag_a, ref_a, "A")
plot_frequencies(common_cookhla_b, common_hibag_b, ref_b, "B")
plot_frequencies(common_cookhla_c, common_hibag_c, ref_c, "C")
plot_hibag_ref(common_hibag_dpb1, ref_dpb1, "DPB1")
plot_frequencies(common_cookhla_dqb1, common_hibag_dqb1, ref_dqb1, "DQB1")
plot_frequencies(common_cookhla_drb1, common_hibag_drb1, ref_drb1, "DRB1")
