library(dplyr)
library(tidyr)

debug <- T

if (debug) {
    args <- c(
        "/home/oystein/hla_imputation_pipeout/2025.08.23/HIBAG/hibag",
        "/home/oystein/hla_imputation_pipeout/2025.08.23/CookHLA/cookhla_output.MHC.HLA_IMPUTATION_OUT.alleles",
        "/home/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fam",
        "/home/oystein/github/HLA-imputation/snakefiles/resources/norwegian_allele_frequencies/common/HLA",
        "/home/oystein/test/2025.08.23/docs",
        "/home/oystein/test/2025.08.23/pipeout",
        "HLA Imputation report TEST"
    )
} else {
    args <- commandArgs(TRUE)
}


hibag_trunk_file <- args[1]
cookhla_file <- args[2]
fam_file <- args[3]
ref_trunk_file <- args[4]
docs_folder <- args[5]
pipeout_folder <- args[6]
title <- args[7]

frequencies_folder <- paste0(docs_folder, "/frequencies/")
hibag_probability_folder <- paste0(docs_folder, "/probability_densities/hibag/")
cookhla_probability_folder <- paste0(docs_folder, "/probability_densities/cookhla/")
mendelian_error_table <- paste0(docs_folder, "/mendelian_error_rates")
md_file <- paste0(docs_folder, "/report.md")

write(
  x = paste("#", title),
  file = md_file,
  append = F
)


# Extract imputed data

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

cookhla_a <- subset(cookhla, hla == "A")
cookhla_b <- subset(cookhla, hla == "B")
cookhla_c <- subset(cookhla, hla == "C")
cookhla_dqb1 <- subset(cookhla, hla == "DQB1")
cookhla_drb1 <- subset(cookhla, hla == "DRB1")


# Extract trios

fam <- read.table(
    file = fam_file,
    header = F,
    col.names = c("fid", "iid", "pat", "mat", "sex", "phen")
)

trios <- subset(fam, !is.na(pat) & pat != "0" & !is.na(mat) & mat != "0")

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

# Check for Mendelian errors

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

# Write mendelian error table to files

write(
    x = "hla\tsoftware\tn_trios\tn_erros\terror_rate",
    file = mendelian_error_table,
    append = F
)


write(
    x = "## Mendelian",
    file = md_file,
    append = T
)

write(
    x = "| HLA | Software | Trios | Errors | Error rate |",
    file = md_file,
    append = T
)

write(
    x = "| --- | --- | --- | --- | --- |",
    file = md_file,
    append = T
)


add_table_row <- function(trios_alleles, software, hla) {
    n_errors <- sum(!trios_alleles$mendel)
    n_rows <- nrow(trios_alleles)
    error_rate <- n_errors / n_rows
    write(
        x = paste0(hla, "\t", software, "\t", n_rows, "\t", n_errors, "\t", error_rate),
        file = mendelian_error_table,
        append = T
    )
    write(
        x = paste("|",hla,"|", software, "|", n_rows, "|", n_errors, "|", error_rate, "|"),
        file = md_file,
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



# Plot allele frequencies and probability densities

extract_common_allele_freqs <- function(allele_table, ref) {
    allele_freqs <- table(c(allele_table$allele1, allele_table$allele2))
    norm_freqs <- allele_freqs / sum(allele_freqs)
    filtered_freqs <- norm_freqs[names(norm_freqs) %in% ref$allele]
    filtered_freqs <- c(filtered_freqs, "rare" = 1 - sum(filtered_freqs, na.rm = TRUE))
    return(filtered_freqs)
}


plot_frequencies <- function(cookhla, hibag, ref, hla, folder) {
    freq_matrix <- rbind(cookhla[ref$allele], hibag[ref$allele], ref$freq)
    png(paste0(folder, "frequencies_HLA_", hla, ".png"), width = 800, height = 600)
    barplot(freq_matrix,
        beside = TRUE,
        main = paste0("Allele frequencies for HLA-", hla),
        xlab = "Allele", ylab = "Frequency", names.arg = ref$allele,
        col = c("lightblue", "pink", "darkgreen"), las = 2
    )
    legend("topright", legend = c("CookHLA", "HIBAG", "Reference"), fill = c("lightblue", "pink", "darkgreen"))
    dev.off()
}

plot_hibag_ref <- function(hibag, ref, hla, folder) {
    freq_matrix <- rbind(hibag[ref$allele], ref$freq)
    png(paste0(folder, "frequencies_HLA_", hla, ".png"), width = 800, height = 600)
    barplot(freq_matrix,
        beside = TRUE,
        main = paste0("Allele frequencies for HLA-", hla),
        xlab = "Allele", ylab = "Frequency", names.arg = ref$allele,
        col = c("pink", "darkgreen"), las = 2
    )
    legend("topright", legend = c("HIBAG", "Reference"), fill = c("pink", "darkgreen"))
    dev.off()
}

plot_prob_density <- function(probs, title, filename, folder) {
    png(paste0(folder, filename))
    plot(density(probs),
        main = title,
        xlab = "Probability",
        ylab = "Density"
    )
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


common_cookhla_a <- extract_common_allele_freqs(cookhla_a, ref_a)
common_cookhla_b <- extract_common_allele_freqs(cookhla_b, ref_b)
common_cookhla_c <- extract_common_allele_freqs(cookhla_c, ref_c)
common_cookhla_dqb1 <- extract_common_allele_freqs(cookhla_dqb1, ref_dqb1)
common_cookhla_drb1 <- extract_common_allele_freqs(cookhla_drb1, ref_drb1)


plot_frequencies(common_cookhla_a, common_hibag_a, ref_a, "A", frequencies_folder)
plot_frequencies(common_cookhla_b, common_hibag_b, ref_b, "B", frequencies_folder)
plot_frequencies(common_cookhla_c, common_hibag_c, ref_c, "C", frequencies_folder)
plot_hibag_ref(common_hibag_dpb1, ref_dpb1, "DPB1", frequencies_folder)
plot_frequencies(common_cookhla_dqb1, common_hibag_dqb1, ref_dqb1, "DQB1", frequencies_folder)
plot_frequencies(common_cookhla_drb1, common_hibag_drb1, ref_drb1, "DRB1", frequencies_folder)

plot_prob_density(hibag_a$prob, "HIBAG HLA-A probability density", "hibag_a_probability_density.png", hibag_probability_folder)
plot_prob_density(hibag_b$prob, "HIBAG HLA-B probability density", "hibag_b_probability_density.png", hibag_probability_folder)
plot_prob_density(hibag_c$prob, "HIBAG HLA-C probability density", "hibag_c_probability_density.png", hibag_probability_folder)
plot_prob_density(hibag_dpb1$prob, "HIBAG HLA-DPB1 probability density", "hibag_dpb1_probability_density.png", hibag_probability_folder)
plot_prob_density(hibag_dqb1$prob, "HIBAG HLA-DQB1 probability density", "hibag_dqb1_probability_density.png", hibag_probability_folder)
plot_prob_density(hibag_drb1$prob, "HIBAG HLA-DRB1 probability density", "hibag_drb1_probability_density.png", hibag_probability_folder)

plot_prob_density(cookhla_a$confidence, "CookHLA HLA-A probability density", "cookhla_a_probability_density.png", cookhla_probability_folder)
plot_prob_density(cookhla_b$confidence, "CookHLA HLA-B probability density", "cookhla_b_probability_density.png", cookhla_probability_folder)
plot_prob_density(cookhla_c$confidence, "CookHLA HLA-C probability density", "cookhla_c_probability_density.png", cookhla_probability_folder)
plot_prob_density(cookhla_dqb1$confidence, "CookHLA HLA-DQB1 probability density", "cookhla_dqb1_probability_density.png", cookhla_probability_folder)
plot_prob_density(cookhla_drb1$confidence, "CookHLA HLA-DRB1 probability density", "cookhla_drb1_probability_density.png", cookhla_probability_folder)