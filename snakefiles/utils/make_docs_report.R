library(dplyr)
library(tidyr)

debug <- T

if (debug) {
    args <- c(
        "/home/oystein/hla_imputation_pipeout/2025.09.16/HIBAG/hibag",
        "/home/oystein/hla_imputation_pipeout/2025.09.16/CookHLA/cookhla_output.MHC.HLA_IMPUTATION_OUT.alleles",
        "/home/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fam",
        "/home/oystein/github/HLA-imputation/snakefiles/resources/norwegian_allele_frequencies/common/HLA",
        "/home/oystein/test/2025.09.16/docs/",
        "/home/oystein/test/2025.09.16/pipeout/",
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

plots_folder <- paste0(docs_folder, "plots/")
mendelian_error_table <- paste0(docs_folder, "mendelian_error_rates")
md_file <- paste0(docs_folder, "report.md")
inconsistencies_pipeout_folder <- paste0(pipeout_folder, "inconsistencies/")
inconsistencies_docs_folder <- paste0(docs_folder, "inconsistencies/")


if(!dir.exists(docs_folder)){
    dir.create(docs_folder)
}


if(!dir.exists(pipeout_folder)){
    dir.create(pipeout_folder)
}

if(!dir.exists(plots_folder)){
    dir.create(plots_folder)
}

if(!dir.exists(inconsistencies_pipeout_folder)){
    dir.create(inconsistencies_pipeout_folder)
}
if(!dir.exists(inconsistencies_docs_folder)){
    dir.create(inconsistencies_docs_folder)
}


if (!file.exists(md_file)) {
  file.create(md_file)
}

write(
  x = paste("#", title),
  file = md_file,
  append = F
)

shared_hla <- c("A", "B", "C", "DQB1", "DRB1")


# Extract imputed data

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
    x = "## Mendelian errors (trios only)",
    file = md_file,
    append = T
)

write(
    x = "| HLA | Software | Trios | Errors | Error rate |\n| --- | --- | --- | --- | --- |",
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

hibag_a <- hibag_a %>% left_join(trios_hibag_a %>% select(iid, mendel), by = "iid")
hibag_b <- hibag_b %>% left_join(trios_hibag_b %>% select(iid, mendel), by = "iid")
hibag_c <- hibag_c %>% left_join(trios_hibag_c %>% select(iid, mendel), by = "iid")
hibag_dpb1 <- hibag_dpb1 %>% left_join(trios_hibag_dpb1 %>% select(iid, mendel), by = "iid")
hibag_dqb1 <- hibag_dqb1 %>% left_join(trios_hibag_dqb1 %>% select(iid, mendel), by = "iid")
hibag_drb1 <- hibag_drb1 %>% left_join(trios_hibag_drb1 %>% select(iid, mendel), by = "iid")


cookhla_a <- cookhla_a %>% left_join(trios_cookhla_a %>% select(iid, mendel), by = "iid")
cookhla_b <- cookhla_b %>% left_join(trios_cookhla_b %>% select(iid, mendel), by = "iid")
cookhla_c <- cookhla_c %>% left_join(trios_cookhla_c %>% select(iid, mendel), by = "iid")
cookhla_dqa1 <- cookhla_dqa1 %>% left_join(trios_cookhla_dqa1 %>% select(iid, mendel), by = "iid")
cookhla_dqb1 <- cookhla_dqb1 %>% left_join(trios_cookhla_dqb1 %>% select(iid, mendel), by = "iid")
cookhla_drb1 <- cookhla_drb1 %>% left_join(trios_cookhla_drb1 %>% select(iid, mendel), by = "iid")

# HIBAG/CookHLA consistency check

consistency_check <- function(row) {
    hibag1 <- row["hibag1"]
    hibag2 <- row["hibag2"]
    cook1 <- row["cook1"]
    cook2 <- row["cook2"]
    if((hibag1 == cook1 && hibag2 == cook2) || (hibag1 == cook2 && hibag2 == cook1)){
        inconsistencies <- 0
    }
    else if((hibag1 == cook1 && hibag2 != cook2) || (hibag1 == cook2 && hibag2 != cook1) || (hibag2 == cook2 && hibag1 != cook1) || (hibag2 == cook1 && hibag1 != cook2) ){
        inconsistencies <- 1
    }
    else{
        inconsistencies <- 2
    }

    return(inconsistencies)
}

find_inconsistencies <- function(row){
    hibag1 <- row["hibag1"]
    hibag2 <- row["hibag2"]
    cook1 <- row["cook1"]
    cook2 <- row["cook2"]
    inconsistent <- ""
    if(hibag1 == cook1 && hibag2 != cook2){
        inconsistent <- paste0(hibag2,",",cook2)#c(hibag2, cook2)
    }
    else if(hibag1 == cook2 && hibag2 != cook1){
        inconsistent <- paste0(hibag2,",",cook1)# c(hibag2, cook1)
    }
    else if(hibag2 == cook2 && hibag1 != cook1){
        inconsistent <- paste0(hibag1, ",", cook1)
    }
    else if(hibag2 == cook1 && hibag1 != cook2){
        inconsistent <- paste0(hibag1,",", cook2)
    }
    return(inconsistent)
}

split_inconsistent_alleles <- function(inconsistent_table){
    inconsistent_table <- inconsistent_table %>%
    separate(inconsistent_alleles, into = c("hibag", "cook"), sep = ",")
    return(inconsistent_table)
}

format_and_compare <- function(hibag, cookhla, hla){
    compare <- cookhla %>% select(iid, hla, cook1 = allele1, cook2 = allele2, cook_post1 = post1, cook_post2 = post2, cook_confidence = confidence, cook_mendel = mendel) %>% full_join(hibag %>% select(iid, hibag1 = allele1, hibag2 = allele2, hibag_prob = prob, hibag_matching = matching, hibag_mendel = mendel), by="iid")
    compare$inconsistencies <- apply(compare, 1, consistency_check)
    single_inconsistencies <- subset(compare, inconsistencies == 1) %>% select(hibag1, hibag2, cook1, cook2)
    single_inconsistencies$inconsistent_alleles <- apply(single_inconsistencies, 1, find_inconsistencies)
    single_inconsistencies <- split_inconsistent_alleles(single_inconsistencies) %>% select(hibag, cook)
    single_inconsistencies <- cbind(hla=hla, single_inconsistencies)
    double_inconsistencies <- subset(compare, inconsistencies == 2) %>% select(hibag1, hibag2, cook1, cook2)
    double_inconsistencies <- cbind(hla=hla, double_inconsistencies)
    inconsistency_table <- table(single_inconsistencies$hibag, single_inconsistencies$cook)
    return(list(compare, single_inconsistencies, double_inconsistencies, inconsistency_table))
}

write(
    x = "## HIBAG/CookHLA consistency check",
    file = md_file,
    append = T
)

write_inconsistency_table_row <- function(compare, hla){
    n_zero <- nrow(subset(compare, inconsistencies == 0))
    n_one <- nrow(subset(compare, inconsistencies == 1))
    n_two <- nrow(subset(compare, inconsistencies == 2))
    write(
        x = paste("|", hla,"|", n_zero, "|", n_one, "|", n_two, "|"),
        file = md_file,
        append = T
    )
}


write(
x = "### Number of samples with 0, 1 and 2 inconsistencies between HIBAG and CookHLA",
file = md_file,
append = T
)
write(
x = "| HLA | 0 | 1 | 2 |\n| --- | --- | --- | --- |",
file = md_file,
append = T
)

results_a <- format_and_compare(hibag_a, cookhla_a, "A")
results_b <- format_and_compare(hibag_b, cookhla_b, "B")
results_c <- format_and_compare(hibag_c, cookhla_c, "C")
results_dqb1 <- format_and_compare(hibag_dqb1, cookhla_dqb1, "DQB1")
results_drb1 <- format_and_compare(hibag_drb1, cookhla_drb1, "DRB1")

write_inconsistency_table_row(results_a[[1]], "A")
write_inconsistency_table_row(results_b[[1]], "B")
write_inconsistency_table_row(results_c[[1]], "C")
write_inconsistency_table_row(results_dqb1[[1]], "DQB1")
write_inconsistency_table_row(results_drb1[[1]], "DRB1")





# Plot allele frequencies and probability densities

write(
x = "## Allele frequencies and probability density plots",
file = md_file,
append = T
)

extract_common_allele_freqs <- function(allele_table, ref) {
    allele_freqs <- table(c(allele_table$allele1, allele_table$allele2))
    norm_freqs <- allele_freqs / sum(allele_freqs)
    filtered_freqs <- norm_freqs[names(norm_freqs) %in% ref$allele]
    filtered_freqs <- c(filtered_freqs, "rare" = 1 - sum(filtered_freqs, na.rm = TRUE))
    return(filtered_freqs)
}


plot_frequencies <- function(cookhla, hibag, ref, hla, absolute_filepath, relative_filepath) {
    freq_matrix <- rbind(cookhla[ref$allele], hibag[ref$allele], ref$freq)
    png(absolute_filepath, width = 800, height = 600)
    barplot(freq_matrix,
        beside = TRUE,
        main = paste0("Allele frequencies for HLA-", hla),
        xlab = "Allele", ylab = "Frequency", names.arg = ref$allele,
        col = c("lightblue", "pink", "darkgreen"), las = 2
    )
    legend("topright", legend = c("CookHLA", "HIBAG", "Reference"), fill = c("lightblue", "pink", "darkgreen"))
    dev.off()
    write(
        x = paste0("![](",relative_filepath,")"),
        file = md_file,
        append = T
    )
}

plot_hibag_ref <- function(hibag, ref, hla, absolute_filepath, relative_filepath) {
    freq_matrix <- rbind(hibag[ref$allele], ref$freq)
    png(absolute_filepath, width = 800, height = 600)
    barplot(freq_matrix,
        beside = TRUE,
        main = paste0("Allele frequencies for HLA-", hla),
        xlab = "Allele", ylab = "Frequency", names.arg = ref$allele,
        col = c("pink", "darkgreen"), las = 2
    )
    legend("topright", legend = c("HIBAG", "Reference"), fill = c("pink", "darkgreen"))
    dev.off()
    write(
        x = paste0("![](",relative_filepath,")"),
        file = md_file,
        append = T
    )
}

plot_prob_density <- function(probs, title, absolute_filepath, relative_filepath) {
    png(filepath)
    plot(density(absolute_filepath),
        main = title,
        xlab = "Probability",
        ylab = "Density"
    )
    dev.off()
    write(
        x = paste0("![](",relative_filepath,")"),
        file = md_file,
        append = T
    )
}

add_rare_row <- function(ref) {
    rare_freq <- 1 - sum(ref$freq)
    rare_row <- data.frame(allele = "rare", freq = rare_freq)
    ref_rare <- rbind(ref, rare_row)
    return(ref_rare)
}


plot_frequencies_and_densities <- function(ref, hibag, cookhla, hla){
    write(
        x = paste0("### HLA-", hla),
        file = md_file,
        append = T
    )
    ref <- add_rare_row(read.table(paste0(ref_trunk_file, ".", hla), header = T))
    common_hibag <- extract_common_allele_freqs(hibag, ref)
    common_cookhla <- extract_common_allele_freqs(cookhla, ref)
    freq_filename <- paste0("frequencies_HLA-", hla, ".png")
    freq_absolute_path <- paste0(plots_folder, freq_filename)
    freq_relative_path <- paste0("plots/", freq_filename)
    hibag_prob_filename <- paste0("hibag_probabilities_HLA-", hla, ".png")
    cook_prob_filename <- paste0("cookhla_probabilities_HLA-", hla, ".png")
    hibag_prob_aboslute_filepath <- paste0(plots_folder, hibag_prob_filename)
    hibag_prob_relative_filepath <- paste0("plots/", hibag_prob_filename)
    cook_prob_absolute_filepath <- paste0(plots_folder, cook_prob_filename)
    cook_prob_relative_filepath <- paste0("plots/", cook_prob_filename)

    plot_frequencies(common_cookhla, common_hibag, ref, "A", freq_absolute_path, freq_relative_path)
    plot_prob_density(hibag$prob, paste0("HIBAG HLA-",hla," probability density"), hibag_prob_aboslute_filepath, hibag_prob_relative_filepath)
    plot_prob_density(cookhla$confidence, paste0("CookHLA HLA-",hla, " probability density"), cook_prob_absolute_filepath, cook_prob_relative_filepath)

}

plot_frequencies_and_densities(ref_trunk_file, hibag_a, cookhla_a, "A")
plot_frequencies_and_densities(ref_trunk_file, hibag_b, cookhla_b, "B")
plot_frequencies_and_densities(ref_trunk_file, hibag_c, cookhla_c, "C")
plot_frequencies_and_densities(ref_trunk_file, hibag_dqb1, cookhla_dqb1, "DQB1")
plot_frequencies_and_densities(ref_trunk_file, hibag_drb1, cookhla_drb1, "DRB1")


dpb1_freq_filename <- "frequencies_HLA-DPB1.png"
dpb1_freq_absolute_filepath <- paste0(plots_folder, dpb1_freq_filename)
dpb1_freq_relative_filepath <- paste0("plots/", dpb1_freq_filename)



dpb1_density_filename <- "hibag_probabilities_HLA-DPB1.png"

dpb1_density_absolute_filepath <- paste0(plots_folder, dpb1_density_filename)
dpb1_density_relative_filepath <- paste0("plots/", dpb1_density_filename)
write(
        x = paste0("### HLA-DPB1"),
        file = md_file,
        append = T
    )
ref_dpb1 <- add_rare_row(read.table(paste0(ref_trunk_file, ".DPB1"), header = T))
common_hibag_dpb1 <- extract_common_allele_freqs(hibag_dpb1, ref_dpb1)
plot_hibag_ref(common_hibag_dpb1, ref_dpb1, "DPB1", dpb1_freq_absolute_filepath, dpb1_density_relative_filepath
plot_prob_density(hibag_dpb1$prob, "HIBAG HLA-DPB1 probability density", dpb1_density_absolute_filepath, dpb1_density_relative_filepath)



dqa1_density_filename <- "cookhla_probabilities_HLA-DQA1.png"
dqa1_density_absolute_filepath <- paste0(plots_folder, dqa1_density_filename)
dqa1_density_relative_filepath <- paste0("plots/", dqa1_density_filename)
write(
        x = paste0("### HLA-DQA1"),
        file = md_file,
        append = T
    )
plot_prob_density(cookhla_dqa1$confidence, "CookHLA HLA-DQA1 probability density", dqa1_density_absolute_filepath, dqa1_density_relative_filepath)


# Merge results

merged_shared <- rbind(results_a[[1]], results_b[[1]], results_c[[1]], results_dqb1[[1]], results_drb1[[1]])
hibag_dpb1_formatted <- cbind(hibag_dpb1[, 1, drop = FALSE], hla = "DPB1", hibag_dpb1[, -1, drop = FALSE]) %>% select(iid, hla, hibag1 = allele1, hibag2 = allele2, hibag_prob = prob, hibag_matching = matching)
cook_dqa1_formatted <- cookhla_dqa1 %>% select(iid, hla, cook1 = allele1, cook2 = allele2, cook_post1 = post1, cook_post2 = post2, cook_confidence = confidence)
merged_all <- bind_rows(merged_shared, hibag_dpb1_formatted, cook_dqa1_formatted)

merged_single_inconsistencies <- rbind(results_a[[2]], results_b[[2]], results_c[[2]], results_dqb1[[2]], results_drb1[[2]])
merged_double_inconsistencies <- rbind(results_a[[3]], results_b[[3]], results_c[[3]], results_dqb1[[3]], results_drb1[[3]])

write.table(x = merged_single_inconsistencies, file = paste0(inconsistencies_pipeout_folder, "single_inconsistencies_samples"), sep="\t", quote=F, row.names=F)
write.table(x = merged_double_inconsistencies, file = paste0(inconsistencies_pipeout_folder, "double_inconsistencies_samples"), sep="\t", quote=F, row.names=F)
write.table(x = merged_all, file = paste0(pipeout_folder,"merged_all"), sep="\t", quote=F, row.names=F)

write_inconsistency_counts <- function(inconsistency_counts, hla){
    write.table(x=inconsistency_counts, file = paste0(inconsistencies_docs_folder, "single_inconsistencies_counts_HLA-", hla), sep="\t", quote=F)
    write.table(x=inconsistency_counts, file = paste0(inconsistencies_pipeout_folder, "single_inconsistencies_counts_HLA-", hla), sep="\t", quote=F)
}

write_inconsistency_counts(results_a[[4]], "A")
write_inconsistency_counts(results_b[[4]], "B")
write_inconsistency_counts(results_c[[4]], "C")
write_inconsistency_counts(results_dqb1[[4]], "DQB1")
write_inconsistency_counts(results_drb1[[4]], "DRB1")