library(dplyr)
library(tidyr)

debug <- T

if (debug) {
    args <- c(
        "/home/oystein/hla_imputation_pipeout/2025.08.23/HIBAG/hibag",
        "/home/oystein/hla_imputation_pipeout/2025.08.23/CookHLA/cookhla_output.MHC.HLA_IMPUTATION_OUT.alleles"
    )
} else {
    args <- commandArgs(TRUE)
}

hibag_trunk_file <- args[1]
cookhla_file <- args[2]

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



compare_a <- hibag_a %>% select(iid, hibag1 = allele1, hibag2 = allele2) %>% left_join(cookhla_a %>% select(iid, cook1 = allele1, cook2 = allele2), by = "iid")
compare_b <- hibag_b %>% select(iid, hibag1 = allele1, hibag2 = allele2) %>% left_join(cookhla_b %>% select(iid, cook1 = allele1, cook2 = allele2), by = "iid")
compare_a$inconsistencies <- apply(compare_a, 1, consistency_check)
compare_b$inconsistencies <- apply(compare_b, 1, consistency_check)

inconsistent_a <- subset(compare_a, inconsistencies == 1)
inconsistent_a$inconsistent_alleles <- apply(inconsistent_a, 1, find_inconsistencies)
inconsistent_a <- split_inconsistent_alleles(inconsistent_a)
double_inconsistent_a <- subset(compare_a, inconsistencies == 2)

inconsistency_table_a <- table(inconsistent_a$hibag, inconsistent_a$cook)