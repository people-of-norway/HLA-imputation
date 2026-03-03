debug <- F

if (debug) {
    in_trunk <- "/home/oystein/hla_imputation_pipeout/2026.02.24/HIBAG/hibag"
    output <- "/home/oystein/test/hibag.chped"
} else {
    args <- commandArgs(TRUE)
    in_trunk <- args[1]
    output <- args[2]
}

library(dplyr)
library(tidyr)  

hibag_a <- read.table(paste0(in_trunk, ".A"), header=T, sep = "\t")
hibag_b <- read.table(paste0(in_trunk, ".B"), header=T, sep = "\t")
hibag_c <- read.table(paste0(in_trunk, ".C"), header=T, sep = "\t")
hibag_drb1 <- read.table(paste0(in_trunk, ".DRB1"), header=T, sep = "\t")
hibag_dqa1 <- read.table(paste0(in_trunk, ".DQA1"), header=T, sep = "\t")
hibag_dqb1 <- read.table(paste0(in_trunk, ".DQB1"), header=T, sep = "\t")
hibag_dpb1 <- read.table(paste0(in_trunk, ".DPB1"), header=T, sep = "\t")   

chped <- hibag_a %>% select(FID = sample.id, IID = sample.id)
chped$PID <- 0
chped$MID <- 0
chped$Sex <- 0
chped$Phe <- -9
hibag_a$A_1 <- paste0("A*", hibag_a$allele1)
hibag_a$A_2 <- paste0("A*", hibag_a$allele2)
hibag_b$B_1 <- paste0("B*", hibag_b$allele1)
hibag_b$B_2 <- paste0("B*", hibag_b$allele2)
hibag_c$C_1 <- paste0("C*", hibag_c$allele1)
hibag_c$C_2 <- paste0("C*", hibag_c$allele2)
hibag_drb1$DRB1_1 <- paste0("DRB1*", hibag_drb1$allele1)
hibag_drb1$DRB1_2 <- paste0("DRB1*", hibag_drb1$allele2)
hibag_dqa1$DQA1_1 <- paste0("DQA1*", hibag_dqa1$allele1)
hibag_dqa1$DQA1_2 <- paste0("DQA1*", hibag_dqa1$allele2)
hibag_dqb1$DQB1_1 <- paste0("DQB1*", hibag_dqb1$allele1)
hibag_dqb1$DQB1_2 <- paste0("DQB1*",    hibag_dqb1$allele2)
hibag_dpb1$DPB1_1 <- paste0("DPB1*", hibag_dpb1$allele1)
hibag_dpb1$DPB1_2 <- paste0("DPB1*", hibag_dpb1$allele2)

chped <- chped %>% left_join(hibag_a %>% select(sample.id, A_1, A_2), by = c("IID" = "sample.id")) %>%
    left_join(hibag_b %>% select(sample.id, B_1, B_2), by = c("IID" = "sample.id")) %>%
    left_join(hibag_c %>% select(sample.id, C_1, C_2), by = c("IID" = "sample.id")) %>%
    left_join(hibag_drb1 %>% select(sample.id, DRB1_1, DRB1_2), by = c("IID" = "sample.id")) %>%
    left_join(hibag_dqa1 %>% select(sample.id, DQA1_1, DQA1_2), by = c("IID" = "sample.id")) %>%
    left_join(hibag_dqb1 %>% select(sample.id, DQB1_1, DQB1_2), by = c("IID" = "sample.id")) %>%
    left_join(hibag_dpb1 %>% select(sample.id, DPB1_1, DPB1_2), by = c("IID" = "sample.id"))

chped$DPA1_1 <- 0
chped$DPA1_2 <- 0

chped <- chped %>% select(FID, IID, PID, MID, Sex, Phe, A_1, A_2, B_1, B_2, C_1, C_2, DPA1_1, DPA1_2, DPB1_1, DPB1_2, DQA1_1, DQA1_2, DQB1_1, DQB1_2,DRB1_1, DRB1_2)
write.table(chped, file=output, col.names=T, row.names=F, sep = "\t", quote=F)