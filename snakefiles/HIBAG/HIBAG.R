
library(HIBAG)

#### INPUT FILE ####
args <- commandArgs(TRUE)
input_trunk <- args[1]
output_trunk <- args[2]
modelfile <- args[3]
threads <- as.integer(args[4])

#### MODEL CREATION ####
mlst<-get(load(modelfile))
geno <- hlaBED2Geno(bed.fn=paste(input_trunk,".bed",sep = ""), 
                        fam.fn=paste(input_trunk,".fam",sep = ""), 
                        bim.fn=paste(input_trunk,".bim",sep = ""))

#### IMPUTATION AND ALLELES FILES FOR EACH GENE ####


create_model_predict_and_write <- function(object, geno, output_full, threads){
    model <- hlaModelFromObj(object)
    hla <- hlaPredict(model, geno, cl=threads)
    write.table(x=hla$value, file=output_full, col.names=T row.names=F, sep = "\t", quote=F)
    return(hla)
}
hla_a <- create_model_predict_and_write(mlst$A, geno, paste0(output_trunk,".A"), threads)
hla_b <- create_model_predict_and_write(mlst$B, geno, paste0(output_trunk,".B"), threads)
hla_c <- create_model_predict_and_write(mlst$C, geno, paste0(output_trunk,".C"), threads)
hla_drb1 <- create_model_predict_and_write(mlst$DRB1, geno, paste0(output_trunk,".DRB1"), threads)
hla_dqb1 <- create_model_predict_and_write(mlst$DQB1, geno, paste0(output_trunk,".DQB1"), threads)
hla_dpb1 <- create_model_predict_and_write(mlst$DPB1, geno, paste0(output_trunk,".DPB1"), threads)

#### CREATING AND WRITING FULL SETS ####
data_merge1 <- merge(hla_a$value, hla_b$value, by = "sample.id")
data_merge2 <- merge(data_merge1, hla_c$value, by = "sample.id")
data_merge3 <- merge(data_merge2, hla_drb1$value, by = "sample.id")
data_merge4 <- merge(data_merge3, hla_dqb1$value, by = "sample.id")
fullset <- merge(data_merge4, hla_dpb1$value, by = "sample.id")

fullset_anon <- fullSet[,-1]

write.table(x=fullset, file=paste0(output_trunk,".fullset"), col.names=T, row.names=F, sep = "\t", quote=F)
write.table(x=fullset_anon, file=paste0(output_trunk,".fullset_anon"), col.names=T, row.names=F, sep = "\t", quote=F)

