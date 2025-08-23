
library(HIBAG)

#### INPUT FILE ####
args <- commandArgs(TRUE)
input_trunk <- args[1]
output_trunk <- args[2]
modelfile <- args[3]
threads <- args[4]

#### MODEL CREATION ####
mlst<-get(load(modelfile))
geno <- hlaBED2Geno(bed.fn=paste(input_trunk,".bed",sep = ""), 
                        fam.fn=paste(input_trunk,".fam",sep = ""), 
                        bim.fn=paste(input_trunk,".bim",sep = ""))

#### IMPUTATION AND ALLELES FILES FOR EACH GENE ####
model_a<-hlaModelFromObj(mlst$A)
hla_a<-hlaPredict(model_a,geno, cl=threads)

write.csv(hla_a$value, paste(output_trunk,"_A.csv",sep = ""), row.names=FALSE)

model_b<-hlaModelFromObj(mlst$B)
hla_b<-hlaPredict(model_b,yourgeno, cl=threads)

write.csv(hla_b$value, paste(output_trunk,"_B.csv",sep = ""), row.names=FALSE)

model_c<-hlaModelFromObj(mlst$C)
hla_c<-hlaPredict(model_c,geno, cl=threads)

write.csv(hla_c$value, paste(output_trunk,"_C.csv",sep = ""), row.names=FALSE)

model_drb1<-hlaModelFromObj(mlst$DRB1)
hla_drb1<-hlaPredict(model_drb1,geno, cl=threads)

write.csv(hla_drb1$value, paste(output_trunk,"_DRB1.csv",sep = ""), row.names=FALSE)

model_dqb1<-hlaModelFromObj(mlst$DQB1)
hla_dqb1<-hlaPredict(model_dqb1,geno, cl=threads)

write.csv(hla_dqb1$value, paste(output_trunk,"_DQB1.csv",sep = ""), row.names=FALSE)

model_dpb1<-hlaModelFromObj(mlst$DPB1)
hla_dpb1<-hlaPredict(model_dpb1,geno, cl=threads)

write.csv(hla_dpb1$value, paste(output_trunk,"_DPB1.csv",sep = ""), row.names=FALSE)

#### CREATING AND WRITING FULL SETS ####
data_merge1 <- merge(hla_a$value, hla_b$value, by = "sample.id")
data_merge2 <- merge(data_merge1, hla_c$value, by = "sample.id")
data_merge3 <- merge(data_merge2, hla_drb1$value, by = "sample.id")
data_merge4 <- merge(data_merge3, hla_dqb1$value, by = "sample.id")
fullSet <- merge(data_merge4, hla_dpb1$value, by = "sample.id")

fullSet_anon <- fullSet[,-1]

write.csv(fullSet, paste(output_trunk,"_fullAlleles",sep = ""), row.names=FALSE)
write.csv(fullSet_anon, paste(output_trunk,"_fullAlleles_anon",sep = ""), row.names=FALSE)

