
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
    write.table(x=hla$value, file=output_full, col.names=T, row.names=F, sep = "\t", quote=F)
}
# create_model_predict_and_write(mlst$A, geno, paste0(output_trunk,".A"), threads)
# create_model_predict_and_write(mlst$B, geno, paste0(output_trunk,".B"), threads)
create_model_predict_and_write(mlst$C, geno, paste0(output_trunk,".C"), threads)
create_model_predict_and_write(mlst$DRB1, geno, paste0(output_trunk,".DRB1"), threads)
create_model_predict_and_write(mlst$DQB1, geno, paste0(output_trunk,".DQB1"), threads)
create_model_predict_and_write(mlst$DPB1, geno, paste0(output_trunk,".DPB1"), threads)

