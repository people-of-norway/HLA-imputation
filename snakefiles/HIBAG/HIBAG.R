
library(HIBAG)
debug <- T

if (debug) {
    input_trunk <- "/home/oystein/hla_imputation_pipeout/2025.09.16/extract_hla"
    output_trunk <- "/home/oystein/test/hibag_test"
    modelfile <- "/home/oystein/github/HLA-imputation/snakefiles/HIBAG/InfiniumOmniExpress-24-European-HLA4-hg19.RData"
    threads <- 32
} else {
    args <- commandArgs(TRUE)
    input_trunk <- args[1]
    output_trunk <- args[2]
    modelfile <- args[3]
    threads <- as.integer(args[4])
}

# #### INPUT FILE ####
# args <- commandArgs(TRUE)
# input_trunk <- args[1]
# output_trunk <- args[2]
# modelfile <- args[3]
# threads <- as.integer(args[4])

#### MODEL CREATION ####
mlst<-get(load(modelfile))
geno <- hlaBED2Geno(bed.fn=paste(input_trunk,".bed",sep = ""), 
                        fam.fn=paste(input_trunk,".fam",sep = ""), 
                        bim.fn=paste(input_trunk,".bim",sep = ""))
|
#### IMPUTATION AND ALLELES FILES FOR EACH GENE ####

create_model_predict_and_write <- function(object, geno, output_hla, output_aa, threads){
    model <- hlaModelFromObj(object)
    hla <- hlaPredict(model, geno, cl=threads)
    hla_conv <- hlaConvSequence(hla = hla, code = "P.code.merge")
    write.table(x=hla$value, file=output_hla, col.names=T, row.names=F, sep = "\t", quote=F)
    write.table(x=hla_conv$value, file=output_aa, col.names=T, row.names=F, sep = "\t", quote=F)
}
create_model_predict_and_write(mlst$A, geno, paste0(output_trunk,".A"), paste0(output_trunk,".aa.A"), threads)
create_model_predict_and_write(mlst$B, geno, paste0(output_trunk,".B"), paste0(output_trunk,".aa.B"), threads)
create_model_predict_and_write(mlst$C, geno, paste0(output_trunk,".C"), paste0(output_trunk,".aa.C"), threads)
create_model_predict_and_write(mlst$DRB1, geno, paste0(output_trunk,".DRB1"), paste0(output_trunk,".aa.DRB1"), threads)
create_model_predict_and_write(mlst$DQA1, geno, paste0(output_trunk,".DQA1"), paste0(output_trunk,".aa.DQA1"), threads)
create_model_predict_and_write(mlst$DQB1, geno, paste0(output_trunk,".DQB1"), paste0(output_trunk,".aa.DQB1"), threads)
create_model_predict_and_write(mlst$DPB1, geno, paste0(output_trunk,".DPB1"), paste0(output_trunk,".aa.DPB1"), threads)  