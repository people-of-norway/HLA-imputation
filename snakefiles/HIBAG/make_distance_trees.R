library(HIBAG)
debug <- T

if (debug) {
    output_trunk <- "/home/oystein/test/distance_tree"
    modelfile <- "/home/oystein/github/HLA-imputation/snakefiles/HIBAG/InfiniumOmniExpress-24-European-HLA4-hg19.RData"
} else {
    args <- commandArgs(TRUE)
    output_trunk <- args[1]
    modelfile <- args[2]
}




#### MODEL CREATION ####
mlst<-get(load(modelfile))


make_tree <- function(object, gene, output_trunk){
    model <- hlaModelFromObj(object)
    d <- hlaDistance(model)
    p <- hclust(as.dist(d))
    png(paste0(output_trunk, ".", gene, ".png"), width=1400, height=800)
    plot(p, main=paste0("HLA-", gene), xlab="", sub="")
    dev.off()
}
make_tree(mlst$A, "A", output_trunk)
make_tree(mlst$B, "B", output_trunk)
make_tree(mlst$C, "C", output_trunk)
make_tree(mlst$DRB1, "DRB1", output_trunk)
make_tree(mlst$DQA1, "DQA1", output_trunk)
make_tree(mlst$DQB1, "DQB1", output_trunk)
make_tree(mlst$DPB1, "DPB1", output_trunk)