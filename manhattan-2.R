setwd("/Users/apple/Documents/03 R文件/drawfile/")
library(gtools)
library(ggplot2)
library(scales)
library(qqman)
source("manhattanplot_function.R")
glm <- read.table('glm_acid.txt', header = T, comment.char = "", stringsAsFactors = F, check.names = F,fill = T)
traitID <- unique(glm$Trait)
for(tr in traitID){
  outfile <- paste0(tr, "_glm")
  sub <- subset(glm, Trait == tr)[,c(2:4,6)]
  write.table(sub, file = paste0(outfile, ".txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  threshold <-c( -log10(0.05/nrow(sub[!is.na(sub$p),])),-log10(0.1/nrow(sub[!is.na(sub$p),])))
  manhattanplot(mydata = paste0(outfile, ".txt"), key = outfile, chr = "all",thresholds = threshold,columns = 2:4,log10 = T, x_tick_labs = "Character")
}
require(Cairo)
CairoPNG("qqplot") 
qq(mlm$p)
dev.off()
