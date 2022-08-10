library(data.table)

get_names <- function(df){
  df$promoter_name <- promoters[df$promoter,]$V4
  df$promoter_chr <- promoters[df$promoter,]$V1
  df$promoter_start <- promoters[df$promoter,]$V2
  df$promoter_end <- promoters[df$promoter,]$V3
  df$pwm_name <- tf.info[df$pwm,]$TF_Name
  return(df)
}

setwd("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/probabilities/")

promoters <- read.table("../../data/gencode/promoter_coordinates_gencode_v19_hs37.bed")
tf.info <- na.omit(as.data.table(read.csv("../tf_info_scored.csv", stringsAsFactors = FALSE)), cols="score_threshold")

all_files <- list.files("pvalues/", full.names = TRUE)

for (i in 1:length(all_files)){
  
  print(i)
  
  # format of the file is x,n,p,promoter,pwm,pvalue
  dfPvalues <- read.csv(all_files[i], header = FALSE)
  colnames(dfPvalues) <- c("x","n","p","promoter","pwm","pvalue")
  
  pval_N <- nrow(tf.info)*nrow(promoters)
  
  to_fill <- pval_N - nrow(dfPvalues)
  
  pvalues <- c(dfPvalues$pvalue,rep(1,to_fill))
  
  fdr <- p.adjust(pvalues,"fdr")
  
  dfPvalues$fdr <- fdr[1:nrow(dfPvalues)]
  
  dfPvalues <- get_names(dfPvalues)
  
  dfPvalues$type <- strsplit(basename(all_files[i]),".txt")[[1]][1]
  
  write.table(dfPvalues,paste0("fdr/",basename(all_files[i])), sep = ",",
              quote = FALSE, row.names = FALSE, col.names = TRUE)
}
