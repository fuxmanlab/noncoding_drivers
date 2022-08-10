library(data.table)

# Rscript alterome_filter.R promoter_alterome tf_info_path out_dir

args = commandArgs(trailingOnly=TRUE)


get_bedOut <- function(bed,tf.info,sampleName){
  bed$wt_score <- pmax(bed$V7,bed$V8)
  bed$A <- pmax(bed$V9,bed$V10)
  bed$C <- pmax(bed$V11,bed$V12)
  bed$G <- pmax(bed$V13,bed$V14)
  bed$T <- pmax(bed$V15,bed$V16)
  
  score <- c()
  for (i in 1:nrow(bed)){
	#print(c(100*i/nrow(bed)))
    score <- c(score,bed[i,grepl(bed$V21[i], colnames(bed))])
  }
  bed$alt_score <- score
  bedOut <- bed[,c("V1","V2","V3","V5","V20","V21","V28","V32","wt_score","alt_score")]
  bedOut$sample <-  gsub(".bed","",sampleName)
  bedOut$tf_name <- tf.info[bed$V5,]$TF_Name
  bedOut$tf_pwm <- tf.info[bed$V5,]$Motif_ID
  bedOut$tf_threshold <- tf.info[bed$V5,]$score_threshold
  bedOut$fis <- abs(bedOut$wt_score-bedOut$alt_score)
  bedOut <- bedOut[bedOut$wt_score>bedOut$tf_threshold | bedOut$alt_score>bedOut$tf_threshold,]
  bedOut <- bedOut[order(bedOut$V1,bedOut$V2),]
  
  bedOut$loss_1 <- ifelse(bedOut$wt_score>bedOut$alt_score & bedOut$alt_score<bedOut$tf_threshold,1,0)
  bedOut$loss_2 <- ifelse(bedOut$fis>2 & bedOut$wt_score>bedOut$alt_score,1,0)
  bedOut$loss_3 <- ifelse(bedOut$fis>3 & bedOut$wt_score>bedOut$alt_score,1,0)
  
  bedOut$gain_1 <- ifelse(bedOut$alt_score>bedOut$wt_score & bedOut$wt_score<bedOut$tf_threshold,1,0)
  bedOut$gain_2 <- ifelse(bedOut$fis>2 & bedOut$alt_score>bedOut$wt_score,1,0)
  bedOut$gain_3 <- ifelse(bedOut$fis>3 & bedOut$alt_score>bedOut$wt_score,1,0)
  
  return(bedOut)
}


bedDir <- args[1]
tf.file <- args[2]
outDir <- args[3]
index <- as.numeric(args[4])

tf.info <- na.omit(as.data.table(read.csv(tf.file,
                              stringsAsFactors = FALSE)), cols="score_threshold")

# bedDir <- "promoter_alterome"
bedFiles <- list.files(bedDir)

bed <- read.csv(paste0(bedDir,"/",bedFiles[index]),
                sep = "\t", header = FALSE,
                stringsAsFactors = FALSE)
bed <- bed[bed$V23=="." & bed$V2==bed$V18,]
if(nrow(bed)>0) {
  bedOut <- get_bedOut(bed,tf.info,bedFiles[index])
  write.table(bedOut, paste0(outDir,"/",bedFiles[index]), quote = FALSE,
              col.names = FALSE, row.names = FALSE, sep = "\t")
}
