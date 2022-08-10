library(data.table)

args = commandArgs(trailingOnly=TRUE)

pFile <- args[1]
nFile <- args[2]
xFile <- args[3]
threshold <- as.numeric(args[4])+2
outDir <- args[5]

# pFile <- "pwm_812.csv"
# nFile <- "number_variants_per_promoter.txt"
# xFile <- "812.txt"
# threshold <- 4+2
# outDir <- args[4]


get_pvalue <- function(index,x,n,p){
  return(binom.test(x[index],n[index],p[index],"greater")$p.value)
}

prob <- read.csv(pFile, header = FALSE)

pwm <- strsplit(strsplit(pFile,".csv")[[1]],"_")[[1]][length(strsplit(strsplit(pFile,".csv")[[1]],"_")[[1]])]

outName <- paste0(outDir,"/",pwm,".out")

n <- fread(nFile, header = FALSE, sep = " ")



if (!basename(xFile) %in% c("700.txt","1162.txt")){
  
  x <- fread(xFile, header = FALSE)
  
  prob <- prob[prob$V2 %in% x$V2,]
  prob <- prob[order(prob$V2),]
  n <- n[n$V2 %in% x$V2,]
  n <- n[order(n$V2),]
  
  x <- x[order(x$V2),]
  
  df <- data.frame("x"=x$V1,"n"=n$V1,"p"=prob[,threshold],"promoter"=x$V2, "pwm"=pwm)
  
  df$pvalue <- sapply(seq(1:nrow(df)), get_pvalue, df$x, df$n, df$p)
  
  write.table(df,outName, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
}

