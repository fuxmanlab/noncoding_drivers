
args = commandArgs(trailingOnly=TRUE)

library(data.table)


pFile <- args[1]
nFile <- args[2]
threshold <- as.numeric(args[3])+2
outDir <- args[4]

prob <- read.csv(pFile, header = FALSE)

pwm <- strsplit(strsplit(pFile,".csv")[[1]],"_")[[1]][2]

outName <- paste0(outDir,"/",pwm,".txt")

n <- fread(nFile, header = FALSE)

prob <- prob[prob$V2 %in% n$V2,]
n <- n[n$V2 %in% prob$V2,]
n <- n[order(n$V2),]

out <- c(pwm,sum(n$V1*prob[,threshold]))
write(out,outName,sep = " ", ncolumns = 2)
