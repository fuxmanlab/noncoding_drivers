
args = commandArgs(trailingOnly=TRUE)


vcf_file <- args[1]
outname <- args[2]

vcf <- read.table(vcf_file, comment.char = "#")

vcf$id <- paste(vcf$V1,vcf$V2,vcf$V5)

vcf <- vcf[!duplicated(vcf$id),]

vcf <- vcf[ , -which(names(vcf) %in% c("id"))]

write.table(vcf,outname,sep = "\t", quote = FALSE,
            row.names = FALSE,col.names = FALSE)
