library(rtracklayer)
library(IRanges)

args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = args[2]
#output_file <- gsub("gencode.v27.bed", "promoter_coordinates_gencode.v27.bed", output_file)


# input_file = "../../data/gencode/gencode.v19.annotation.gtf"
# output_file = "../../data/gencode/cds_coordinates_gencode.v27.bed"

# Read GTF file
#gtf <- readGFF("gencode_protein_coding_basic.gtf")
command <- paste("grep -w 'CDS'",input_file, ">",'coding')
system(command)
gtf <- readGFF('coding')
command <- "rm coding"
system(command)

# Load chr lengths
# load("scr/chr_lengths.rda")
# chr_lengths <- data.frame(chr_lengths, stringsAsFactors = FALSE)

# Filter for transcript and gene
gtf <- gtf[gtf$type %in% c("CDS") & !(gtf$seqid %in% c('chrM')),!(colnames(gtf) %in% c("score","phase"))]

# Create gene id and gene name dictionary
gene_dict <- unique(gtf[,colnames(gtf) %in% c("gene_id","gene_name")])
gene_dict <- gene_dict[order(gene_dict$gene_id),]


# Calculate promoters per gene
coordinates <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), 
                        c("chr", "start", "end", "gene_name", "level","strand","gene_id"))

for (i in 1:length(gene_dict$gene_id)) {
  tmp <- gtf[gtf$gene_id %in% c(gene_dict$gene_id[i]),]
  chromosome <- as.character(unique(tmp$seqid))
  gene_id <- as.character(unique(tmp$gene_id))
  strand <- as.character(unique(tmp$strand))
  gene_name <- as.character(unique(tmp$gene_name))
  query <- reduce(IRanges(tmp$start,tmp$end))
  for (j in 1:length(query)){
    coordinates[nrow(coordinates)+1,] <- c(chromosome,
                                           start(query)[j],
                                           end(query)[j],
                                           gene_name,
                                           "1",
                                           strand,
                                           gene_id)
  }
  if (i%%100==0){
    print(i)
  }
}

write.table(coordinates, output_file,quote = F,sep = "\t",row.names = F, col.names = F)
