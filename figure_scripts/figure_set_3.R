library(data.table)
library(tidyr)
library(stringr)
library(dunn.test)
library(ggplot2)
library(stringr)
library(readxl)
library(survival)
library(survminer)
library(limma)
library(ggplot2)
library(matrixStats)


get_pval <- function(x){
  
  return(wilcox.test(x, mu = 0, alternative = "two.sided")$p.value)
  
}

cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

process_df <- function(df, ids){
  thr <- 0.05
  df <- df[df$A.logPadj_BF > -log10(thr) | df$B.logPadj_BF > -log10(thr) ,]
  tmp <- df[df$project!="literature_snv" & df$SNP %in% ids,]
  tmp$project <- "literature_snv"
  df <- rbind(df,tmp)
  return(df)
  
}

donors <- fread("~/Documents/Fuxman/noncoding_cancer/data/pcawg/donor.all_projects.tsv")

all_snvs <- fread("~/Documents/Fuxman/noncoding_cancer/data/SNVs_in_promoters/all_snvs_promoters_uniq.txt")
all_snvs <- all_snvs[!duplicated(paste(all_snvs$V1,all_snvs$V6)),]
all_snvs$V1 <- str_split(all_snvs$V1,"[.]",simplify = TRUE)[,1]
colnames(all_snvs) <- c("sample","chr","position","ref","alt","gene")
all_snvs$value <- 1


# all_snvs <- spread(all_snvs, V1, value, fill = 0)

pcawg <- read.csv("~/Documents/Fuxman/noncoding_cancer/data/pcawg/pcawg_sample.csv", stringsAsFactors = FALSE)
pcawg <-  dplyr::select(pcawg,c(aliquot_id,dcc_project_code,icgc_donor_id))
projects <- fread("~/Documents/Fuxman/noncoding_cancer/data/pcawg/pcawg_projects.csv", header = FALSE, stringsAsFactors = FALSE)
projects <-  dplyr::select(projects,c(V1,V3))


all_snvs <- merge(all_snvs,pcawg, by.x = "sample", by.y = "aliquot_id")
all_snvs <- merge(all_snvs,projects, by.x = "dcc_project_code", by.y = "V1")
colnames(all_snvs)[10] <- "cancer_type"

# load("~/Documents/Fuxman/noncoding_cancer/data/rnaseq/gene_quant/transcrip2gene.RData") 
# rna_files <- list.files("~/Documents/Fuxman/noncoding_cancer/data/rnaseq/gene_quant", full.names = TRUE)
# rna_files <- grep("sf",rna_files, value = TRUE)
# for (i in 1:length(rna_files)){
#   print(i)
#   rna <- fread(rna_files[i])
#   rna <- merge(rna,transcrip2gene, by.x = "Name", by.y = "ensembl_transcript_id_version")
#   rna <- aggregate(rna$TPM, by=list(Category=rna$hgnc_symbol), FUN=sum)
#   rna <- rna[rna$Category!="",]
#   colnames(rna) <- c("gene","sample")
#   if(i==1){
#     out <- rna
#   } else {
#     out <- merge(out,rna, by = "gene", all.x = TRUE)
#   }
#   colnames(out)[i+1] <- str_split(basename(rna_files[i]),"[.]")[[1]][1]
#   
# }
# 
# out[is.na(out)] <- 0
# 
# fwrite(out,"~/Documents/Fuxman/noncoding_cancer/data/rnaseq/gene_quant/gene_tpm.tsv",
#        row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


####################################
# RNA-seq analysis for TPM in mpra+, predicted driver, mutated promoter & no mutated promoter

tpm <- fread("~/Documents/Fuxman/noncoding_cancer/data/rnaseq/gene_quant/gene_tpm.tsv")

sig_snv <- fread("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/probabilities/filtered_fdr/mutation_total_set_sig_only.csv")
sig_snv[sig_snv$promoter_name=="1-Mar"]$promoter_name <- "MARCH1"


s <- fread("~/Documents/Fuxman/noncoding_cancer/data/for_ryan/snvs_mpra_total.csv")
s$neg <- ifelse(s$alt=="A","T",ifelse(s$alt=="T","A",ifelse(s$alt=="G","C","G")))
s$neg_ref <- ifelse(s$ref=="A","T",ifelse(s$ref=="T","A",ifelse(s$ref=="G","C","G")))
s$id <- paste(s$chr,s$start,s$ref,s$ref, sep = ":")
s$id_neg <- paste(s$chr,s$start,s$neg_ref,s$neg, sep = ":")

ids <- c(s[s$type=="literature_snv",]$id,s[s$type=="literature_snv",]$id_neg)


ht29 <- process_df(fread("~/Documents/Fuxman/noncoding_cancer/mpra/analysis_pilot_20200217/results/OL14_pilot_HT29_emVAR_20200217.out"),ids)
jurkat <- process_df(fread("~/Documents/Fuxman/noncoding_cancer/mpra/analysis_pilot_20200217/results/OL14_pilot_JURKAT_emVAR_20200217.out"),ids)
skmel28 <- process_df(fread("~/Documents/Fuxman/noncoding_cancer/mpra/analysis_pilot_20200217/results/OL14_pilot_SKMEL28_emVAR_20200217.out"),ids)
df <- rbind(ht29,jurkat,skmel28)
df <- df[df$project %in% c("cancer_snv_sig"),]
df$effect <- ifelse(df$LogSkew > 0, "overexpression","underexpression")
snv <- sig_snv[!duplicated(sig_snv$mut_id),]
snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
df <- merge(snv,df, by = "SNP")
thr <- 0.05
sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR) & df$project=="cancer_snv_sig",]
counts_gene <- as.data.frame(table(sig_df$promoter_name))
out <- data.frame()
thr_samples <- 2
for (gene in counts_gene$Var1){
  c_over <- paste0(unique(sig_df[sig_df$LogSkew > 0 & sig_df$promoter_name==gene,]$icgc_donor_id), collapse = ",")
  c_under <- paste0(unique(sig_df[sig_df$LogSkew < 0 & sig_df$promoter_name==gene,]$icgc_donor_id), collapse = ",")
  flag <- ifelse((nchar(c_over) > 1 & nchar(c_under) < 1)  & str_count(c_over,",")+1 >= thr_samples | (nchar(c_over) < 1 & nchar(c_under) > 1 & str_count(c_under,",")+1 >= thr_samples),1,0)
  c_out <- data.frame(gene=gene,over_ids=c_over,under_ids=c_under,flag=flag)
  out <- rbind(out,c_out)
}

mpra_positive <- out[out$flag==1,]


all_snvs <- data.frame(all_snvs)
tpm_cancer <- dplyr::select(tpm,c("gene",grep("_c",colnames(tpm),value = TRUE)))
rna_validated <- data.frame()
for( cancer_type in unique(all_snvs$cancer_type)){
  c_all_snvs <- all_snvs[all_snvs$cancer_type==cancer_type,]
  for( gene in as.character(mpra_positive$gene)){
    
    c_donors <- unique(c_all_snvs$icgc_donor_id)
    mpra_positive_donors <- unique(sig_df[sig_df$icgc_donor_id %in% c_donors & sig_df$promoter_name==gene,]$icgc_donor_id)
    predicted_driver <- unique(sig_snv[sig_snv$icgc_donor_id %in% c_donors & sig_snv$promoter_name==gene,]$icgc_donor_id)
    mutated_promoter <- unique(c_all_snvs[c_all_snvs$gene==gene,]$icgc_donor_id)
    no_mutated_promoter <- c_donors[!c_donors %in% mutated_promoter]
    
    if(length(mpra_positive_donors) > 0 & length(grep(paste(mpra_positive_donors,collapse="|"),colnames(tpm_cancer),value = TRUE))>1){
      tmp <- as.data.frame(dplyr::select(tpm,c("gene",grep(paste(mpra_positive_donors,collapse="|"),colnames(tpm_cancer),value = TRUE))))
      tmp <- tmp[tmp$gene==gene,]
      tmp <- gather(tmp,id,tpm,-gene)
      tmp <- cbind(group="mpra_positive_donors",cancer_type,tmp)
      rna_validated <- rbind(rna_validated,tmp)
      
      tmp <- as.data.frame(dplyr::select(tpm,c("gene",grep(paste(predicted_driver,collapse="|"),colnames(tpm_cancer),value = TRUE))))
      tmp <- tmp[tmp$gene==gene,]
      tmp <- gather(tmp,id,tpm,-gene)
      tmp <- cbind(group="predicted_driver",cancer_type,tmp)
      rna_validated <- rbind(rna_validated,tmp)
      
      tmp <- as.data.frame(dplyr::select(tpm,c("gene",grep(paste(mutated_promoter,collapse="|"),colnames(tpm_cancer),value = TRUE))))
      tmp <- tmp[tmp$gene==gene,]
      tmp <- gather(tmp,id,tpm,-gene)
      tmp <- cbind(group="mutated_promoter",cancer_type,tmp)
      rna_validated <- rbind(rna_validated,tmp)
      
      tmp <- as.data.frame(dplyr::select(tpm,c("gene",grep(paste(no_mutated_promoter,collapse="|"),colnames(tpm_cancer),value = TRUE))))
      tmp <- tmp[tmp$gene==gene,]
      tmp <- gather(tmp,id,tpm,-gene)
      tmp <- cbind(group="no_mutated_promoter",cancer_type,tmp)
      rna_validated <- rbind(rna_validated,tmp)
      
    }
  }
}

pdf("~/Documents/Fuxman/noncoding_cancer/data/rnaseq/expression_plots_by_group/expression_plots_by_group.pdf")
for( gene in unique(rna_validated$gene)){
  
  tmp <- rna_validated[rna_validated$gene==gene,]
  k <- kruskal.test(tmp$tpm ~ tmp$group)
  dunn.test(tmp$tpm, tmp$group)
  p <- ggplot(tmp,aes(x=tmp$group,y=log10(tmp$tpm+1))) + geom_boxplot() + geom_point() + ggtitle(paste(gene, tmp$cancer_type[1], sig_df[sig_df$promoter_name==gene]$effect[1],round(k$p.value,4)))
  plot(p)
}
dev.off()




########################
# TF Gain/Loss and +/-/no validation


tpm <- fread("~/Documents/Fuxman/noncoding_cancer/data/rnaseq/gene_quant/gene_tpm.tsv")

sig_snv <- fread("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/probabilities/filtered_fdr/mutation_total_set_sig_only.csv")
sig_snv[sig_snv$promoter_name=="1-Mar"]$promoter_name <- "MARCH1"


s <- fread("~/Documents/Fuxman/noncoding_cancer/data/for_ryan/snvs_mpra_total.csv")
s$neg <- ifelse(s$alt=="A","T",ifelse(s$alt=="T","A",ifelse(s$alt=="G","C","G")))
s$neg_ref <- ifelse(s$ref=="A","T",ifelse(s$ref=="T","A",ifelse(s$ref=="G","C","G")))
s$id <- paste(s$chr,s$start,s$ref,s$ref, sep = ":")
s$id_neg <- paste(s$chr,s$start,s$neg_ref,s$neg, sep = ":")

ids <- c(s[s$type=="literature_snv",]$id,s[s$type=="literature_snv",]$id_neg)


ht29 <- process_df(fread("~/Documents/Fuxman/noncoding_cancer/mpra/analysis_pilot_20200217/results/OL14_pilot_HT29_emVAR_20200217.out"),ids)
jurkat <- process_df(fread("~/Documents/Fuxman/noncoding_cancer/mpra/analysis_pilot_20200217/results/OL14_pilot_JURKAT_emVAR_20200217.out"),ids)
skmel28 <- process_df(fread("~/Documents/Fuxman/noncoding_cancer/mpra/analysis_pilot_20200217/results/OL14_pilot_SKMEL28_emVAR_20200217.out"),ids)
ht29$cell_line <- "ht29"
jurkat$cell_line <- "jurkat"
skmel28$cell_line <- "skmel28"
df <- rbind(ht29,jurkat,skmel28)
df <- df[df$project %in% c("cancer_snv_sig"),]
df$effect <- ifelse(df$LogSkew > 0, "overexpression","underexpression")

sig_snv$id_tf <- paste(sig_snv$mut_id,sig_snv$tf_name)
sig_snv$tf_effect <- ifelse(sig_snv$gain_1==1 | sig_snv$gain_2==1 | sig_snv$gain_3==1,"gain","loss")
snv <- sig_snv[!duplicated(sig_snv$id_tf),]
snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
df <- merge(snv,df, by = "SNP", allow.cartesian=TRUE)
thr <- 0.05
sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR) & df$project=="cancer_snv_sig",]
df <- as.data.frame(df)
sig_df <- as.data.frame(sig_df)
out <- data.frame()
header <- c("cell_line","tf_name","gain_over","gain_under","gain_noDiff","loss_over","loss_under","loss_noDiff")
for( cell_line in unique(df$cell_line)){
  
  for ( tf in unique(sig_df$tf_name)){
    
    c_sig_df <- sig_df[sig_df$cell_line==cell_line & sig_df$tf_name==tf ,]
    c_df <- df[df$cell_line==cell_line & df$tf_name==tf ,]
    c_row <- data.frame(cell_line,tf,
               nrow(c_sig_df[c_sig_df$tf_effect=="gain" & c_sig_df$effect=="overexpression",]),
               nrow(c_sig_df[c_sig_df$tf_effect=="gain" & c_sig_df$effect=="underexpression",]),
               nrow(c_df[c_df$tf_effect=="gain",])-nrow(c_sig_df[c_sig_df$tf_effect=="gain" ,]),
               nrow(c_sig_df[c_sig_df$tf_effect=="loss" & c_sig_df$effect=="overexpression",]),
               nrow(c_sig_df[c_sig_df$tf_effect=="loss" & c_sig_df$effect=="underexpression",]),
               nrow(c_df[c_df$tf_effect=="loss",])-nrow(c_sig_df[c_sig_df$tf_effect=="loss" ,])
              )
    colnames(c_row) <- header
    out <- rbind(out,c_row)
  }
  
}

out$tf_name <- as.character(out$tf_name)

######  ACTIVATOR
out$ratio <- (out$gain_over + out$loss_under + 1)/(out$gain_under + out$loss_over + 1)

fwrite(out,"~/Documents/Fuxman/noncoding_cancer/data/rnaseq/expression_plots_by_group/tf_effect.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

f_out <- out[log10(out$ratio) > 0.5,]
tfs_sig <- as.data.frame(table(f_out$tf_name))
tfs_sig <- as.character(tfs_sig[tfs_sig$Freq > 2,]$Var1)

f_out <- f_out[f_out$tf_name %in% tfs_sig,]
f_out$tf_type <- ifelse(f_out$gain_over+f_out$loss_under > f_out$gain_under+f_out$loss_over,"activator","repressor")
f_out <- f_out[!duplicated(f_out$tf_name),]

normalized_tmp <- data.frame()
tf_sig_snv <- sig_snv[sig_snv$tf_name %in% tfs_sig,]
for( cancer_type in unique(all_snvs$cancer_type)){
  c_all_snvs <- all_snvs[all_snvs$cancer_type==cancer_type,]
  for( gene in as.character(unique(tf_sig_snv$promoter_name))){
    
    c_donors <- unique(c_all_snvs$icgc_donor_id)
    predicted_driver <- unique(sig_snv[sig_snv$icgc_donor_id %in% c_donors & sig_snv$promoter_name==gene,]$icgc_donor_id)
    mutated_promoter <- unique(c_all_snvs[c_all_snvs$gene==gene,]$icgc_donor_id)
    no_mutated_promoter <- c_donors[!c_donors %in% mutated_promoter]
    
    if(gene %in% tpm$gene & length(predicted_driver) > 1 & length(grep(paste(predicted_driver,collapse="|"),colnames(tpm_cancer),value = TRUE))>1){
      
      tmp_pd <- as.data.frame(dplyr::select(tpm,c("gene",grep(paste(predicted_driver,collapse="|"),colnames(tpm_cancer),value = TRUE))))
      tmp_pd <- tmp_pd[tmp_pd$gene==gene,]
      tmp_pd <- gather(tmp_pd,id,tpm,-gene)
      tmp_pd <- cbind(group="predicted_driver",cancer_type,tmp_pd)
      
      
      tmp <- as.data.frame(dplyr::select(tpm,c("gene",grep(paste(no_mutated_promoter,collapse="|"),colnames(tpm_cancer),value = TRUE))))
      tmp <- tmp[tmp$gene==gene,]
      tmp <- gather(tmp,id,tpm,-gene)
      tmp <- cbind(group="no_mutated_promoter",cancer_type,tmp)
      # tmp_pd$mean_normalized <- tmp_pd$tpm/mean(tmp$tpm)
      tmp_pd$median_normalized <- log((tmp_pd$tpm+1)/median(tmp$tpm+1), base = 2)
      
      normalized_tmp <- rbind(normalized_tmp,tmp_pd)
      
    }
  }
}

sig_snv <- as.data.frame(sig_snv)
sig_snv$tf_effect <- ifelse(sig_snv$gain_1==1 | sig_snv$gain_2==1 | sig_snv$gain_3==1,"gain","loss")

out_binomial <- data.frame()
header <- c("tf_name","gained_overexpressed","gained_total","over_pval","loss_underexpressed",
            "loss_total","under_pval")
out_boxplot_gain <- data.frame()
out_boxplot_loss <- data.frame()
header_boxplot <- c("tf","tpm_normalized","tfbs_effect")
for( tf in tfs_sig) {
  
  tf_gain_genes <- unique(sig_snv[sig_snv$tf_name==tf & sig_snv$tf_effect=="gain",]$promoter_name)
  tf_loss_genes <- unique(sig_snv[sig_snv$tf_name==tf & sig_snv$tf_effect=="loss",]$promoter_name)
  tf_gain_tpm <- normalized_tmp[normalized_tmp$gene %in% tf_gain_genes,]
  tf_loss_tpm <- normalized_tmp[normalized_tmp$gene %in% tf_loss_genes,]
  pval_gain <- ifelse(nrow(tf_gain_tpm)>0,binom.test(sum(tf_gain_tpm$median_normalized>0),nrow(tf_gain_tpm),p=0.5,"greater")$p.value,NA)
  pval_loss <- ifelse(nrow(tf_loss_tpm)>0,binom.test(sum(tf_loss_tpm$median_normalized<0),nrow(tf_loss_tpm),p=0.5,"greater")$p.value,NA)
  
  c_row <- data.frame(tf,sum(tf_gain_tpm$median_normalized>0),nrow(tf_gain_tpm),pval_gain,
             sum(tf_loss_tpm$median_normalized<0),nrow(tf_loss_tpm),pval_loss
             )
  colnames(c_row) <- header
  
  if (length(tf_gain_tpm$median_normalized)>0){
    c_row_boxplot_gain <- data.frame(tf=tf_gain_tpm$median_normalized)
    colnames(c_row_boxplot_gain)[1] <- tf
    out_boxplot_gain <- cbind.fill(out_boxplot_gain,c_row_boxplot_gain)
  }
  
  if (length(tf_loss_tpm$median_normalized)>0){
    c_row_boxplot_loss <- data.frame(tf=tf_loss_tpm$median_normalized)
    colnames(c_row_boxplot_loss)[1] <- tf
    out_boxplot_loss <- cbind.fill(out_boxplot_loss,c_row_boxplot_loss)
  }
  
  out_binomial <- rbind(out_binomial,c_row)
  
}

out_binomial$gain_ratio <- out_binomial$gained_overexpressed/out_binomial$gained_total
out_binomial$loss_ratio <- out_binomial$loss_underexpressed/out_binomial$loss_total

fwrite(out_binomial,"~/Documents/Fuxman/noncoding_cancer/data/rnaseq/expression_plots_by_group/tf_effect_pval.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

apply(out_boxplot_gain, 2, get_pval)
apply(out_boxplot_loss, 2, get_pval)


fwrite(out_boxplot_gain,"~/Documents/Fuxman/noncoding_cancer/data/rnaseq/expression_plots_by_group/tf_effect_gene_tpm_gain.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(out_boxplot_loss,"~/Documents/Fuxman/noncoding_cancer/data/rnaseq/expression_plots_by_group/tf_effect_gene_tpm_loss.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
out_boxplot_gain <- as.data.frame(out_boxplot_gain)
out_boxplot_loss <- as.data.frame(out_boxplot_loss)
out_boxplot_loss <- out_boxplot_loss[,colnames(out_boxplot_loss) %in% colnames(out_boxplot_gain)]
out_boxplot_total <- data.frame()
for ( i in 1:ncol(out_boxplot_loss)){
  out_boxplot_total <- cbind.fill(out_boxplot_total,out_boxplot_gain[,i],out_boxplot_loss[,i])
}
colnames(out_boxplot_total) <- rep(colnames(out_boxplot_gain),each=2)

fwrite(out_boxplot_total,"~/Documents/Fuxman/noncoding_cancer/data/rnaseq/expression_plots_by_group/tf_effect_gene_tpm_total.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

############################
# fraction total snvs by tf family

sig_snv <- fread("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/predicted_drivers.tsv")
sig_snv[sig_snv$promoter_name=="1-Mar"]$promoter_name <- "MARCH1"

tf_info <- fread("~/Documents/Fuxman/noncoding_cancer/motif_benchmark/tf_info_processed.csv")
tf_info <- dplyr::select(tf_info,c("TF_Name","Family_Name"))
tf_info <- tf_info[!duplicated(tf_info$TF_Name),]

sig_snv <- sig_snv[!duplicated(paste(sig_snv$mut_id,sig_snv$tf_name)),]
sig_snv <- merge(sig_snv, tf_info, by.x = "tf_name", by.y="TF_Name")

tf_families <- c("bHLH","bZIP","C2H2 ZF","Ets","Forkhead","Homeodomain","Nuclear receptor","Rel","Sox","Other")
sig_snv$Family_Name <- ifelse(sig_snv$Family_Name %in% tf_families,sig_snv$Family_Name,"Other")
out <- matrix(NA,ncol = length(tf_families),nrow = length(unique(sig_snv$cancer_type)))
for(i in 1:length(sort(unique(sig_snv$cancer_type)))){
  c_sig_snv <- sig_snv[sig_snv$cancer_type==sort(unique(sig_snv$cancer_type))[i],]
  c_row <- c()
  for(j in 1:length(tf_families)){
    print(c(i,j))
    out[i,j] <- sum(str_count(c_sig_snv$Family_Name, tf_families[j]))/nrow(c_sig_snv)
  }
  
}

out <- as.data.frame(out)
colnames(out) <- tf_families
rownames(out) <- sort(unique(sig_snv$cancer_type))
test <- pheatmap::pheatmap(out, cluster_cols = TRUE, cluster_rows = TRUE,fontsize_row = 8, 
                           treeheight_row = 0, treeheight_col = 0)

write.table(out,"~/Documents/Fuxman/noncoding_cancer/sig_tfs/tf_family_freq_sig_snv.tsv",
            col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

###############
# Fitness and essential

tf_info <- fread("~/Documents/Fuxman/noncoding_cancer/motif_benchmark/tf_info_processed.csv")
cgc <- fread("~/Documents/Fuxman/noncoding_cancer/data/cosmic/Census_allThu Aug  2 18_21_58 2018.csv")
cgc <- cgc[cgc$`Gene Symbol` %in% tf_info$TF_Name,]
promoters <- fread("~/Documents/Fuxman/noncoding_cancer/data/gencode/no_cds_promoter_v19_hs37.bed")
sig_promoters <- unique(promoters[promoters$V8 %in% unique(sig$promoter_id),]$V4)
n_promoters <- unique(promoters$V4[!promoters$V4 %in% sig_promoters])
sig_snv <- fread("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/probabilities/filtered_fdr/mutation_total_set_sig_only.csv")
test <- unique(sig_snv[sig_snv$tf_name %in% tf_info$TF_Name,]$mut_id)
sig_tfs <- unique(sig_snv$tf_name)[unique(sig_snv$tf_name) %in% tf_info$TF_Name]
no_sig_tfs <- unique(tf_info$TF_Name)
no_sig_tfs <- no_sig_tfs[!no_sig_tfs %in% sig_tfs]

fitness <- fread("~/Downloads/fitness.csv")
fitness <- fitness[fitness$`#NumTKOHits` >=3,]

s_f <- sum(sig_tfs %in% fitness$`#GENE`)
t_f <- sum(no_sig_tfs %in% fitness$`#GENE`)
# c_f <- sum(cgc$`Gene Symbol` %in% fitness$`#GENE`)
# t_c <-  unique(promoters$V4[!promoters$V4 %in% cgc$`Gene Symbol`])

res <- prop.test(x = c(s_f, t_f), n = c(length(sig_tfs), length(no_sig_tfs)))
# res <- prop.test(x = c(c_f, t_f), n = c(length(cgc$`Gene Symbol`), length(t_c)))
# sum(unique(promoters$V4) %in% fitness$`#GENE`)/length(unique(promoters$V4))

fitness <- fread("~/Downloads/Achilles_common_essentials.csv")

s_f <- sum(sig_tfs %in% fitness$V1)
t_f <- sum(no_sig_tfs %in% fitness$V1)
# c_f <- sum(cgc$`Gene Symbol` %in% fitness$V1)
# t_c <-  unique(promoters$V4[!promoters$V4 %in% cgc$`Gene Symbol`])

res <- prop.test(x = c(s_f, t_f), n = c(length(sig_tfs), length(no_sig_tfs)))
# res <- prop.test(x = c(c_f, t_f), n = c(length(cgc$`Gene Symbol`), length(t_c)))
# sum(unique(promoters$V4) %in% fitness$V1)/length(unique(promoters$V4))

N <- 20229
p <- 0.1
sqrt(p*(1-p)/N)

############
# Prognosis

prognostics <- fread("~/Documents/Fuxman/noncoding_cancer/data/human_protein_atlas/counts.csv")

sig_fav <- prop.test(c(sum(sig_tfs %in% prognostics[prognostics$fcount>0,]$gene),
                       sum(no_sig_tfs %in% prognostics[prognostics$fcount>0,]$gene)),
                     c(length(sig_tfs),
                       length(no_sig_tfs)))
sig_unf <- prop.test(c(sum(sig_tfs %in% prognostics[prognostics$ucount>0,]$gene),
                       sum(no_sig_tfs %in% prognostics[prognostics$ucount>0,]$gene)),
                     c(length(sig_tfs),
                       length(no_sig_tfs)))
sig_eit <- prop.test(c(sum(sig_tfs %in% prognostics[prognostics$total>0,]$gene),
                       sum(no_sig_tfs %in% prognostics[prognostics$total>0,]$gene)),
                     c(length(sig_tfs),
                       length(no_sig_tfs)))


# cgc_fav <- prop.test(c(sum(cgc$`Gene Symbol` %in% prognostics[prognostics$fcount>0,]$gene),
#                        sum(t_c %in% prognostics[prognostics$fcount>0,]$gene)),
#                      c(length(cgc$`Gene Symbol`),
#                        length(t_c)))
# cgc_unf <- prop.test(c(sum(cgc$`Gene Symbol` %in% prognostics[prognostics$ucount>0,]$gene),
#                        sum(t_c %in% prognostics[prognostics$ucount>0,]$gene)),
#                      c(length(sig_promoters),
#                        length(t_c)))
# cgc_eit <- prop.test(c(sum(cgc$`Gene Symbol` %in% prognostics[prognostics$total>0,]$gene),
#                        sum(t_c %in% prognostics[prognostics$total>0,]$gene)),
#                      c(length(cgc$`Gene Symbol`),
#                        length(t_c)))
# 
# sum(unique(promoters$V4) %in% prognostics[prognostics$fcount>0,]$gene)/length(unique(promoters$V4))
# sum(unique(promoters$V4) %in% prognostics[prognostics$ucount>0,]$gene)/length(unique(promoters$V4))
# sum(unique(promoters$V4) %in% prognostics[prognostics$total>0,]$gene)/length(unique(promoters$V4))


##############
# Monti data all

monti <- read_excel("~/Downloads/pancancer_cis_genes_by_freq.xlsx")
sig_snv <- fread("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/probabilities/filtered_fdr/mutation_total_set_sig_only.csv")
sig_snv[sig_snv$promoter_name=="1-Mar"]$promoter_name <- "MARCH1"

promoters <- fread("~/Documents/Fuxman/noncoding_cancer/data/gencode/no_cds_promoter_v19_hs37.bed")
sig_promoters <- unique(promoters[promoters$V8 %in% unique(sig_snv$promoter_id),]$V4)
n_promoters <- unique(promoters$V4[!promoters$V4 %in% sig_promoters])

out_amp <- matrix(NA, ncol = 8, nrow = 0)
out_del <- matrix(NA, ncol = 8, nrow = 0)
out_all <- matrix(NA, ncol = 8, nrow = 0)
for ( i  in c(1:6)){
  
  if(i<6){
    c_monti_amp <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types == i & monti$SCNA.Type=="Amplification",]$Gene)
    c_monti_del <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types == i & monti$SCNA.Type=="Deletion",]$Gene)
    
  } else{
    c_monti_amp <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types >= i & monti$SCNA.Type=="Amplification",]$Gene)
    c_monti_del <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types >= i & monti$SCNA.Type=="Deletion",]$Gene)
    
  }
 
  s_amp <- sum(sig_promoters %in% c_monti_amp)
  t_amp <- sum(n_promoters %in% c_monti_amp)
  s_del <- sum(sig_promoters %in% c_monti_del)
  t_del <- sum(n_promoters %in% c_monti_del)
  s_all <- sum(sig_promoters %in% unique(c(c_monti_amp,c_monti_del)))
  t_all <- sum(n_promoters %in% unique(c(c_monti_amp,c_monti_del)))
  
  amp_or <- (s_amp/(length(sig_promoters)-s_amp))/(t_amp/(length(n_promoters)-t_amp))
  del_or <- (s_del/(length(sig_promoters)-s_del))/(t_del/(length(n_promoters)-t_del))
  total_or <- (s_all/(length(sig_promoters)-s_all))/(t_all/(length(n_promoters)-t_all))
  
  res_amp <- prop.test(x = c(s_amp, t_amp), n = c(length(sig_promoters), length(n_promoters)))
  res_del <- prop.test(x = c(s_del, t_del), n = c(length(sig_promoters), length(n_promoters)))
  res_total <- prop.test(x = c(s_all, t_all), n = c(length(sig_promoters), length(n_promoters)))

  c_row_amp <- unname(c(i,length(c_monti_amp),res_amp$estimate[1],res_amp$estimate[2],res_amp$p.value,
             sqrt(res_amp$estimate[1]*(1-res_amp$estimate[1])/length(sig_promoters)),
            sqrt(res_amp$estimate[1]*(1-res_amp$estimate[1])/length(sig_promoters)),amp_or))
  
  c_row_del <- unname(c(i,length(c_monti_del),res_del$estimate[1],res_del$estimate[2],res_del$p.value,
                    sqrt(res_del$estimate[1]*(1-res_del$estimate[1])/length(sig_promoters)),
                    sqrt(res_del$estimate[1]*(1-res_del$estimate[1])/length(sig_promoters)),del_or))
  
  c_row_total <- unname(c(i,length(unique(c(c_monti_amp,c_monti_del))),res_total$estimate[1],res_total$estimate[2],res_total$p.value,
                        sqrt(res_total$estimate[1]*(1-res_total$estimate[1])/length(sig_promoters)),
                        sqrt(res_total$estimate[1]*(1-res_total$estimate[1])/length(sig_promoters)),total_or))
  
  out_amp <- rbind(out_amp,c_row_amp)
  out_del <- rbind(out_del,c_row_del)
  out_all <- rbind(out_all,c_row_total)
  
}

out_amp <- as.data.frame(out_amp)
colnames(out_amp) <- c("N_cell_lines","N_genes","p1","p2","pval","e1","e2","OR")

out_del <- as.data.frame(out_del)
colnames(out_del) <- c("N_cell_lines","N_genes","p1","p2","pval","e1","e2","OR")

out_all <- as.data.frame(out_all)
colnames(out_all) <- c("N_cell_lines","N_genes","p1","p2","pval","e1","e2","OR")

fwrite(out_amp,"~/Documents/Fuxman/noncoding_cancer/monti/amp_predicted_drivers.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(out_del,"~/Documents/Fuxman/noncoding_cancer/monti/del_predicted_drivers.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(out_all,"~/Documents/Fuxman/noncoding_cancer/monti/all_predicted_drivers.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

##########
# monti CGC

monti <- read_excel("~/Downloads/pancancer_cis_genes_by_freq.xlsx")
sig_snv <- fread("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/probabilities/filtered_fdr/mutation_total_set_sig_only.csv")
sig_snv[sig_snv$promoter_name=="1-Mar"]$promoter_name <- "MARCH1"

cgc <- fread("~/Documents/Fuxman/noncoding_cancer/data/cosmic/Census_allThu Aug  2 18_21_58 2018.csv")
cgc <- unique(cgc$`Gene Symbol`)
promoters <- fread("~/Documents/Fuxman/noncoding_cancer/data/gencode/no_cds_promoter_v19_hs37.bed")
no_cgc <- unique(promoters[!promoters$V4 %in% cgc_genes,]$V4)


out_amp <- matrix(NA, ncol = 8, nrow = 0)
out_del <- matrix(NA, ncol = 8, nrow = 0)
out_all <- matrix(NA, ncol = 8, nrow = 0)
for ( i  in c(1:6)){
  
  if(i<6){
    c_monti_amp <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types == i & monti$SCNA.Type=="Amplification",]$Gene)
    c_monti_del <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types == i & monti$SCNA.Type=="Deletion",]$Gene)
    
  } else{
    c_monti_amp <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types >= i & monti$SCNA.Type=="Amplification",]$Gene)
    c_monti_del <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types >= i & monti$SCNA.Type=="Deletion",]$Gene)
    
  }
  
  s_amp <- sum(cgc %in% c_monti_amp)
  t_amp <- sum(no_cgc %in% c_monti_amp)
  s_del <- sum(cgc %in% c_monti_del)
  t_del <- sum(no_cgc %in% c_monti_del)
  s_all <- sum(cgc %in% unique(c(c_monti_amp,c_monti_del)))
  t_all <- sum(no_cgc %in% unique(c(c_monti_amp,c_monti_del)))
  
  amp_or <- (s_amp/(length(cgc)-s_amp))/(t_amp/(length(no_cgc)-t_amp))
  del_or <- (s_del/(length(cgc)-s_del))/(t_del/(length(no_cgc)-t_del))
  total_or <- (s_all/(length(cgc)-s_all))/(t_all/(length(no_cgc)-t_all))
  
  res_amp <- prop.test(x = c(s_amp, t_amp), n = c(length(cgc), length(no_cgc)))
  res_del <- prop.test(x = c(s_del, t_del), n = c(length(cgc), length(no_cgc)))
  res_total <- prop.test(x = c(s_all, t_all), n = c(length(cgc), length(no_cgc)))
  
  c_row_amp <- unname(c(i,length(c_monti_amp),res_amp$estimate[1],res_amp$estimate[2],res_amp$p.value,
                        sqrt(res_amp$estimate[1]*(1-res_amp$estimate[1])/length(cgc)),
                        sqrt(res_amp$estimate[1]*(1-res_amp$estimate[1])/length(cgc)),amp_or))
  
  c_row_del <- unname(c(i,length(c_monti_del),res_del$estimate[1],res_del$estimate[2],res_del$p.value,
                        sqrt(res_del$estimate[1]*(1-res_del$estimate[1])/length(cgc)),
                        sqrt(res_del$estimate[1]*(1-res_del$estimate[1])/length(cgc)),del_or))
  
  c_row_total <- unname(c(i,length(unique(c(c_monti_amp,c_monti_del))),res_total$estimate[1],res_total$estimate[2],res_total$p.value,
                          sqrt(res_total$estimate[1]*(1-res_total$estimate[1])/length(cgc)),
                          sqrt(res_total$estimate[1]*(1-res_total$estimate[1])/length(cgc)),total_or))
  
  out_amp <- rbind(out_amp,c_row_amp)
  out_del <- rbind(out_del,c_row_del)
  out_all <- rbind(out_all,c_row_total)
  
}

out_amp <- as.data.frame(out_amp)
colnames(out_amp) <- c("N_cell_lines","N_genes","p1","p2","pval","e1","e2","OR")

out_del <- as.data.frame(out_del)
colnames(out_del) <- c("N_cell_lines","N_genes","p1","p2","pval","e1","e2","OR")

out_all <- as.data.frame(out_all)
colnames(out_all) <- c("N_cell_lines","N_genes","p1","p2","pval","e1","e2","OR")

fwrite(out_amp,"~/Documents/Fuxman/noncoding_cancer/monti/amp_cgc.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(out_del,"~/Documents/Fuxman/noncoding_cancer/monti/del_cgc.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(out_all,"~/Documents/Fuxman/noncoding_cancer/monti/all_cgc.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

########


df  <- rbind(jurkat,ht29,skmel28)
df <- df[df$project == "cancer_snv_sig",]
snv <- sig_snv[!duplicated(sig_snv$mut_id),]
snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
df <- merge(snv,df, by = "SNP")
thr <- 0.05
sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR),]
sig_df$effect <- ifelse(sig_df$LogSkew>0,"up","down")
total <- aggregate(sig_df$effect, list(sig_df$SNP), paste, collapse=",")
total$cell_count <- str_count(total$x,",")+1
total$up <- ifelse(str_count(total$x,"up")>0,1,0)
total$down <- ifelse(str_count(total$x,"down")>0,1,0)
total$direction <- total$up + total$down


promoters <- unique(fread("~/Documents/Fuxman/noncoding_cancer/data/gencode/no_cds_promoter_v19_hs37.bed")$V4)

sig_up <- unique(sig_df[sig_df$SNP %in% total[total$direction==1 & total$cell_count>1 & total$up>0,]$Group.1,]$promoter_name)
sig_down <- unique(sig_df[sig_df$SNP %in% total[total$direction==1 & total$cell_count>1 & total$down>0,]$Group.1,]$promoter_name)

# promoters <- unique(df[!df$promoter_name %in% c(sig_up,sig_down),]$promoter_name)


out_amp <- matrix(NA, ncol = 10, nrow = 0)
out_del <- matrix(NA, ncol = 10, nrow = 0)
for ( i  in c(1:6)){
  
  if(i<6){
    c_monti_amp <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types == i & monti$SCNA.Type=="Amplification",]$Gene)
    c_monti_del <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types == i & monti$SCNA.Type=="Deletion",]$Gene)
    
  } else{
    c_monti_amp <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types > i & monti$SCNA.Type=="Amplification",]$Gene)
    c_monti_del <- unique(monti[monti$Number.of.Occurrence.as.Cis.Gene..Out.of.19.Cancer.Types == i & monti$SCNA.Type=="Deletion",]$Gene)
    
  }
  
  s_amp <- sum(sig_up %in% c_monti_amp)
  t_amp <- sum(n_promoters %in% c_monti_amp)
  s_del <- sum(sig_down %in% c_monti_del)
  t_del <- sum(n_promoters %in% c_monti_del)
  
  amp_or <- (s_amp/(length(sig_up)-s_amp))/(t_amp/(length(n_promoters)-t_amp))
  del_or <- (s_del/(length(sig_down)-s_del))/(t_del/(length(n_promoters)-t_del))
  
  res_amp <- prop.test(x = c(s_amp, t_amp), n = c(length(sig_up), length(n_promoters)))
  res_del <- prop.test(x = c(s_del, t_del), n = c(length(sig_down), length(n_promoters)))
  
  c_row_amp <- unname(c(i,length(c_monti_amp),res_amp$estimate[1],res_amp$estimate[2],res_amp$p.value,
                        sqrt(res_amp$estimate[1]*(1-res_amp$estimate[1])/length(sig_up)),
                        sqrt(res_amp$estimate[1]*(1-res_amp$estimate[1])/length(sig_up)),
                        s_amp,length(sig_up),amp_or))
  
  c_row_del <- unname(c(i,length(c_monti_del),res_del$estimate[1],res_del$estimate[2],res_del$p.value,
                        sqrt(res_del$estimate[1]*(1-res_del$estimate[1])/length(sig_down)),
                        sqrt(res_del$estimate[1]*(1-res_del$estimate[1])/length(sig_down)),
                        s_del,length(sig_down),del_or))
  
  out_amp <- rbind(out_amp,c_row_amp)
  out_del <- rbind(out_del,c_row_del)
  
}

out_amp <- as.data.frame(out_amp)
colnames(out_amp) <- c("N_cell_lines","N_genes","p1","p2","pval","e1","e2","sig_N","total_N","OR")

out_del <- as.data.frame(out_del)
colnames(out_del) <- c("N_cell_lines","N_genes","p1","p2","pval","e1","e2","sig_N","total_N","OR")


#################
# Survival analysis

df <- rbind(ht29,jurkat,skmel28)
df <- df[df$project %in% c("cancer_snv_sig"),]
df$effect <- ifelse(df$LogSkew > 0, "overexpression","underexpression")
snv <- sig_snv[!duplicated(sig_snv$mut_id),]
snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
df <- merge(snv,df, by = "SNP")
thr <- 0.05
sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR) & df$project=="cancer_snv_sig",]
counts_gene <- as.data.frame(table(sig_df$promoter_name))
out <- data.frame()
thr_fraction <- 1
for (gene in counts_gene$Var1){
  c_over <- paste0(unique(sig_df[sig_df$LogSkew > 0 & sig_df$promoter_name==gene,]$icgc_donor_id), collapse = ",")
  c_under <- paste0(unique(sig_df[sig_df$LogSkew < 0 & sig_df$promoter_name==gene,]$icgc_donor_id), collapse = ",")
  c_over <- ifelse(c_over=="",0,str_count(c_over,",")+1)
  c_under <- ifelse(c_under=="",0,str_count(c_under,",")+1)
  if((c_over>1 | c_under>1)){
    flag <- ifelse( c_over/(c_over+c_under) >= thr_fraction,1, ifelse(c_under/(c_over+c_under) >= thr_fraction,1,0))
    c_out <- data.frame(gene=gene,over_ids=c_over,under_ids=c_under,flag=flag)
    out <- rbind(out,c_out)
    
  }
}

mpra_positive <- out[out$flag==1,]
out <- data.frame()
all_snvs <- as.data.frame(all_snvs)
for( cancer_type in unique(all_snvs$cancer_type)){
  c_all_snvs <- all_snvs[all_snvs$cancer_type==cancer_type,]
  for( gene in as.character(mpra_positive$gene)){
  # for ( gene in unique(sig_snv$promoter_name)){
    
    c_donors <- unique(c_all_snvs$icgc_donor_id)
    predicted_driver <- unique(sig_snv[sig_snv$icgc_donor_id %in% c_donors & sig_snv$promoter_name==gene,]$icgc_donor_id)
    mutated_promoter <- unique(c_all_snvs[c_all_snvs$gene==gene,]$icgc_donor_id)
    no_mutated_promoter <- c_donors[!c_donors %in% mutated_promoter]
    
    c_pd <- donors[donors$icgc_donor_id %in% predicted_driver & !is.na(donors$donor_survival_time),]
    c_np <- donors[donors$icgc_donor_id %in% no_mutated_promoter & !is.na(donors$donor_survival_time),]
    
    if(nrow(c_pd) > 1 & nrow(c_np)>1 ){
      
      os <- rbind(c_pd,c_np)
      os$type <- c(rep("predicted_driver",nrow(c_pd)),rep("no_mutated_promoter",nrow(c_np)))
      os$donor_vital_status <- ifelse(os$donor_vital_status=="deceased",2,1)
      res.cox <- coxph(Surv(donor_survival_time, donor_vital_status) ~ type, data = os)
      hazard_ratio <- exp(res.cox$coefficients)
      pw <- summary(res.cox)$waldtest[3]
      pl <- summary(res.cox)$logtest[3]
      ps <- summary(res.cox)$sctest[3]
      # fit <- survfit(Surv(donor_survival_time, donor_vital_status) ~ type,data = os)
      # ggsurvplot(fit, data = os, risk.table = TRUE, pval = FALSE)
      out <- rbind(out,data.frame(cancer=cancer_type,gene=gene,n_pd=nrow(c_pd),
                                  n_nm=nrow(c_np),hazard_ratio=unname(hazard_ratio),p_wald=pw,
                                  p_log=pl,p_score=ps
                                  ))
      
    }
  }
}

sum(out$p_wald<0.05 & out$p_log<0.05 & out$p_score<0.05)
sum(all$p_wald<0.05 & all$p_log<0.05 & all$p_score<0.05)


###################
# venn diagram

ht29_sig <- ht29[ht29$Skew.logFDR > -log10(0.05) & ht29$project %in% c("cancer_snv_sig"),]$SNP
jurkat_sig <- jurkat[jurkat$Skew.logFDR > -log10(0.05) &  jurkat$project %in% c("cancer_snv_sig"),]$SNP
skmel28_sig <- skmel28[skmel28$Skew.logFDR > -log10(0.05) &  skmel28$project %in% c("cancer_snv_sig"),]$SNP

df  <- list(ht29=ht29_sig,jurkat=jurkat_sig,skmel28=skmel28_sig)

universe <- sort(unique(c(ht29_sig, jurkat_sig, skmel28_sig)))
Counts <- matrix(0, nrow=length(universe), ncol=3)
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% ht29_sig
  Counts[i,2] <- universe[i] %in% jurkat_sig
  Counts[i,3] <- universe[i] %in% skmel28_sig
}

colnames(Counts) <- c("ht29","jurkat","skmel28")

cols<-c("Red", "Green", "Blue")
vennDiagram(vennCounts(Counts), circle.col=cols)


##############

df <- rbind(ht29,skmel28,jurkat)
df <- df[df$Skew.logFDR > -log10(0.05) & df$project %in% c("cancer_snv_sig"),]
df <- df[!duplicated(df$SNP),]
snv <- sig[!duplicated(sig$id_samples),]
snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
snv <- snv[snv$SNP %in% df$SNP,]
df <- merge(snv,df, by = "SNP", all.y=TRUE)
df <- select(snv,c("SNP","icgc_donor_id"))
donor_mpra_snvs <- as.data.frame(table(df$icgc_donor_id))


pcawg <- fread("~/Documents/Fuxman/noncoding_cancer/data/pcawg/pcawg_sample.csv")
pcawg <- merge(pcawg,projects,by.x="dcc_project_code",by.y="V1", all.x=TRUE)
pdc2pcawg <- fread("~/Documents/Fuxman/noncoding_cancer/data/pcawg/pdc2pcawg.csv")

donor_predicted_driver_snvs <- as.data.frame(table(sig[!duplicated(sig$id_samples)]$icgc_donor_id))
donor_promoter_snvs <- as.data.frame(table(all_snvs$icgc_donor_id))

freq <- fread("~/Documents/Fuxman/noncoding_cancer/data/all_snvs/freq_samples.txt")
freq$sample <- str_split(freq$V2,"[.]",simplify = TRUE)[,1]
freq <- merge(freq, select(pcawg,c("icgc_donor_id","aliquot_id")), by.x="sample", by.y="aliquot_id", all.x=TRUE)

freq <- merge(freq, select(pdc2pcawg,c("Object ID","ICGC Donor")), by.x="sample", by.y="Object ID", all.x=TRUE)
freq[is.na(freq)] <- ""
freq$donor <- paste0(freq$icgc_donor_id,freq$`ICGC Donor`)
freq <- merge(freq,select(pcawg,c("V3","icgc_donor_id")),by.x="donor",by.y="icgc_donor_id", all.x=TRUE)
freq <- select(freq,c("donor","V1","V3"))
colnames(freq) <- c("donor","all_snvs","cancer_type")
freq <- merge(freq,donor_predicted_driver_snvs,by.x="donor",by.y="Var1",all.x=TRUE)
freq <- merge(freq,donor_promoter_snvs,by.x="donor",by.y="Var1",all.x=TRUE)
freq <- merge(freq,donor_mpra_snvs,by.x="donor",by.y="Var1",all.x=TRUE)
freq[is.na(freq)] <- 0
freq <- freq[!duplicated(freq),]
colnames(freq) <- c("donor","all_snvs","cancer_type","predicted_drivers","promoter_snvs","mpra_snvs")
freq <- as.data.frame(freq)
nums <- unlist(lapply(freq, is.numeric))  
# freq[,nums] <- log10(freq[,nums]+1)
freq <- freq[order(freq$cancer_type),]


pd <- select(freq,c(donor,cancer_type,predicted_drivers))
pd <- spread(pd,cancer_type,predicted_drivers)[,-1]

g <- select(freq,c(donor,cancer_type,all_snvs))
g <- spread(g,cancer_type,all_snvs)[,-1]

out <- data.frame()
cor_df <- matrix(NA,nrow = 0,ncol = 3) 
for ( i in 1:ncol(g)){
  pos2delete <- is.infinite(pd[,i]) | is.na(pd[,i]) | is.infinite(g[,i]) | is.na(g[,i])
  r <- cor.test(pd[,i][!pos2delete],g[,i][!pos2delete])
  cor_df <- rbind(cor_df,c(colnames(g)[i],unname(r$estimate),r$p.value))
  out <- cbind.fill(out,pd[,i],g[,i])
}
colnames(out) <- rep(colnames(g),each=2)
cor_df <- as.data.frame(cor_df)
colnames(cor_df) <- c("cancer_type","r","p")

cancer_median_snvs <- data.frame(cancer_type=colnames(g),
                                 median_svns=colMedians(as.matrix(g),na.rm = TRUE))

fwrite(freq,"~/Documents/Fuxman/noncoding_cancer/data/all_snvs/snvs_table.tsv",
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


pos2delete <- is.infinite(freq$all_snvs)  | is.infinite(freq$predicted_drivers)
r <- cor.test(freq$all_snvs[!pos2delete],freq$predicted_drivers[!pos2delete])

test <- fread("~/Downloads/snvs_table.tsv")

freq$ratio <- freq$predicted_drivers/freq$promoter_snvs
median(freq[freq$cancer_type=="head_neck",]$ratio, na.rm = TRUE)
test <- spread(freq,cancer_type, ratio)
fwrite(test,"~/Documents/Fuxman/noncoding_cancer/data/all_snvs/ratio_promoters.tsv",
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


