library(data.table)
library(stringr)
library(tidyr)
library(pheatmap)
library(viridis)
library(ReactomePA)
library(biomaRt)
library(data.table)
library(readxl)

# sig <- fread("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/probabilities/filtered_fdr/mutation_total_set_sig_only.csv")
# sig <- sig[sig$tf_name %in% tf_info$TF_Name,]
# sig <- sig[sig$pwm_cisbp_id %in% tf_info$Motif_ID,]
# sig <- sig[!duplicated(sig),]
# fwrite(sig,"~/Documents/Fuxman/noncoding_cancer/tfbs_profile/predicted_drivers.tsv",
#        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# sig$cancer_type <- ifelse(sig$all_samples_significant==1,"all_samples",sig$cancer_type)
sig <- fread("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/predicted_drivers.tsv")
length(unique(sig$mut_id))
length(unique(sig$promoter_name))
length(unique(sig$tf_name))
sig$effect <- ifelse((sig$gain_1==1 | sig$gain_2==1 | sig$gain_3==1), "gain","loss")


s <- fread("~/Documents/Fuxman/noncoding_cancer/data/for_ryan/snvs_mpra_total.csv")
s$neg <- ifelse(s$alt=="A","T",ifelse(s$alt=="T","A",ifelse(s$alt=="G","C","G")))
s$id <- paste(s$chr,s$start,s$alt, sep = "_")
s$id_neg <- paste(s$chr,s$start,s$neg, sep = "_")

ids <- c(s[s$type=="literature_snv",]$id,s[s$type=="literature_snv",]$id_neg)

literature_snv <- sig[(sig$mut_id %in% ids ) & !duplicated(sig$mut_id),]

sig_gain_id <- unique(paste(sig$cancer_type,sig$mut_id,sig$effect))
sig_gain_id <- str_split(sig_gain_id," ", simplify = TRUE)
sig_gain_id <- setnames(as.data.frame(sig_gain_id),c("cancer","mut","effect"))
sig_gain_id$id <- paste(sig_gain_id$cancer, sig_gain_id$mut)
test <- aggregate(effect ~ id, data = sig_gain_id, paste, collapse = ",")
test$effect <- ifelse(grepl(",",test$effect),"two",test$effect)
test$id <- str_split(test$id," ", simplify = TRUE)[,1]
View(table(paste(test$id,test$effect)))

pan_cancer <- sig[sig$all_samples_significant==1,]
pan_cancer$cancer_type <- "all_samples"
sig_gain_id <- unique(paste(pan_cancer$cancer_type,pan_cancer$mut_id,pan_cancer$effect))
sig_gain_id <- str_split(sig_gain_id," ", simplify = TRUE)
sig_gain_id <- setnames(as.data.frame(sig_gain_id),c("cancer","mut","effect"))
sig_gain_id$id <- paste(sig_gain_id$cancer, sig_gain_id$mut)
test <- aggregate(effect ~ id, data = sig_gain_id, paste, collapse = ",")
test$effect <- ifelse(grepl(",",test$effect),"two",test$effect)
test$id <- str_split(test$id," ", simplify = TRUE)[,1]
View(table(paste(test$id,test$effect)))

####################

samples <- fread("~/Documents/Fuxman/noncoding_cancer/data/pcawg/cancer_freq.csv")
sig$id_samples <- paste(sig$aliquot_id,sig$cancer_type,sig$promoter_name,sep = "_")
sig$cacer_id <-paste(sig$cancer_type,sig$id_samples)
test <- sig[!duplicated(sig$cacer_id),]
test$cacer_id <-paste(test$cancer_type,test$promoter_name)
test <- as.data.frame(table(test$cacer_id))
test$Var1 <- as.character(test$Var1)
test <- separate(test, Var1, into = c("cancer", "gene"), sep = " ")
test <- merge(test,samples, by.x = "cancer", by.y = "V1")
test$Freq <- 100*test$Freq/test$V2
test <- test[,-c(4)]
genes <- unique(test[test$Freq>=20,]$gene)
test <- test[test$gene %in% genes,]
tmp <- spread(test, cancer, Freq, fill = 0)
colnames(tmp)[1] <- "a_genes"
missing <- samples$V1[!samples$V1 %in% colnames(tmp)]
for (m in missing){
  tmp[m] <- 0
}
tmp <- tmp[order(rowSums(tmp[,-1]),decreasing=TRUE),order(colnames(tmp))]
rownames(tmp) <- tmp$a_genes
# tmp[tmp>0] <- 1

pheatmap(tmp[,-1], cluster_rows = FALSE, cluster_cols = FALSE, 
         show_rownames = TRUE, color = viridis(100))

# tmp[tmp==0] <- 0.000001
# hist.data = hist(unlist(tmp[,-1]), plot=FALSE, breaks = 20)
# hist.data$counts = ifelse(is.finite(log10(hist.data$counts)),log10(hist.data$counts),0)
# plot(hist.data, ylab='log10(Frequency)')

write.table(tmp,"~/Documents/Fuxman/noncoding_cancer/figures/heatmap.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


pan_cancer <- sig[sig$all_samples_significant==1,]
pan_cancer$cancer_type <- "all_samples"
sig_gain_id <- unique(paste(pan_cancer$cancer_type,pan_cancer$mut_id,pan_cancer$effect))
sig_gain_id <- str_split(sig_gain_id," ", simplify = TRUE)
sig_gain_id <- setnames(as.data.frame(sig_gain_id),c("cancer","mut","effect"))
sig_gain_id$id <- paste(sig_gain_id$cancer, sig_gain_id$mut)
test <- aggregate(effect ~ id, data = sig_gain_id, paste, collapse = ",")
test$effect <- ifelse(grepl(",",test$effect),"two",test$effect)
test$id <- str_split(test$id," ", simplify = TRUE)[,1]
View(table(paste(test$id,test$effect)))


###############
# Fitness and essential

cgc <- fread("~/Documents/Fuxman/noncoding_cancer/data/cosmic/Census_allThu Aug  2 18_21_58 2018.csv")
promoters <- fread("~/Documents/Fuxman/noncoding_cancer/data/gencode/no_cds_promoter_v19_hs37.bed")
sig_promoters <- unique(promoters[promoters$V8 %in% unique(sig$promoter_id),]$V4)
n_promoters <- unique(promoters$V4[!promoters$V4 %in% sig_promoters])
sig_tfs <- unique(sig$tf_name)
no_sig_tfs <- unique(tf_info$TF_Name[!tf_info$TF_Name %in% sig_tfs])



fitness <- fread("~/Downloads/fitness.csv")
fitness <- fitness[fitness$`#NumTKOHits` >=3,]

s_f <- sum(sig_promoters %in% fitness$`#GENE`)
t_f <- sum(n_promoters %in% fitness$`#GENE`)
c_f <- sum(cgc$`Gene Symbol` %in% fitness$`#GENE`)
t_c <-  unique(promoters$V4[!promoters$V4 %in% cgc$`Gene Symbol`])

res <- prop.test(x = c(s_f, t_f), n = c(length(sig_promoters), length(n_promoters)))
res <- prop.test(x = c(c_f, t_f), n = c(length(cgc$`Gene Symbol`), length(t_c)))
sum(unique(promoters$V4) %in% fitness$`#GENE`)/length(unique(promoters$V4))

fitness <- fread("~/Downloads/Achilles_common_essentials.csv")

s_f <- sum(sig_promoters %in% fitness$V1)
t_f <- sum(n_promoters %in% fitness$V1)
c_f <- sum(cgc$`Gene Symbol` %in% fitness$V1)
t_c <-  unique(promoters$V4[!promoters$V4 %in% cgc$`Gene Symbol`])

res <- prop.test(x = c(s_f, t_f), n = c(length(sig_promoters), length(n_promoters)))
res <- prop.test(x = c(c_f, t_f), n = c(length(cgc$`Gene Symbol`), length(t_c)))
sum(unique(promoters$V4) %in% fitness$V1)/length(unique(promoters$V4))

N <- 813
p <- 0.14
sqrt(p*(1-p)/N)

############
# Prognosis

prognostics <- fread("~/Documents/Fuxman/noncoding_cancer/data/human_protein_atlas/counts.csv")

sig_fav <- prop.test(c(sum(sig_promoters %in% prognostics[prognostics$fcount>0,]$gene),
                       sum(n_promoters %in% prognostics[prognostics$fcount>0,]$gene)),
                       c(length(sig_promoters),
                       length(n_promoters)))
sig_unf <- prop.test(c(sum(sig_promoters %in% prognostics[prognostics$ucount>0,]$gene),
                       sum(n_promoters %in% prognostics[prognostics$ucount>0,]$gene)),
                     c(length(sig_promoters),
                       length(n_promoters)))
sig_eit <- prop.test(c(sum(sig_promoters %in% prognostics[prognostics$total>0,]$gene),
                       sum(n_promoters %in% prognostics[prognostics$total>0,]$gene)),
                     c(length(sig_promoters),
                       length(n_promoters)))


cgc_fav <- prop.test(c(sum(cgc$`Gene Symbol` %in% prognostics[prognostics$fcount>0,]$gene),
                       sum(t_c %in% prognostics[prognostics$fcount>0,]$gene)),
                     c(length(cgc$`Gene Symbol`),
                       length(t_c)))
cgc_unf <- prop.test(c(sum(cgc$`Gene Symbol` %in% prognostics[prognostics$ucount>0,]$gene),
                       sum(t_c %in% prognostics[prognostics$ucount>0,]$gene)),
                     c(length(sig_promoters),
                       length(t_c)))
cgc_eit <- prop.test(c(sum(cgc$`Gene Symbol` %in% prognostics[prognostics$total>0,]$gene),
                       sum(t_c %in% prognostics[prognostics$total>0,]$gene)),
                     c(length(cgc$`Gene Symbol`),
                       length(t_c)))

sum(unique(promoters$V4) %in% prognostics[prognostics$fcount>0,]$gene)/length(unique(promoters$V4))
sum(unique(promoters$V4) %in% prognostics[prognostics$ucount>0,]$gene)/length(unique(promoters$V4))
sum(unique(promoters$V4) %in% prognostics[prognostics$total>0,]$gene)/length(unique(promoters$V4))


N <- 813
p <- 0.3517835
sqrt(p*(1-p)/N)













