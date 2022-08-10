library(stringr)
library(data.table)
require(ggplot2)



# promoter_set <- read.table("../gencode/promoter_gencode_v19_sequences_hs37.txt", header = TRUE, sep = "\t", 
#                            stringsAsFactors = FALSE)

create_pwm <- function(pwm){
  pwm <- as.matrix(t(pwm))
  if(any(pwm==0)){
    pwm <- (pwm+0.001)/colSums((pwm+0.001))
  }
  return(pwm)
}

score_pwm <- function(pwm) {
  npos<-ncol(pwm)
  ic<-numeric(length=npos)
  for (i in 1:npos) {
    ic[i]<-2 + sum(sapply(pwm[, i], function(x) { 
      if (x > 0) { x*log2(x) } else { 0 }
    }))
  }    
  return(ic)
}

read_vcf <- function(vcf_path){
  vcf <- read.table(vcf_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  header <- c("chr","pos","ref","alt","merged","total","id","r_ref","r_alt")
  vcf$id <- paste(vcf$CHR,vcf$POS,sep = "_")
  vcf <- vcf[!grepl(",",vcf$ALT),]
  vcf <- cbind(vcf,str_split_fixed(vcf$X20.AD, ",", 2), stringsAsFactors = FALSE)
  vcf[8:9] <- lapply(vcf[8:9], as.integer)
  colnames(vcf) <- header
  vcf <- vcf[nchar(vcf$alt)==1,]
  vcf <- vcf[nchar(vcf$ref)==1,]
  vcf <- unique(vcf)
  # vcf <- vcf[vcf$total >=20,]
  # vcf <- vcf[vcf$r_ref>0 & vcf$r_alt>0,]
  # vcf$pval_ref <- mapply(bt_ref, vcf$r_ref, vcf$total)
  # vcf$pval_alt <- mapply(bt_alt, vcf$r_ref, vcf$total)
  return(vcf)
}

bt_ref <- function(a, b, p = 0.5) {
  binom.test(a, b, 0.5, alternative=c("greater"), conf.level = 0.95)$p.value
}

bt_alt <- function(a, b, p = 0.5) {
  binom.test(a, b, 0.5, alternative=c("less"), conf.level = 0.95)$p.value
}

get_out <-function(gabpa,vcf,chip_exp,tf_threshold){
  gabpa$id <- paste(gabpa$V1,gabpa$V2,sep = "_")
  
  vcf <- vcf[vcf$id %in% gabpa$id,]
  gabpa <- gabpa[gabpa$id %in% vcf$id,]
  if(nrow(vcf)>0){
    gabpa <- gabpa[!duplicated(gabpa[,c('id')]),]
    
    gabpa$wt_score <- pmax(gabpa$V7,gabpa$V8)
    gabpa$A <- pmax(gabpa$V9,gabpa$V10)
    gabpa$C <- pmax(gabpa$V11,gabpa$V12)
    gabpa$G <- pmax(gabpa$V13,gabpa$V14)
    gabpa$T <- pmax(gabpa$V15,gabpa$V16)
    
    vcf <- vcf[order(vcf$id),]
    gabpa <- gabpa[order(gabpa$id),]
    
    tmp <- cbind(vcf,gabpa)
    score <- c()
    for (i in 1:nrow(tmp)){
      i_score <- tmp[i,grepl(tmp$alt[i], colnames(tmp))]
      if(i_score>0){
        score <- c(score,i_score)
      } else {
        score <- c(score,NA)
      }
      
    }
    tmp$alt_score <- score
    tmp <- tmp[,!colnames(tmp) %in% c("merged","V1","V2","V3","V4","id","id.1","V5","V6","V7",
                                      "V8","V9","V10","V11","V12","V13","V14","V15","V16",
                                      "A","C","G","T")]
    tmp <- tmp[complete.cases(tmp), ]
    if (nrow(tmp)>0){
      tmp <- tmp[(tmp$alt_score>=tf_threshold | tmp$wt>=tf_threshold),]
      # out <- tmp[tmp$pval_ref<0.05 | tmp$pval_alt<0.05,]
      tmp$tf_threshold <- tf_threshold
      tmp$pwm <- chip_exp[chip_exp$encode_id %in% strsplit(basename(vcf_path_1),"_")[[1]][1],]$Motif_ID
      tmp$tf <- chip_exp[chip_exp$encode_id %in% strsplit(basename(vcf_path_1),"_")[[1]][1],]$TF_Name
      # out$label <- ifelse(out$pval_ref<0.05,"wt","a")
      tmp$motif_pred <- ifelse(tmp$wt_score>tmp$alt_score,1,0)
    } else{
      tmp <- c()
    }
    
  } else {
    tmp <- c()
  }
  return(tmp)
}

get_vcf <- function(vcf1,vcf2){
  test <- rbind(vcf1,vcf2)
  vcf <- test[!duplicated(test[,c('id')]),]
  test <- data.table(test)
  test <- test[,list(total_t = sum(total), total_r = sum(r_ref),total_a = sum(r_alt)), by = "id"]
  vcf$total <- test$total_t
  vcf$r_ref <- test$total_r
  vcf$r_alt <- test$total_a
  vcf <- get_allelic_imbalance(vcf)
  return(vcf)
}

get_allelic_imbalance <- function(vcf){
  vcf <- vcf[vcf$total >=20,]
  vcf <- vcf[vcf$r_ref>0 & vcf$r_alt>0,]
  vcf$pval_ref <- mapply(bt_ref, vcf$r_ref, vcf$total)
  vcf$pval_alt <- mapply(bt_alt, vcf$r_ref, vcf$total)
  vcf$q_ref <- p.adjust(vcf$pval_ref, method = "fdr")
  vcf$q_alt <- p.adjust(vcf$pval_alt, method = "fdr")
  return(vcf)
}

get_pr <- function(final,wt){
  flag <- ifelse(wt=="wt",1,0)
  final$diff <- abs(final$wt_score - final$alt_score)
  precision <- c()
  recall <- c()
  for (i in seq(0.5,7,0.5)){
    p <- sum(final$allelic_imbalance==1 &
               final$allele_count==flag)
    tp <- sum(final$allelic_imbalance==1 & 
                final$allele_count==flag &
                final$motif_pred==flag &
                final$diff>=i)
    pp <- sum( final$motif_pred==flag &
                 final$diff>=i)
    precision <- c(precision,100*tp/pp)
    # print(c(i,p,tp,pp))
    recall <- c(recall,100*tp/p)
  }
  
  precision[is.na(precision)] <- 0
  recall[is.na(recall)] <- 0
  
  df <- data.frame(i=seq(0.5,7,0.5),
                   precision=precision,
                   recall=recall)
  return(df)
}

###### CREATE CONTINGENCY TABLE

get_values_tf_threshold <- function(final){
  ll <- sum(final$allelic_imbalance==1 & final$motif_pred==1 & 
              final$allele_count==1 & final$alt_score<final$tf_threshold
  )
  
  gg <- sum(final$allelic_imbalance==1 & final$motif_pred==0 & 
              final$allele_count==0 & final$wt_score<final$tf_threshold
  )
  
  lg <- sum(final$allelic_imbalance==1 & final$motif_pred==0 & 
              final$allele_count==1 & final$wt_score<final$tf_threshold
  )
  
  gl <- sum(final$allelic_imbalance==1 & final$motif_pred==1 & 
              final$allele_count==0 & final$alt_score<final$tf_threshold
  )
  
  return(c(ll,lg,gl,gg,
           sum(final$motif_pred==1 & final$alt_score<final$tf_threshold),
           sum(final$motif_pred==0 & final$wt_score<final$tf_threshold),
           sum(final$allelic_imbalance==1 & final$allele_count==1),
           sum(final$allelic_imbalance==1 & final$allele_count==0),
           nrow(final)
  ))
}


#### FIS

get_values_fis <- function(fis_threshold,final){
  ll <- sum(final$allelic_imbalance==1 & final$motif_pred==1 &
              final$allele_count==1 & final$fis>=fis_threshold
  )
  
  gg <- sum(final$allelic_imbalance==1 & final$motif_pred==0 & 
              final$allele_count==0 & final$fis>=fis_threshold
  )
  
  lg <- sum(final$allelic_imbalance==1 & final$motif_pred==0 & 
              final$allele_count==1 & final$fis>=fis_threshold
  )
  
  gl <- sum(final$allelic_imbalance==1 & final$motif_pred==1 & 
              final$allele_count==0 & final$fis>=fis_threshold
  )
  
  return(c(ll,lg,gl,gg,
           sum(final$motif_pred==1 & final$fis>=fis_threshold),
           sum(final$motif_pred==0 &final$fis>=fis_threshold),
           sum(final$allelic_imbalance==1 & final$allele_count==1),
           sum(final$allelic_imbalance==1 & final$allele_count==0),
           nrow(final)
  ))
}


get_statistics_total <- function(values){
  true_p <- values[1] + values[4]
  pred_p <- values[5] + values[6]
  pos <- values[7] + values[8]
  precision <- true_p/pred_p
  recall <- true_p/pos
  f_val <- 2*(precision*recall/(precision+recall))
  f_05 <- 1.25*(precision*recall/(0.25*precision+recall))
  acc <- (true_p)/(true_p+values[2]+values[3])
  return(c(precision,recall,f_val,f_05,acc))
}

get_statistics_gl <- function(values,flag){
  if(flag=="loss"){
    true_p <- values[1]
    pred_p <- values[5]
    pos <- values[7]
    div <- values[3]
  } else{
    true_p <- values[4]
    pred_p <- values[6]
    pos <- values[8]
    div <- values[2]
  }
  precision <- true_p/pred_p
  recall <- true_p/pos
  f_val <- 2*(precision*recall/(precision+recall))
  f_05 <- 1.25*(precision*recall/(0.25*precision+recall))
  acc <- (true_p)/(true_p+div)
  return(c(precision,recall,f_val,f_05,acc))
}

get_stats <- function(final){
  stats_total <- as.data.frame(matrix(NA, nrow = 0, ncol = 5))
  stats_loss <- as.data.frame(matrix(NA, nrow = 0, ncol = 5))
  stats_gain <- as.data.frame(matrix(NA, nrow = 0, ncol = 5))
  for (fis_threshold in seq(0,6)){
    print(fis_threshold)
    values <- get_values_fis(fis_threshold,final)
    stats_total <- rbind(stats_total,get_statistics_total(values))
    stats_loss <- rbind(stats_loss,get_statistics_gl(values,"loss"))
    stats_gain <- rbind(stats_gain,get_statistics_gl(values,"gain"))
  }
  
  stats_total <- rbind(stats_total,get_statistics_total(get_values_tf_threshold(final)))
  stats_loss <- rbind(stats_loss,get_statistics_gl(get_values_tf_threshold(final),"loss"))
  stats_gain <- rbind(stats_gain,get_statistics_gl(get_values_tf_threshold(final),"gain"))
  
  stats_total <- cbind(0:7,stats_total)
  stats_loss <- cbind(0:7,stats_loss)
  stats_gain <- cbind(0:7,stats_gain)
  
  colnames(stats_total) <- colnames(stats_loss) <- colnames(stats_gain) <- c("threshold","precision","recall","F","F_05","acc")
  
  return(list(total=stats_total,loss=stats_loss,gain=stats_gain))
}

tf_info_ic <- function(tf.info){
  pwm_length <- c()
  pwm_score <- c()
  for (i in 1:nrow(tf.info)){
    print(i)
    pwm.file <- paste("../../motif_benchmark/homo_sapiens_cisbp/pwms_all_motifs/",
                      tf.info$file[i],sep = "")
    pwm <- create_pwm(read.table(pwm.file, header = TRUE, row.names = 1))
    pwm_score <- c(pwm_score,sum(score_pwm(pwm)))
    pwm_length <- c(pwm_length,ncol(pwm))
  }
  
  tf.info$pwm_length <- pwm_length
  tf.info$pwm_score <- pwm_score
  return(tf.info)
}

process_final <- function(final,qval){
  final$allelic_imbalance <- ifelse(final$q_ref<qval | final$q_alt<qval, 1,0)
  final$allele_count <- ifelse(final$r_ref>final$r_alt, 1,0)
  final$fis <- abs(final$wt_score-final$alt_score)
  return(final)
}


# final <- as.data.frame(matrix(NA, nrow = 0, ncol = 17))
# count_files <- list.files("allele_counts_all/")


#for (i in chip_exp$encode_id[47:length(chip_exp$encode_id)]){
#for (i in chip_exp$encode_id[1:22]){
# for (i in chip_exp$encode_id){
#   print(i)
#   vcf_files <- count_files[grepl(i,count_files)]
# 
#   vcf_path_1 <- paste0("allele_counts_all/",vcf_files[1])
#   vcf_path_2 <- paste0("allele_counts_all/",vcf_files[2])
# 
#   tf_pwm <- chip_exp[chip_exp$encode_id %in% strsplit(basename(vcf_path_1),"_")[[1]][1],]$index
#   tf_threshold <- tf.info[tf_pwm,]$score_threshold
#   gabpa <- tryCatch(read.table(paste0("alterome_all/tf_",tf_pwm,"_peak.bed")),error=function(e) NULL)
# 
#   if(!is.null(gabpa)){
#     vcf1 <- read_vcf(vcf_path_1)
#     vcf2 <- read_vcf(vcf_path_2)
#     vcf <- get_vcf(vcf1,vcf2)
# 
#     #print(c(i,nrow(get_out(gabpa,vcf,chip_exp,tf_threshold))))
#     final <- rbind(final,get_out(gabpa,vcf,chip_exp,tf_threshold))
#   }
# }


# save(final, file = "../../tfbs_profile/chip_seq_analysis/allelic_imbalance_all_peak.RData")

setwd("~/Documents/Fuxman/noncoding_cancer/data/encode_fastq/")

chip_exp <- read.csv("../../tfbs_profile/chip_seq_analysis/chip_exp_old.csv", stringsAsFactors = FALSE)
tf.info <- na.omit(as.data.table(read.csv("../../tfbs_profile/tf_info_scored.csv",stringsAsFactors = FALSE)), 
                   cols="score_threshold")
tf.info <- tf_info_ic(tf.info)
tf.info.chip <- tf.info[tf.info$Motif_ID %in% chip_exp$Motif_ID,]
tf.info.chip <- tf.info.chip[order(match(tf.info.chip$Motif_ID,chip_exp$Motif_ID))]
chip_exp$pwm_length <- tf.info.chip$pwm_length
chip_exp$pwm_ic <- tf.info.chip$pwm_score

chip_exp <- chip_exp[order(chip_exp$pwm_ic),]

load("../../tfbs_profile/chip_seq_analysis/allelic_imbalance_all_peak.RData")

final <- process_final(final,0.05)
# top <- process_final(top,0.05)
# bottom <- process_final(bottom,0.05)

print(paste("Total allelic imbalance:",sum(final$allelic_imbalance==1)))

# out <- final[final$allelic_imbalance==1,]
# 
# write.table(out,"../for_ryan/chipseq_allelic_imbalance.csv", sep = ",",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)

# print(paste("PRED WT: ",sum(final$motif_pred==1 &
#                               final$alt_score<final$tf_threshold
# )))

# print(paste("PRED ALT: ",sum(final$motif_pred==0 &
#                               final$wt_score<final$tf_threshold
# )))

out_all <- get_stats(final)
out_bottom <- get_stats(bottom)
out_top <- get_stats(top)


melt_total <- melt(stats_total, id="threshold")
ggplot(data=melt_total,
       aes(x=threshold, y=value, colour=variable)) +
  geom_line() + 
  ggtitle("Total")

melt_loss <- melt(stats_loss, id="threshold")
ggplot(data=melt_loss,
       aes(x=threshold, y=value, colour=variable)) +
  geom_line() + 
  ggtitle("Loss")

melt_gain <- melt(stats_gain, id="threshold")
ggplot(data=melt_gain,
       aes(x=threshold, y=value, colour=variable)) +
  geom_line() + 
  ggtitle("Gain")





