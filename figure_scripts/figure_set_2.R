library(data.table)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(rtracklayer)

process_df <- function(df, ids){
  thr <- 0.05
  df <- df[df$A.logPadj_BF > -log10(thr) | df$B.logPadj_BF > -log10(thr) ,]
  tmp <- df[df$project!="literature_snv" & df$SNP %in% ids,]
  tmp$project <- "literature_snv"
  df <- rbind(df,tmp)
  return(df)
  
}


qval_thresholding <- function(df,title){
  df <- as.data.frame(df)
  out <- data.frame()  
  for (thr in seq(0.01,0.2,0.01)){
    sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR),]
    sig_df <- sig_df[abs(sig_df$LogSkew) > 0,]
    cases <- table(sig_df$project)
    total <- table(df$project)
    total <- total[names(total) %in% names(cases)]
    tmp <- data.frame(g=names(cases),v=as.numeric(cases/total),t=rep(thr,length(names(cases))))
    out <- rbind(out,tmp)
    # colnames(out) <- names(cases)
  }
  
  p <- ggplot(out, aes(x=t, y=v, col=g)) + geom_line() + ggtitle(title) + xlab("qvalue") + ylab("validated fraction")
  out <- spread(out, g, v)
  return(out)
  
}


qval_thresholding_tf_family <- function(df,tf_info,sig_snv,title){
  
  df <- df[df$project %in% c("cancer_snv_sig"),]
  
  tf_info <- dplyr::select(tf_info,c("TF_Name","Family_Name"))
  tf_info <- tf_info[!duplicated(tf_info$TF_Name),]
  
  sig_snv$mut_id <- paste(sig_snv$mut_id,sig_snv$tf_name,sep = "_")
  
  snv <- sig_snv[!duplicated(sig_snv$mut_id),]
  snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
  
  snv <- merge(snv,tf_info,by.x = "tf_name",by.y = "TF_Name", all.x = TRUE)
  snv$mut_fam_id <- paste(snv$chr,snv$start,snv$REF,snv$ALT,snv$Family_Name,sep = ":")
  snv <- snv[!duplicated(snv$mut_fam_id),]
  snv <- snv[!is.na(snv$Family_Name)]
  
  df <- merge(snv,df,by = "SNP")
  
  df <- as.data.frame(df)
  out <- data.frame()
  
  tf_families <- c("bHLH","bZIP","C2H2 ZF","Ets","Forkhead","Homeodomain","Nuclear receptor","Rel","Sox","T-box")
  
  for (fam in tf_families){
    for (thr in seq(0.01,0.2,0.01)){
      sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR) & df$Family_Name==fam & df$project=="cancer_snv_sig",]
      sig_df <- sig_df[abs(sig_df$LogSkew) > 0,]
      cases <- nrow(sig_df)
      total <- nrow(df[df$Family_Name==fam & df$project=="cancer_snv_sig",])
      tmp <- data.frame(g=fam,v=as.numeric(cases/total),t=thr)
      out <- rbind(out,tmp)
      # colnames(out) <- names(cases)
    }
  }
  
  #out[is.na(out)] <- 0
  
  p <- ggplot(out, aes(x=t, y=v, col=g)) + geom_line() + ggtitle(title) + xlab("qvalue") + ylab("validated fraction")
  
  out <- spread(out, g, v)
  
  return(out)
  
}

gain_loss_VR <- function(df,sig_snv){
  
  df <- df[df$project %in% c("cancer_snv_sig"),]
  
  sig_snv$effect <- ifelse((sig_snv$gain_1==1 | sig_snv$gain_2==1 | sig_snv$gain_3==1), "gain","loss")
  
  sig_snv$mut_id_effect <- paste(sig_snv$mut_id,sig_snv$effect,sep = "_")
  
  snv <- sig_snv[!duplicated(sig_snv$mut_id_effect),]
  both <- names(table(snv$mut_id)[table(snv$mut_id)>1])
  snv$effect <- ifelse(snv$mut_id %in% both, "both",snv$effect)
  snv  <- snv[!duplicated(snv$mut_id),]
  snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
  
  df <- merge(snv,df,by = "SNP")
  
  df <- as.data.frame(df)
  out <- data.frame()
  
  effect_types <- unique(df$effect)[order(unique(df$effect))]
  
  for (e in effect_types){
    for (thr in seq(0.01,0.2,0.01)){
      sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR) & df$effect==e & df$project=="cancer_snv_sig",]
      sig_df <- sig_df[abs(sig_df$LogSkew) > 0,]
      cases <- nrow(sig_df)
      total <- nrow(df[df$effect==e & df$project=="cancer_snv_sig",])
      tmp <- data.frame(g=e,v=as.numeric(cases/total),t=thr)
      out <- rbind(out,tmp)
      # colnames(out) <- names(cases)
    }
  }
  out <- spread(out, g, v)
  return(out)
  
}

vr_by_mutFreq <- function(df,sig_snv){
  
  thr <- 0.05
  
  df <- df[df$project %in% c("cancer_snv_sig"),]
  
  snv <- sig_snv[!duplicated(sig_snv$id_samples),]
  snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
  counts <- as.data.frame(table(snv$SNP))
  
  df <- merge(counts,df,by.x = "Var1", by.y = "SNP")
  
  df <- as.data.frame(df)
  out <- data.frame()
  
  df <- df[order(df$Freq),]
  
  for (i in seq(1:5)){
    if(i < 5){
      sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR) & df$Freq==i & df$project=="cancer_snv_sig",]
      sig_df <- sig_df[abs(sig_df$LogSkew) > 0,]
      cases <- nrow(sig_df)
      total <- nrow(df[df$Freq==i & df$project=="cancer_snv_sig",])
    } else{
      sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR) & df$Freq>=i & df$project=="cancer_snv_sig",]
      sig_df <- sig_df[abs(sig_df$LogSkew) > 0,]
      cases <- nrow(sig_df)
      total <- nrow(df[df$Freq>=i & df$project=="cancer_snv_sig",])
    }
    tmp <- data.frame(g=i,v=as.numeric(cases/total),cases=cases,total=total)
    out <- rbind(out,tmp)
  }
  return(out)
  
}

fraction_over_under_expr <- function(df,sig_snv,freq_thr){
  
  df <- df[df$project %in% c("cancer_snv_sig"),]
  
  df$effect <- ifelse(df$LogSkew > 0, "overexpression","underexpression")
  
  snv <- sig_snv[!duplicated(sig_snv$mut_id),]
  snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
  df <- merge(snv,df, by = "SNP")
  
  counts_gene <- as.data.frame(table(df$promoter_name))
  counts_gene <- counts_gene[counts_gene$Freq>freq_thr,]
  
  thr <- 0.05
  sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR) & df$project=="cancer_snv_sig",]
  
  out <- data.frame()
  for (gene in counts_gene$Var1){
    c_over <- nrow(sig_df[sig_df$LogSkew > 0 & sig_df$promoter_name==gene,])/counts_gene[counts_gene==gene,]$Freq
    c_under <- nrow(sig_df[sig_df$LogSkew < 0 & sig_df$promoter_name==gene,])/counts_gene[counts_gene==gene,]$Freq
    
    c_out <- data.frame(gene=gene,over_fraction=c_over,under_fraction=c_under,total=counts_gene[counts_gene==gene,]$Freq)
    out <- rbind(out,c_out)
  }
  
  return(out)
  
}

tss_analysis <- function(df,tss,cell_line,all_tss_distante){
  
  tmp <- all_tss_distante
  df <- df[df$project %in% c("cancer_snv_sig"),]
  thr <- 0.05
  
  snv <- sig_snv[!duplicated(sig_snv$id_promoter),]
  snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
  snv <- dplyr::select(snv,c("promoter_name","SNP","start"))
  
  df <- merge(snv,df, by = "SNP")
  
  tss_distante <- data.frame()
  for (gene in unique(snv$promoter_name)){
    print(which(unique(snv$promoter_name)==gene))
    c_tmp <- tmp[tmp$V7==gene,]
    c_gene <- snv[snv$promoter_name==gene,]
    c_tss <- tss[tss$gene_name==gene,]
    
    if("-" %in% c_tss$strand){
      
      for( s in c_gene$SNP){
        
        c_tss <- tss[tss$gene_name==gene,]
        c_start <- c_gene[c_gene$SNP==s & c_gene$promoter_name==gene,]$start
        c_tss$dist <- -(c_tss$end - c_start)
        c_tss <- c_tss[order(abs(c_tss$dist)),]
        c_tss$dist <- -c_tss$dist
        dist_tss <- c_tss[c_tss$dist <= 250,]$dist[1]
        tss_distante <- rbind(tss_distante,data.frame(SNP=s,gene=gene,tss_dist=dist_tss))
      } 
      # for( k in 1:nrow(c_tmp)){
      #   
      #   c_tss <- tss[tss$gene_name==gene,]
      #   c_start <- c_tmp$V3[k]
      #   c_tss$dist <- -(c_tss$end - c_start)
      #   c_tss <- c_tss[order(abs(c_tss$dist)),]
      #   c_tss$dist <- -c_tss$dist
      #   dist_tss <- c_tss[c_tss$dist <= 250,]$dist[1]
      #   all_tss_distante <- rbind(all_tss_distante,data.frame(SNP=paste(c_tmp$V1[k],c_tmp$V3[k],sep = "_"),gene=gene,tss_dist=dist_tss))
      # }
      
    } else {
      for( s in c_gene$SNP){
        
        c_tss <- tss[tss$gene_name==gene,]
        c_start <- c_gene[c_gene$SNP==s & c_gene$promoter_name==gene,]$start
        c_tss$dist <-  c_start - c_tss$start
        c_tss <- c_tss[order(abs(c_tss$dist)),]
        dist_tss <- c_tss[c_tss$dist <= 250,]$dist[1]
        tss_distante <- rbind(tss_distante,data.frame(SNP=s,gene=gene,tss_dist=dist_tss))
      } 
      # for( k in 1:nrow(c_tmp)){
      #   
      #   c_tss <- tss[tss$gene_name==gene,]
      #   c_start <- c_tmp$V3[k]
      #   c_tss$dist <-  c_start - c_tss$start
      #   c_tss <- c_tss[order(abs(c_tss$dist)),]
      #   dist_tss <- c_tss[c_tss$dist <= 250,]$dist[1]
      #   all_tss_distante <- rbind(all_tss_distante,data.frame(SNP=paste(c_tmp$V1[k],c_tmp$V3[k],sep = "_"),gene=gene,tss_dist=dist_tss))
      # } 
      
    }
  }
  
  tss_distante <- tss_distante[tss_distante$tss_dist < 250 & tss_distante$tss_dist > -2000,]
  tss_distante <- tss_distante[order(tss_distante$tss_dist),]
  
  df$id_promoter <- paste(df$SNP,df$promoter_name)
  tss_distante$id_promoter <- paste(tss_distante$SNP,tss_distante$gene)
  df <- merge(df,tss_distante, by="id_promoter")
  df$tss_dist <- ifelse(df$tss_dist < -2000,-2000,ifelse(df$tss_dist > 250, 250,df$tss_dist))
  df$tss_dist[which(is.na(df$tss_dist))] <- 250
  df <- df[order(df$tss_dist),]
  
  out <- data.frame()
  for(start in seq(-2000,151,by = 1)){
    
    toFilter <- seq(start,start+99)
    c_df <- df[df$tss_dist %in% toFilter,]
    f_snv <- tss_distante[tss_distante$tss_dist %in% toFilter,]
    c_tmp <- tmp[tmp$tss_dist %in% toFilter,]
    c_df_sig <- c_df[(c_df$Skew.logFDR > -log10(thr)) & !is.na(c_df$Skew.logFDR) & c_df$project=="cancer_snv_sig",]
    
    out <- rbind(out,data.frame(range=toFilter[length(toFilter)/2],predicted=nrow(f_snv)/nrow(tss_distante),active=nrow(c_df)/nrow(df),all_snvs=nrow(c_tmp)/nrow(tmp),validation_rate = nrow(c_df_sig)/nrow(c_df)))
    
  }
  return(list(out=out,df=df))
}

get_signatures <- function(df,signatures){
  
  df <- df[df$project=="cancer_snv_sig",]
  thr <- 0.05
  df_jurkat <- merge(signatures, df, by = "SNP")
  df_jurkat$significant <- ifelse(df_jurkat$Skew.logFDR > -log10(thr),"differential","not")
  df_jurkat$significant[is.na(df_jurkat$significant)] <- "not"
  df_jurkat$significant <- as.factor(df_jurkat$significant)
  
  df_jurkat$apobec <- df_jurkat$BI_COMPOSITE_SNV_SBS2_P + df_jurkat$BI_COMPOSITE_SNV_SBS13_P + 
    df_jurkat$BI_COMPOSITE_SNV_SBS69_P
  df_jurkat$uv <- df_jurkat$BI_COMPOSITE_SNV_SBS7a_S + df_jurkat$BI_COMPOSITE_SNV_SBS7b_S +
    df_jurkat$BI_COMPOSITE_SNV_SBS7c_S + df_jurkat$BI_COMPOSITE_SNV_SBS3_P +
    df_jurkat$BI_COMPOSITE_SNV_SBS55_S + df_jurkat$BI_COMPOSITE_SNV_SBS67_S +
    df_jurkat$BI_COMPOSITE_SNV_SBS75_S
  
  df_jurkat$age <- df_jurkat$BI_COMPOSITE_SNV_SBS1_P + df_jurkat$BI_COMPOSITE_SNV_SBS5_P
  df_jurkat <- dplyr::select(df_jurkat,-c(grep("COMPOSITE",colnames(df_jurkat))))
  df_jurkat <- as.data.frame(df_jurkat[!duplicated(df_jurkat$SNP),])
  
  out <- data.frame()
  for(start in seq(0,0.9,by = 0.01)){
    
    c_df_apobec <- df_jurkat[(df_jurkat$apobec >= start) & (df_jurkat$apobec < start+0.1),]
    c_df_apobec_sig <- c_df_apobec[(c_df_apobec$Skew.logFDR > -log10(thr)) & !is.na(c_df_apobec$Skew.logFDR),]
    
    c_df_uv <- df_jurkat[(df_jurkat$uv >= start) & (df_jurkat$uv < start+0.1),]
    c_df_uv_sig <- c_df_uv[(c_df_uv$Skew.logFDR > -log10(thr)) & !is.na(c_df_uv$Skew.logFDR),]
    
    out <- rbind(out,data.frame(probability=(start+start+0.1)/2,
                                uv_validated = nrow(c_df_uv_sig),
                                uv_total = nrow(c_df_uv),
                                uv_vr = nrow(c_df_uv_sig)/nrow(c_df_uv),
                                apobec_validated =  nrow(c_df_apobec_sig),
                                apobec_total = nrow(c_df_apobec),
                                apobec_vr = nrow(c_df_apobec_sig)/nrow(c_df_apobec)))
    
  }
  
  out[is.na(out)]  <- 0
  
  
  return(out)
  
}

get_prognostics <- function(df,sig_snv,cell_line){
  
  promoters <- fread("~/Documents/Fuxman/noncoding_cancer/data/gencode/no_cds_promoter_v19_hs37.bed")
  promoters <- unique(promoters$V4)
  out <- data.frame()
  df <- df[df$project %in% c("cancer_snv_sig"),]
  thr <- 0.05
  snv <- sig_snv[!duplicated(sig_snv$mut_id),]
  snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
  snv <- dplyr::select(snv,c("promoter_name","SNP"))
  df <- merge(snv,df, by = "SNP")
  sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR),]
  
  genes_df <- data.frame()
  for(g in unique(sig_df$promoter_name)){
    c_log <- sig_df[sig_df$promoter_name==g,]$LogSkew
    gene_state <- ifelse(sum(c_log>0)==length(c_log),"up",ifelse(sum(c_log<0)==length(c_log),"down","mixed"))
    genes_df <- rbind(genes_df,data.frame(gene=g,gene_state=gene_state,stringsAsFactors = FALSE))
  }
  
  
  down_genes <- unique(genes_df[genes_df$gene_state=="down",]$gene)
  up_genes <- unique(genes_df[genes_df$gene_state=="up",]$gene)
  mixed_genes <- unique(genes_df[genes_df$gene_state=="mixed",]$gene)
  
  prognostics <- fread("~/Documents/Fuxman/noncoding_cancer/data/human_protein_atlas/counts.csv")
  
  s_d <- sum(down_genes %in% prognostics[prognostics$fcount>0,]$gene)
  s_u <- sum(up_genes %in% prognostics[prognostics$fcount>0,]$gene)
  s_m <- sum(mixed_genes %in% prognostics[prognostics$fcount>0,]$gene)
  s_t <- sum(promoters %in% prognostics[prognostics$fcount>0,]$gene)
  
  # res <- prop.test(x = c(s_f, t_f), n = c(length(sig_promoters), length(n_promoters)))
  # prop1 <- as.numeric(res$estimate[1])
  # prop2 <- as.numeric(res$estimate[2])
  out <- rbind(out, data.frame(cell_line=cell_line,case="prognostics",prognostics="favorable",down_prop=s_d/length(down_genes),
                               up_prop=s_u/length(up_genes),mixed_prop=s_m/length(mixed_genes),
                               all_genes_prop=s_t/length(promoters)
                               ))
  
  s_d <- sum(down_genes %in% prognostics[prognostics$ucount>0,]$gene)
  s_u <- sum(up_genes %in% prognostics[prognostics$ucount>0,]$gene)
  s_m <- sum(mixed_genes %in% prognostics[prognostics$ucount>0,]$gene)
  s_t <- sum(promoters %in% prognostics[prognostics$ucount>0,]$gene)
  
  # res <- prop.test(x = c(s_f, t_f), n = c(length(sig_promoters), length(n_promoters)))
  # prop1 <- as.numeric(res$estimate[1])
  # prop2 <- as.numeric(res$estimate[2])
  out <- rbind(out, data.frame(cell_line=cell_line,case="prognostics",prognostics="unfavorable",down_prop=s_d/length(down_genes),
                               up_prop=s_u/length(up_genes),mixed_prop=s_m/length(mixed_genes),
                               all_genes_prop=s_t/length(promoters)
  ))
  
  s_d <- sum(down_genes %in% prognostics[prognostics$total>0,]$gene)
  s_u <- sum(up_genes %in% prognostics[prognostics$total>0,]$gene)
  s_m <- sum(mixed_genes %in% prognostics[prognostics$total>0,]$gene)
  s_t <- sum(promoters %in% prognostics[prognostics$total>0,]$gene)
  
  # res <- prop.test(x = c(s_f, t_f), n = c(length(sig_promoters), length(n_promoters)))
  # prop1 <- as.numeric(res$estimate[1])
  # prop2 <- as.numeric(res$estimate[2])
  out <- rbind(out, data.frame(cell_line=cell_line,case="prognostics",prognostics="either",down_prop=s_d/length(down_genes),
                               up_prop=s_u/length(up_genes),mixed_prop=s_m/length(mixed_genes),
                               all_genes_prop=s_t/length(promoters)
  ))
  
  
  return(out)
}


get_signatures_probability <- function(df,signatures){
  
  df <- df[df$project=="cancer_snv_sig",]
  thr <- 0.05
  df_jurkat <- merge(signatures, df, by = "SNP")
  df_jurkat$significant <- ifelse(df_jurkat$Skew.logFDR > -log10(thr),"differential","not")
  df_jurkat$significant[is.na(df_jurkat$significant)] <- "not"
  df_jurkat$significant <- as.factor(df_jurkat$significant)
  
  df_jurkat$apobec <- df_jurkat$BI_COMPOSITE_SNV_SBS2_P + df_jurkat$BI_COMPOSITE_SNV_SBS13_P + 
    df_jurkat$BI_COMPOSITE_SNV_SBS69_P
  df_jurkat$uv <- df_jurkat$BI_COMPOSITE_SNV_SBS7a_S + df_jurkat$BI_COMPOSITE_SNV_SBS7b_S +
    df_jurkat$BI_COMPOSITE_SNV_SBS7c_S + df_jurkat$BI_COMPOSITE_SNV_SBS3_P +
    df_jurkat$BI_COMPOSITE_SNV_SBS55_S + df_jurkat$BI_COMPOSITE_SNV_SBS67_S +
    df_jurkat$BI_COMPOSITE_SNV_SBS75_S
  
  df_jurkat$age <- df_jurkat$BI_COMPOSITE_SNV_SBS1_P + df_jurkat$BI_COMPOSITE_SNV_SBS5_P
  df_jurkat <- dplyr::select(df_jurkat,-c(grep("COMPOSITE",colnames(df_jurkat))))
  df_jurkat <- as.data.frame(df_jurkat[!duplicated(df_jurkat$SNP),])
  
  apobec_positive <- sum(df_jurkat$apobec>0.5 & df_jurkat$significant=="differential")/sum(df_jurkat$apobec>0.5)
  apobec_negative <- sum(df_jurkat$apobec <= 0.5 & df_jurkat$significant=="differential")/sum(df_jurkat$apobec <= 0.5)
  uv_positive <- sum(df_jurkat$uv>0.5 & df_jurkat$significant=="differential")/sum(df_jurkat$uv>0.5)
  uv_negative <- sum(df_jurkat$uv <= 0.5 & df_jurkat$significant=="differential")/sum(df_jurkat$uv <= 0.5)
  
  out <- data.frame(apobec_positive = apobec_positive,
                    apobec_negative = apobec_negative,
                    uv_positive = uv_positive,
                    uv_negative = uv_negative)
  
  out$apobec_pos_error <- sqrt(apobec_positive*(1-apobec_positive)/sum(df_jurkat$apobec>0.5))
  out$apobec_neg_error <- sqrt(apobec_negative*(1-apobec_negative)/sum(df_jurkat$apobec <= 0.5))
  out$uv_pos_error <- sqrt(uv_positive*(1-uv_positive)/sum(df_jurkat$uv>0.5))
  out$uv_neg_error <- sqrt(uv_negative*(1-uv_negative)/sum(df_jurkat$u <= 0.5))
  out$N_apobec_pos <- sum(df_jurkat$apobec>0.5)
  out$N_apobec_neg <- sum(df_jurkat$apobec <= 0.5)
  out$N_uv_pos <- sum(df_jurkat$uv > 0.5)
  out$N_uv_neg <- sum(df_jurkat$uv <= 0.5)
  
  
  return(out)
  
}


#####################
# validation rate by qvalue

sig_snv <- fread("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/predicted_drivers.tsv")
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
tf_info <- fread("~/Documents/Fuxman/noncoding_cancer/motif_benchmark/tf_info_processed.csv")

ht29_qval_thresholding <- qval_thresholding(ht29,"HT29")
jurkat_qval_thresholding <-  qval_thresholding(jurkat,"Jurkat")
skmel28_qval_thresholding <- qval_thresholding(skmel28,"SKMel28")


fwrite(ht29_qval_thresholding,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/ht29.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(jurkat_qval_thresholding,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/jurkat.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(skmel28_qval_thresholding,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/skmel28.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

##############
# validation rate by TF family

ht29_qval_thresholding_tf_family <- qval_thresholding_tf_family(ht29,tf_info,sig_snv,"HT29 by TF family")
jurkat_qval_thresholding_tf_family <- qval_thresholding_tf_family(jurkat,tf_info,sig_snv,"Jurkat by TF family")
skmel28_qval_thresholding_tf_family <- qval_thresholding_tf_family(skmel28,tf_info,sig_snv,"SKMel by TF family")

fwrite(ht29_qval_thresholding_tf_family,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/ht29_tf.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(jurkat_qval_thresholding_tf_family,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/jurkat_tf.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(skmel28_qval_thresholding_tf_family,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/skmel28_tf.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


################
# gain, loss, both validation rate

ht29_effect <- gain_loss_VR(ht29,sig_snv)
jurkat_effect <- gain_loss_VR(jurkat,sig_snv)
skmel28_effect <- gain_loss_VR(skmel28,sig_snv)


fwrite(ht29_effect,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/ht29_effect.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(jurkat_effect,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/jurkat_effect.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(skmel28_effect,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/skmel28_effect.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


################
# validation rate by mutation frequency

ht29_freq <- vr_by_mutFreq(ht29,sig_snv)
jurkat_freq <- vr_by_mutFreq(jurkat,sig_snv)
skmel28_freq <- vr_by_mutFreq(skmel28,sig_snv)

out <- cbind(ht29_freq,jurkat_freq[,2],skmel28_freq[,2])
colnames(out) <- c("SNV frequency,","HT-29","Jurkat","SK-MEL-28")

fwrite(out,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/heatmap_freq.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



#########################
# fraction of SNVs by +/-

ht29_expr <- fraction_over_under_expr(ht29,sig_snv,4)
jurkat_expr <- fraction_over_under_expr(jurkat,sig_snv,4)
skmel28_expr <- fraction_over_under_expr(skmel28,sig_snv,4)

out_expr <- merge(merge(ht29_expr[,-4],jurkat_expr[,-4], by = "gene"),skmel28_expr[,-4], by = "gene")
colnames(out_expr) <- c("gene","ht29+","ht29-","jurkat+","jurkat-","skmel28+","skmel28-")
rownames(out_expr) <- out_expr$gene
test <- pheatmap::pheatmap(out_expr[,-1], cluster_cols = FALSE,fontsize_row = 8, 
                   treeheight_row = 0, treeheight_col = 0)
test <- out_expr[test$tree_row$order,]
fwrite(test,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/heatmap_jurkat.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

  
#########################
# validation rate by TSS distance

# all_snvs <- fread("~/Documents/Fuxman/noncoding_cancer/data/SNVs_in_promoters/snvs_promoters_uniq.txt")
# all_snvs <- add_column(all_snvs, start=all_snvs$V2, .after="V1")
# all_snvs$start <- all_snvs$start - 1
# all_snvs <- dplyr::select(all_snvs,-c(V3,V4))
# write.table(all_snvs,"~/Documents/Fuxman/noncoding_cancer/data/SNVs_in_promoters/snvs_promoters_uniq.bed",
#             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

load("~/Documents/Fuxman/noncoding_cancer/data/SNVs_in_promoters/all_tss_distance.RData")
genes <- fread("~/Documents/Fuxman/noncoding_cancer/data/gencode/gene_coordinates_gencode.v19.bed")
promoters <- fread("~/Documents/Fuxman/noncoding_cancer/data/gencode/promoter_coordinates_gencode.v19.bed")
tss <- readGFF("~/Documents/Fuxman/noncoding_cancer/data/gencode/tss.gtf")
tss <- tss[tss$type %in% c("gene","transcript") & tss$gene_type %in% c("protein_coding") & !(tss$seqid %in% c('chrM')),!(colnames(tss) %in% c("score","phase"))]
tss$id <- paste(tss$start, tss$end, tss$gene_name)
tss <- tss[!duplicated(tss$id),]

ht29_tss <- tss_analysis(ht29,tss,"ht29",all_tss_distante)
jurkat_tss <- tss_analysis(jurkat,tss,"jurkat",all_tss_distante)
skmel28_tss <- tss_analysis(skmel28,tss,"skmel28",all_tss_distante)

# out <- cbind(ht29_tss,jurkat_tss[,2],skmel28_tss[,2])
# colnames(out) <- c("position","ht29","jurkat","skmel28")
# 
# fwrite(out,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/tss_vr.tsv",
#        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(c$out,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/ht29_tss.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(jurkat_tss$out,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/jurkat_tss.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(skmel28_tss$out,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/skmel28_tss.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

nrow(skmel28_tss$df[skmel28_tss$df$Skew.logFDR > -log10(0.05) & skmel28_tss$df$tss_dist < -250,])/nrow(skmel28_tss$df[skmel28_tss$df$Skew.logFDR > -log10(0.05),])

#########################
# mutational signatures

signatures <- fread("~/Documents/Fuxman/noncoding_cancer/signatures/u_cancer_mutations_id_seq_last_prob_signatures.bed")
signatures$SNP <- paste(signatures$Chr,signatures$Position,signatures$Ref,signatures$Alt, sep = ":")

ht29_signature <- get_signatures(ht29,signatures)
jurkat_signature <- get_signatures(jurkat,signatures)
skmel28_signature <- get_signatures(skmel28,signatures)


fwrite(ht29_signature,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/ht29_signature.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(jurkat_signature,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/jurkat_signature.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(skmel28_signature,"~/Documents/Fuxman/noncoding_cancer/mpra/figure_tables/qval_threshold/skmel28_signature.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

ht29_signature <- get_signatures_probability(ht29,signatures)
jurkat_signature <- get_signatures_probability(jurkat,signatures)
skmel28_signature <- get_signatures_probability(skmel28,signatures)

total <- rbind(ht29_signature,jurkat_signature,skmel28_signature)
total$cell_line <- c("ht29","jurkat","skmel28")

fwrite(total,"~/Documents/Fuxman/noncoding_cancer/sig_tfs/signatures.tsv",
       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



#########################
# prognostics TFs


get_prognostics <- function(df,sig_snv,cell_line){
  
  sig_snv$id_tf <- paste(sig_snv$mut_id,sig_snv$tf_name)
  out <- data.frame()
  df <- df[df$project %in% c("cancer_snv_sig"),]
  thr <- 0.05
  snv <- sig_snv[!duplicated(sig_snv$id_tf),]
  snv$SNP <- paste(snv$chr,snv$start,snv$REF,snv$ALT,sep = ":")
  snv$effect <- ifelse(snv$gain_1==1 | snv$bone_soft_tissue_significant==1 | snv$gain_3==1,"gain","loss")
  snv <- dplyr::select(snv,c("tf_name","effect","SNP"))
  df <- merge(snv,df, by = "SNP")
  df$regulation <- NA
  df$regulation <- ifelse(df$LogSkew > 0, "up",ifelse(df$LogSkew < 0,"down",df$regulation))
  sig_df <- df[(df$Skew.logFDR > -log10(thr)) & !is.na(df$Skew.logFDR),]
  gain_up <- unique(sig_df[sig_df$effect=="gain" & sig_df$regulation=="up"]$tf_name)
  gain_down <- unique(sig_df[sig_df$effect=="gain" & sig_df$regulation=="down"]$tf_name)
  loss_up <- unique(sig_df[sig_df$effect=="loss" & sig_df$regulation=="up"]$tf_name)
  loss_down <- unique(sig_df[sig_df$effect=="loss" & sig_df$regulation=="down"]$tf_name)
  
  log10( (length(gain_up)/length(loss_down))/ (length(gain_down)/length(loss_up)) )
  
  genes_df <- data.frame()
  for(g in unique(sig_df$promoter_name)){
    c_log <- sig_df[sig_df$promoter_name==g,]$LogSkew
    gene_state <- ifelse(sum(c_log>0)==length(c_log),"up",ifelse(sum(c_log<0)==length(c_log),"down","mixed"))
    genes_df <- rbind(genes_df,data.frame(gene=g,gene_state=gene_state,stringsAsFactors = FALSE))
  }
  
  
  down_genes <- unique(genes_df[genes_df$gene_state=="down",]$gene)
  up_genes <- unique(genes_df[genes_df$gene_state=="up",]$gene)
  mixed_genes <- unique(genes_df[genes_df$gene_state=="mixed",]$gene)
  
  prognostics <- fread("~/Documents/Fuxman/noncoding_cancer/data/human_protein_atlas/counts.csv")
  
  s_d <- sum(down_genes %in% prognostics[prognostics$fcount>0,]$gene)
  s_u <- sum(up_genes %in% prognostics[prognostics$fcount>0,]$gene)
  s_m <- sum(mixed_genes %in% prognostics[prognostics$fcount>0,]$gene)
  s_t <- sum(promoters %in% prognostics[prognostics$fcount>0,]$gene)
  
  # res <- prop.test(x = c(s_f, t_f), n = c(length(sig_promoters), length(n_promoters)))
  # prop1 <- as.numeric(res$estimate[1])
  # prop2 <- as.numeric(res$estimate[2])
  out <- rbind(out, data.frame(cell_line=cell_line,case="prognostics",prognostics="favorable",down_prop=s_d/length(down_genes),
                               up_prop=s_u/length(up_genes),mixed_prop=s_m/length(mixed_genes),
                               all_genes_prop=s_t/length(promoters)
  ))
  
  s_d <- sum(down_genes %in% prognostics[prognostics$ucount>0,]$gene)
  s_u <- sum(up_genes %in% prognostics[prognostics$ucount>0,]$gene)
  s_m <- sum(mixed_genes %in% prognostics[prognostics$ucount>0,]$gene)
  s_t <- sum(promoters %in% prognostics[prognostics$ucount>0,]$gene)
  
  # res <- prop.test(x = c(s_f, t_f), n = c(length(sig_promoters), length(n_promoters)))
  # prop1 <- as.numeric(res$estimate[1])
  # prop2 <- as.numeric(res$estimate[2])
  out <- rbind(out, data.frame(cell_line=cell_line,case="prognostics",prognostics="unfavorable",down_prop=s_d/length(down_genes),
                               up_prop=s_u/length(up_genes),mixed_prop=s_m/length(mixed_genes),
                               all_genes_prop=s_t/length(promoters)
  ))
  
  s_d <- sum(down_genes %in% prognostics[prognostics$total>0,]$gene)
  s_u <- sum(up_genes %in% prognostics[prognostics$total>0,]$gene)
  s_m <- sum(mixed_genes %in% prognostics[prognostics$total>0,]$gene)
  s_t <- sum(promoters %in% prognostics[prognostics$total>0,]$gene)
  
  # res <- prop.test(x = c(s_f, t_f), n = c(length(sig_promoters), length(n_promoters)))
  # prop1 <- as.numeric(res$estimate[1])
  # prop2 <- as.numeric(res$estimate[2])
  out <- rbind(out, data.frame(cell_line=cell_line,case="prognostics",prognostics="either",down_prop=s_d/length(down_genes),
                               up_prop=s_u/length(up_genes),mixed_prop=s_m/length(mixed_genes),
                               all_genes_prop=s_t/length(promoters)
  ))
  
  
  return(out)
}





