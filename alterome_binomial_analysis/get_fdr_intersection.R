
setwd("~/Documents/Fuxman/noncoding_cancer/tfbs_profile/probabilities/")

# promoters <- read.table("../../data/gencode/promoter_coordinates_gencode_v19_hs37.bed", stringsAsFactors = FALSE)
# tf.info <- na.omit(as.data.table(read.csv("../tf_info_scored.csv", stringsAsFactors = FALSE)), cols="score_threshold")

# gain1 <- read.csv("fdr/gain_1.txt")
# gain2 <- read.csv("fdr/gain_2.txt")
# gain3 <- read.csv("fdr/gain_3.txt")
# loss1 <- read.csv("fdr/loss_1.txt")
# loss2 <- read.csv("fdr/loss_2.txt")
# loss3 <- read.csv("fdr/loss_3.txt")

total_set <- fread("filtered_fdr/mutation_total_set190426.csv")
total_ids <- unique(paste(total_set$promoter_id, total_set$pwm_id, sep = "_"))

fdr_threshold <- 0.01
cancer <- "uterus"

gain1 <- read.csv(paste0("fdr/",cancer,"_gain_1.txt"))
gain2 <- read.csv(paste0("fdr/",cancer,"_gain_2.txt"))
gain3 <- read.csv(paste0("fdr/",cancer,"_gain_3.txt"))
loss1 <- read.csv(paste0("fdr/",cancer,"_loss_1.txt"))
loss2 <- read.csv(paste0("fdr/",cancer,"_loss_2.txt"))
loss3 <- read.csv(paste0("fdr/",cancer,"_loss_3.txt"))

gain1$id <- paste(gain1$promoter, gain1$pwm, sep = "_")
gain2$id <- paste(gain2$promoter, gain2$pwm, sep = "_")
gain3$id <- paste(gain3$promoter, gain3$pwm, sep = "_")

loss1$id <- paste(loss1$promoter, loss1$pwm, sep = "_")
loss2$id <- paste(loss2$promoter, loss2$pwm, sep = "_")
loss3$id <- paste(loss3$promoter, loss3$pwm, sep = "_")

# ig <- intersect(intersect(gain1[gain1$fdr < fdr_threshold, ]$id, gain2[gain2$fdr < fdr_threshold, ]$id),gain3[gain3$fdr < fdr_threshold, ]$id)
#  
# il <- intersect(intersect(loss1[loss1$fdr < fdr_threshold, ]$id, loss2[loss2$fdr < fdr_threshold, ]$id),loss3[loss3$fdr < fdr_threshold, ]$id)
# 
# ig_df <- gain1[ gain1$id %in% ig, ]
# il_df <- loss1[ loss1$id %in% il, ]

ig_1 <- gain1[ gain1$id %in% total_ids, ]
ig_2 <- gain2[ gain2$id %in% total_ids, ]
ig_3 <- gain3[ gain3$id %in% total_ids, ]

il_1 <- loss1[ loss1$id %in% total_ids, ]
il_2 <- loss2[ loss2$id %in% total_ids, ]
il_3 <- loss3[ loss3$id %in% total_ids, ]

ig_df <- rbind(ig_1,ig_2,ig_3)
ig_df <- ig_df[!duplicated(ig_df$id),]
il_df <- rbind(il_1,il_2,il_3)
il_df <- il_df[!duplicated(il_df$id),]

ig_df <- ig_df[ , -which(names(ig_df) %in% c("type","id"))]
il_df <- il_df[ , -which(names(il_df) %in% c("type","id"))]

# gain1 <- gain1[ gain1$id %in% ig,]
# gain2 <- gain2[ gain2$id %in% ig,]
# gain3 <- gain3[ gain3$id %in% ig,]
# loss1 <- loss1[ loss1$id %in% il,]
# loss2 <- loss2[ loss2$id %in% il,]
# loss3 <- loss3[ loss3$id %in% il,]

# gain <- rbind(gain1,gain2,gain3)
# loss <- rbind(loss1,loss2,loss3)

# write.table(gain,paste0("filtered_fdr/xnp/",cancer,"_gain.csv"), sep = ",", row.names = FALSE, quote = FALSE)
# write.table(loss,paste0("filtered_fdr/xnp/",cancer,"_loss.csv"), sep = ",", row.names = FALSE, quote = FALSE)


# write.table(loss1,paste0("filtered_fdr/xnp/",cancer,"_loss_1.csv"), sep = ",", row.names = FALSE, quote = FALSE)
# write.table(loss2,paste0("filtered_fdr/xnp/",cancer,"_loss_2.csv"), sep = ",", row.names = FALSE, quote = FALSE)
# write.table(loss3,paste0("filtered_fdr/xnp/",cancer,"_loss_3.csv"), sep = ",", row.names = FALSE, quote = FALSE)
# write.table(gain1,paste0("filtered_fdr/xnp/",cancer,"_gain_1.csv"), sep = ",", row.names = FALSE, quote = FALSE)
# write.table(gain2,paste0("filtered_fdr/xnp/",cancer,"_gain_2.csv"), sep = ",", row.names = FALSE, quote = FALSE)
# write.table(gain3,paste0("filtered_fdr/xnp/",cancer,"_gain_3.csv"), sep = ",", row.names = FALSE, quote = FALSE)

write.table(ig_df, paste0("filtered_fdr/intersection_cancers/all_thresholds/",cancer,"_gain_intersection.txt"), 
            quote = FALSE, row.names = FALSE, sep = " ", col.names = TRUE)

write.table(il_df, paste0("filtered_fdr/intersection_cancers/all_thresholds/",cancer,"_loss_intersection.txt"), 
            quote = FALSE, row.names = FALSE, sep = " ", col.names = TRUE)

