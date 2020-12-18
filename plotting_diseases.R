library(ggplot2)
library(viridis)
library(cowplot)
library(reshape2)
theme_set(theme_cowplot())

PCA <- read.table("PCA_components.mahajan.to_plot.gz")
weights <- read.table("weights.mahajan.to_plot.gz")
means <- read.table("means.mahajan.to_plot.gz")
probs <- read.table("probabilities.mahajan.to_plot.gz")
rsids <- read.table("pca_rsids.mahajan.to_plot.gz", stringsAsFactors = F)
#supersnps <- read.table("superclump.snps.mahajan.to_plot", stringsAsFactors = F)
annot <- read.table("annot.mahajan.to_plot.gz")
annot_names <- read.table("annot_names.mahajan.to_plot.gz", stringsAsFactors = F)

for_plot <- data.frame(pc1 = PCA$V1, pc2 = PCA$V2, pc3 = PCA$V3, pc4 = PCA$V4,
                       prob = probs[,1], super_snp = 1)
#for_plot$super_snp[rsids[,1] %in% supersnps[,1]] <- 0
#for_plot$super_snp <- as.factor(for_plot$super_snp)

#PCA plot
#ggplot(for_plot, aes(pc1, pc2, color = super_snp)) + geom_point()


#plot PCAs
#for(i in 1:nrow(annot_names)){
#  for_plot <- data.frame(pc1 = PCA$V1, pc2 = PCA$V2, pc3 = PCA$V3, pc4 = PCA$V4,
#                         annot = annot[,i])
#  the_plot <- ggplot(for_plot, aes(pc1, pc2, color = annot)) + geom_point() +
#    scale_color_viridis() + labs(caption = annot_names[i,1])
#  ggsave(paste0("pca.", annot_names[i,1], ".png"), the_plot, "png",
#         path = "plots/", width = 6, height = 6)
#}

#plot PCA1 vs PCA2
#ggplot(for_plot, aes(pc1, pc2, color = prob)) + geom_point() +
#  scale_color_viridis()


#########################################
#difference plot
#final_weights <- read.table("total_weights.to_plot.gz", stringsAsFactors = F)
#PGS_annot <- read.table("PGS_annot.to_plot.gz", stringsAsFactors = F)

final_weights <- read.table("total_weights.mahajan.to_plot.gz", stringsAsFactors = F)
PGS_annot <- read.table("PGS_annot_sums.mahajan.to_plot.gz", stringsAsFactors = F)
superclump_annot <- read.table("super_sums.mahajan.to_plot.gz", stringsAsFactors = F)

#PGS_annot_second <- read.table("PGS_annot_sums.to_plot.gz", stringsAsFactors = F)

PGS_diff <- read.table("PGS_diff.mahajan.to_plot.gz", stringsAsFactors = F)
gmm_diff <- read.table("gmm_diff.mahajan.to_plot.gz", stringsAsFactors = F)
super_diff <- read.table("super_diff.mahajan.to_plot.gz", stringsAsFactors = F)

final_weights <- as.data.frame(t(final_weights))
colnames(final_weights) <- c("Annotation", "Nonrisky", "Risky")
final_weights$Nonrisky <- as.numeric(final_weights$Nonrisky)
final_weights$Risky <- as.numeric(final_weights$Risky)

PGS_annot_T <- as.data.frame(t(PGS_annot))
final_weights$PRS <- as.numeric(PGS_annot_T[,2])

superclump_annot_T <- as.data.frame(t(superclump_annot))
final_weights$super <- as.numeric(superclump_annot_T[,2])

final_weights <- melt(final_weights, id.vars = "Annotation")


diff_mtx <- data.frame(Annotations = PGS_annot_T[,1], PRS = PGS_diff, super = super_diff,
                       gmm = gmm_diff)
#diff_mtx$Annotations <- PGS_annot_T[,1]
#diff_mtx$PRS <- PGS_diff
#diff_mtx$super <- super_diff
#diff_mtx$gmm <- gmm_diff
colnames(diff_mtx) <- c("Annotations", "PRS", "Super", "GMM")

diff_mtx <- melt(diff_mtx, id.vars = "Annotations")



#ggplot(final_weights, aes(y = Annotation, x = value, fill = variable)) +
# geom_bar(stat = "identity", position = position_dodge())

#ggplot(final_weights, aes(y = Annotation, x = value, color = variable)) +
#  geom_point(size = 3) + geom_hline(aes(yintercept = (1:nrow(final_weights))+0.5)) +
#  xlab("Functional Annotation Score")

diff_mtx$nice_anno <- gsub("_", " ", lapply(strsplit(diff_mtx$Annotations,
                                                     ".", fixed = T), function(x) x[1]))

#plot to make 3 methods, 24 annotations differences
ggplot(diff_mtx, aes(y = nice_anno, x = value, color = variable)) +
  geom_point(size = 3) + geom_hline(aes(yintercept = (1:nrow(diff_mtx))+0.5),  color = "grey80") +
  geom_vline(aes(xintercept = 0)) +
  labs(x = "Functional Annotation Score Difference",
       y = "Annotations",
       title = "Systemic Sclerosis - Lopez-Isac et al",
       color = "Method") +
  scale_color_manual(values = c("skyblue1", "palegreen1", "palevioletred1"))




#############################################################################
#this tells whether nonrisky vs risky gmm annots are different from each other
nonrisky_gmm <- read.table('gmm_nonrisky.ttest.mahajan.to_plot.gz', stringsAsFactors = F)
risky_gmm <- read.table('gmm_risky.ttest.mahajan.to_plot.gz', stringsAsFactors = F)
gmm_tstats <- rep(0, ncol(risky_gmm))
for(i in 1:ncol(risky_gmm)){
  gmm_tstats[i] <- t.test(nonrisky_gmm[nonrisky_gmm[,i] != 0,i], risky_gmm[risky_gmm[,i] != 0, i])$stat
}

nonrisky_superclump <- read.table('super_nonrisky.ttest.mahajan.to_plot.gz', stringsAsFactors = F)
risky_superclump <- read.table('super_risky.ttest.mahajan.to_plot.gz', stringsAsFactors = F)
super_tstats <- rep(0, ncol(risky_superclump))
for(i in 1:ncol(risky_superclump)){
  if((length(unique(risky_superclump[,i])) > 1)  & (length(unique(nonrisky_superclump[,i]))>1)){
    super_tstats[i] <- t.test(nonrisky_superclump[nonrisky_superclump[,i] != 0,i], risky_superclump[risky_superclump[,i] != 0, i])$stat
  }
}

nonrisky_PGS <- read.table('PGS_nonrisky.ttest.mahajan.to_plot.gz', stringsAsFactors = F)
risky_PGS <- read.table('PGS_risky.ttest.mahajan.to_plot.gz', stringsAsFactors = F)
PGS_tstats <- rep(0, ncol(risky_PGS))
for(i in 1:ncol(risky_PGS)){
  if((length(unique(risky_PGS[,i])) > 1)  & (length(unique(nonrisky_PGS[,i]))>1)){
    PGS_tstats[i] <- t.test(nonrisky_PGS[nonrisky_PGS[,i] != 0,i], risky_PGS[risky_PGS[,i] != 0, i])$stat
  }
}

write.csv(PGS_tstats, "mahajan_PGStstat.csv")
write.csv(super_tstats, "mahajan_supertstat.csv")
write.csv(gmm_tstats, "mahajan_gmmtstat.csv")
