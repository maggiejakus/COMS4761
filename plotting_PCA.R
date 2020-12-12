library(ggplot2)
library(viridis)
library(cowplot)
library(reshape2)
theme_set(theme_cowplot())

gmm_PCA <- read.table("PCA_components.GMM.to_plot.gz")

gmm_plot <- data.frame(pc1 = gmm_PCA$V1, pc2 = gmm_PCA$V2, pc3 = gmm_PCA$V3, pc4 = gmm_PCA$V4)

rownames(gmm_plot) <- c("Coronary Artery Disease", "Atrial Fibrillation", "Type 2 Diabetes", "Breast Cancer",
                        "IBS", "Alzheimers", "Systemic Sclerosis", "Stroke")

#PCA plot
ggplot(gmm_plot, aes(pc1, pc2)) + geom_point(size =3) + geom_text(
  label=rownames(gmm_plot), 
  nudge_x = c(0.012, .009, .009, .008, .003, .006, .01, .004), 
  size = 5,
  check_overlap = T) +
  labs(x = "PC1",
       y = "PC2",
       title = "PCA Decomposition - GMM Method")

