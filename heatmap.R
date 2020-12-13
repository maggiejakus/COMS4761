library(gplots)
library(stringr)

all_scores <- list(read.table("gmm_diff.nikpay.to_plot.gz"), 
read.table("gmm_diff.lambert.to_plot.gz"),
read.table("gmm_diff.liu.to_plot.gz"),
read.table("gmm_diff.michailidou.to_plot.gz"),
read.table("gmm_diff.lopez-isac.to_plot.gz"),
read.table("gmm_diff.mahajan.to_plot.gz"),
read.table("gmm_diff.malik.to_plot.gz"),
read.table("gmm_diff.christophersen.to_plot.gz"))

PGS_annot <- read.table("PGS_annot_sums.lopez-isac.to_plot.gz", stringsAsFactors = F)
annots <- str_split(PGS_annot[1,], fixed("."), simplify = T)[,1]

all_scores <- as.matrix(do.call("cbind", all_scores))
author_names <- c("Coronary Artery Disease", "Alzheimers", "IBS", "Breast Cancer",
                  "Systemic Sclerosis", "Type 2 Diabetes", "Stroke", 
                  "Atrial Fibrillation")

#colnames(all_scores) <- author_names
par(mar=c(7,4,8,2)+.1)
png("test.png", width = 800, height = 800, bg = "white")
heatmap.2(all_scores, scale="row", 
          key=TRUE, symkey=FALSE, density.info="none", labCol = author_names, labRow = annots,
          cexRow=1,cexCol=1,margins=c(12,8),trace="none",srtCol=45)
dev.off()


