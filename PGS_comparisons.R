#import all GMM diff data

nikpay <- read.table("PGS_diff.nikpay.to_plot.gz")
lambert <- read.table("PGS_diff.lambert.to_plot.gz")
liu <- read.table("PGS_diff.liu.to_plot.gz")
mich <- read.table("PGS_diff.michailidou.to_plot.gz")
lopez <- read.table("PGS_diff.lopez-isac.to_plot.gz")
mahajan <- read.table("PGS_diff.mahajan.to_plot.gz")
malik <- read.table("PGS_diff.malik.to_plot.gz")
chris <- read.table("PGS_diff.christophersen.to_plot.gz")

PGS_annot <- read.table("PGS_annot_sums.lopez-isac.to_plot.gz", stringsAsFactors = F)
annots <- str_split(PGS_annot[1,], fixed("."), simplify = T)[,1]

diff_mtx <- data.frame(Annotations = annots, CAD = nikpay, Alzheimers = lambert,
                       IBS = liu, BC = mich, Sclerosis = lopez, Diabetes = mahajan, 
                       Stroke = malik, AtrialFib = chris)
colnames(diff_mtx) <- c("Annotations", "CAD", "Alzheimers", "IBS", "Breast Cancer",
                        "Sclerosis", "Diabetes", "Stroke", "Atrial Fib")


diff_mtx <- melt(diff_mtx, id.vars = "Annotations")

ggplot(diff_mtx, aes(y = Annotations, x = value, color = variable)) +
  geom_point(size = 3) + geom_hline(aes(yintercept = (1:nrow(diff_mtx))+0.5),  color = "grey80") +
  geom_vline(aes(xintercept = 0)) +
  labs(x = "Functional Annotation Score Difference",
       y = "Annotations",
       title = "PGS",
       color = "Disease") 



