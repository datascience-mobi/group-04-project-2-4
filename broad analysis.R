
library(RColorBrewer)
library(plyr)

fold_changes <- NCI_TPW_gep_treated - NCI_TPW_gep_untreated
pca <- prcomp(fold_changes)
plot(pca$rotation[,1], pca$rotation[,2])

cb <- matrix(, nrow = 70, ncol = 0)
as.data.frame(drug_annotation)

for (i in 1:15){
  columns_of_one_drug <- c(grep (drug_annotation$Drug[i], colnames(fold_changes), value = TRUE))
  columns_of_one_drug <- c(columns_of_one_drug, rep(NA, 70-length(columns_of_one_drug)))
  cb <- cbind (cb, columns_of_one_drug)
}

colnames(cb) <- drug_annotation$Drug

colors <- brewer.pal(n = 15, name = "Dark2")
plot(pca$rotation[,1], pca$rotation[,2], col = colors)