
#PCA of all celllines and drugs

#calculate fold change due to drug treatment
fold_changes <- NCI_TPW_gep_treated - NCI_TPW_gep_untreated

#matrix annotation_of_celllines_per_drug contains all names 
#of celllines treated with each drug in a column (15 drugs and columns)
annotation_of_celllines_per_drug <- matrix(, nrow = 70, ncol = 0)
as.data.frame(drug_annotation)

for (i in 1:15){
  columns_of_one_drug <- c(grep (drug_annotation$Drug[i], colnames(fold_changes), value = TRUE))
  columns_of_one_drug <- c(columns_of_one_drug, rep(NA, 70-length(columns_of_one_drug)))
  annotation_of_celllines_per_drug <- cbind (annotation_of_celllines_per_drug, columns_of_one_drug)
}

colnames(annotation_of_celllines_per_drug) <- drug_annotation$Drug

#PCA
pca <- prcomp(fold_changes)

#color all erlotinib cell lines in red, rest in green
color_palette <- rainbow(15)
color_erlotinib <- ifelse (colnames(annotation_of_celllines_per_drug) == "erlotinib", color_palette[1], color_palette[8])
plot(pca$rotation[,1], pca$rotation[,2], col = color_erlotinib, pch = 19, xlab = "PC1", ylab = "PC2")

#define vector with 15 colors, one for each drug
colors <- cbind(color_palette, colnames(annotation_of_celllines_per_drug))
colnames(colors) <- c("Color", "Drug")

plot(pca$rotation[,1], pca$rotation[,2], col = colors, pch = 19, xlab = "PC1", ylab = "PC2")