
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

library(RColorBrewer)
data.frame(subset(annotation_of_celllines_per_drug, select = ,
           matchRetVal = match(colors$Drug, colors$Color))
for (i in 1:15){
  if (colnames(annotation_of_celllines_per_drug) == colors[i,2]) {
    color_all_drugs <- colors[i,1]
  }}


 
#of celllines treated with each drug in a column (61 cell lines)
annotation_sorted_by_cell_lines <- matrix(, nrow = 15, ncol = 0)
as.data.frame(cellline_annotation)

names (cellline_annotation$Cell_Line_Name[33]) <- "SK-MEL_2_"
for (i in 1:61){
  columns_of_one_cell_line <- c(grep (cellline_annotation$Cell_Line_Name[i], colnames(fold_changes), value = TRUE))
  columns_of_one_cell_line <- c(columns_of_one_cell_line, rep(NA, 15-length(columns_of_one_cell_line)))
  annotation_sorted_by_cell_lines <- cbind (annotation_sorted_by_cell_lines, columns_of_one_cell_line)
}

colnames(annotation_sorted_by_cell_lines) <- cellline_annotation$Cell_Line_Name
