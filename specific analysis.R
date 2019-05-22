#new dataframe with only erlotinib columns: 249 - 307 e=erlotinib
#e_treated <- NCI_TPW_gep_treated[,which(colnames(NCI_TPW_gep_treated) == grep ("erlotinib", colnames(NCI_TPW_gep_treated))]
e_treated <- NCI_TPW_gep_treated[, 249:307]
e_untreated <- NCI_TPW_gep_untreated[, 249:307]

#save as dataframe
e_treated <- as.data.frame(e_treated)
e_untreated <- as.data.frame(e_untreated)

library(GeneNet)
z.transform(e_foldchange)
#e_foldchange_z_transformed: treated-untreated + z-Transformation + PCA 
e_foldchange <- e_treated - e_untreated
median_gene <- apply (e_foldchange, 1, function (x){
  + median(x) })
sd_gene <- apply(e_foldchange, 1, function (x){
  + sd(x) })
e_foldchange_z_transformed <- (e_folzdchange-median(median_gene)) / sd(sd_gene)

pca <- prcomp(e_foldchange_z_transformed)
plot(pca$rotation[,1], pca$rotation[,2], main = "z-transformed")
pca <- prcomp(e_foldchange)
plot(pca$rotation[,1], pca$rotation[,2], main = "foldchange")

#select cell lines with highest variance
var_cell_line <- apply(e_foldchange_z_transformed, 2, function (x) {
  + var(x)})
cell_line_var_greater_75quantile <- var_cell_line [which (var_cell_line > quantile(var_cell_line,0.75))]
cell_line_var_decreasing_top15 <- sort(cell_line_var_greater_75quantile, decreasing = TRUE)

table_cell_lines_var_top15 <- cbind(names(cell_line_var_decreasing_top15), cell_line_var_decreasing_top15)
rownames(table_cell_lines_var_top15) <- c(1:15)
colnames(table_cell_lines_var_top15) <- c("cell line", "variance")
#kable(): simple table formatting function
kable(table_cell_lines_var_top15)


#PCA visualization with factoextra package: 
fviz_eig(pca) #shows percentage explained by each PC
fviz_pca_ind(pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

#PCA with transformed matrix (each dot represents a gene):
transformed_e_foldchange_z_transformed <- t(e_foldchange_z_transformed)
pca <- prcomp(transformed_erlotinib_normalized_z_transformed)
plot(pca$rotation[,1], pca$rotation[,2])
text(pca$rotation, labels = rownames(erlotinib_normalized_z_transformed), cex = 0.4, pos = 3)

#PCA visualization with factoextra package (according to www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-t-prcomp-vs-princomp/) : 
fviz_eig(pca) #shows percentage explained by each PC

#plot of cell lines with gradient colors by the quality of representation
fviz_pca_ind(pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

#plot of genes showing their contribution
fviz_pca_var(pca, col.var ="contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

#get contribution of genes to the PCs
results.genes <- get_pca_var(pca)
genes_pca_highest_contribution <- results.genes$contrib
