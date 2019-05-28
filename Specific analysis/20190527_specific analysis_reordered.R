#new matrix only with samples/columns treated with erlotinib  (e=erlotinib)
e_treated <- NCI_TPW_gep_treated[,grep ("erlotinib", colnames(NCI_TPW_gep_treated))]
e_untreated <- NCI_TPW_gep_untreated[,grep ("erlotinib", colnames(NCI_TPW_gep_untreated))]
e_foldchange <- e_treated - e_untreated

#e_foldchange_z_transformed: z-Transformation 
library(GeneNet)
z.transform(e_foldchange)
median_gene <- apply (e_foldchange, 1, function (x){
  + median(x) })
sd_gene <- apply(e_foldchange, 1, function (x){
  + sd(x) })
e_foldchange_z_transformed <- (e_foldchange-median(median_gene)) / sd(sd_gene)


# table of 15 cell lines with highest variance to show most regulated cell lines
#select 15 cell lines with highest variance (greater than 75% quantile, sorted by decreasing value)
var_cell_line <- apply(e_foldchange_z_transformed, 2, function (x) {
  + var(x)})
cell_line_var_greater_75quantile <- sort(var_cell_line [which (var_cell_line > quantile(var_cell_line,0.75))], decreasing = TRUE)
rm(var_cell_line)

#create a table containing the name and the variance of the 15 cell lines with highest variance
table_cell_lines_var_top15 <- cbind(names(cell_line_var_greater_75quantile), cell_line_var_greater_75quantile)
rownames(table_cell_lines_var_top15) <- c(1:nrow(table_cell_lines_var_top15))
colnames(table_cell_lines_var_top15) <- c("cell line", "variance")
#kable(): simple table formatting function
kable(table_cell_lines_var_top15)


#PCA to show most regulated cell lines
#PCA with transformed matrix (each point represents a sample):
e_foldchange_matrix_transformed <- t(e_foldchange_z_transformed)
pca <- prcomp(e_foldchange_matrix_transformed)
plot(pca$rotation[,1], pca$rotation[,2])
text(pca$rotation, labels = rownames(e_foldchange_z_transformed), cex = 0.4, pos = 3)

#PCA visualization with factoextra package (according to www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-t-prcomp-vs-princomp/) : 
#plot of cell lines with gradient colors by their quality of representation
fviz_pca_ind(pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) #cos2 values correspond to their quality of representation


#bar plot showing the percentage explained by each PC
'fviz_eig(pca) 
#plot of genes showing their contribution to the two PCs 
fviz_pca_var(pca, col.var ="contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
#get contribution of genes to the PCs
results.genes <- get_pca_var(pca)
genes_pca_highest_contribution <- results.genes$contrib '


#mean of gene expression of each gene over all cell lines
'e_treated_mean_over_cell_lines <- rowMeans(e_treated)
e_untreated_mean_over_cell_lines <- rowMeans(e_untreated)'
e_foldchange_mean_over_cell_lines <- rowMeans(e_foldchange) #equal to e_treated_mean_over_cell_lines - e_untreated_mean_over_cell_line

#determine the p-value for a paired two-sample t-test 
p_values <- sapply(1:nrow(e_treated), function(x) {
  t.test(e_treated[x,], e_untreated[x,],paired= T)$p.value}) # perform t-test and save p-values of each gene in p_vales-vector
FDR_values <- p.adjust(p_values, method = "BH", n = length(p_values)) #calculate FDR with benjamini-hochberg (BH)
statistics_values <- cbind(e_foldchange_mean_over_cell_lines,p_values, FDR_values) #combine mean, p_values and FDR in one matrix


cb1 <-if(FDR_values<0.05) col = "forestgreen"
cb2 <-if (abs(e_foldchange_mean_over_cell_lines)>1 ) col = "cyan"
cb3 <-if (FDR_values < 0.05 & abs(e_foldchange_mean_over_cell_lines) < 1) col = "firebrick"
cb <- cbind("FDR_values < 0.05" = cb1, "e_foldchange > 1 " = cb2, "FDR_values < 0.05 and e_foldchange > 1" =cb3)
plot(e_foldchange_mean_over_cell_lines, -log10(FDR_values), pch=20, main="Volcano plot", xlim=c(-1,2), col= cb)