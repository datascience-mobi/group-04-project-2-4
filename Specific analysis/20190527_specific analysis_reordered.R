#new matrix only with samples/columns treated with erlotinib  (e=erlotinib)
e_treated <- NCI_TPW_gep_treated[,grep ("erlotinib", colnames(NCI_TPW_gep_treated))]
e_untreated <- NCI_TPW_gep_untreated[,grep ("erlotinib", colnames(NCI_TPW_gep_untreated))]
e_foldchange <- e_treated - e_untreated

#e_foldchange_normalized: z-Transformation to get mean=0 and sd=1
e_foldchange_normalized <- apply(e_foldchange, 2, function(x){
  (x - mean(x)) / sd(x)
})


# table of 15 cell lines with highest variance to show most regulated cell lines
#select 15 cell lines with highest variance (greater than 75% quantile, sorted by decreasing value)
var_cell_line <- apply(e_foldchange, 2, var)
cell_line_var_greater_75quantile <- sort(var_cell_line [which (var_cell_line > quantile(var_cell_line,0.75))], decreasing = TRUE)
rm(var_cell_line)

#create a table containing the name and the variance of the 15 cell lines with highest variance
table_cell_lines_var_top15 <- cbind(names(cell_line_var_greater_75quantile), cell_line_var_greater_75quantile)
rownames(table_cell_lines_var_top15) <- c(1:nrow(table_cell_lines_var_top15))
colnames(table_cell_lines_var_top15) <- c("cell line", "variance")


#PCA to show most regulated cell lines
#PCA with transformed matrix (each point represents a sample):
e_foldchange_matrix_transformed <- t(e_foldchange_normalized)
pca <- prcomp(e_foldchange_matrix_transformed)
plot(pca$rotation[,1], pca$rotation[,2])
text(pca$rotation, labels = rownames(e_foldchange_normalized), cex = 0.4, pos = 3)

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


#Volcano plot
#mean of gene expression of each gene over all cell lines
e_foldchange_mean_over_cell_lines <- rowMeans(e_foldchange) #equal to e_treated_mean_over_cell_lines - e_untreated_mean_over_cell_line

#determine the p-value for a paired two-sample t-test 
p_values <- sapply(rownames(e_treated), function(x) {
  t.test(e_treated[x,], e_untreated[x,],paired= T)$p.value}) # perform t-test and save p-values of each gene in p_vales-vector
FDR_values <- p.adjust(p_values, method = "BH", n = length(p_values))#calculate FDR with benjamini-hochberg (BH)


#table of results 
statistics_values <- cbind(e_foldchange_mean_over_cell_lines, FDR_values)
#coloring with package enhanced volcano
#install package EnhancedVolcano (needs ggplot2, ggrepel)
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

EnhancedVolcano(statistics_values, 
                lab = rownames(statistics_values),
                x = "e_foldchange_mean_over_cell_lines", #colname of FC values in this table (statistics_values)
                y = "FDR_values", #colname of FDR (statistics_values)
                title = "Volcano plot of all genes",
                pCutoff = 10e-15, #threshold for coloring significant ones
                FCcutoff = 1, #threshold for coloring high FC
                transcriptPointSize = 2,
                transcriptLabSize = 4.0)

#further questions: what "do" the red genes? are they involved in certain pathways?
#maybe we could color the genes in the volcano plot according to pathways (after milestone 4)

#save the "red" genes seen in the volcano plot in a vector for further analysis
volcano_genes <- rownames(statistics_values)[which(abs(statistics_values[, 1]) > 1 
                                                   & statistics_values[, 2] < 10e-15)]

#Density plot with these genes (untreated vs. treated)
plot(density(e_treated[volcano_genes, ]), "Density plot of gene expression", col = "red")
lines(density(e_untreated[volcano_genes, ]), col = "black")
legend("topright", legend = c("untreated", "treated"), col = c("black", "red"), pch = 15)


#MA-Plot
M <- e_foldchange # M= log2(treated) - log2 (untreated)
A <- 1/2*(e_treated+ e_untreated) # average log2-expression value A = 1/2 (log2(treated)+log2(untreated))
MA <- cbind("M"= rowMeans(M), "A" = rowMeans(A), FDR_values)
MA <- as.data.frame(MA)
ggplot(data=MA, aes(x=A, y=M), fdr=0.05, fc=1.5, genenames=NULL, detection_call=NULL, size=NULL) #plot with ggplot

install.packages("BiocManager") #plot with limma
BiocManager::install("limma")
limma::plotMA()