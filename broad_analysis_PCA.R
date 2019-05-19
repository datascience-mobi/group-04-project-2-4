#Broad analysis
#install package plyr
library(plyr)

#calculate fold change due to drug treatment
fold_changes <- NCI_TPW_gep_treated - NCI_TPW_gep_untreated
fold_changes <- as.data.frame(fold_changes)

#matrix annotation_of_celllines_per_drug contains all names 
#of celllines treated with each drug in a column (15 drugs, each in a seperated column containing all celllines treated with that drug)
annotation_of_celllines_per_drug <- matrix(, nrow = 60, ncol = 0)

#to cbind all celllines of one drug in the annotation matrix, they need to have the same number of rows
#all columns are filled to 60 rows with "NA"s and before that the actual number of rows gets stored in length_without_NAs
length_without_NAs <- c()

for (i in 1:15){
  columns_of_one_drug <- c(grep (drug_annotation$Drug[i], colnames(fold_changes), value = TRUE))
  length_without_NAs <- c(length_without_NAs, length(columns_of_one_drug))
  columns_of_one_drug <- c(columns_of_one_drug, rep(NA, 60-length(columns_of_one_drug)))
  annotation_of_celllines_per_drug <- cbind (annotation_of_celllines_per_drug, columns_of_one_drug)
}

#sort annotation alphabetically
colnames(annotation_of_celllines_per_drug) <- drug_annotation$Drug
annotation_of_celllines_per_drug <- annotation_of_celllines_per_drug [ , c(8, 1, 15, 7, 4, 14, 3, 5, 11, 6, 2, 12, 10, 13, 9)]
length_without_NAs <- length_without_NAs[c(8, 1, 15, 7, 4, 14, 3, 5, 11, 6, 2, 12, 10, 13, 9)]

#define a color palette with 15 chosen colors
color_palette <- c("aquamarine", "brown", "forestgreen", "slategrey", "chartreuse", "darkgoldenrod1", "cadetblue","purple", "firebrick1", "deepskyblue", "gold", "violetred4", "deeppink", "plum2", "blue" )

#define vector with 15 colors each assigned to one drug (needed for legend)
colors <- cbind(color_palette, colnames(annotation_of_celllines_per_drug))
colnames(colors) <- c("Color", "Drug")

#create vector containing each color name so often as the treated celllines with one drug
color_vector_all_drugs <- c()
for (i in 1:15){
  for (j in 1:length_without_NAs[i] )
  color_vector_all_drugs <- c(color_vector_all_drugs, colors [i, 1])
}

#Boxplot
#par makes spaces outside the plot larger, xaxt: removes labels on x-axis
#title() used to move xlab nearer to the axis
par(oma = c(1, 1, 1, 8))
boxplot(NCI_TPW_gep_untreated, xaxt = "n", ylab = "Gene expression profile",vertical =  T, main = "Boxplot: gene expression profile of untreated NCI60 celllines", boxcol = color_vector_all_drugs)
title(xlab = "Celllines treated with different drugs", line = 1.0)
legend(x = 860, y = 14.5, legend = colors[, 2], col = colors[, 1], pch = 19, xpd = "TRUE")

#PCA
pca <- prcomp(fold_changes)

#color all celllines in PCA according to drug treatment
par(oma = c(1, 1, 1, 8))
#PC1 and PC2
plot(pca$rotation[,1], pca$rotation[,2], col = color_vector_all_drugs, pch = 19, xlab = "PC1", ylab = "PC2", main = "PCA with FC of all celllines treated with different drugs")
legend(x = 0.08, y = 0.143, legend = colors[, 2], col = colors[, 1], pch = 19, xpd = "TRUE")
#PC2 and PC3
plot(pca$rotation[,2], pca$rotation[,3], col = color_vector_all_drugs, pch = 19, xlab = "PC2", ylab = "PC3", main = "PCA with FC of all celllines treated with different drugs")
legend(x = 0.16, y = 0.096, legend = colors[, 2], col = colors[, 1], pch = 19, xpd = "TRUE")




#Coloring according to cancer type (9 types, see cellline_annotation$Cancer_Type)

#annotation of celllines (in columns 61 cell lines) treated with drugs (in rows 15 drugs)
annotation_sorted_by_cell_lines <- matrix(, nrow = 15, ncol = 0)

#Problem: name of cellline SK-MEL_2 is part of cellline SK-MEL_28
#Solution: rename SK-MEL_2 to SK-MEL_2_ (first define it as new factor level)
levels(cellline_annotation$Cell_Line_Name) <- c(levels(cellline_annotation$Cell_Line_Name), "SK-MEL_2_")
cellline_annotation[33, 1] <- "SK-MEL_2_"

#problem: grep does not find SK-MEL_2_ in colnames --> creates column of NAs instead of values in SK-MEL_2 column
length_without_NAs <- c()
for (i in 1:61){
  columns_of_one_cell_line <- c(grep (cellline_annotation$Cell_Line_Name[i], colnames(fold_changes), value = TRUE))
  length_without_NAs <- c(length_without_NAs, length(columns_of_one_cell_line))
  columns_of_one_cell_line <- c(columns_of_one_cell_line, rep(NA, 15-length(columns_of_one_cell_line)))
  annotation_sorted_by_cell_lines <- cbind (annotation_sorted_by_cell_lines, columns_of_one_cell_line)
}
colnames(annotation_sorted_by_cell_lines) <- cellline_annotation$Cell_Line_Name

#define vector with 9 colors each assigned to one cancer type
cancertypes <- c()
for (i in 1:61){
  cancertypes <- union(cancertypes, cellline_annotation$Cancer_type[i])
}
colors_cancertype <- cbind(color_palette[1:9], cancertypes)
colnames(colors) <- c("Color", "Cancertypes")

#create vector containing color names according to cancer type 
color_vector_cancertype <- c()
for (i in 1:61){ #for each cell line
  for (j in 1:9){ #for each cancer type
    if (cellline_annotation$Cancer_type[i] == cancertypes[j]){
      for (k in 1:length_without_NAs[i]){ #add color so often as stored in legth (number of drugs tested on this cellline)
        color_vector_cancertype <- c(color_vector_cancertype, colors_cancertype [j, 1])
      }
    }
  }
}

#Color PCA according to cancertype
pca <- prcomp(fold_changes)
par(oma = c(1, 1, 1, 10))
plot(pca$rotation[,3], pca$rotation[,4], col = color_vector_cancertype, pch = 19, xlab = "PC3", ylab = "PC4", main = "PCA with FC of all celllines treated with different drugs")
legend(x = 0.11, y = 0.06, legend = colors_cancertype[, 2], col = colors_cancertype[, 1], pch = 19, xpd = "TRUE")



#Density plot of all celllines and drugs, in black treated, red untreated
plot(density(NCI_TPW_gep_untreated), "Density plot of gene expression")
lines(density(NCI_TPW_gep_treated), col = "red")
legend("topright", legend = c("untreated", "treated"), col = c("black", "red"), pch = 15)



#Find biomarker: 100 genes that have highest fold change for each cellline
genes_highest_foldchange_100 <- data.frame(matrix(nrow = 100, ncol = 819)) #100 rows, for most upregulated genes, 819 columns for celllines
gene_names_highest_foldchange_100 <- data.frame(matrix(nrow = 100, ncol = 819))

for (i in 1:819){ #for all celllines
  genes_highest_foldchange_100[, i] <- fold_changes[1:100, i]
  gene_names_highest_foldchange_100[, i] <- rownames(fold_changes)[1:100]
  for(j in 101:13299){ #for all genes
    min <- min(genes_highest_foldchange_100[,i])
    if (fold_changes[j, i] > min){
      gene_names_highest_foldchange_100[genes_highest_foldchange_100[ , i] == min , i] <- rownames(fold_changes)[j]
      genes_highest_foldchange_100[genes_highest_foldchange_100[ , i] == min , i] <- fold_changes[j, i]
    }
  }
}

counts_gene_names_highest_FC_100 <- c()
for (i in 1:13299){#for each gene
  counter = 0
  for(j in 1:100){ #for each row in gene_names_highest_foldchange_100
    for (k in 1:819){#for each column in gene_names_highest_foldchange_100
      if (rownames(fold_changes)[i] == gene_names_highest_foldchange_100[j, k] ){
        counter = counter + 1
      }
    }
  }
  counts_gene_names_highest_FC_100[i] <- counter
}
rm(counter, i, j, k)

#Name the counts with the genenames, sort them by starting with highest and delete all genes that did not occur
names(counts_gene_names_highest_FC_100) <- rownames(fold_changes)
counts_gene_names_highest_FC_100 <- sort(counts_gene_names_highest_FC_100, decreasing = TRUE)
counts_gene_names_highest_FC_100 <- counts_gene_names_highest_FC_100[-which(counts_gene_names_highest_FC_100 == 0)]

#Barplot of genes showing how often their fold change is in the top 100
barplot(counts_gene_names_highest_FC_100[1:20], main = "Most upregulated genes", xlab = "Genes", xaxt = "n", ylab = "Counts (top100 FC)")
text(x = barplot[,1], y = -50, labels = names(counts_gene_names_highest_FC_100)[1:20], xpd = TRUE, cex = 1, srt = 60)


#Biomarkers for each cancer type seperately
mean_FC_cancertype <- matrix(nrow = 13299, ncol = 9)
cancertypes_sorted_to_FC_columns <- data.frame(matrix(nrow = 105, ncol = 9))

for (j in 1:9){ #for each cancer type
   for (i in 1:61){ #for each cell line
      if (cellline_annotation$Cancer_type[i] == cancertypes[j]){
        cancertypes_sorted_to_FC_columns[, j] <- rbind.data.frame(annotation_sorted_by_cell_lines[ , i])
    }
  }
}

for (i in 1:9){ #for each cancer type
  for (j in 1:13299){ #for each gene
    mean_FC_cancertype[j, i] <- mean(fold_changes[j, ])
  }
}
