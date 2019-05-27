#Broad analysis

#calculate fold change
fold_changes <- NCI_TPW_gep_treated - NCI_TPW_gep_untreated
fold_change_numbers <- fold_changes #in fold_changes annotation is added, original is needed for calculations
fold_changes <- as.data.frame(fold_changes)


#1. annotation of samples by drug (annotation_of_celllines_per_drug) 
#contains samples treated with one drug in a column (15 drugs = columns)
annotation_of_celllines_per_drug <- matrix(, nrow = 60, ncol = 0)

#to cbind all celllines of one drug in the annotation matrix, they need to have the same number of rows
#all columns are filled to 60 rows with "NA"s 

for (i in 1:nrow(drug_annotation)){
  columns_of_one_drug <- c(grep (drug_annotation$Drug[i], colnames(fold_changes), value = TRUE))
  #fill with NAs to 60
  columns_of_one_drug <- c(columns_of_one_drug, rep(NA, 60-length(columns_of_one_drug)))
  annotation_of_celllines_per_drug <- cbind (annotation_of_celllines_per_drug, columns_of_one_drug)
}
rm(i, columns_of_one_drug)

#name and sort annotation alphabetically
colnames(annotation_of_celllines_per_drug) <- drug_annotation$Drug
annotation_of_celllines_per_drug <- annotation_of_celllines_per_drug[, order(colnames(annotation_of_celllines_per_drug))]


#2. annotation of samples by cellline
#61 columns = cell lines
annotation_sorted_by_cell_lines <- matrix(, nrow = 15, ncol = 0)

#Problem: name of cellline SK-MEL_2 is part of cellline SK-MEL_28
#Solution: rename SK-MEL-2 to SK-MEL-2_ (first define it as new factor level)
levels(cellline_annotation$Cell_Line_Name) <- c(levels(cellline_annotation$Cell_Line_Name), "SK-MEL-2_")
cellline_annotation[33, 1] <- "SK-MEL-2_"

for (i in 1:61){
  columns_of_one_cell_line <- c(grep (cellline_annotation$Cell_Line_Name[i], colnames(fold_changes), value = TRUE))
  #fill with NAs to 15
  columns_of_one_cell_line <- c(columns_of_one_cell_line, rep(NA, 15-length(columns_of_one_cell_line)))
  annotation_sorted_by_cell_lines <- cbind (annotation_sorted_by_cell_lines, columns_of_one_cell_line)
}
rm(i, columns_of_one_cell_line)
colnames(annotation_sorted_by_cell_lines) <- cellline_annotation$Cell_Line_Name



#add drug, cellline and cancertype annotation in FC
drugs <- sapply(colnames(fold_changes), function(x){
  #arr.ind array indices are returned(row and col), [2] because we want to know in which col
  colnames(annotation_of_celllines_per_drug)[which(x == annotation_of_celllines_per_drug, arr.ind = TRUE)[2]]
})

cellline <- sapply(colnames(fold_changes), function(x){
  colnames(annotation_sorted_by_cell_lines)[which(x == annotation_sorted_by_cell_lines, arr.ind = TRUE)[2]]
})

fold_changes <- rbind("Drug" = drugs, "Cellline" = cellline, fold_changes)

cancertype <- sapply(fold_changes[2, ], function(x){ #2nd row contains cellline annotation of samples
  cellline_annotation$Cancer_type[cellline_annotation$Cell_Line_Name == x]
})
cancertype <- as.vector(unlist(cancertype))

fold_changes <- rbind("Cancertype" = cancertype, fold_changes)
rm(drugs, cellline, cancertype)



#Coloring: 1. drug (color_vector_all_drugs), 2. cancertype (color_vector_cancertype)
#1. define a color palette with 15 chosen colors
color_palette_drug <- c("aquamarine", "brown", "forestgreen", "slategrey", "chartreuse", "darkgoldenrod1", "cadetblue","purple", "firebrick1", "deepskyblue", "gold", "violetred4", "deeppink", "plum2", "blue" )
names(color_palette_drug) <- colnames(annotation_of_celllines_per_drug)
  
#create vector containing a color name for each sample according to drug
color_vector_drug <- sapply(colnames(fold_changes), function(x){
  unname(color_palette_drug[fold_changes[2, x]]) #unname: only color without drug is stored, 2nd row of FC contains drug
})

#2. define a color palette with 9 chosen colors
color_palette_cancertype <- c("aquamarine", "brown", "forestgreen", "chartreuse", "darkgoldenrod1", "cadetblue","purple", "firebrick1", "deepskyblue")
names(color_palette_cancertype) <- levels(cellline_annotation$Cancer_type)

#create vector containing a color name for each sample according to cancertype
color_vector_cancertype <- sapply(colnames(fold_changes), function(x){
  unname(color_palette_cancertype[fold_changes[1, x]]) #1st row of FC contains cancertype 
})



#Boxplot
#par makes spaces outside the plot larger, xaxt: removes labels on x-axis
#title() used to move xlab nearer to the axis
par(oma = c(1, 1, 1, 8))
boxplot(NCI_TPW_gep_untreated, xaxt = "n", ylab = "Gene expression profile", vertical =  T, 
        main = "Boxplot: gene expression profile of untreated NCI60 celllines", 
        boxcol = color_vector_drug)
title(xlab = "Celllines treated with different drugs", line = 1.0)
legend(x = 860, y = 14.5, legend = names(color_palette_drug), col = color_palette_drug, pch = 19, xpd = "TRUE")



#PCA
pca <- prcomp(fold_change_numbers)

#color PCA according to drug 
par(oma = c(1, 1, 1, 8))
#PC1 and PC2
plot(pca$rotation[,1], pca$rotation[,2], col = color_vector_drug, pch = 19, xlab = "PC1", ylab = "PC2", main = "PCA with FC of all samples")
legend(x = 0.08, y = 0.143, legend = names(color_palette_drug), col = color_palette_drug, pch = 19, xpd = "TRUE")
#PC2 and PC3
plot(pca$rotation[,2], pca$rotation[,3], col = color_vector_drug, pch = 19, xlab = "PC2", ylab = "PC3", main = "PCA with FC of all samples")
legend(x = 0.16, y = 0.096, legend = names(color_palette_drug), col = color_palette_drug, pch = 19, xpd = "TRUE")

#Color PCA according to cancertype
par(oma = c(1, 1, 1, 10))
#PC3 and PC4
plot(pca$rotation[,3], pca$rotation[,4], col = color_vector_cancertype, pch = 19, xlab = "PC3", ylab = "PC4", main = "PCA with FC of all samples")
legend(x = 0.11, y = 0.06, legend = names(color_palette_cancertype), col = color_palette_cancertype, pch = 19, xpd = "TRUE")

rm(pca)


#Density plot of all celllines and drugs, in black treated, red untreated
plot(density(NCI_TPW_gep_untreated), "Density plot of gene expression")
lines(density(NCI_TPW_gep_treated), col = "red")
legend("topright", legend = c("untreated", "treated"), col = c("black", "red"), pch = 15)



#Barplot of genes with highest mean FC over all samples
mean_FC <- apply(fold_change_numbers, 1, mean)

#sort: starting with highest mean FC, abs() makes all negative values positive
mean_FC <- sort(abs(mean_FC), decreasing = TRUE)
par(oma = c(10, 1, 1, 1))
barplot(mean_FC[1:20], main = "Genes with highest mean FC",  ylab = "mean FC", las = 2) #las = 2: vertical x labels

#alternative: calculating the mean FC over positive FC values
mean_FC <- apply(fold_change_numbers, 1, function(x){
  mean(abs(x)) 
})
par(oma = c(10, 1, 1, 1))
barplot(sort(mean_FC, decreasing = TRUE) [1:20], main = "Genes with highest mean FC",  ylab = "mean FC", las = 2)
#--> Result: nearly the same genes (little change in order) with directly calculating the mean or calculating 
#the mean of positive FC values



#TO DO:
"#Biomarkers for each cancer type seperately
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
"

#Boxplot showed batch effect --> we have to normalize the data
# each sample should have mean 0 and sd 1
FC_normalized <- apply(fold_change_numbers, 2, function(x){
  (x - mean(x)) / sd(x)
})

boxplot(FC_normalized, xaxt = "n", ylab = "Fold changes", vertical =  T, 
        main = "Boxplot: fold change distribution after normalization")

