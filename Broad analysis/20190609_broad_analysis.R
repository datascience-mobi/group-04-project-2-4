#Broad analysis

#Problem: name of cellline SK-MEL-2 is part of cellline SK-MEL_28
#Solution: rename SK-MEL-2 to SK-MEL-2_ (first define it as new factor level)
levels(cellline_annotation$Cell_Line_Name) <- c(levels(cellline_annotation$Cell_Line_Name), "SK-MEL-2_")
cellline_annotation[33, 1] <- "SK-MEL-2_"
#delete level SK-MEL-2 (otherwise we would have 62, instead of 61 levels)
cellline_annotation$Cell_Line_Name <- factor(as.character(cellline_annotation$Cell_Line_Name))



#calculate fold change
fold_changes <- NCI_TPW_gep_treated - NCI_TPW_gep_untreated
fold_changes <- as.data.frame(fold_changes)



#Create annotation: for each sample name(rows) the columns contain the drug, cellline and cancertype
#1. Drug
sample_drug <- as.data.frame(sapply(levels(drug_annotation$Drug), grepl, colnames(fold_changes), ignore.case = TRUE))
  #creates table with TRUE and FALSE for each sample and drug
rownames(sample_drug) <- colnames(fold_changes)
drugs <- as.vector(apply(sample_drug, 1, function(x){
  colnames(sample_drug[which(x)])
}))

#2. Cellline
sample_cellline <- as.data.frame(sapply(levels(cellline_annotation$Cell_Line_Name), grepl, colnames(fold_changes), ignore.case = TRUE)) 
  #creates table with TRUE and FALSE for each sample and cellline
rownames(sample_cellline) <- colnames(fold_changes)
cellline <- as.vector(unlist(apply(sample_cellline, 1, function(x){
  colnames(sample_cellline[which(x)])
})))

annotation <- cbind("Drug" = drugs, "Cellline" = cellline)
rownames(annotation) <- colnames(fold_changes)

#3. Cancertype
cancertype <- sapply(annotation[, 2], function(x){ #2nd column contains cellline annotation of samples
  cellline_annotation$Cancer_type[cellline_annotation$Cell_Line_Name == x]
})
cancertype <- as.vector(unlist(cancertype))

annotation <- cbind(annotation, "Cancertype" = cancertype)
rm(drugs, sample_drug, cellline, sample_cellline, cancertype)



#Coloring: 1. drug (color_vector_all_drugs), 2. cancertype (color_vector_cancertype)
#1. define a color palette with 15 chosen colors
color_palette_drug <- c("aquamarine", "brown", "forestgreen", "slategrey", "chartreuse", "darkgoldenrod1", "cadetblue","purple", "firebrick1", "deepskyblue", "gold", "violetred4", "deeppink", "plum2", "blue" )
names(color_palette_drug) <- levels(drug_annotation$Drug)
  
#create vector containing a color name for each sample according to drug
color_vector_drug <- sapply(rownames(annotation), function(x){
  unname(color_palette_drug[annotation[x, 1]]) #unname: only color without drug is stored, first column of annotation contains drug
})

#2. define a color palette with 9 chosen colors
color_palette_cancertype <- c("aquamarine", "brown", "forestgreen", "chartreuse", "darkgoldenrod1", "cadetblue","purple", "firebrick1", "deepskyblue")
names(color_palette_cancertype) <- levels(cellline_annotation$Cancer_type)

#create vector containing a color name for each sample according to cancertype
color_vector_cancertype <- sapply(rownames(annotation), function(x){
  unname(color_palette_cancertype[annotation[x, 3]]) #3rd columns of annotation contains cancertype 
})



#Boxplot
#par makes spaces outside the plot larger, xaxt: removes labels on x-axis
#title() used to move xlab nearer to the axis
par(oma = c(1, 1, 1, 8), xpd = "TRUE")
boxplot(NCI_TPW_gep_untreated, 
        xaxt = "n", 
        ylab = "Gene expression profile", 
        vertical =  T, 
        main = "Gene expression profile of untreated NCI60 celllines")
title(xlab = "Samples", line = 1.0)
  #batch effect was seen --> corresponding to drugs?

#Color plot according to drugs
boxplot(NCI_TPW_gep_untreated, 
        xaxt = "n", 
        ylab = "Gene expression profile", 
        vertical =  T, 
        main = "Gene expression profile of untreated NCI60 celllines", 
        boxcol = color_vector_drug)
title(xlab = "Samples", line = 1.0)
legend(x = 860, 
       y = 14.5, 
       legend = names(color_palette_drug), 
       col = color_palette_drug, 
       pch = 19)

#Normalization of data is necessary
#each sample should have mean 0 and sd 1
untreated_normalized <- apply(NCI_TPW_gep_untreated, 2, function(x){
  (x - mean(x)) / sd(x)
})
FC_normalized <- apply(fold_changes, 2, function(x){
  (x - mean(x)) / sd(x)
})

#boxplot of normalized untreated values
par(oma = c(1, 1, 1, 8), xpd = "TRUE")
boxplot(untreated_normalized, xaxt = "n", ylab = "Gene expression profile", vertical =  T, 
        main = "Normalized gene expression profile of untreated NCI60 celllines", 
        boxcol = color_vector_drug)
title(xlab = "Samples", line = 1.0)
legend(x = 860, y = 3.4, legend = names(color_palette_drug), col = color_palette_drug, pch = 19)



#Density plot of all celllines and drugs, in black treated, red untreated
plot(density(NCI_TPW_gep_untreated), "Density plot of gene expression")
lines(density(NCI_TPW_gep_treated), col = "red")
legend("topright", legend = c("untreated", "treated"), col = c("black", "red"), pch = 15)



#PCA
pca <- prcomp(FC_normalized)

#color PCA according to drug 
par(oma = c(1, 1, 1, 8), mfrow = c(2, 2)) #mfrow to create multiple plots
#PC1 and PC2
plot(pca$rotation[,1], 
     pca$rotation[,2], 
     col = color_vector_drug, 
     pch = 19, 
     xlab = "PC1", 
     ylab = "PC2")
#PC2 and PC3
plot(pca$rotation[,2], 
     pca$rotation[,3], 
     col = color_vector_drug, 
     pch = 19, xlab = "PC2", 
     ylab = "PC3")
#create legend on the right side
legend(x = 0.07, 
       y = 0.096, 
       legend = names(color_palette_drug), 
       col = color_palette_drug, 
       pch = 19, 
       xpd = "TRUE",
       cex = 0.9)
#Title: mtext = margin text, side = 3 (upside)
mtext("Coloring according to drug", side = 3, line = -2, outer = TRUE)

#Color PCA according to cancertype
#PC1 and PC2
plot(pca$rotation[,1], 
     pca$rotation[,2], 
     col = color_vector_cancertype, 
     pch = 19, xlab = "PC1", 
     ylab = "PC2")
#PC2 and PC3
plot(pca$rotation[,2], 
     pca$rotation[,3], 
     col = color_vector_cancertype, 
     pch = 19, xlab = "PC2", 
     ylab = "PC3")
legend(x = 0.07, 
       y = 0.096, 
       legend = names(color_palette_cancertype), 
       col = color_palette_cancertype, 
       pch = 19, 
       xpd = "TRUE",
       cex = 0.9)
mtext("Coloring according to cancertype", side = 3, line = -18, outer = TRUE)

rm(pca)



#Barplot of genes with highest mean FC over all samples
mean_FC <- apply(fold_changes, 1, mean)

#sort: starting with highest mean FC, abs() makes all negative values positive
mean_FC <- sort(abs(mean_FC), decreasing = TRUE)
par(oma = c(10, 1, 1, 1), mfrow = c(1,1))
barplot(mean_FC[1:20], main = "Genes with highest mean FC",  ylab = "mean FC", las = 2) #las = 2: vertical x labels

#alternative: calculating the mean FC over positive FC values
mean_FC <- apply(fold_changes, 1, function(x){
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




