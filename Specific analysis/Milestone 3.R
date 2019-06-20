#Barplot of genes with highest mean in FC over Erlotinib
genes_FC_erlotinb <- apply(e_foldchange_normalized, 1, mean)
genes_FC_erlotinb <- sort(abs(genes_FC_erlotinb), decreasing = TRUE)
par(oma =c(10,1,1,1))
barplot(genes_FC_erlotinb [1:20], main = "Genes with highest FC after erloinib treatment", ylab = "FC", las = 2)

genes_FC_erlotinb[1:20]

#Find the row number with the name
which(rownames(e_foldchange_normalized) == "CHAC1") 

# vector which only include celllines which were used in e_foldchange_normalized
NegLogGI50_59_celllines <- NegLogGI50 ["erlotinib", -c(8,29)]
#@Anna: should we change the name from NegLogGI50_59_celllines_neg to simply LogGI50?
NegLogGI50_59_celllines_neg <- NegLogGI50_59_celllines * (-1)

#Coloring according to cancertype
e_color_cancertype <- color_vector_cancertype[grep("erlotinib", names(color_vector_cancertype), value = TRUE)]

#Scatter plot
par(oma = c(1,1,1,10), xpd = "TRUE") #size of outer margins: bottom, top, left, right
plot(NegLogGI50_59_celllines_neg, e_foldchange_normalized ["EGR1",], 
     col = e_color_cancertype, 
     pch = 19, 
     xlab = "logGI50", 
     ylab = "EGR1 Expression (log2, relative to control)",        
     main = "Erlotinib 24h")

legend(x = -3.8, y = 2, 
       legend = names(color_palette_cancertype), 
       col = color_palette_cancertype, 
       pch = 19)

#label only points in the left bottom quarter
labeled_celllines <- names(NegLogGI50_59_celllines_neg)[NegLogGI50_59_celllines_neg < - 5.5
                                                        & e_foldchange_normalized["EGR1", ] < 0]

text(NegLogGI50_59_celllines_neg[labeled_celllines], e_foldchange_normalized ["EGR1", labeled_celllines], 
     labels = labeled_celllines,
     cex = 0.7,
     pos = 3) #position of text at the top of the point


#Pearson correlation
res <- cor.test(NegLogGI50_59_celllines_neg, e_foldchange_normalized ["EGR1",], 
     method = "pearson")
res

lm()




# Expressionsdaten von EGFR Expressione relativ to controle 
#--> Einfluss der therapie auf die Erpression des Her 1 genes

plot(NegLogGI50_59_celllines_neg, e_foldchange_normalized ["EGFR",], 
      col = e_color_cancertype, 
      pch = 19, 
      xlab = "logGI50", 
      ylab = "EGFR Expression (log2, relative to control)",        
      main = "Erlotinib 24h")

legend(x = -3.8, y = 2, 
      legend = names(color_palette_cancertype), 
      col = color_palette_cancertype, 
      pch = 19)

#label only points which have a dicreased Her 1 expression after erlotinib treatment
labeled_celllines <- names(NegLogGI50_59_celllines_neg)[ e_foldchange_normalized["EGFR", ] < 0]

text(NegLogGI50_59_celllines_neg[labeled_celllines], e_foldchange_normalized ["EGFR", labeled_celllines], 
      labels = labeled_celllines,
      cex = 0.7,
      pos = 3) #position of text at the top of the point




   
# Her 1 gene expression in the untreated celllines gegen die GI50     
plot(NegLogGI50_59_celllines_neg, e_untreated ["EGFR",], 
      col = e_color_cancertype, 
      pch = 19, 
     xlab = "logGI50", 
     ylab = "EGFR Expression (untreated)",        
     main = "Expression of the epidermal growth factor receptor (Her 1)")
   
legend(x = -3.8, y = 11, 
    legend = names(color_palette_cancertype), 
    col = color_palette_cancertype, 
    pch = 19)

text(NegLogGI50_59_celllines_neg, e_untreated ["EGFR",], labels  = NegLogGI50_59_cellline_names,cex= 0.7, pos=3)   


#ERROR      
# lable celline which belong to cluster 3 in the heatmap in milestone 4
labeled_cellline_clusters <- names(NegLogGI50_59_celllines_neg)[ cellline_clusters [,1]  > 2 ]
text(NegLogGI50_59_celllines_neg[labeled_cellline_clusters], e_untreated[ "EGFR",],
    labels = labeled_cellline_clusters,
    cex = 0.7,
    pos = 3)  
   



   
'to do:
linear regression'