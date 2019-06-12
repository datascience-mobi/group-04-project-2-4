#Barplot of genes with highest mean in FC over Erlotinib
genes_FC_erlotinb <- apply(e_foldchange_normalized, 1, mean)
genes_FC_erlotinb <- sort(abs(genes_FC_erlotinb), decreasing = TRUE)
par(oma =c(10,1,1,1))
barplot(genes_FC_erlotinb [1:20], main = "Genes with highest FC after erloinib treatment", ylab = "FC", las = 2)

genes_FC_erlotinb[1:20]

#Find the row number with the name
which(rownames(e_foldchange_normalized) == "CHAC1") 

# vector which only include celllines which were used in e_foldchange_normalized
NegLogGI50_59_celllines <- NegLogGI50 [6, -c(8,29)]
#should we change the name from NegLogGI50_59_celllines_neg to LogGI50?
NegLogGI50_59_celllines_neg <- NegLogGI50_59_celllines*(-1)

#plot(NegLogGI50_59_celllines_neg, e_treated [3308,], col = color_vector_cancertype, pch = 19, xlab = "logGI50", ylab = "EGR1 Expression (log2, relative to controle",        main = "Erlotinib 24h")
#legend(x = 0.11, y = 0.06, legend = names(color_palette_cancertype), col = color_palette_cancertype, pch = 19,      xpd = "TRUE")
#foldchange muss genommen werden, da sonst keine pos uns negativen werte herauskommen. 

par(oma = c(1,1,1,10), xpd = "TRUE") #size of outer margins: bottom, top, left, right
plot(NegLogGI50_59_celllines_neg, e_foldchange_normalized ["EGR1",], 
     col = color_vector_cancertype, 
     pch = 19, 
     xlab = "logGI50", 
     ylab = "EGR1 Expression (log2, relative to control)",        
     main = "Erlotinib 24h")

NegLogGI50_59_cellline_names <- colnames(NegLogGI50 [,-c (8,29)])
 text(NegLogGI50_59_celllines_neg, e_foldchange_normalized ["EGR1",], labels=NegLogGI50_59_cellline_names, pos =4)
# text(e_foldchange_normalized ["EGR1",], label=ifelse (e_foldchange_normalized ["EGR1",4] < -4,
#as.character (NegLogGI50_59_cellline_names),''), pos =4)
 #läuft durch macht abe nichts 
 
legend(x = -3.8, y = 2, 
       legend = names(color_palette_cancertype), 
       col = color_palette_cancertype, 
       pch = 19)

'to do:
nur spezifische Cellinien benennen
pearson correlation
linear regression'