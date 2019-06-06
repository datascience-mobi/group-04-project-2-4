#Barplot of genes with highest mean in FC over Erlotinib
genes_FC_erlotinb <- apply(e_foldchange_z_transformed, 1, mean)
genes_FC_erlotinb <- sort(abs(genes_FC_erlotinb), decreasing = TRUE)
par(oma =c(10,1,1,1))
barplot(genes_FC_erlotinb [1:20], main = "Genes with highest FC after erloinib treatment", ylab = "FC", las = 2)

genes_FC_erlotinb[1:20]

#Finde the row number with the name
which(rownames(e_foldchange_z_transformed) == "CHAC1") 

# vector whichonly include celllines which were used in e_foldchange_z_transformed
NegLogGI50_59_celllines <- NegLogGI50 [, -c(8,29)]
NegLogGI50_59_celllines <- NegLogGI50_59_celllines [-c(1:5, 7:15),]
NegLogGI50_59_celllines_neg <- NegLogGI50_59_celllines*(-1)

#plot(NegLogGI50_59_celllines_neg, e_treated [3308,], col = color_vector_cancertype, pch = 19, xlab = "logGI50", ylab = "EGR1 Expression (log2, relative to controle",        main = "Erlotinib 24h")
#legend(x = 0.11, y = 0.06, legend = names(color_palette_cancertype), col = color_palette_cancertype, pch = 19,      xpd = "TRUE")
#foldchange muss genommen werden, da sonst keine pos uns negativen werte herauskommen. 

plot(NegLogGI50_59_celllines_neg, e_foldchange_z_transformed [3308,], col = color_vector_cancertype, pch = 19, xlab = "logGI50", ylab = "EGR1 Expression (log2, relative to controle",        main = "Erlotinib 24h")
legend(x = -3, y = 7, legend = names(color_palette_cancertype), col = color_palette_cancertype, pch = 19)

'to do:
legende hinzufügen
pearson correlation
linear regression'