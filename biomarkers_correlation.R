
par(mar = c(4, 4, 4, 0), xpd = "TRUE")

correlation_b <- c(e_foldchange[biomarkers_final [1],],e_foldchange[biomarkers_final [2],],e_foldchange[biomarkers_final [3],],e_foldchange[biomarkers_final [4],],e_foldchange[biomarkers_final [5],])
correlation_g <-rep(LogGI50_59_celllines,5)
correlation_c <- c(rep("green",59),rep("yellow",59),rep("red",59),rep("blue",59),rep("black",59))
correlation_bn <- c(rep("ATAD2",59),rep("E2F8",59),rep("ERCC6L",59),rep("MIR15A//DLEU2L//DLEU2",59),rep("MKI67",59))
correlation_biomarkers <- cbind (correlation_b, correlation_g, correlation_c, correlation_bn )

color_palette_biomarkers <- c("green","yellow","red" , "blue", "black" )
names(color_palette_biomarkers) <- biomarkers_final

plot(correlation_biomarkers[,"correlation_g"],correlation_biomarkers[,"correlation_b"], 
     col = correlation_c, 
     pch = 19,
     xlab = "logGI50", 
     ylab = "foldchange of biomarkers",        
     main = " Correlation between biomarkers and the GI50 values of the different celllines",
     cex.main = 1.2)

legend("topright", 
       legend = names( color_palette_biomarkers),
       col = color_palette_biomarkers,
       pch = 19,
       inset = c(-0.3, 0),
       xpd = TRUE)
#HALLO