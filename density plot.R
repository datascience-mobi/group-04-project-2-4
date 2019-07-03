biomarkers <- rownames(statistics_values)[which(abs(statistics_values[, 1]) > 0.25)]
biomarkers_MA_up <- rownames(statistics_values)[which(statistics_values[, 1] > 0.25)]
biomarkers_MA_down <- rownames(statistics_values)[which(statistics_values[, 1] < -0.25)]

#Density plot with these genes (untreated vs. treated_upregulated)
plot(density(e_untreated[biomarkers_MA_up, ]), "Density plot of gene expression of upregulated biomarkers", cex = 1.2, col = "green", ylim= c(0.0, 0.5))

lines(density(e_treated[biomarkers_MA_up, ]), col = "red")
legend("topright", legend = c("untreated", "treated"), col = c("green", "red"), pch = 15)

#Density plot with these genes (untreated vs. treated_downregulated)
plot(density(e_untreated[biomarkers_MA_down, ]), "Density plot of gene expression of downregulated biomarkers", cex = 1.2, col = "green", ylim= c(0.0, 0.7))

lines(density(e_treated[biomarkers_MA_down, ]), col = "red")
legend("topright", legend = c("untreated", "treated"), col = c("green", "red"), pch = 15)