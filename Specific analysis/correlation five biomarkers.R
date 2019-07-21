scatter plots GI50/fold change with our five biomarkers

cor_plot_biomarkers_final <- sapply (biomarkers_final, function(x){
  e_foldchange[x,]
})
cor_plot_biomarkers_final <- cbind("LogGI50" = LogGI50_59_celllines, annotation_cancertype_df, cor_plot_biomarkers_final)

cor_plot_biomarkers_final$Cancertypes <- as.factor(cor_plot_biomarkers_final$Cancertypes)
cor_plot_biomarkers_final <- as.data.frame(cor_plot_biomarkers_final)


plot <- lapply(colnames(cor_plot_biomarkers_final[3:7]), function(gene){
  ggplot_title <- paste("Foldchange of" , gene)
  ggplot(cor_plot_biomarkers_final, aes(x = LogGI50, y = cor_plot_biomarkers_final[ ,gene], colour = Cancertypes)) +
    geom_point()+
    scale_colour_manual(values = color_palette_cancertype)+
    theme(axis.title.y = element_blank())+
    
    ggtitle(paste0 (ggplot_title, "\n"))
})

# load ggpubr
library(ggpubr)
# draw all five scatter plots next to each other
ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[5]], 
          ncol = 2, nrow = 3, common.legend = TRUE, legend = "top")

```