#892
celllines_final <- unlist(sapply(celllines_overlap, function(x){
  x[which (x == LogGI50_top10)]
}))

### Fold change of biomarkers 
**Barplot** of foldchange of our five biomarkers in our sensitive cell lines
```{r}
df_barplot_results <- as.data.frame(sapply ((biomarkers_final), function(x) {
  e_foldchange[x, labeled_celllines_PCA]
}))

df_barplot_results <- as.data.frame(cbind(cellline = labeled_celllines_PCA, annotation_cancertype[labeled_celllines_PCA], df_barplot_results))

plot <- lapply(colnames(df_barplot_results)[3:7], function(gene){
  ggplot_title <- paste("Foldchange of" , gene)
  ggplot(data=df_barplot_results, aes(x=cellline, y= df_barplot_results[ ,gene], fill=cancertype)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values=c("aquamarine", "forestgreen", "chartreuse", 
                               "cadetblue","purple", 
                               "firebrick1", "deepskyblue"))+
    ggtitle(paste0 (ggplot_title, "\n"))+
    theme(axis.title.x = element_blank())+
    coord_flip()
})


library(ggpubr)
ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[5]] + rremove("x.text"), 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2,
          common.legend = TRUE)