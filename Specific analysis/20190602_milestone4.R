
#assign genes to pathways (as in supplementary material)
MAPK <- c("AURKA", "AURKB", "SON", "CENPA", "KIF11", "DUSP6", "TPX2", "EGR1")
AKT_PI3K <- c("CFLAR", "XIAP", "BIRC5", "BIRC3", "GADD45A", "MCL1", "BCL2", "BCL2L1", "HIF1A", "AP1", "AR", "STAT3", "IL6", "VEGF")
ER_Stress_Apoptotic <- c("ATF3", "DDIT3", "TNFRSF10B", "TRIB3", "BCL2L11", "BBC3", "CASP2", "GADD34")
pathways <- list(MAPK = MAPK, AKT_PI3K = AKT_PI3K, ER_Stress_Apoptotic = ER_Stress_Apoptotic) #create named list


heatmap_data <- sapply(1:length(pathways), function(x){ #for each pathway
  pathway <- sapply(colnames(e_foldchange_normalized), function(y){ #for each sample
    mean(e_foldchange_normalized[unlist(pathways[x]), y]) #calculate mean FC of genes of one pathway
  })
  return(pathway)
})
colnames(heatmap_data) <- c(names(pathways))

heatmap(heatmap_data, col = cm.colors(256), scale="column", margins=c(5,10))