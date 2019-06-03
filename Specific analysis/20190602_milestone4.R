
#assign genes to pathways (as in supplementary material)
MAPK <- c("AURKA", "AURKB", "SON", "CENPA", "KIF11", "DUSP6", "TPX2", "EGR1")
AKT_PI3K <- c("CFLAR", "XIAP", "BIRC5", "BIRC3", "GADD45A", "MCL1", "BCL2", "BCL2L1", "HIF1A", "AR", "STAT3", "IL6", "JUN", "FOS")
  #AP1 for AKT_PI3K excluded (genename does not exist in FC)
  #do they mean with AP1 these 5 genes? AP1AR, AP1B1, AP1G1, AP1G2, AP1M2, AP1S1, AP1S2
  #same problem for VEGF --> VEGFA, VEGFB, VEGFB
Cell_Cycle_Checkpoint <- c("CHEK1", "CHEK2", "TP53", "MDM2", "CDKN1A", "CCNE1", "CDK2", "CDC25A", "SMC1A", "CCNB1", "CDK1", "CDC25B", "CDC25C", "PLK1", "WEE1", "CCND1")
JAK_STAT <- c("SOCS1", "NMI", "BCL2L1", "CDKN1A", "MYC")
Immunity_Related <- c("BCL2L1", "XIAP", "BIRC5", "CFLAR", "IL6", "IL1B", "TNF")
  #PTGS does not exist
DNA_repair <- c("XRCC6", "PRKDC", "DCLRE1C", "XRCC4", "LIG4", "BRCA1", "RAD52", "BRCA2", "RAD54L", "ATM", "ATR", "PARP1")
  #NHEJ1 does not exist
DNA_damage <- c("TP53", "FAS", "BAX", "PMAIP1")
  #BBC3 & ARG do not exist
ER_Stress_Survival <- c("HSPA5", "HSP90B1", "XBP1", "P4HB", "ATF4", "GADD45A")
ER_Stress_Apoptotic_Response <- c("ATF3", "DDIT3", "TNFRSF10B", "TRIB3", "BCL2L11", "CASP2")
  #BBC3 & GADD34 do not exist
Apoptosis_extrinsic_activation <- c("BID", "FAS", "CASP8", "CASP9", "APAF1")
Apoptosis_intrinsic_activation <- c("BAX", "BAK1", "BID", "PMAIP1", "CASP6", "CASP9", "APAF1")
  #BBC3 does not exist
Autophagy_recycling_starvation <- c("MAP1LC3B", "ULK1", "ATG13", "BECN1", "SQSTM1", "ATG5", "ATG12", "CTSB")
Autophagy_toxic <- c("SQSTM1", "ATG7", "RIPK1", "CTSB")

#create named list of all pathways
pathways <- list(MAPK = MAPK, 
                 AKT_PI3K = AKT_PI3K, 
                 Cell_Cycle_Checkpoint = Cell_Cycle_Checkpoint, 
                 JAK_STAT = JAK_STAT,
                 Immunity_Related = Immunity_Related,
                 DNA_repair = DNA_repair,
                 DNA_damage = DNA_damage,
                 ER_Stress_Survival = ER_Stress_Survival,
                 ER_Stress_Apoptotic_Response = ER_Stress_Apoptotic_Response,
                 Apoptosis_extrinsic_activation = Apoptosis_extrinsic_activation, 
                 Apoptosis_intrinsic_activation = Apoptosis_intrinsic_activation, 
                 Autophagy_recycling_starvation = Autophagy_recycling_starvation,
                 Autophagy_toxic = Autophagy_toxic) 
rm(MAPK, AKT_PI3K, Cell_Cycle_Checkpoint, JAK_STAT, Immunity_Related, DNA_repair, DNA_damage, ER_Stress_Survival,
   ER_Stress_Apoptotic_Response, Apoptosis_extrinsic_activation, Apoptosis_intrinsic_activation, 
   Autophagy_recycling_starvation, Autophagy_toxic)

#calculate mean FC over genes of one pathway for each sample
heatmap_data <- sapply(1:length(pathways), function(x){ #for each pathway
  sapply(colnames(e_foldchange_normalized), function(y){ #for each sample
    abs(mean(e_foldchange_normalized[pathways[[x]], y])) #calculate mean FC of genes of one pathway
  })
})

colnames(heatmap_data) <- c(names(pathways))
rownames(heatmap_data) <- colnames(e_foldchange_normalized)

#Generate heatmap with package gplots
library(gplots)
par(oma = c(1, 1, 1, 10))
heatmap.2(t(heatmap_data),
        scale = "column",
        main = "Mean FC of genes per pathway for 59 celllines",
        cex.main = 0.3,
        xlab = "Celllines",
        col = bluered(100),
        trace = "none", density.info = "none") #trace shows blue lines, in color key density would be shown

#TO DO
#Alternative: package pheatmap (pretty heatmap)
