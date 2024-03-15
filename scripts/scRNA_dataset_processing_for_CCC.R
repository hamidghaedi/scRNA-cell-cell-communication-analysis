## Preparing scRNA dataset: revised cell labels-----------------
library(Seurat)
library(tidyverse)

##loading and processing Seurat object
harmonized_seurat_04032024 <- readRDS("G:/Documents/chen2020_scRNA/harmonized_seurat_04032024.RDS")
# 
allEpi <- subset(harmonized_seurat_04032024, subset=clusters=="Epithelial cells")

# integration and clsuter identiufication
merged_seurat <- allEpi %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData()

merged_seurat <- SCTransform(merged_seurat, vars.to.regress = c("mitoRatio", "orig.ident"))

# Calculate PCs using variable features determined by SCTransform (3000 by default)
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)
# Integration
harmonized_seurat <- RunHarmony(merged_seurat, 
                                group.by.vars = c("orig.ident", "Surgery_Type"), 
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:40)
# to set reduction to harmony and finding the clusters
harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = c(0.1, 0.2, 0.3, 0.4))

# visualization
Idents(harmonized_seurat) <- harmonized_seurat@meta.data$SCT_snn_res.0.2

#  Visualization of clustering

png(filename = "./images/allEpi_cellDimPlot.png", width = 24, height = 8.135, units = "in", res = 300)
CellDimPlot(harmonized_seurat,
            group.by = c("Invasiveness","SCT_snn_res.0.1","SCT_snn_res.0.2"),
            reduction = "UMAP",
            theme_use = "theme_blank", 
            label = TRUE)
dev.off()

# find markers
allEpi_markers <- FindAllMarkers(object = harmonized_seurat, 
                                 only.pos = TRUE,
                                 logfc.threshold = 0.25)   

top10 <- allEpi_markers %>%
  mutate(delta_pct = (pct.1 - pct.2)) %>%
  #filter(avg_log2FC > 1.5) %>%  # only keep rows where avg_log2FC > 1.5
  group_by(cluster) %>%
  top_n(n = 10, wt = delta_pct)

data.table::fwrite(top10, "allEpi_markers_top10_all_markers.csv")
# for res 0.2
data.table::fwrite(top10, "allEpi_markers_snn_res.0.2_top10_all_markers.csv")

#top10 <- read.csv("allEpi_markers_top10_all_markers.csv")

# Assuming your data is stored in a dataframe called 'top10'
grouped_genes <- top10 %>%
  group_by(cluster) %>%
  summarize(genes = paste(gene, collapse = ", "))



# using markers to name the clusters

# 0       S100A9, IRS2, UCA1, FCRLB, BPGM, MYO16, YEATS4, FRS2, SLC1A6, AC025159.1                
# 1       NRG1, LAMC2, KRT5, BCAM, ITGA2, IGFBP7, LAMB3, PALLD, LINC00513, IL18                   
# 2       HSPA6, YEATS4, FRS2, MYO16, ZFAND2A, SLC35E3, LTO1, AC025159.1, IRS2, MYEOV             
# 3       H19, DLGAP1, TNNT3, BOK-AS1, LINC01980, MYCL, TMEM178B, SEZ6L2, IGF2BP2, AC068587.4     
# 4       MUC4, OLFM4, TRIM31, PLAT, SGMS2, AGR3, CPA6, SGPP2, DSP, SELENOP                       
# 5       LCN15, PLA2G2A, LINC01088, PLAC8, TPM2, SLC7A11-AS1, TRBC2, ANXA10, SLC8A1-AS1, TACC3   
# 6       UPK1B, IDH1, SNX31, VSIG2, UPK1A, UPK2, SPAG16, DMAC1, SMIM30, SCHLAP1                  
# 7       NDUFA4L2, LINC02163, CRH, MSMB, DLGAP1, KRT20, LINC01980, AP005230.1, TESC, TSHZ2       
# 8       LCN15, SPARC, PLA2G2A, LINC01088, ROBO2, IGKC, IGLC1, FABP4, CRTAC1, LGALS1             
# 9       MSMB (basal>> luminal), KCND2 (neuronal), CLIC6, CLDN23, MYCL, PLBD1, HS6ST2(neuronal), C11orf96, PPP1R9A, TMEM178B            
# 10      CDH13, BCAM, IGFBP2, SEMA4A, IGFBP7, LINC00513, STK17A, KRT15, CLU, COL4A5              
# 11      SKAP1, SRGN, AP005262.2, TBC1D5, KRT20, TNRC6A, CRH, CDC14A, SCHLAP1, TSHZ2             
# 12      HLA-DQA1, HLA-DQB1, HLA-DPA1, HLA-DPB1, HLA-DRB1, PSMB9, IFI44L, HLA-DRA, HLA-DMA, HLA-F

# clsuer 0: luminal_basal
# cluster 1: basal_cells
# cluster 2: basal_luminal
# cluster 3: luminal_basal (means more luminal than basal),
# cluster 4: squamous_cells
# cluster5 : basal_luminal
# cluster 6: luminal_differentiated
# cluster 7: basal_luminal_KRT20
# cluster 8: basal_luminal_metastatic: because of expression of metastasis related genes SPARC, PLA2G2A, and ROBO2
# cluster 9: basal_luminal_neuronal
# cluster 10: basal_luminal_KRT15
# Cluster 11: luminal_immunomodulators: luminal because of expressing KRT20and immunomodulator because of SKAP1, SRGN, TNRC6A, CRH, and SCHLAP1
# cluster 12: These cells are APCs becuase of expression HLA classII

# 

# visualize epithelial markers
sel_markers <- c("EPCAM", "CDH1", "KRT20", "UPK1A","UPK2", "ERBB3", "CD44","CDH3","KRT14","KRT16","KRT5")

png(filename = "epithelial_markers_plots.png", width = 16, height = 8.135, units = "in", res = 300)
FeatureDimPlot(
  srt = harmonized_seurat,
  features = sel_markers,
  reduction = "UMAP",
  theme_use = "theme_blank")
dev.off()

# BLCA markers 
blca_markers <- data.frame(
  stringsAsFactors = FALSE,
  CellType = c("Luminal",
               "EMT and smooth muscle","EMT and claudin","Basal","P53-like","Squamous",
               "Immune","Neuroendocrine","Neuronal differentiation",
               "Downregulated CIS","Upregulated CIS","Cancer stem cell"),
  Genes = c("CYP2J2,ERBB2,ERBB3,FGFR3,FOXA1,GATA3,GPX2,KRT18,KRT20,PPARG,XBP1,UPK1A,UPK2","PGM5,DES,C7,SRFP4,COMP,SGCD",
            "ZEB1,ZEB2,VIM,SNAI1,TWIST1,FOXC2,CDH2,CLDN3,CLDN7,CLDN4,CDH1,SNAI2",
            "CD44,CDH3,KRT1,KRT14,KRT16,KRT5,KRT6A,KRT6B,KRT6C",
            "ACTG2,CNN1,MYH11,MFAP4,PGM5,FLNC,ACTC1,DES,PCP4",
            "DSC1,DSC2,DSC3,DSG1,DSG2,DSG3,S100A7,S100A8",
            "CD274,PDCD1LG2,IDO1,CXCL11,L1CAM,SAA1","CHGA,CHGB,SCG2,ENO2,SYP,NCAM1",
            "MSI1,PLEKHG4B,GNG4,PEG10,RND2,APLP1,SOX2,TUBB2B",
            "CRTAC1,CTSE,PADI3","MSN,NR3C1","CD44,KRT5,RPSA,ALDH1A1")
)

# create plot for all markers
for(i in 1:nrow(blca_markers)){
  CellType <- blca_markers$CellType[i]
  genes <- blca_markers$Genes[i]
  png(filename = paste0("./images/",CellType,"_markers_plots.png"), width = 16, height = 8.135, units = "in", res = 300)
  FeatureDimPlot(
    srt = harmonized_seurat,
    features = strsplit(genes, ",")[[1]],
    reduction = "UMAP",
    theme_use = "theme_blank",
    title= paste0(CellType, " markers")
  )
  dev.off()
}


# Rename all identities
harmonized_seurat <- RenameIdents(object = harmonized_seurat, 
                                  "0" = "luminal_basal",
                                  "1" = "basal_cells", 
                                  "2" = "basal_luminal",
                                  "3" = "luminal_basal",
                                  "4" = "squamous_cells",
                                  "5" = "basal_luminal",
                                  "6" = "luminal_differentiated",
                                  "7" = "basal_luminal_KRT20",
                                  "8" = "basal_luminal_metastatic",
                                  "9" = "basal_luminal_neuronal",
                                  "10" = "basal_luminal_KRT15",
                                  "11" = "luminal_immunomodulators",
                                  "12" = "APCs")
# adding to metadata
harmonized_seurat$allEpiCluster <- Idents(harmonized_seurat)

# Add old clsuters names
normEpi_seurat <- readRDS("G:/Documents/normEpi_seurat.RDS")
epi_seurat <- readRDS("G:/Documents/updated_epi_seurat.RDS")

old_annot <- data.frame(cell_id = c(colnames(normEpi_seurat), colnames(epi_seurat)),
                        old_annot <- c(normEpi_seurat$cluster, epi_seurat$clusters))

tmpMet <- harmonized_seurat@meta.data
tmpMet$cell_id <- rownames(tmpMet)
# joing 
tmpMet <- dplyr::left_join(tmpMet, old_annot)
rownames(tmpMet) <- tmpMet$cell_id
# 
harmonized_seurat@meta.data <- tmpMet


# Adding new epithelial cell annotation to the main dataset
# 
allMetDat <- harmonized_seurat_04032024@meta.data
cell_orders <- rownames(allMetDat)

epiMetDat <- harmonized_seurat@meta.data
colnames(epiMetDat)[24] <- "oldEpiClusters"
#
epiMetDat$oldEpiClusters <- as.character(epiMetDat$oldEpiClusters)
epiMetDat$allEpiCluster <- as.character(epiMetDat$allEpiCluster)

allMetDat$cellSubtype2 <-as.character(allMetDat$cellSubtype)
#
# Loop through rows of allMetDat
for (i in 1:nrow(allMetDat)) {
  cell <- allMetDat$cells[i]
  
  # Find the corresponding allEpiCluster value from epiMetDat
  matching_row <- epiMetDat$cells == cell
  
  if (any(matching_row)) {
    allMetDat$cellSubtype2[i] <- epiMetDat$allEpiCluster[matching_row]
  }
}

harmonized_seurat_04032024@meta.data <- allMetDat

# sacing for later use
#saveRDS(harmonized_seurat_04032024, "G:/Documents/chen2020_scRNA/harmonized_seurat_13032024.RDS")