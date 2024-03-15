
# Data Input & Processing, Initialization of CellChat Object and CCC analysis-------------

## Loading libraries and helper functions---------------
library(Seurat)
library(CellChat)
library(patchwork)
library(harmony)
library(SCP)
library(NMF)
library(ggalluvial)

options(stringsAsFactors = FALSE)





#' A function to perform cell-cell communication analysis using CellChat
#' 
#' @param object A Seurat object containing single-cell RNA-seq data
#' @param cellTypeColumn Column in Seurat metadata containing cell type information
#' @param cellTypes Cell type(s) of interest for analysis
#' @param groupColumn Column in Seurat metadata for further grouping
#' @param group Group of samples (not cells)
#' @param nWorker Number of cores to be used for parallel processing
#' @param min.cells Minimum number of cells in a cell type to be considered for analysis
#' @param saveAs File name to save the resulting CellChat object as an RDS file in the working directory
#' @param db.use Interaction type to be considered (if NULL, all except non-protein signaling will be considered)
#' 
#' @return A CellChat object containing the results of cell-cell communication analysis
#'
cccAnalyzer <- function(object, cellTypeColumn, cellTypes, groupColumn, group, nWorker = 12, min.cells = 20, saveAs = NULL, db.use = NULL) {
  
  # Print progress message
  cat("Subsetting Seurat object...\n")
  
  # Define group column in Seurat object
  object$groupColumn <- object[[groupColumn]]
  
  # Selected cells based on cell types
  cell_ids <- rownames(object@meta.data)[object@meta.data[[cellTypeColumn]] %in% cellTypes]
  
  # Subset based on cell ids and group
  subObj <- subset(object, cell = cell_ids)
  subObj <- subset(subObj, subset = groupColumn == group)
  
  # Define a sample column in Seurat object (required for CellChat)
  subObj$samples <- subObj$orig.ident
  
  # Convert to CellChat object
  cellchat <- createCellChat(object = subObj, group.by = cellTypeColumn, assay = "RNA")
  
  # Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human 
  
  if (is.null(db.use)) {
    CellChatDB.use <- subsetDB(CellChatDB) # Select all except non-protein signaling
  } else {
    # Use a subset of CellChatDB for cell-cell communication analysis
    CellChatDB.use <- subsetDB(CellChatDB, search = db.use, key = "annotation")
  }
  
  # Preprocess the expression data for cell-cell communication analysis
  # Set the used database in the CellChat object
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multisession", workers = nWorker) # Set up parallel processing
  
  # Print progress message
  cat("Identifying overexpressed genes and interactions...\n")
  
  # Identify overexpressed genes and interactions
  cellchat <- identifyOverExpressedGenes(cellchat)
  print( # to print messages 
    cellchat <- identifyOverExpressedInteractions(cellchat)
  )
  
  # Compute communication probabilities
  print( # to print messages from inner function
    cellchat <- computeCommunProb(cellchat,
                                  population.size = TRUE, # Control abundant cell populations
                                  type = "triMean")
  )
  # Filtration
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  
  # Compute CCC at pathway-level
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Compute network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  # Save the resulting CellChat object as an RDS file
  if (!is.null(saveAs)) {
    saveRDS(cellchat, paste0(saveAs, ".rds"))
    cat("CellChat object saved as", paste0(saveAs, ".rds"), "\n")
  }
  
  # Print progress message
  cat("Cell-cell communication analysis completed.\n")
  
  # Return the CellChat object
  return(cellchat)
}



# Function to extract the top pathways based on their strength in the network
# 
# @param object A CellChat object containing signaling network data with centrality measures computed
# @param n The number of top pathways to extract
# @param pattern A character vector specifying the type of signaling pattern to consider: 
#                "outgoing" for outgoing signaling pathways,
#                "incoming" for incoming signaling pathways, 
#                or "all" for both outgoing and incoming pathways.
#                Defaults to "outgoing".
# @return A vector containing the names of the top pathways
# 
topPathways <- function(object, n = 20, pattern = c("outgoing", "incoming", "all")) {
  # Check if the object is a CellChat object
  if (!inherits(object, "CellChat")) {
    stop("The object must be a CellChat object.")
  }
  
  # Extract centrality measures from the CellChat object
  centr <- slot(object, "netP")$centr
  
  # Initialize matrices to store outgoing and incoming signaling strengths
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  
  # Populate the matrices with outgoing and incoming signaling strengths
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  
  # Determine the matrix based on the specified pattern
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  } else {
    stop("Invalid pattern. Choose from 'outgoing', 'incoming', or 'all'.")
  }
  
  # Compute the total signaling strength for each pathway
  pSum <- rowSums(mat)
  
  # Return the top 'n' pathways based on their total signaling strength
  if (n <= nrow(mat)) {
    return(names(head(sort(pSum, decreasing = TRUE), n)))
  } else {
    return(names(sort(pSum, decreasing = TRUE)))
  }
}


## cell chat interaction database --------

png(filename = "./images/CellChatDBCategories.png", width = 16, height = 8.135, units = "in", res = 300)
showDatabaseCategory(CellChatDB.human)
dev.off()


## Preparing scRNA dataset: revised cell labels-----------------


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
harmonized_seurat_13032024 <- readRDS("G:/Documents/chen2020_scRNA/harmonized_seurat_13032024.RDS")

##Convert seurat to cell chat objec-----
# for normal samples we need to convert iCAF and myoCAF to fibroblast/ or progenitors 
nSeu <- subset(harmonized_seurat_13032024, subset=Invasiveness=="normal")

tmpMet <- nSeu@meta.data
tmpMet$clusters <- as.character(tmpMet$clusters)
# convert labeles
tmpMet$clusters <- ifelse(tmpMet$clusters == "i-CAF" & tmpMet$Invasiveness == "normal", "iCAF-progens",
                          ifelse(tmpMet$clusters == "Myo-CAF" & tmpMet$Invasiveness == "normal", "myoCAF-progens",
                                 as.character(tmpMet$clusters)))
# convert back to factor
tmpMet$clusters <- as.factor(tmpMet$clusters)
# set the factor levels comparable to other objects
levs <- c("Epithelial cells","T-cells", "Endothelial cells", "iCAF-progens", "Myeloid cells","myoCAF-progens","APCs","B-cells","Mast cells")
# Relevel the factor variable to match desired levels
tmpMet$clusters <- factor(tmpMet$clusters, levels = levs)

# put it back
nSeu@meta.data <- tmpMet
Idents(nSeu) <- nSeu$clusters
rm(tmpMet)

# 
normal_superCluster <- cccAnalyzer(object = nSeu,
                                   cellTypeColumn="clusters",
                                   cellTypes=unique(nSeu$clusters),
                                   groupColumn="Invasiveness",
                                   group="normal",
                                   saveAs = "normal_superClusters")
# 
# NMIBC
nmibc_superCluster <- cccAnalyzer(object = harmonized_seurat_13032024,
                                   cellTypeColumn="clusters",
                                   cellTypes=unique(harmonized_seurat_13032024$clusters),
                                   groupColumn="Invasiveness",
                                   group="Noninvasive",
                                  saveAs = "nmibc_superClusters")
# MIBC
mibc_superCluster <- cccAnalyzer(object = harmonized_seurat_13032024,
                                   cellTypeColumn="clusters",
                                   cellTypes=unique(harmonized_seurat_13032024$clusters),
                                   groupColumn="Invasiveness",
                                   group="Invasive",
                                   saveAs = "mibc_superClusters")


# Part II: visualization of cell-cell communication -------------

## Aggregated cell-cell communication network visualization-------------------
normal_superCluster <- aggregateNet(normal_superCluster)
nmibc_superCluster <- aggregateNet(nmibc_superCluster)
mibc_superCluster <- aggregateNet(mibc_superCluster)

# Visualization the number of interactions or the total interaction strength (weights)
#png(filename = "./images/int_nCount_Strength.png", width = 20, height = 13.135, units = "in", res = 300)
par(mfrow = c(2,3), xpd=TRUE)
par(mar = c(1, 1, 3, 1)) 
netVisual_circle(normal_superCluster@net$count, vertex.weight = as.numeric(table(normal_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "Benign\nnCount_interactions")
netVisual_circle(nmibc_superCluster@net$count, vertex.weight = as.numeric(table(nmibc_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "NMIBC\nnCount_interactions")
netVisual_circle(mibc_superCluster@net$count, vertex.weight = as.numeric(table(mibc_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "MIBC\nnCount_interactions")
# Weights
netVisual_circle(normal_superCluster@net$weight, vertex.weight = as.numeric(table(normal_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(nmibc_superCluster@net$weight, vertex.weight = as.numeric(table(nmibc_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(mibc_superCluster@net$weight, vertex.weight = as.numeric(table(mibc_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#dev.off()



# Visualize outgoing signals from each cell type: COUNTS
normMat_count <- normal_superCluster@net$count
nmibcMat_count <- nmibc_superCluster@net$count
mibcMat_count <- mibc_superCluster@net$count
#

# select epi, t-cells, endothelials and fibroblasts
normMat_count <- normMat_count[c(1,2,3,4,6), c(1,2,3,4,6)]
nmibcMat_count <- nmibcMat_count[c(1,2,3,4,6), c(1,2,3,4,6)]
mibcMat_count <- mibcMat_count[c(1,2,3,4,6), c(1,2,3,4,6)]



# Visualize outgoing signals from each cell type" WEIGTHS
normMat_weight <- normal_superCluster@net$weight
nmibcMat_weight <- nmibc_superCluster@net$weight
mibcMat_weight <- mibc_superCluster@net$weight
#
# select epi, t-cells, endothelials and fibroblasts
normMat_weight <- normMat_weight[c(1,2,3,4,6), c(1,2,3,4,6)]
nmibcMat_weight <- nmibcMat_weight[c(1,2,3,4,6), c(1,2,3,4,6)]
mibcMat_weight <- mibcMat_weight[c(1,2,3,4,6), c(1,2,3,4,6)]


for (i in 1:nrow(normMat_count)) {
  filename <- paste0("./images/", rownames(nmibcMat_weight)[i], "_sent_out_sigs.png")
  
  png(filename = filename, width = 15, height = 8, units = "in", res = 300)
  par(mfrow = c(2, 3), xpd = TRUE)
  # Extract data for normal, NMIBC, and MIBC
  cmat_nor <- matrix(0, nrow = nrow(normMat_count), ncol = ncol(normMat_count), dimnames = dimnames(normMat_count))
  cmat_nor[i, ] <- normMat_count[i, ]
  
  cmat_nmibc <- matrix(0, nrow = nrow(normMat_count), ncol = ncol(normMat_count), dimnames = dimnames(nmibcMat_count))
  cmat_nmibc[i, ] <- nmibcMat_count[i, ]
  
  cmat_mibc <- matrix(0, nrow = nrow(normMat_count), ncol = ncol(normMat_count), dimnames = dimnames(nmibcMat_count))
  cmat_mibc[i, ] <- mibcMat_count[i, ]
  
  wmat_nor <- matrix(0, nrow = nrow(normMat_weight), ncol = ncol(normMat_weight), dimnames = dimnames(normMat_weight))
  wmat_nor[i, ] <- normMat_weight[i, ]
  
  wmat_nmibc <- matrix(0, nrow = nrow(normMat_weight), ncol = ncol(normMat_weight), dimnames = dimnames(nmibcMat_count))
  wmat_nmibc[i, ] <- nmibcMat_weight[i, ]
  
  wmat_mibc <- matrix(0, nrow = nrow(normMat_weight), ncol = ncol(normMat_weight), dimnames = dimnames(nmibcMat_count))
  wmat_mibc[i, ] <- mibcMat_weight[i, ]
  
  # Plotting
  netVisual_circle(cmat_nor, vertex.weight = as.numeric(table(normal_superCluster@idents))[c(4, 9, 3, 5, 8)], 
                   weight.scale = TRUE, edge.weight.max = max(normMat_count), title.name = paste0("Benign\n", rownames(normMat_count)[i], " nCount_interactions"))
  
  netVisual_circle(cmat_nmibc, vertex.weight = as.numeric(table(nmibc_superCluster@idents))[c(1, 2, 3, 4, 6)], 
                   weight.scale = TRUE, edge.weight.max = max(nmibcMat_count), title.name = paste0("NMIBC\n", rownames(nmibcMat_count)[i], " nCount_interactions"))
  
  netVisual_circle(cmat_mibc, vertex.weight = as.numeric(table(mibc_superCluster@idents))[c(1, 2, 3, 4, 6)], 
                   weight.scale = TRUE, edge.weight.max = max(mibcMat_count), title.name = paste0("MIBC\n", rownames(mibcMat_count)[i], " nCount_interactions"))
  
  netVisual_circle(wmat_nor, vertex.weight = as.numeric(table(normal_superCluster@idents))[c(4, 9, 3, 5, 8)], 
                   weight.scale = TRUE, edge.weight.max = max(normMat_weight), title.name = "Benign Interaction weights/strength")
  
  netVisual_circle(wmat_nmibc, vertex.weight = as.numeric(table(nmibc_superCluster@idents))[c(1, 2, 3, 4, 6)], 
                   weight.scale = TRUE, edge.weight.max = max(nmibcMat_weight), title.name = "NMIBC Interaction weights/strength")
  
  netVisual_circle(wmat_mibc, vertex.weight = as.numeric(table(mibc_superCluster@idents))[c(1, 2, 3, 4, 6)], 
                   weight.scale = TRUE, edge.weight.max = max(mibcMat_weight), title.name = "MIBC Interaction weights/strength")
  
  dev.off()
}



## Signaling pathways showing significant communications--------

# obtain list of significant pathways
normal_pat <- normal_superCluster@netP$pathways
nmibc_pat <- nmibc_superCluster@netP$pathways
mibc_pat <- mibc_superCluster@netP$pathways




# Part III: Systems analysis of cell-cell communication network----

## Identify dominant participants in CCC ----

# Compute the network centrality scores

### the slot 'netP' means the inferred intercellular communication network of signaling pathways
normal_superCluster <- netAnalysis_computeCentrality(normal_superCluster, slot.name = "netP") 
nmibc_superCluster <- netAnalysis_computeCentrality(nmibc_superCluster, slot.name = "netP") 
mibc_superCluster <- netAnalysis_computeCentrality(mibc_superCluster, slot.name = "netP") 

# Find elements present in just one group and not in the other groups
normal_only <- setdiff(normal_pat, union(nmibc_pat, mibc_pat))
nmibc_only <- setdiff(nmibc_pat, union(normal_pat, mibc_pat))
mibc_only <- setdiff(mibc_pat, union(normal_pat, nmibc_pat))

# Find elements present in normal_pat but not in nmibc_pat or mibc_pat
normal_exclusive <- setdiff(normal_pat, union(nmibc_pat, mibc_pat))

# Find elements present in nmibc_pat or mibc_pat but not in normal_pat
cancer_exclsuive <- setdiff(union(nmibc_pat, mibc_pat), normal_pat)

# Find elements present in both nmibc_pat and mibc_pat
cancer_common <- intersect(nmibc_pat, mibc_pat)


# Visualization using scatter plot
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(normal_superCluster) + ggtitle("Benign")
gg2 <- netAnalysis_signalingRole_scatter(nmibc_superCluster)+ ggtitle("NMIBC")
gg3 <- netAnalysis_signalingRole_scatter(mibc_superCluster)+ ggtitle("MIBC")

png(filename = "./images/dominant_sender_receiver_cell_groups.png", width = 14, height = 5.5, units = "in", res = 300)
gg1+gg2+gg3
dev.off()



##Identify signals contributing the most to outgoing or incoming signaling of certain cell groups-----

# select top20 pathways, when pattern is all 

# Define the list of sample groups
sample_groups <- list(normal_superCluster, nmibc_superCluster, mibc_superCluster)
names(sample_groups) <- c("Benign", "NMIBC", "MIBC")

# Loop through each sample group
for (i in 1:length(sample_groups)) {
  # Get top pathways for the current sample group
  top_path <- topPathways(sample_groups[[i]], n = 20, pattern = "all")
  
  # Generate heatmap plots for outgoing, incoming, and all patterns
  ht1 <- netAnalysis_signalingRole_heatmap(sample_groups[[i]], signaling = top_path, width = 5, height = 15, pattern = "outgoing", title = names(sample_groups)[i])
  ht2 <- netAnalysis_signalingRole_heatmap(sample_groups[[i]], signaling = top_path, width = 5, height = 15, pattern = "incoming", title = names(sample_groups)[i])
  ht3 <- netAnalysis_signalingRole_heatmap(sample_groups[[i]], signaling = top_path, width = 5, height = 15, pattern = "all", title = names(sample_groups)[i])
  # Plot the heatmaps
  png(filename = paste0("./images/", names(sample_groups)[i], "_signalingRole_heatmap.png"), width = 9, height = 9, units = "in", res = 300)
  ht1 + ht2 + ht3
  dev.off()
}



## Identify global communication patterns-----

##Identify and visualize outgoing communication pattern
# 
# selectK(sample_groups[[1]], pattern = "outgoing") #3
# selectK(sample_groups[[2]], pattern = "outgoing") #4
# selectK(sample_groups[[3]], pattern = "outgoing") #7
# #
# out_k <- c(3,4,7)
# 
# nPatterns = 3
# 
# png(filename ="./images/test.png", width = 8, height = 12, units = "in", res = 300)
# sample_groups[[1]] <- identifyCommunicationPatterns(sample_groups[[1]], pattern = "outgoing",width = 5, height = 20, k = nPatterns)
# dev.off()

## Manifold and classification learning analysis of signaling networks----
###Identify signaling groups based on their functional similarity


# Initialize an empty plot
plotL <- list()
# Loop through each sample group
for (i in seq_along(sample_groups)) {
  cellchat <- sample_groups[[i]]
  nm <- names(sample_groups)[i]
  # Perform computations
  cellchat <- computeNetSimilarity(cellchat, type = "functional")
  cellchat <- netEmbedding(cellchat, type = "functional")
  cellchat <- netClustering(cellchat, type = "functional")
  # Visualization in 2D-space
  gg <- netVisual_embedding(cellchat, type = "functional", label.size = 3.5) + ggtitle(nm)
  # Add each subplot to the combined plot
  plotL[[i]] <- gg
}

png(filename ="./images/functional_similarities.png", width = 15, height = 6, units = "in", res = 300)
ggarrange(plotL[[1]], plotL[[2]], plotL[[3]],
          ncol = 3, nrow = 1)
dev.off()
