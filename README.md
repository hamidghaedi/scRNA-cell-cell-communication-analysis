# scRNA cell-cell communication (CCC) analysis

To figure out how cells communicate all over the body, we need to accurately show the links between them and analyze those links on a big scale.To infer interactions between cells in a single-cell RNA dataset, CellChat comes in handy. It's an R package designed for analyzing and visualizing cell-cell communication. The goal is to make it easy for researchers to identify and understand these interactions through clear and visually appealing representations.

CellChatDB is an additional resource that complements CellChat. It's a curated database that includes literature-supported information about ligand-receptor interactions in various species.

In this repo, I'll be using bladder cancer scRNA datasets that I have used for other repositories on scRNA data analysis. Briefly the dataset consisted of eight primary bladder tumor tissues (2 low-grade bladder urothelial tumors, six high-grade bladder urothelial tumors) along with 3 adjacent normal mucosae (PMID: 33033240)


<!--1.  **Part I: Data input & processing and initialization of CellChat object**: this includes load and preprocess single-cell RNA-seq data,
    creating a CellChat object and set ligand-receptor interaction
    database and finally preprocess expression data for cell-cell
    communication analysis

2.  **Part II: Inference of cell-cell communication network** : In this step we compute
    communication probabilities, filter out weak interactions, and will extract
    and visualize the inferred cellular communication network

3.  **Part III: Visualization of cell-cell communication network**: 
    This includes visualizing communication at various levels (e.g., ligand-receptor
    pairs, signaling pathways) and exploring communication patterns using
    different visualization techniques

4.  **Part IV: Systems analysis of cell-cell communication network** :
    In this, we will Compute network centrality scores, identify dominant senders and
    receivers, analyze signaling roles and contributions, and finally explore
    global communication patterns using manifold and classification
    learning

5.  **Part V: Comparison analysis of multiple datasets using CellChat** :
    In this part, we merge cellChat objects from different datasets, compare and
    visualize interactions, strengths, and major sources/targets, 
    predict general principles of cell-cell communication across
    datasets

6.  **Part VI: Predict General principles of cell-cell communication**: In this part we will
    compare total interactions and interaction strength between datasets
    , visualize differential interactions and strengths among different
    cell population, identify specific signaling changes between cell
    types and conditions

7.  **Part VII: Identify conserved and context-Specific signaling pathways**: In this step, we will compare the overall information flow of each signaling pathway, try to identify conserved and context-specific pathways, and finally Visualize pathway
    distances in joint manifold
-->

## Data Input & Processing, Initialization of CellChat Object and CCC analysis

We have three sample groups present in the dataset: normal, non-muscle invasive bladder cancer (NMIBC), and muscle-invasive (MIBC) bladder cancer cells. Each group will be extracted from the harmonized Seurat object (prepared as outlined here) and converted into a CellChat object. To prepare the Seurat object for CCC analysis, the cell types needed to be revised, and additional columns have been added to the metadata(see `./scripts/scRNA_dataset_processing_for_CCC.R` for more details)

I wrote a function for cell-cell communication analysis (`cccAnalyzer()`) to streamline the analysis steps. This function performs cell-cell communication analysis using the CellChat tool on single-cell RNA-seq data stored in a Seurat object. Here's an overview of the arguments:

`object` (Seurat object): The Seurat object containing single-cell RNA-seq data.

`cellTypeColumn` (character): The name of the column in the Seurat metadata containing cell type information.

`cellTypes` (character vector): A vector specifying the cell types of interest.

`groupColumn` (character): The name of the column in the Seurat metadata used for further grouping.

`group` (character): The specific group within the groupColumn to analyze.

`nWorker` (numeric, default: 12): The number of cores to be used for parallel processing.

`min.cells` (numeric, default: 20): The minimum number of cells in a cell type to be considered for analysis.

`saveAs` (character, default: NULL): A filename to save the CellChat object as an RDS file in WD.

`db.use` (character, default: NULL): The type of interaction to be considered.

The result is a CellChat object that can be used for visualization and further processing.

``` r
# Part I: Data input & processing and initialization of CellChat object-------------

# Loading libraries---------------
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)


# A function to take Seurat obje, make subsets, create cellchat obj and process 

cccAnalyzer <- function(object, # Seurat object
                        cellTypeColumn, # Column in Seurat metadata containing cell type information
                        cellTypes, # Cell type(s) of interest
                        groupColumn, # Column in Seurat metadata for further grouping
                        group, # Group of samples (not cells)
                        nWorker = 12, # Number of cores to be used for cellchat analysis
                        min.cells = 20, # Minimum number of cells in a cell type to be considered for cell-cell communication analysis
                        saveAs = NULL, # File name to save the resulting CellChat object as an RDS file in WD
                        db.use = NULL) { # Interaction type to be considered (if NULL, all except non-protein signaling will be considered)
  
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
```

### CellChat interaction database

Let see what type of interaction is included in CellChat tool.

```r
#Set the ligand-receptor interaction database----------------------
png(filename = "./images/CellChatDBCategories.png.png", width = 16, height = 8.135, units = "in", res = 300)
showDatabaseCategory(CellChatDB.human)
dev.off()
```

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/CellChatDBCategories.png" width="90%"/>


For human CellChatDB v2 contains ~3,300 validated molecular interactions, including ~40% of secrete autocrine/paracrine signaling interactions, ~17% of extracellular matrix (ECM)-receptor interactions, ~13% of cell-cell contact interactions and ~30% non-protein signaling(i.e., metabolic and synaptic signaling).
CellChatDB v2 contains ~3,300 validated molecular interactions, including ~40% of secrete autocrine/paracrine signaling interactions, ~17% of extracellular matrix (ECM)-receptor interactions, ~13% of cell-cell contact interactions and ~30% non-protein signaling(i.e., metabolic and synaptic signaling).



### CCC between cell types in normal, NMIBC and MIBC samples

The cell types are Epithelial cells, T-cells, Endothelial cells, i-CAF, Myeloid cells, Myo-CAF, APCs, B-cells, and Mast cells. Our aim here is to explore the differences in signaling pathways between these groups of cells, considering the sample types (normal, NMIBC, and MIBC).Later, we may follow up on the obtained results from this step by further subtyping cell super clusters, such as epithelial cells into basal, luminal, etc.

``` r
##Convert seurat to cell chat objec-----
# # for normal samples we need to convert iCAF and myoCAF to fibroblast/ or progenitors 
nSeu <- subset(harmonized_seurat_13032024, subset=Invasiveness=="normal")

tmpMet <- nSeu@meta.data
tmpMet$clusters <- as.character(tmpMet$clusters)
# convert labeles
tmpMet$clusters <- ifelse(tmpMet$clusters == "i-CAF" & tmpMet$Invasiveness == "normal", "iCAF-progens",
                          ifelse(tmpMet$clusters == "Myo-CAF" & tmpMet$Invasiveness == "normal", "myoCAF-progens",
                                 as.character(tmpMet$clusters)))
# convert back to factor
tmpMet$clusters <- as.factor(tmpMet$clusters)
# put it back
nSeu@meta.data <- tmpMet
Idents(nSeu) <- nSeu$clusters
rm(tmpMet)


normal_superCluster <- cccAnalyzer(object = nSeu,
                                   cellTypeColumn="clusters",
                                   cellTypes=unique(nSeu$clusters),
                                   groupColumn="Invasiveness",
                                   group="normal",
                                   saveAs = "normal_superClusters")
                                   
# Subsetting Seurat object...
# [1] "Create a CellChat object from a Seurat object"
# The `meta.data` slot in the Seurat object is used as cell meta information 
# Set cell identities for the new CellChat object 
# The cell groups used for CellChat analysis are  APCs, B-cells, Endothelial cells, Epithelial cells, iCAF-progens, Mast cells, Myeloid cells, myoCAF-progens, T-cells 
# Identifying overexpressed genes and interactions...
# The number of highly variable ligand-receptor pairs used for signaling inference is 1789 
# An object of class CellChat created from a single dataset 
#  23190 genes.
#  20787 cells. 
# CellChat analysis of single cell RNA-seq data! 
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-03-14 16:01:12.033785]"
#   |================================================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-03-14 16:34:43.947512]"
# An object of class CellChat created from a single dataset 
#  23190 genes.
#  20787 cells. 
# CellChat analysis of single cell RNA-seq data! 
# CellChat object saved as normal_superClusters.rds 
# Cell-cell communication analysis completed.  

# NMIBC
nmibc_superCluster <- cccAnalyzer(object = harmonized_seurat_13032024,
                                   cellTypeColumn="clusters",
                                   cellTypes=unique(harmonized_seurat_13032024$clusters),
                                   groupColumn="Invasiveness",
                                   group="Noninvasive",
                                   saveAs = "nmibc_superClusters")
                                   
# Subsetting Seurat object...
# [1] "Create a CellChat object from a Seurat object"
# The `meta.data` slot in the Seurat object is used as cell meta information 
# Set cell identities for the new CellChat object 
# The cell groups used for CellChat analysis are  Epithelial cells, T-cells, Endothelial cells, i-CAF, Myeloid cells, Myo-CAF, APCs, B-cells, Mast cells 
# Identifying overexpressed genes and interactions...
# The number of highly variable ligand-receptor pairs used for signaling inference is 1789 
# An object of class CellChat created from a single dataset 
#  23190 genes.
#  20787 cells. 
# CellChat analysis of single cell RNA-seq data! 
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-03-13 23:42:11.814997]"
#   |================================================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-03-14 00:16:23.161631]"
# An object of class CellChat created from a single dataset 
#  23190 genes.
#  20787 cells. 
# CellChat analysis of single cell RNA-seq data! 
# CellChat object saved as normal_superClusters.rds 
# Cell-cell communication analysis completed.
                                   
# MIBC
mibc_superCluster <- cccAnalyzer(object = harmonized_seurat_13032024,
                                   cellTypeColumn="clusters",
                                   cellTypes=unique(harmonized_seurat_13032024$clusters),
                                   groupColumn="Invasiveness",
                                   group="Invasive",
                                   saveAs = "mibc_superClusters")
                                   
# Subsetting Seurat object...
# [1] "Create a CellChat object from a Seurat object"
# The `meta.data` slot in the Seurat object is used as cell meta information 
# Set cell identities for the new CellChat object 
# The cell groups used for CellChat analysis are  Epithelial cells, T-cells, Endothelial cells, i-CAF, Myeloid cells, Myo-CAF, APCs, B-cells, Mast cells 
# Identifying overexpressed genes and interactions...
# The number of highly variable ligand-receptor pairs used for signaling inference is 1718 
# An object of class CellChat created from a single dataset 
#  23190 genes.
#  36270 cells. 
# CellChat analysis of single cell RNA-seq data! 
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-03-14 09:07:25.499606]"
#   |================================================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-03-14 09:24:56.80803]"
# An object of class CellChat created from a single dataset 
#  23190 genes.
#  36270 cells. 
# CellChat analysis of single cell RNA-seq data! 
# CellChat object saved as mibc_superClusters.rds 
# Cell-cell communication analysis completed.
```

## visualization of cell-cell communication

Earlier, we saved a CellChat object for each group of cells. Now, we will use these objects to visualize active signaling pathways in each group of samples. Let's take a look at aggregated number of signaling pathways and their relative weights for each cell type in different sample groups. 

``` r
# Part II: visualization of cell-cell communication -------------

## Aggregated cell-cell communication network visualization-------------------
normal_superCluster <- aggregateNet(normal_superCluster)
nmibc_superCluster <- aggregateNet(nmibc_superCluster)
mibc_superCluster <- aggregateNet(mibc_superCluster)

# Visualization the number of interactions or the total interaction strength (weights)
png(filename = "./images/int_nCount_Strength.png", width = 20, height = 13.135, units = "in", res = 300)
par(mfrow = c(2,3), xpd=TRUE)
netVisual_circle(normal_superCluster@net$count, vertex.weight = as.numeric(table(normal_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "Benign\nnCount_interactions")
netVisual_circle(nmibc_superCluster@net$count, vertex.weight = as.numeric(table(nmibc_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "NMIBC\nnCount_interactions")
netVisual_circle(mibc_superCluster@net$count, vertex.weight = as.numeric(table(mibc_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "MIBC\nnCount_interactions")
# Weights
netVisual_circle(normal_superCluster@net$weight, vertex.weight = as.numeric(table(normal_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(nmibc_superCluster@net$weight, vertex.weight = as.numeric(table(nmibc_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(mibc_superCluster@net$weight, vertex.weight = as.numeric(table(mibc_superCluster@idents)), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

```


<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/int_nCount_Strength.png" width="100%"/>

Based on the above figures, the CCC numbers between epithelial cells, fibroblasts, endothelial cells, and T-cells are significant. In terms of CCC weights/strength, autocrine signaling in endothelial cells and crosstalk between i-CAF and T-cells are noteworthy in benign bladder cells. In NMIBC, autocrine signaling in T-cells and endothelial cells, and crosstalk between these two, are noteworthy. In the case of MIBC, the most significant aspect is autocrine signaling in epithelial cells.

Lets take a look at signaling sent out from these cells in more details:

```R
# Visualize outgoing signals from each cell type: COUNTS
normMat_count <- normal_superCluster@net$count
nmibcMat_count <- nmibc_superCluster@net$count
mibcMat_count <- mibc_superCluster@net$count
#
# reorder normal matrix
normMat_count <- normMat_count[c(4,9,3,5,7,8,1,2,6), c(c(4,9,3,5,7,8,1,2,6))]

# select epi, t-cells, endothelials and fibroblasts
normMat_count <- normMat_count[c(1,2,3,4,6), c(1,2,3,4,6)]
nmibcMat_count <- nmibcMat_count[c(1,2,3,4,6), c(1,2,3,4,6)]
mibcMat_count <- mibcMat_count[c(1,2,3,4,6), c(1,2,3,4,6)]



# Visualize outgoing signals from each cell type" WEIGTHS
normMat_weight <- normal_superCluster@net$weight
nmibcMat_weight <- nmibc_superCluster@net$weight
mibcMat_weight <- mibc_superCluster@net$weight
#
# reorder normal matrix
normMat_weight <- normMat_weight[c(4,9,3,5,7,8,1,2,6), c(c(4,9,3,5,7,8,1,2,6))]

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
  
  cmat_nmibc <- matrix(0, nrow = nrow(normMat_count), ncol = ncol(normMat_count), dimnames = dimnames(normMat_count))
  cmat_nmibc[i, ] <- nmibcMat_count[i, ]
  
  cmat_mibc <- matrix(0, nrow = nrow(normMat_count), ncol = ncol(normMat_count), dimnames = dimnames(normMat_count))
  cmat_mibc[i, ] <- mibcMat_count[i, ]
  
  wmat_nor <- matrix(0, nrow = nrow(normMat_weight), ncol = ncol(normMat_weight), dimnames = dimnames(normMat_weight))
  wmat_nor[i, ] <- normMat_weight[i, ]
  
  wmat_nmibc <- matrix(0, nrow = nrow(normMat_weight), ncol = ncol(normMat_weight), dimnames = dimnames(normMat_weight))
  wmat_nmibc[i, ] <- nmibcMat_weight[i, ]
  
  wmat_mibc <- matrix(0, nrow = nrow(normMat_weight), ncol = ncol(normMat_weight), dimnames = dimnames(normMat_weight))
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
```

**Epithelial cells**

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/Epithelialcells_sent_out_sigs.png" width="100%"/>

**T-cells**

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/T-cells_sent_out_sigs.png" width="100%"/>

**Endothelial cells**

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/Endothelialcells_sent_out_sigs.png" width="100%"/>

**Myo-CAF**

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/Myo-CAF_sent_out_sigs.png" width="100%"/>

**i-CAF**

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/i-CAF_sent_out_sigs.png" width="100%"/>


To identify all the signaling pathways showing significant communications in different sample group, we can access the pathway list for each group, and then visualize those that are more interesting 

