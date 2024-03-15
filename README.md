# scRNA cell-cell communication (CCC) analysis

To figure out how cells communicate all over the body, we need to accurately show the links between them and analyze those links on a big scale.To infer interactions between cells in a single-cell RNA dataset, CellChat comes in handy. It's an R package designed for analyzing and visualizing cell-cell communication. The goal is to make it easy for researchers to identify and understand these interactions through clear and visually appealing representations.

CellChatDB is an additional resource that complements CellChat. It's a curated database that includes literature-supported information about ligand-receptor interactions in various species.

In this repo, I'll be using bladder cancer scRNA datasets that I have used for other repositories on scRNA data analysis. Briefly the dataset consisted of eight primary bladder tumor tissues (2 low-grade bladder urothelial tumors, six high-grade bladder urothelial tumors) along with 3 adjacent normal mucosae (PMID: 33033240)


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

```

### CellChat interaction database

Let see what type of interaction is included in CellChat tool.

```r
#Set the ligand-receptor interaction database----------------------
png(filename = "./images/CellChatDBCategories.png", width = 16, height = 8.135, units = "in", res = 300)
showDatabaseCategory(CellChatDB.human)
dev.off()
```

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/CellChatDBCategories.png" width="100%"/>


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
# set the factor levels comparable to other objects
levs <- c("Epithelial cells","T-cells", "Endothelial cells", "iCAF-progens", "Myeloid cells","myoCAF-progens","APCs","B-cells","Mast cells")
# Relevel the factor variable to match desired levels
tmpMet$clusters <- factor(tmpMet$clusters, levels = levs)

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
# The cell groups used for CellChat analysis are  Epithelial cells, T-cells, Endothelial cells, iCAF-progens, Myeloid cells, myoCAF-progens, APCs, B-cells, Mast cells 
# Identifying overexpressed genes and interactions...
# The number of highly variable ligand-receptor pairs used for signaling inference is 1789 
# An object of class CellChat created from a single dataset 
#  23190 genes.
#  20787 cells. 
# CellChat analysis of single cell RNA-seq data! 
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-03-15 10:05:26.91681]"
#   |=============================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-03-15 10:34:01.148103]"
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


<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/int_nCount_Strength.png" width="110%"/>

Based on the above figures, the CCC numbers between epithelial cells, fibroblasts, endothelial cells, and T-cells are significant. In terms of CCC weights/strength, autocrine signaling in endothelial cells and crosstalk between i-CAF and T-cells are noteworthy in benign bladder cells. In NMIBC, autocrine signaling in T-cells and endothelial cells, and crosstalk between these two, are noteworthy. In the case of MIBC, the most significant aspect is autocrine signaling in epithelial cells.

Lets take a look at signaling sent out from these cells in more details:

```R

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


To identify all the signaling pathways showing significant communications in different sample group, we can access the pathway list for each group

```R
## Signaling pathways showing significant communications--------

# Find elements present in mibc but not in nmibc and normal
mibc_only <- setdiff(mibc_superCluster@netP$pathways, c(nmibc_superCluster@netP$pathways, normal_superCluster@netP$pathways))

# Find elements present in nmibc or mibc but not in normal
nmibc_mibc_only <- setdiff(union(nmibc_superCluster@netP$pathways, mibc_superCluster@netP$pathways), normal_superCluster@netP$pathways)

# Find elements present in nmibc only
nmibc_only <- setdiff(nmibc_superCluster@netP$pathways, c(mibc_superCluster@netP$pathways, normal_superCluster@netP$pathways))

# Find elements present in all three groups
common <- intersect(intersect(nmibc_superCluster@netP$pathways, mibc_superCluster@netP$pathways), normal_superCluster@netP$pathways)

# 
# > mibc_only
# character(0)
# > nmibc_mibc_only
# [1] "IFN-II"    "ApoE"      "ncWNT"     "CD39"      "CSPG4"     "APRIL"     "LAIR1"     "HGF"       "SLITRK"    "DESMOSOME"
# > nmibc_only
# [1] "IFN-II"    "ApoE"      "ncWNT"     "CD39"      "APRIL"     "LAIR1"     "HGF"       "SLITRK"    "DESMOSOME"
# > common
# [1] "COLLAGEN" "MIF"      "LAMININ"  "APP"      "CD99"     "CypA"     "FN1"      "MK"       "PTPRM"    "MHC-II"   "ADGRE"    "VISFATIN"
# [13] "GALECTIN" "PARs"     "PECAM1"   "VEGF"     "JAM"      "IGFBP"    "SEMA3"    "ANGPT"    "CALCR"    "ESAM"     "GAP"      "NOTCH"   
# [25] "NECTIN"   "ADGRG"    "CD46"     "ANGPTL"   "GDF"      "TENASCIN" "EPHA"     "EGF"      "CDH5"     "VCAM"     "CDH1"     "PLAU"    
# [37] "GAS"      "EPHB"     "PECAM2"   "MPZ"      "CDH"      "SEMA4"    "ADGRA"    "GRN"      "SEMA6"    "PDGF"     "SELL"     "NRXN"    
# [49] "CLDN"     "TWEAK"    "BMP"      "CADM"     "OCLN"    
```

## Systems analysis of cell-cell communication network

CellChat finds the most important cells in communication networks by computing several network centrality measures for each cell group. It looks at different measures to figure out which cells send, receive, mediate, or influence signals the most.

These measures include:

*Out-Degree* : How much a cell sends signals.

*In-Degree* : How much a cell receives signals.

*Flow Betweenness* : How much a cell helps signals flow between other cells.

*Information Centrality* : How important a cell is for spreading information.

These measures work best in networks where the strength of connections between cells is known (weighteddirected network). In these networks, Out-Degree helps find major senders, In-Degree finds major receivers, Flow Betweenness identifies mediators, and Information Centrality shows influential cells.


We can visualize dominant sender and receiver per cell type in different sample group using the following 

```R
# Visualization using scatter plot
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(normal_superCluster) + ggtitle("Benign")
gg2 <- netAnalysis_signalingRole_scatter(nmibc_superCluster)+ ggtitle("NMIBC")
gg3 <- netAnalysis_signalingRole_scatter(mibc_superCluster)+ ggtitle("MIBC")

png(filename = "./images/dominant_sender_receiver_cell_groups.png", width = 14, height = 5.5, units = "in", res = 300)
gg1+gg2+gg3
dev.off()
```
<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/dominant_sender_receiver_cell_groups.png" width="100%"/>

As it can be seen in the above figures, there is a sort of balanced CCC in benign tissue where endothelial cells and T-cells emerge as receivers (high value in y-axis), while endothelial cells and iCAF-progenitors show high value along the x-axis, indicating that they are dominant senders.

In the case of NMIBC, T-cells appear to be dominant receivers and endothelial cells appear to be dominant senders. In NMIBC, epithelial cells are both dominant senders and receivers (autocrine activity). These patterns, especially the epithelial cell subtype in MIBC, warrant further investigation.

Now let identify signals contributing the most to outgoing or incoming signaling of certain cell groups.We can  answer the question on which signals contributing most to outgoing or incoming signaling of certain cell groups by looking at the below heatmaps

```R
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

```

**Dominant signaling in benign cells**

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/Benign_signalingRole_heatmap.png" width="100%"/>

**Dominant signaling in NMIBC cells**

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/NMIBC_signalingRole_heatmap.png" width="100%"/>


**Dominant signaling in MIBC cells**

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/MIBC_signalingRole_heatmap.png" width="100%"/>


Manifold and classification learning analysis of signaling networks

CellChat employs manifold and classification learning analyses of signaling networks to assess the similarity between significant signaling pathways and group them accordingly. Manifold learning techniques help in understanding the intrinsic structure or geometry of high-dimensional data like signaling networks. These methods aim to represent the data in a lower-dimensional space while retaining essential properties (PCA, tSNE, UMAP).

Functional similarity analysis in CellChat identifies pathways or ligand-receptor pairs that serve similar roles by assessing their degree of similarity. This analysis requires consistency in cell population composition between datasets. On the other hand, structural similarity analysis in CellChat compares the structure of signaling networks regardless of the similarity of senders and receivers. This allows for a broader understanding of network organization beyond specific cell interactions.

```R
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
```
<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/functional_similarities.png" width="100%"/>


