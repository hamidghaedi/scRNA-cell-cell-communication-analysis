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

We have three groups of cells, present in the dataset, normal, non-muscle invasive bladder cancer(NMIBC) and muscle-invasive bladder cancer cells. Each will be extracted from the harmonized Seurat object (prepared as outlined here) and they need to be converted to CellChat object

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
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
```

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/CellChatDBCategories.png" width="90%"/>

CellChatDB v2 contains ~3,300 validated molecular interactions, including ~40% of secrete autocrine/paracrine signaling interactions, ~17% of extracellular matrix (ECM)-receptor interactions, ~13% of cell-cell contact interactions and ~30% non-protein signaling(i.e., metabolic and synaptic signaling).


### CCC between cell types in normal, NMIBC and MIBC samples

The cell types are Epithelial cells, T-cells, Endothelial cells, i-CAF, Myeloid cells, Myo-CAF, APCs, B-cells, and Mast cells. Our aim here is to explore the differences in signaling pathways between these groups of cells, considering the sample types (normal, NMIBC, and MIBC).Later, we may follow up on the obtained results from this step by further subtyping cell super clusters, such as epithelial cells into basal, luminal, etc.

``` r
##Convert seurat to cell chat objec-----
normal_superCluster <- cccAnalyzer(object = harmonized_seurat_13032024,
                                   cellTypeColumn="clusters",
                                   cellTypes=unique(harmonized_seurat_13032024$clusters),
                                   groupColumn="Invasiveness",
                                   group="normal",
                                   saveAs = "normal_superClusters")
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


nmibc_superCluster <- cccAnalyzer(object = harmonized_seurat_13032024,
                                   cellTypeColumn="clusters",
                                   cellTypes=unique(harmonized_seurat_13032024$clusters),
                                   groupColumn="Invasiveness",
                                   group="Noninvasive",
                                   saveAs = "nmibc_superClusters")
                                   
mibc_superCluster <- cccAnalyzer(object = harmonized_seurat_13032024,
                                   cellTypeColumn="clusters",
                                   cellTypes=unique(harmonized_seurat_13032024$clusters),
                                   groupColumn="Invasiveness",
                                   group="Invasive",
                                   saveAs = "mibc_superClusters")
```

## Part II: Inference of Cell-Cell Communication Network

Earlier, we saved a CellChat object for each group of cells. Now, we will use these objects to add more information, and subsequently, they will be saved with the same names. As we progress through our analysis, the CellChat objects will be updated iteratively. At the end of this part we visualize the number of interactions as well as the strength/weights of interactions for normal, NMIBC, and MIBC cells. Please note that each of the objects `normal_cellchat.rds`, `nmibc_cellChat.rds` and `mibc_cellChat.rds`, need to be read in R at a time using `cellchat <- readRDS("SAVED_CELLCHAT_OBJ")`.

``` r
# Part II: Inference of cell-cell communication network-------------

# Compute the communication probability/strength between any interacting cell groups -----

## Two options based on the methods for computing the average gene expression per cell group

### 1) "triMean", producing fewer but stronger interactions
#cellchat <- computeCommunProb(cellchat, 
#                             population.size = TRUE, 
#                              raw.use = FALSE,
#                              type = "triMean") # defult.It approximates 25% truncated mean, implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%

### 2) truncatedMean, producing more interactions. A value should be assigned to ’trim'
# If the above returned not much data, run the following:
cellchat <- computeCommunProb(cellchat, 
                              population.size = TRUE, 
                              raw.use = FALSE, 
                              type = "truncatedMean",
                              trim = 0.05) # 

# Filter out the cell-cell communications ------------
## if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Extract the inferred cellular communication network as a data frame

df.net <- subsetCommunication(cellchat) # returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors.

df.net <- subsetCommunication(cellchat, slot.name = "netP") # to access the the inferred communications at the level of signaling pathways

df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) # gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) # gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

# Infer the cell-cell communication at a signaling pathway level

cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network-----

#cellchat <- aggregateNet(cellchat,
#                         sources.use = c("basal_cell", "cancer_associated_luminal_cell", "differentiated_luminal_cell", 
#                                         "early_basal_cell", "immunomodulatory_luminal_cell","unique_luminal_cell"),
#                         targets.use = c("B cells", "T-cells", "myeloid cells"))


cellchat <- aggregateNet(cellchat)

# Visualization the number of interactions or the total interaction strength (weights)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


# This snipet generate subplots of the signaling sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
```

For normal cells:

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/normal_interaction_count_strength.png" width="60%"/>

For NMIBC cells:

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/nmibc_interaction_count_strength.png" width="60%"/>

For MIBC cells:

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/mibc_interaction_count_strength.png" width="60%"/>
