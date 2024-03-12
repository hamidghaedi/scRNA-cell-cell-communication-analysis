---
editor_options: 
  markdown: 
    wrap: 72
---

# scRNA cell-cell communication analysis

To figure out how cells communicate all over the body, we need to
accurately show the links between them and analyze those links on a big
scale.To infer interactions between cells in a single-cell RNA dataset,
CellChat comes in handy. It's an R package designed for analyzing and
visualizing cell-cell communication. The goal is to make it easy for
researchers to identify and understand these interactions through clear
and visually appealing representations.

CellChatDB is an additional resource that complements CellChat. It's a
curated database that includes literature-supported information about
ligand-receptor interactions in various species.

In this repo, I'll be using bladder cancer scRNA datasets that I have
used for other repositories on scRNA data analysis. Briefly the dataset
consisted of eight primary bladder tumor tissues (2 low-grade bladder
urothelial tumors, six high-grade bladder urothelial tumors) along with
3 adjacent normal mucosae (PMID: 33033240)

## Outlines:

1.  **Part I: Data Input & Processing and Initialization of CellChat Object**: this includes load and preprocess single-cell RNA-seq data,
    creating a CellChat object and set ligand-receptor interaction
    database and finally preprocess expression data for cell-cell
    communication analysis

2.  **Part II: Inference of Cell-Cell Communication Network** : In this step we compute
    communication probabilities, filter out weak interactions and will extract
    and visualize the inferred cellular communication network

3.  **Part III: Visualization of Cell-Cell Communication Network**: 
    This include visualizing of communication at various levels (e.g., ligand-receptor
    pairs, signaling pathways) and exploring communication patterns using
    different visualization techniques

4.  **Part IV: Systems Analysis of Cell-Cell Communication Network** :
    In this we will Compute network centrality scores, identify dominant senders and
    receivers, analyze signaling roles and contributions and finally explore
    global communication patterns using manifold and classification
    learning

5.  **Part V: Comparison Analysis of Multiple Datasets using CellChat** :
    In this part we Merge CellChat objects from different datasets, compare and
    visualize interactions, strengths, and major sources/targets, 
    predict general principles of cell-cell communication across
    datasets

6.  **Part VI: Predict General Principles of Cell-Cell Communication**: In this part we will
    compare total interactions and interaction strength between datasets
    , visualize differential interactions and strengths among different
    cell population, identify specific signaling changes between cell
    types and conditions

7.  **Part VII: Identify Conserved and Context-Specific Signaling Pathways**: In this step we will compare overall information flow of each signaling pathway, try to identify conserved and context-specific pathways and finally Visualize pathway
    distances in joint manifold

## Part I: Data Input & Processing and Initialization of CellChat Object

We have three group of cells, present in the dataset, normal, non muscle
invasive bladder cancer(NMIBC) and mucel invasive bladder cancer cells.
Each will be extracted from the harmonized Seurat object (prepared as
outlined here) and they need to be converted to CellChat object

``` r
# Part I: Data input & processing and initialization of CellChat object-------------

# Loading libraries---------------
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)



#Set the ligand-receptor interaction database----------------------
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

```
<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/CellChatDBCategories.png" width="90%">

```
## Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)


# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# OR use a subset of CellChatDB for cell-cell communication analysis
##CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling



# Convert seurat to cell chat object -----
## Normal cells---------------------------
# Obtaining the cells
 su <- readRDS("G:/Documents/chen2020_scRNA/harmonized_seurat.RDS")
 su <- su[, rownames(su@meta.data)[su@meta.data$Invasiveness == "normal"]]
 
 # Obtain normalized data as CellChat needs it
 data.input <- GetAssayData(su, assay = "RNA", slot = "data") # normalized data matrix
 # define a new meta obj
 meta <- data.frame(group = su$newClust, row.names = rownames(su@meta.data)) # create a dataframe of the cell labels
 ### Create a CellChat object-----

 cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
 
# set the used database in the object
 cellchat@DB <- CellChatDB.use
 
 ### Preprocessing the expression data for cell-cell communication analysis----

 # subset the expression data of signaling genes for saving computation cost
 cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
 #
 future::plan("multicore", workers =20) # do parallel
 
 cellchat <- identifyOverExpressedGenes(cellchat)
 cellchat <- identifyOverExpressedInteractions(cellchat)
 # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
 cellchat <- projectData(cellchat, PPI.human)
 
 # saving objects for later use
saveRDS(cellchat, "G:/Documents/SCP_env/renv/normal_cellchat.rds")

## NMIBC cells ---------------------

su <- readRDS("C:/Users/User1/Documents/SCP_env/renv/suCancer.RDS")
# drop low quality epithelial cells
su <- su[, rownames(su@meta.data)[su@meta.data$newClust != 'Epithelial cells']]
su <- su[, rownames(su@meta.data)[su@meta.data$newClust != 'adhesion_signaling_luminal_cell']]
su <- su[, rownames(su@meta.data)[su@meta.data$Invasiveness == "Noninvasive"]]

# Defining a new clustering based on luminal and basal markers:
metDat <- su@meta.data

metDat$newClust <- gsub("cancer_associated_", "", metDat$newClust)
metDat$newClust <- gsub("differentiated_", "", metDat$newClust)
metDat$newClust <- gsub("early_", "", metDat$newClust)
metDat$newClust <- gsub("immunomodulatory_", "", metDat$newClust)
metDat$newClust <- gsub("unique_", "", metDat$newClust)
# Re-assign
su@meta.data <- metDat

# converting to CellChat object
# Obtain normalized data as CellChat needs it
data.input <- GetAssayData(su, assay = "RNA", slot = "data") # normalized data matrix
# define a new meta obj
meta <- data.frame(group = su$newClust, row.names = rownames(su@meta.data)) # create a dataframe of the cell labels
### Create a CellChat object-----
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
# set the used database in the object
cellchat@DB <- CellChatDB.use
### Preprocessing the expression data for cell-cell communication analysis----

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#
future::plan("multicore", workers =20) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)
# saving objects for later use
saveRDS(cellchat, "G:/Documents/SCP_env/renv/nmibc_cellChat.rds")

## MIBC cells----------------
su <- readRDS("C:/Users/User1/Documents/SCP_env/renv/suCancer.RDS")
# drop low quality epithelial cells
su <- su[, rownames(su@meta.data)[su@meta.data$newClust != 'Epithelial cells']]
su <- su[, rownames(su@meta.data)[su@meta.data$newClust != 'adhesion_signaling_luminal_cell']]
su <- su[, rownames(su@meta.data)[su@meta.data$Invasiveness == "Invasive"]]
# Defining a new claustering based on luminal and basal markers:
metDat <- su@meta.data
metDat$newClust <- gsub("cancer_associated_", "intCancer_", metDat$newClust)
metDat$newClust <- gsub("differentiated_", "", metDat$newClust)
metDat$newClust <- gsub("early_", "", metDat$newClust)
metDat$newClust <- gsub("immunomodulatory_", "intCancer_", metDat$newClust)
metDat$newClust <- gsub("unique_", "", metDat$newClust)
# Re-assign
su@meta.data <- metDat
# Obtain normalized data as CellChat needs it
data.input <- GetAssayData(su, assay = "RNA", slot = "data") # normalized data matrix
meta <- data.frame(group = su$newClust, row.names = rownames(su@meta.data)) # create a dataframe of the cell labels

### Create a CellChat object-----

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# set the used database in the object
cellchat@DB <- CellChatDB.use

### Preprocessing the expression data for cell-cell communication analysis----

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#
future::plan("multicore", workers =20) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)
# saving objects for later use
saveRDS(cellchat, "G:/Documents/SCP_env/renv/mibc_cellChat.rds")
```

## Part II: Inference of Cell-Cell Communication Network

Earlier, we saved a CellChat object for each group of cells. Now, we will use these objects to add more information, and subsequently, they will be saved with the same names. As we progress through our analysis, the CellChat objects will be updated iteratively. At the end of this part we  visualize the number of interactions as well as the strength/weights of interactions for normal, NMIBC, and MIBC cells.
Please note that each of the objects `normal_cellchat.rds`, `nmibc_cellChat.rds` and `mibc_cellChat.rds`, need to be read in R at a time using `cellchat <- readRDS("SAVED_CELLCHAT_OBJ")`. 


```R
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
 
<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/normal_interaction_count_strength.png" width="60%">

 For NMIBC cells:
 
<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/nmibc_interaction_count_strength.png" width="60%">

For MIBC cells:

<img src="https://github.com/hamidghaedi/scRNA-cell-cell-communication-analysis/blob/main/images/mibc_interaction_count_strength.png" width="60%">

