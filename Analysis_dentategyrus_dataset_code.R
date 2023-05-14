# Here, were going to analyse the dentate_gyrus dataset using trajectory inference methods
load("data/DentateGyrus_dataset_further_analysed.RData")

# Libraries needed
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(SeuratObject)
library(loomR)
library(dyno)
library(tidyverse)
library(cowplot)

options(dynwrap_backend = c('container'))







############################################### SEURAT ########################################################

# Create the expression matrix by loading the necessary files
expression_matrix_dentategyrus <- ReadMtx(
  mtx = "./data/dentategyrus_manno_sparse_matrix.mtx", features = "./data/dentategyrus_manno_genes.csv",
  cells = "./data/dentategyrus_manno_barcodes.csv"
)

# Creating a Seurat object with this counts matrix
dentate_gyrus_seurat_object <- CreateSeuratObject(counts = expression_matrix_dentategyrus)

# Both matrices are the same, which they shouldn't be
dentate_gyrus_seurat_object@assays$RNA@data[1:10,1:10]
dentate_gyrus_seurat_object@assays$RNA@counts[1:10,1:10]

# Normalize the counts matrix to obtain the expression matrix
dentate_gyrus_seurat_object <- NormalizeData(dentate_gyrus_seurat_object)

# Now, data has the normalized counts matrix (which corresponds to the expression matrix) and counts has the counts matrix
dentate_gyrus_seurat_object@assays$RNA@data[1:10,1:10]
dentate_gyrus_seurat_object@assays$RNA@counts[1:10,1:10]

















############################################### DYNVERSE ########################################################

# Once we have the Seurat object with the counts and expression matrices, we're going to start with dynverse 
# in order to use the other inference methods

# It's important that we remember that Seurat has the genes (features) as rows and the cells as columns, whereas
# dynverse follows the opposite rule (genes=columns; cells=rows), so that's why we have to transpose the matrices

dentate_gyrus_dataset <- wrap_expression(
  counts = t(as.matrix(dentate_gyrus_seurat_object@assays$RNA@counts)),
  expression = t(as.matrix(dentate_gyrus_seurat_object@assays$RNA@data))
)


# Selecting the best methods for our dataset

guidelines <- guidelines_shiny(dentate_gyrus_dataset)
methods_selected <- guidelines$methods_selected  

# We're also going to be adding our clusters

clusters <- read.csv("./data/dentategyrus_manno_clusters.csv", sep = ",", header = FALSE)
colnames(clusters) <- c("cell_id", "group_id")
dentate_gyrus_dataset <- add_grouping(dentate_gyrus_dataset, clusters)

# Color palette

dentategyrus_color_palette <- c('CA'= '#1f77b4', 'CA1-Sub'= '#aec7e8', 'CA2-3-4'= '#9edae5', 
                                'GlialProg'= '#b24502', 'Granule'= '#f7b6d2', 'ImmAstro'= '#ff9896',
                                'ImmGranule1'= '#c5b0d5', 'ImmGranule2'= '#e377c2', 'Nbl1'= '#bcbd22', 
                                'Nbl2'= '#17becf', 'OPC'= '#c7c7c7', 'RadialGlia'= '#FFCB47',
                                'RadialGlia2'= '#FF6D00', 'nIPC'= '#2ca02c')


# We're going to use our own dimensionality reduction 
# First, we read the csv file with the coordinates for each cell

tsne_plot <- read.csv("./data/dentategyrus_manno_tsne.csv", sep = "\t", header = FALSE)

# Use this function we've created to transform it to a matrix with the cell names as row names
to.matrix <- function(x) {
  m <- as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m 
}

tsne_matrix <- to.matrix(tsne_plot)


# Since we are going to be using paga_tree (and PAGA), we need to add a start cell. 
# Also add grouping and dimred as prior information. 

dentate_gyrus_dataset <- add_prior_information(
  dentate_gyrus_dataset,
  start_id = "10X83_2:AACTCCCTCAAAGTAGx",
  groups_id = dentate_gyrus_dataset$grouping
)
















############################################### Running the methods ############################################### 

### 1.-Slingshot

# Infer trajectory (we're going to do three different approaches, one where we don't give
# any previous information, one when we add the root cell and one where we add the root cell and the clustering)

model_slingshot <- infer_trajectory(dentate_gyrus_dataset, methods_selected[1])

model_slingshot_root <- infer_trajectory(dentate_gyrus_dataset, methods_selected[1], give_priors = c("start_id"))

model_slingshot_root_cluster <- infer_trajectory(dentate_gyrus_dataset, methods_selected[1], 
                                                 give_priors = c("start_id", "groups_id"))


## Plotting the trajectory
# Plot the normal graph

plot_graph(model_slingshot, expression_source = dentate_gyrus_dataset$expression,
           grouping = dentate_gyrus_dataset$grouping, color_cells = "grouping")

plot_graph(model_slingshot_root, expression_source = dentate_gyrus_dataset$expression,
           grouping = dentate_gyrus_dataset$grouping, color_cells = "grouping")

plot_graph(model_slingshot_root_cluster, expression_source = dentate_gyrus_dataset$expression,
           grouping = dentate_gyrus_dataset$grouping, color_cells = "grouping")


# Add the matrix as a dimred

model_slingshot_tsne <- add_dimred(model_slingshot, dimred = tsne_matrix, 
                                   expression_source = dentate_gyrus_dataset$expression)

model_slingshot_root_tsne <- add_dimred(model_slingshot_root, dimred = tsne_matrix, 
                                        expression_source = dentate_gyrus_dataset$expression)

model_slingshot_root_cluster_tsne <- add_dimred(model_slingshot_root_cluster, dimred = tsne_matrix, 
                                                expression_source = dentate_gyrus_dataset$expression)

# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_slingshot_tsne, 
  expression_source = dentate_gyrus_dataset$expression, 
  grouping = dentate_gyrus_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=dentategyrus_color_palette)

plot_dimred(
  model_slingshot_root_tsne, 
  expression_source = dentate_gyrus_dataset$expression, 
  grouping = dentate_gyrus_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=dentategyrus_color_palette)

plot_dimred(
  model_slingshot_root_cluster_tsne, 
  expression_source = dentate_gyrus_dataset$expression, 
  grouping = dentate_gyrus_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=dentategyrus_color_palette)


## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_slingshot_tsne <- add_root(model_slingshot_tsne,
                                 root_cell_id = dentate_gyrus_dataset$prior_information$start_id)

model_slingshot_root_tsne <- add_root(model_slingshot_root_tsne,
                                      root_cell_id = dentate_gyrus_dataset$prior_information$start_id)

model_slingshot_root_cluster_tsne <- add_root(model_slingshot_root_cluster_tsne,
                                              root_cell_id = dentate_gyrus_dataset$prior_information$start_id)

## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_slingshot_tsne,
  expression_source = dentate_gyrus_dataset$expression,
  grouping = dentate_gyrus_dataset$grouping,
  features_oi = 50)

plot_heatmap(
  model_slingshot_root_tsne,
  expression_source = dentate_gyrus_dataset$expression,
  grouping = dentate_gyrus_dataset$grouping,
  features_oi = 50)

plot_heatmap(
  model_slingshot_root_cluster_tsne,
  expression_source = dentate_gyrus_dataset$expression,
  grouping = dentate_gyrus_dataset$grouping,
  features_oi = 50)

## Calculate pseudotime

model_slingshot_tsne <- add_pseudotime(model_slingshot_tsne,pseudotime = calculate_pseudotime(model_slingshot_tsne))

model_slingshot_root_tsne <- add_pseudotime(model_slingshot_root_tsne,
                                            pseudotime = calculate_pseudotime(model_slingshot_root_tsne))

model_slingshot_root_cluster_tsne <- add_pseudotime(model_slingshot_root_cluster_tsne,
                                                    pseudotime = calculate_pseudotime(model_slingshot_root_cluster_tsne))

# Plot the pseudotime

plot_dimred(
  model_slingshot_tsne, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_slingshot_root_tsne, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_slingshot_root_cluster_tsne, "pseudotime",
  label_milestones = FALSE)

















### 2.-Paga_tree

# Infer trajectory

model_pagatree <- infer_trajectory(dentate_gyrus_dataset, methods_selected[2], give_priors = c("start_id"))

model_pagatree_cluster <- infer_trajectory(dentate_gyrus_dataset, methods_selected[2], 
                                           give_priors = c("start_id", "groups_id"))


## Plotting the trajectory
# Plot the normal graph

plot_graph(model_pagatree, expression_source = dentate_gyrus_dataset$expression,
           grouping = dentate_gyrus_dataset$grouping, color_cells = "grouping")

plot_graph(model_pagatree_cluster, expression_source = dentate_gyrus_dataset$expression,
           grouping = dentate_gyrus_dataset$grouping, color_cells = "grouping")


# Add the matrix as a dimred

model_pagatree_tsne <- add_dimred(model_pagatree, dimred = tsne_matrix, 
                                  expression_source = dentate_gyrus_dataset$expression)

model_pagatree_cluster_tsne <- add_dimred(model_pagatree_cluster, 
                                          dimred = tsne_matrix, expression_source = dentate_gyrus_dataset$expression)

# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_pagatree_tsne, 
  expression_source = dentate_gyrus_dataset$expression, 
  grouping = dentate_gyrus_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=dentategyrus_color_palette)

plot_dimred(
  model_pagatree_cluster_tsne, 
  expression_source = dentate_gyrus_dataset$expression, 
  grouping = dentate_gyrus_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=dentategyrus_color_palette)


## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_pagatree_tsne <- add_root(model_pagatree_tsne,
                                root_cell_id = dentate_gyrus_dataset$prior_information$start_id)

model_pagatree_cluster_tsne <- add_root(model_pagatree_cluster_tsne,
                                        root_cell_id = dentate_gyrus_dataset$prior_information$start_id)

## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_pagatree_tsne,
  expression_source = dentate_gyrus_dataset$expression,
  grouping = dentate_gyrus_dataset$grouping,
  features_oi = 50
)

plot_heatmap(
  model_pagatree_cluster_tsne,
  expression_source = dentate_gyrus_dataset$expression,
  grouping = dentate_gyrus_dataset$grouping,
  features_oi = 50
)


## Calculate pseudotime 

model_pagatree_tsne <- add_pseudotime(model_pagatree_tsne,pseudotime = calculate_pseudotime(model_pagatree_tsne))

model_pagatree_cluster_tsne <- add_pseudotime(model_pagatree_cluster_tsne,
                                              pseudotime = calculate_pseudotime(model_pagatree_cluster_tsne))


# Plot the pseudotime

plot_dimred(
  model_pagatree_tsne, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_pagatree_cluster_tsne, "pseudotime",
  label_milestones = FALSE)













### 3.-Scorpius

# Infer trajectory

model_scorpius <- infer_trajectory(dentate_gyrus_dataset, methods_selected[3])

model_scorpius_root <- infer_trajectory(dentate_gyrus_dataset, methods_selected[3], give_priors = c("start_id"))

model_scorpius_root_cluster <- infer_trajectory(dentate_gyrus_dataset, methods_selected[3], 
                                                give_priors = c("start_id", "groups_id"))

## Plotting the trajectory
# Plot the normal graph

plot_graph(model_scorpius, expression_source = dentate_gyrus_dataset$expression,
           grouping = dentate_gyrus_dataset$grouping, color_cells = "grouping")

plot_graph(model_scorpius_root, expression_source = dentate_gyrus_dataset$expression,
           grouping = dentate_gyrus_dataset$grouping, color_cells = "grouping")

plot_graph(model_scorpius_root_cluster, expression_source = dentate_gyrus_dataset$expression,
           grouping = dentate_gyrus_dataset$grouping, color_cells = "grouping")

# Add the matrix as a dimred

model_scorpius_tsne <- add_dimred(model_scorpius, dimred = tsne_matrix, 
                                  expression_source = dentate_gyrus_dataset$expression)

model_scorpius_root_tsne <- add_dimred(model_scorpius_root, dimred = tsne_matrix, 
                                       expression_source = dentate_gyrus_dataset$expression)

model_scorpius_root_cluster_tsne <- add_dimred(model_scorpius_root_cluster, dimred = tsne_matrix, 
                                               expression_source = dentate_gyrus_dataset$expression)

# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_scorpius_tsne, 
  expression_source = dentate_gyrus_dataset$expression, 
  grouping = dentate_gyrus_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=dentategyrus_color_palette)

plot_dimred(
  model_scorpius_root_tsne, 
  expression_source = dentate_gyrus_dataset$expression, 
  grouping = dentate_gyrus_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=dentategyrus_color_palette)

plot_dimred(
  model_scorpius_root_cluster_tsne, 
  expression_source = dentate_gyrus_dataset$expression, 
  grouping = dentate_gyrus_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=dentategyrus_color_palette)

## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_scorpius_tsne <- add_root(model_scorpius_tsne,
                                root_cell_id = dentate_gyrus_dataset$prior_information$start_id)

model_scorpius_root_tsne <- add_root(model_scorpius_root_tsne,
                                     root_cell_id = dentate_gyrus_dataset$prior_information$start_id)

model_scorpius_root_cluster_tsne <- add_root(model_scorpius_root_cluster_tsne,
                                             root_cell_id = dentate_gyrus_dataset$prior_information$start_id)


## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_scorpius_tsne,
  expression_source = dentate_gyrus_dataset$expression,
  grouping = dentate_gyrus_dataset$grouping,
  features_oi = 50
)

plot_heatmap(
  model_scorpius_root_tsne,
  expression_source = dentate_gyrus_dataset$expression,
  grouping = dentate_gyrus_dataset$grouping,
  features_oi = 50
)
plot_heatmap(
  model_scorpius_root_cluster_tsne,
  expression_source = dentate_gyrus_dataset$expression,
  grouping = dentate_gyrus_dataset$grouping,
  features_oi = 50
)

## Calculate pseudotime --> uses milestone_begin as root

model_scorpius_tsne <- add_pseudotime(model_scorpius_tsne,pseudotime = calculate_pseudotime(model_scorpius_tsne))

model_scorpius_root_tsne <- add_pseudotime(model_scorpius_root_tsne,
                                           pseudotime = calculate_pseudotime(model_scorpius_root_tsne))

model_scorpius_root_cluster_tsne <- add_pseudotime(model_scorpius_root_cluster_tsne,
                                                   pseudotime = calculate_pseudotime(model_scorpius_root_cluster_tsne))

# Plot the pseudotime 

plot_dimred(
  model_scorpius_tsne, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_scorpius_root_tsne, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_scorpius_root_cluster_tsne, "pseudotime",
  label_milestones = FALSE)













### 4.-PAGA

# Infer trajectory

model_paga <- infer_trajectory(dentate_gyrus_dataset, "paga", give_priors = c("start_id"))

model_paga_cluster <- infer_trajectory(dentate_gyrus_dataset, "paga", give_priors = c("start_id", "groups_id"))


## Plotting the trajectory
# Plot the normal graph

plot_graph(model_paga, expression_source = dentate_gyrus_dataset$expression,
           grouping = dentate_gyrus_dataset$grouping, color_cells = "grouping")

plot_graph(model_paga_cluster, expression_source = dentate_gyrus_dataset$expression,
           grouping = dentate_gyrus_dataset$grouping, color_cells = "grouping")


# Add the matrix as a dimred

model_paga_tsne <- add_dimred(model_paga, dimred = tsne_matrix, expression_source = dentate_gyrus_dataset$expression)

model_paga_cluster_tsne <- add_dimred(model_paga_cluster, dimred = tsne_matrix, 
                                      expression_source = dentate_gyrus_dataset$expression)

# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_paga_tsne, 
  expression_source = dentate_gyrus_dataset$expression, 
  grouping = dentate_gyrus_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=dentategyrus_color_palette)

plot_dimred(
  model_paga_cluster_tsne, 
  expression_source = dentate_gyrus_dataset$expression, 
  grouping = dentate_gyrus_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=dentategyrus_color_palette)

## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_paga_tsne <- add_root(model_paga_tsne, root_cell_id = dentate_gyrus_dataset$prior_information$start_id)

model_paga_cluster_tsne <- add_root(model_paga_cluster_tsne, 
                                    root_cell_id = dentate_gyrus_dataset$prior_information$start_id)


## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_paga_tsne,
  expression_source = dentate_gyrus_dataset$expression,
  grouping = dentate_gyrus_dataset$grouping,
  features_oi = 50
)

plot_heatmap(
  model_paga_cluster_tsne,
  expression_source = dentate_gyrus_dataset$expression,
  grouping = dentate_gyrus_dataset$grouping,
  features_oi = 50
)

## Calculate pseudotime --> uses milestone_begin as root

model_paga_tsne <- add_pseudotime(model_paga_tsne, pseudotime = calculate_pseudotime(model_paga_tsne))

model_paga_cluster_tsne <- add_pseudotime(model_paga_cluster_tsne, 
                                          pseudotime = calculate_pseudotime(model_paga_cluster_tsne))

# Plot the pseudotime 

plot_dimred(
  model_paga_tsne, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_paga_cluster_tsne, "pseudotime",
  label_milestones = FALSE)


