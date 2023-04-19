# Here, were going to analyse the pancreas dataset using trajectory inference methods
load("data/pancreas_dataset_further_analysed.RData")

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
expression_matrix_pancreas <- ReadMtx(
  mtx = "./data/pancreas_sparse_matrix.mtx", features = "./data/pancreas_genes.csv",
  cells = "./data/pancreas_barcodes.csv"
)

# Creating a Seurat object with this counts matrix
pancreas_seurat_object <- CreateSeuratObject(counts = expression_matrix_pancreas)

# Both matrices are the same, which they shouldn't be
pancreas_seurat_object@assays$RNA@data[1:10,1:10]
pancreas_seurat_object@assays$RNA@counts[1:10,1:10]

# Normalize the counts matrix to obtain the expression matrix
pancreas_seurat_object <- NormalizeData(pancreas_seurat_object)

# Now, data has the normalized counts matrix (which corresponds to the expression matrix) and counts has the counts matrix
pancreas_seurat_object@assays$RNA@data[1:10,1:10]
pancreas_seurat_object@assays$RNA@counts[1:10,1:10]













############################################### DYNVERSE ########################################################

# Once we have the Seurat object with the counts and expression matrices, we're going to start with dynverse 
# in order to use the other inference methods

# It's important that we remember that Seurat has the genes (features) as rows and the cells as columns, whereas
# dynverse follows the opposite rule (genes=columns; cells=rows), so that's why we have to transpose the matrices

pancreas_dataset <- wrap_expression(
  counts = t(as.matrix(pancreas_seurat_object@assays$RNA@counts)),
  expression = t(as.matrix(pancreas_seurat_object@assays$RNA@data))
)

# Selecting the best methods for our dataset

guidelines <- guidelines_shiny(pancreas_dataset)
methods_selected <- guidelines$methods_selected 


# We're going to be adding our clusters

clusters <- read.csv("./data/pancreas_clusters.csv", sep = ",", header = TRUE)
colnames(clusters) <- c("cell_id", "group_id")
pancreas_dataset <- add_grouping(pancreas_dataset, clusters)


# Color palette

pancreas_color_palette <- c('Ductal'= '#6F4375', 'Ngn3 low EP'= '#ED070B', 'Ngn3 high EP'= '#F8961E',
                            "Pre-endocrine"= "#F9C74F","Beta"= "#90BE6D", "Alpha"= "#43AA8B", 
                            "Epsilon"= "#577590", "Delta"= "#B2ABF2")


# We're going to use our own dimensionality reduction 
# First, we read the csv file with the coordinates for each cell

umap_plot <- read.csv("./data/pancreas_umap.csv", sep = "\t", header = TRUE)

# Use this function we've created to transform it to a matrix with the cell names as row names
to.matrix <- function(x) {
  m <- as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m 
}

umap_matrix <- to.matrix(umap_plot)


# Since we are going to be using paga_tree (and PAGA), we need to add a start cell. 
# Also add grouping and dimred as prior information. 

pancreas_dataset <- add_prior_information(
  pancreas_dataset, 
  start_id = "AAACCTGCATCATCCC",
  groups_id = pancreas_dataset$grouping,
  dimred = umap_matrix
)














############################################### Running the methods ############################################### 


### 1.-Slingshot

# Infer trajectory (we're going to do three different approaches, one where we don't give
# any previous information, one when we add the root cell and one where we add the root cell and the clustering)


model_slingshot <- infer_trajectory(pancreas_dataset, methods_selected[1])

model_slingshot_root <- infer_trajectory(pancreas_dataset, methods_selected[1], give_priors = c("start_id"))

model_slingshot_root_cluster <- infer_trajectory(pancreas_dataset, methods_selected[1], 
                                                 give_priors = c("start_id", "groups_id"))


## Plotting the trajectory
# Plot the normal graph

plot_graph(model_slingshot, expression_source = pancreas_dataset$expression,
           grouping = pancreas_dataset$grouping, color_cells = "grouping")

plot_graph(model_slingshot_root, expression_source = pancreas_dataset$expression,
           grouping = pancreas_dataset$grouping, color_cells = "grouping")

plot_graph(model_slingshot_root_cluster, expression_source = pancreas_dataset$expression,
           grouping = pancreas_dataset$grouping, color_cells = "grouping")


# Add the matrix as a dimred

model_slingshot_umap <- add_dimred(model_slingshot, dimred = umap_matrix, 
                                   expression_source = pancreas_dataset$expression)

model_slingshot_root_umap <- add_dimred(model_slingshot_root, dimred = umap_matrix, 
                                        expression_source = pancreas_dataset$expression)

model_slingshot_root_cluster_umap <- add_dimred(model_slingshot_root_cluster, dimred = umap_matrix, 
                                                expression_source = pancreas_dataset$expression)


# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_slingshot_umap, 
  expression_source = pancreas_dataset$expression, 
  grouping = pancreas_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=pancreas_color_palette)

plot_dimred(
  model_slingshot_root_umap, 
  expression_source = pancreas_dataset$expression, 
  grouping = pancreas_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=pancreas_color_palette)

plot_dimred(
  model_slingshot_root_cluster_umap, 
  expression_source = pancreas_dataset$expression, 
  grouping = pancreas_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=pancreas_color_palette)


## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_slingshot_umap <- add_root(model_slingshot_umap,
                                 root_cell_id = pancreas_dataset$prior_information$start_id)

model_slingshot_root_umap <- add_root(model_slingshot_root_umap,
                                      root_cell_id = pancreas_dataset$prior_information$start_id)

model_slingshot_root_cluster_umap <- add_root(model_slingshot_root_cluster_umap,
                                              root_cell_id = pancreas_dataset$prior_information$start_id)

## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_slingshot_umap,
  expression_source = pancreas_dataset$expression,
  grouping = pancreas_dataset$grouping,
  features_oi = 50)

plot_heatmap(
  model_slingshot_root_umap,
  expression_source = pancreas_dataset$expression,
  grouping = pancreas_dataset$grouping,
  features_oi = 50)

plot_heatmap(
  model_slingshot_root_cluster_umap,
  expression_source = pancreas_dataset$expression,
  grouping = pancreas_dataset$grouping,
  features_oi = 50)

## Calculate pseudotime

model_slingshot_umap <- add_pseudotime(model_slingshot_umap,pseudotime = calculate_pseudotime(model_slingshot_umap))

model_slingshot_root_umap <- add_pseudotime(model_slingshot_root_umap,
                                            pseudotime = calculate_pseudotime(model_slingshot_root_umap))

model_slingshot_root_cluster_umap <- add_pseudotime(model_slingshot_root_cluster_umap,
                                                    pseudotime = calculate_pseudotime(model_slingshot_root_cluster_umap))

# Plot the pseudotime

plot_dimred(
  model_slingshot_umap, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_slingshot_root_umap, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_slingshot_root_cluster_umap, "pseudotime",
  label_milestones = FALSE)

























### 2.-Paga_tree

# Infer trajectory

model_pagatree <- infer_trajectory(pancreas_dataset, methods_selected[2], give_priors = c("start_id"))

model_pagatree_cluster <- infer_trajectory(pancreas_dataset, methods_selected[2], 
                                           give_priors = c("start_id", "groups_id"))


## Plotting the trajectory
# Plot the normal graph

plot_graph(model_pagatree, expression_source = pancreas_dataset$expression,
           grouping = pancreas_dataset$grouping, color_cells = "grouping")

plot_graph(model_pagatree_cluster, expression_source = pancreas_dataset$expression,
           grouping = pancreas_dataset$grouping, color_cells = "grouping")

# Add the matrix as a dimred

model_pagatree_umap <- add_dimred(model_pagatree, dimred = umap_matrix, 
                                  expression_source = pancreas_dataset$expression)

model_pagatree_cluster_umap <- add_dimred(model_pagatree_cluster, dimred = umap_matrix, 
                                          expression_source = pancreas_dataset$expression)

# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_pagatree_umap, 
  expression_source = pancreas_dataset$expression, 
  grouping = pancreas_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=pancreas_color_palette)

plot_dimred(
  model_pagatree_cluster_umap, 
  expression_source = pancreas_dataset$expression, 
  grouping = pancreas_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=pancreas_color_palette)


## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_pagatree_umap <- add_root(model_pagatree_umap,
                                root_cell_id = pancreas_dataset$prior_information$start_id)

model_pagatree_cluster_umap <- add_root(model_pagatree_cluster_umap,
                                        root_cell_id = pancreas_dataset$prior_information$start_id)

## Predicting and visualizing genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_pagatree_umap,
  expression_source = pancreas_dataset$expression,
  grouping = pancreas_dataset$grouping,
  features_oi = 50
)

plot_heatmap(
  model_pagatree_cluster_umap,
  expression_source = pancreas_dataset$expression,
  grouping = pancreas_dataset$grouping,
  features_oi = 50
)



## Calculate pseudotime

model_pagatree_umap <- add_pseudotime(model_pagatree_umap,pseudotime = calculate_pseudotime(model_pagatree_umap))

model_pagatree_cluster_umap <- add_pseudotime(model_pagatree_cluster_umap,
                                              pseudotime = calculate_pseudotime(model_pagatree_cluster_umap))


# Plot the pseudotime

plot_dimred(
  model_pagatree_umap, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_pagatree_cluster_umap, "pseudotime",
  label_milestones = FALSE)





















### 3.-Scorpius

# Infer trajectory

model_scorpius <- infer_trajectory(pancreas_dataset, methods_selected[3])

model_scorpius_root <- infer_trajectory(pancreas_dataset, methods_selected[3], give_priors = c("start_id"))

model_scorpius_root_cluster <- infer_trajectory(pancreas_dataset, methods_selected[3], 
                                                give_priors = c("start_id", "groups_id"))

## Plotting the trajectory
# Plot the normal graph

plot_graph(model_scorpius, expression_source = pancreas_dataset$expression,
           grouping = pancreas_dataset$grouping, color_cells = "grouping")

plot_graph(model_scorpius_root, expression_source = pancreas_dataset$expression,
           grouping = pancreas_dataset$grouping, color_cells = "grouping")

plot_graph(model_scorpius_root_cluster, expression_source = pancreas_dataset$expression,
           grouping = pancreas_dataset$grouping, color_cells = "grouping")

# Add the matrix as a dimred

model_scorpius_umap <- add_dimred(model_scorpius, dimred = umap_matrix, 
                                  expression_source = pancreas_dataset$expression)

model_scorpius_root_umap <- add_dimred(model_scorpius_root, dimred = umap_matrix, 
                                       expression_source = pancreas_dataset$expression)

model_scorpius_root_cluster_umap <- add_dimred(model_scorpius_root_cluster, dimred = umap_matrix, 
                                               expression_source = pancreas_dataset$expression)

# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_scorpius_umap, 
  expression_source = pancreas_dataset$expression, 
  grouping = pancreas_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=pancreas_color_palette)

plot_dimred(
  model_scorpius_root_umap, 
  expression_source = pancreas_dataset$expression, 
  grouping = pancreas_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=pancreas_color_palette)

plot_dimred(
  model_scorpius_root_cluster_umap, 
  expression_source = pancreas_dataset$expression, 
  grouping = pancreas_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=pancreas_color_palette)

## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_scorpius_umap <- add_root(model_scorpius_umap,
                                root_cell_id = pancreas_dataset$prior_information$start_id)

model_scorpius_root_umap <- add_root(model_scorpius_root_umap,
                                     root_cell_id = pancreas_dataset$prior_information$start_id)

model_scorpius_root_cluster_umap <- add_root(model_scorpius_root_cluster_umap,
                                             root_cell_id = pancreas_dataset$prior_information$start_id)


## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_scorpius_umap,
  expression_source = pancreas_dataset$expression,
  grouping = pancreas_dataset$grouping,
  features_oi = 50
)

plot_heatmap(
  model_scorpius_root_umap,
  expression_source = pancreas_dataset$expression,
  grouping = pancreas_dataset$grouping,
  features_oi = 50
)

plot_heatmap(
  model_scorpius_root_cluster_umap,
  expression_source = pancreas_dataset$expression,
  grouping = pancreas_dataset$grouping,
  features_oi = 50
)

## Calculate pseudotime --> uses milestone_begin as root

model_scorpius_umap <- add_pseudotime(model_scorpius_umap,pseudotime = calculate_pseudotime(model_scorpius_umap))

model_scorpius_root_umap <- add_pseudotime(model_scorpius_root_umap,
                                           pseudotime = calculate_pseudotime(model_scorpius_root_umap))

model_scorpius_root_cluster_umap <- add_pseudotime(model_scorpius_root_cluster_umap,
                                                   pseudotime = calculate_pseudotime(model_scorpius_root_cluster_umap))

# Plot the pseudotime 

plot_dimred(
  model_scorpius_umap, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_scorpius_root_umap, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_scorpius_root_cluster_umap, "pseudotime",
  label_milestones = FALSE)

















### 4.-PAGA

# Infer trajectory

model_paga <- infer_trajectory(pancreas_dataset, "paga", give_priors = c("start_id"))

model_paga_cluster <- infer_trajectory(pancreas_dataset, "paga", give_priors = c("start_id", "groups_id"))

## Plotting the trajectory
# Plot the normal graph

plot_graph(model_paga, expression_source = pancreas_dataset$expression,
           grouping = pancreas_dataset$grouping, color_cells = "grouping")

plot_graph(model_paga_cluster, expression_source = pancreas_dataset$expression,
           grouping = pancreas_dataset$grouping, color_cells = "grouping")

# Add the matrix as a dimred

model_paga_umap <- add_dimred(model_paga, dimred = umap_matrix, expression_source = pancreas_dataset$expression)

model_paga_cluster_umap <- add_dimred(model_paga_cluster, dimred = umap_matrix, 
                                      expression_source = pancreas_dataset$expression)

# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_paga_umap, 
  expression_source = pancreas_dataset$expression, 
  grouping = pancreas_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=pancreas_color_palette)

plot_dimred(
  model_paga_cluster_umap, 
  expression_source = pancreas_dataset$expression, 
  grouping = pancreas_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=pancreas_color_palette)

## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_paga_umap <- add_root(model_paga_umap, root_cell_id = pancreas_dataset$prior_information$start_id)

model_paga_cluster_umap <- add_root(model_paga_cluster_umap, 
                                    root_cell_id = pancreas_dataset$prior_information$start_id)


## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_paga_umap,
  expression_source = pancreas_dataset$expression,
  grouping = pancreas_dataset$grouping,
  features_oi = 50
)

plot_heatmap(
  model_paga_cluster_umap,
  expression_source = pancreas_dataset$expression,
  grouping = pancreas_dataset$grouping,
  features_oi = 50
)

## Calculate pseudotime

model_paga_umap <- add_pseudotime(model_paga_umap, pseudotime = calculate_pseudotime(model_paga_umap))

model_paga_cluster_umap <- add_pseudotime(model_paga_cluster_umap, 
                                          pseudotime = calculate_pseudotime(model_paga_cluster_umap))

# Plot the pseudotime 

plot_dimred(
  model_paga_umap, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_paga_cluster_umap, "pseudotime",
  label_milestones = FALSE)
