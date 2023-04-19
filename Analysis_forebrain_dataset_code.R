# Here, were going to analyse the pancreas dataset using trajectory inference methods
load("data/forebrain_dataset_further_analysed.RData")

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
expression_matrix_forebrain <- ReadMtx(
  mtx = "./data/forebrain_sparse_matrix.mtx", features = "./data/forebrain_genes.csv",
  cells = "./data/forebrain_barcodes.csv"
)

# Creating a Seurat object with this counts matrix
forebrain_seurat_object <- CreateSeuratObject(counts = expression_matrix_forebrain)

# Both matrices are the same, which they shouldn't be
forebrain_seurat_object@assays$RNA@data[1:20,1:20]
forebrain_seurat_object@assays$RNA@counts[1:20,1:20]

# Normalize the counts matrix to obtain the expression matrix
forebrain_seurat_object <- NormalizeData(forebrain_seurat_object)

# Now, data has the normalized counts matrix (which corresponds to the expression matrix) and counts has the counts matrix
forebrain_seurat_object@assays$RNA@data[1:20,1:20]
forebrain_seurat_object@assays$RNA@counts[1:20,1:20]





















############################################### DYNVERSE ########################################################

# Once we have the Seurat object with the counts and expression matrices, we're going to start with dynverse 
# in order to use the other inference methods

# It's important that we remember that Seurat has the genes (features) as rows and the cells as columns, whereas
# dynverse follows the opposite rule (genes=columns; cells=rows), so that's why we have to transpose the matrices

forebrain_dataset <- wrap_expression(
  counts = t(as.matrix(forebrain_seurat_object@assays$RNA@counts)),
  expression = t(as.matrix(forebrain_seurat_object@assays$RNA@data))
)


# Selecting the best methods for our dataset

guidelines <- guidelines_shiny(forebrain_dataset)
methods_selected <- guidelines$methods_selected 

# We're going to be adding our clusters

clusters <- read.csv("./data/forebrain_clusters.csv", sep = ",", header = FALSE)
colnames(clusters) <- c("cell_id", "group_id")
forebrain_dataset <- add_grouping(forebrain_dataset, clusters)

# Color palette

forebrain_color_palette <- c("0"= '#7A0000', "1"= '#FF6D00', "2"= '#FF9E00', "3"= "#20003D", "4"= "#7C33BC", 
                             "5"= "#5686D9", "6"= "#18AF9D")


# We're going to use our own dimensionality reduction 
# First, we read the csv file with the coordinates for each cell

pca_plot <- read.csv("./data/forebrain_pca_6.csv", sep = "\t", header = FALSE)

# Use this function we've created to transform it to a matrix with the cell names as row names
to.matrix <- function(x) {
  m <- as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m 
}

# Our dimred matrix
pca_matrix <- to.matrix(pca_plot)


# Since we are going to be using paga_tree (and PAGA), we need to add a start cell. 
# Also add grouping and dimred as prior information.   

forebrain_dataset <- add_prior_information(
  forebrain_dataset,
  start_id = "10X_17_028:AAACGGGTCTGCGTAAx",
  groups_id = forebrain_dataset$grouping
)

















############################################### Running the methods ############################################### 

### 1.-Slingshot

# Infer trajectory (we're going to do three different approaches, one where we don't give
# any previous information, one when we add the root cell and one where we add the root cell and the clustering)

model_slingshot <- infer_trajectory(forebrain_dataset, methods_selected[1])

model_slingshot_root <- infer_trajectory(forebrain_dataset, methods_selected[1], give_priors = c("start_id"))

model_slingshot_root_cluster <- infer_trajectory(forebrain_dataset, methods_selected[1], 
                                                 give_priors = c("start_id", "groups_id"))


## Plotting the trajectory
# Plot the normal graph

plot_graph(model_slingshot, expression_source = forebrain_dataset$expression,
           grouping = forebrain_dataset$grouping, color_cells = "grouping")

plot_graph(model_slingshot_root, expression_source = forebrain_dataset$expression,
           grouping = forebrain_dataset$grouping, color_cells = "grouping")

plot_graph(model_slingshot_root_cluster, expression_source = forebrain_dataset$expression,
           grouping = forebrain_dataset$grouping, color_cells = "grouping")


# Add the matrix as a dimred

model_slingshot_pca <- add_dimred(model_slingshot, dimred = pca_matrix, 
                                  expression_source = forebrain_dataset$expression)

model_slingshot_root_pca <- add_dimred(model_slingshot_root, dimred = pca_matrix, 
                                       expression_source = forebrain_dataset$expression)

model_slingshot_root_cluster_pca <- add_dimred(model_slingshot_root_cluster, dimred = pca_matrix, 
                                               expression_source = forebrain_dataset$expression)


# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_slingshot_pca, 
  expression_source = forebrain_dataset$expression, 
  grouping = forebrain_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=forebrain_color_palette)

plot_dimred(
  model_slingshot_root_pca, 
  expression_source = forebrain_dataset$expression, 
  grouping = forebrain_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=forebrain_color_palette)

plot_dimred(
  model_slingshot_root_cluster_pca, 
  expression_source = forebrain_dataset$expression, 
  grouping = forebrain_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=forebrain_color_palette)



## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_slingshot_pca <- add_root(model_slingshot_pca,
                                root_cell_id = forebrain_dataset$prior_information$start_id)

model_slingshot_root_pca <- add_root(model_slingshot_root_pca,
                                     root_cell_id = forebrain_dataset$prior_information$start_id)

model_slingshot_root_cluster_pca <- add_root(model_slingshot_root_cluster_pca,
                                             root_cell_id = forebrain_dataset$prior_information$start_id)

## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_slingshot_pca,
  expression_source = forebrain_dataset$expression,
  grouping = forebrain_dataset$grouping,
  features_oi = 50)

plot_heatmap(
  model_slingshot_root_pca,
  expression_source = forebrain_dataset$expression,
  grouping = forebrain_dataset$grouping,
  features_oi = 50)

plot_heatmap(
  model_slingshot_root_cluster_pca,
  expression_source = forebrain_dataset$expression,
  grouping = forebrain_dataset$grouping,
  features_oi = 50)

## Calculate pseudotime --> used 2 as root, which is correct

model_slingshot_pca <- add_pseudotime(model_slingshot_pca,pseudotime = calculate_pseudotime(model_slingshot_pca))

model_slingshot_root_pca <- add_pseudotime(model_slingshot_root_pca,
                                           pseudotime = calculate_pseudotime(model_slingshot_root_pca))

model_slingshot_root_cluster_pca <- add_pseudotime(model_slingshot_root_cluster_pca,
                                                   pseudotime = calculate_pseudotime(model_slingshot_root_cluster_pca))

# Plot the pseudotime

plot_dimred(
  model_slingshot_pca, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_slingshot_root_pca, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_slingshot_root_cluster_pca, "pseudotime",
  label_milestones = FALSE)

















### 2.-Paga_tree

# Infer trajectory

model_pagatree <- infer_trajectory(forebrain_dataset, methods_selected[2], give_priors = c("start_id"))

model_pagatree_cluster <- infer_trajectory(forebrain_dataset, methods_selected[2], 
                                           give_priors = c("start_id", "groups_id"))


## Plotting the trajectory
# Plot the normal graph

plot_graph(model_pagatree, expression_source = forebrain_dataset$expression,
           grouping = forebrain_dataset$grouping, color_cells = "grouping")

plot_graph(model_pagatree_cluster, expression_source = forebrain_dataset$expression,
           grouping = forebrain_dataset$grouping, color_cells = "grouping")


# Add the matrix as a dimred

model_pagatree_pca <- add_dimred(model_pagatree, dimred = pca_matrix, 
                                 expression_source = forebrain_dataset$expression)

model_pagatree_cluster_pca <- add_dimred(model_pagatree_cluster, 
                                         dimred = pca_matrix, expression_source = forebrain_dataset$expression)

# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_pagatree_pca, 
  expression_source = forebrain_dataset$expression, 
  grouping = forebrain_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=forebrain_color_palette)

plot_dimred(
  model_pagatree_cluster_pca, 
  expression_source = forebrain_dataset$expression, 
  grouping = forebrain_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=forebrain_color_palette)


## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_pagatree_pca <- add_root(model_pagatree_pca,
                               root_cell_id = forebrain_dataset$prior_information$start_id)

model_pagatree_cluster_pca <- add_root(model_pagatree_cluster_pca,
                                       root_cell_id = forebrain_dataset$prior_information$start_id)

## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_pagatree_pca,
  expression_source = forebrain_dataset$expression,
  grouping = forebrain_dataset$grouping,
  features_oi = 50
)

plot_heatmap(
  model_pagatree_cluster_pca,
  expression_source = forebrain_dataset$expression,
  grouping = forebrain_dataset$grouping,
  features_oi = 50
)


## Calculate pseudotime --> Uses 8 as root (idk why, because in the previous case it used a good one)

model_pagatree_pca <- add_pseudotime(model_pagatree_pca,pseudotime = calculate_pseudotime(model_pagatree_pca))

model_pagatree_cluster_pca <- add_pseudotime(model_pagatree_cluster_pca,
                                             pseudotime = calculate_pseudotime(model_pagatree_cluster_pca))


# Plot the pseudotime

plot_dimred(
  model_pagatree_pca, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_pagatree_cluster_pca, "pseudotime",
  label_milestones = FALSE)














### 3.-Scorpius

# Infer trajectory

model_scorpius <- infer_trajectory(forebrain_dataset, methods_selected[3])

model_scorpius_root <- infer_trajectory(forebrain_dataset, methods_selected[3], give_priors = c("start_id"))

model_scorpius_root_cluster <- infer_trajectory(forebrain_dataset, methods_selected[3], 
                                                give_priors = c("start_id", "groups_id"))

## Plotting the trajectory
# Plot the normal graph

plot_graph(model_scorpius, expression_source = forebrain_dataset$expression,
           grouping = forebrain_dataset$grouping, color_cells = "grouping")

plot_graph(model_scorpius_root, expression_source = forebrain_dataset$expression,
           grouping = forebrain_dataset$grouping, color_cells = "grouping")

plot_graph(model_scorpius_root_cluster, expression_source = forebrain_dataset$expression,
           grouping = forebrain_dataset$grouping, color_cells = "grouping")

# Add the matrix as a dimred

model_scorpius_pca <- add_dimred(model_scorpius, dimred = pca_matrix, 
                                 expression_source = forebrain_dataset$expression)

model_scorpius_root_pca <- add_dimred(model_scorpius_root, dimred = pca_matrix, 
                                      expression_source = forebrain_dataset$expression)

model_scorpius_root_cluster_pca <- add_dimred(model_scorpius_root_cluster, dimred = pca_matrix, 
                                              expression_source = forebrain_dataset$expression)

# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_scorpius_pca, 
  expression_source = forebrain_dataset$expression, 
  grouping = forebrain_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=forebrain_color_palette)

plot_dimred(
  model_scorpius_root_pca, 
  expression_source = forebrain_dataset$expression, 
  grouping = forebrain_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=forebrain_color_palette)

plot_dimred(
  model_scorpius_root_cluster_pca, 
  expression_source = forebrain_dataset$expression, 
  grouping = forebrain_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=forebrain_color_palette)

## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_scorpius_pca <- add_root(model_scorpius_pca,
                               root_cell_id = forebrain_dataset$prior_information$start_id)

model_scorpius_root_pca <- add_root(model_scorpius_root_pca,
                                    root_cell_id = forebrain_dataset$prior_information$start_id)

model_scorpius_root_cluster_pca <- add_root(model_scorpius_root_cluster_pca,
                                            root_cell_id = forebrain_dataset$prior_information$start_id)


## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_scorpius_pca,
  expression_source = forebrain_dataset$expression,
  grouping = forebrain_dataset$grouping,
  features_oi = 50
)

plot_heatmap(
  model_scorpius_root_pca,
  expression_source = forebrain_dataset$expression,
  grouping = forebrain_dataset$grouping,
  features_oi = 50
)

plot_heatmap(
  model_scorpius_root_cluster_pca,
  expression_source = forebrain_dataset$expression,
  grouping = forebrain_dataset$grouping,
  features_oi = 50
)


## Calculate pseudotime --> uses milestone_begin as root

model_scorpius_pca <- add_pseudotime(model_scorpius_pca,pseudotime = calculate_pseudotime(model_scorpius_pca))

model_scorpius_root_pca <- add_pseudotime(model_scorpius_root_pca,
                                          pseudotime = calculate_pseudotime(model_scorpius_root_pca))

model_scorpius_root_cluster_pca <- add_pseudotime(model_scorpius_root_cluster_pca,
                                                  pseudotime = calculate_pseudotime(model_scorpius_root_cluster_pca))

# Plot the pseudotime 

plot_dimred(
  model_scorpius_pca, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_scorpius_root_pca, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_scorpius_root_cluster_pca, "pseudotime",
  label_milestones = FALSE)

















### 4.-PAGA

# Infer trajectory

model_paga <- infer_trajectory(forebrain_dataset, "paga", give_priors = c("start_id"))

model_paga_cluster <- infer_trajectory(forebrain_dataset, "paga", give_priors = c("start_id", "groups_id"))


## Plotting the trajectory
# Plot the normal graph

plot_graph(model_paga, expression_source = forebrain_dataset$expression,
           grouping = forebrain_dataset$grouping, color_cells = "grouping")

plot_graph(model_paga_cluster, expression_source = forebrain_dataset$expression,
           grouping = forebrain_dataset$grouping, color_cells = "grouping")


# Add the matrix as a dimred

model_paga_pca <- add_dimred(model_paga, dimred = pca_matrix, expression_source = forebrain_dataset$expression)

model_paga_cluster_pca <- add_dimred(model_paga_cluster, dimred = pca_matrix, 
                                     expression_source = forebrain_dataset$expression)

# Combining a dimensionality reduction, a trajectory model and a cell clustering

plot_dimred(
  model_paga_pca, 
  expression_source = forebrain_dataset$expression, 
  grouping = forebrain_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=forebrain_color_palette)

plot_dimred(
  model_paga_cluster_pca, 
  expression_source = forebrain_dataset$expression, 
  grouping = forebrain_dataset$grouping, 
  label_milestones = FALSE)+ scale_color_manual(values=forebrain_color_palette)

## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_paga_pca <- add_root(model_paga_pca, root_cell_id = forebrain_dataset$prior_information$start_id)

model_paga_cluster_pca <- add_root(model_paga_cluster_pca, 
                                   root_cell_id = forebrain_dataset$prior_information$start_id)


## Predicting and visualising genes of interest
# A global overview of the most predictive genes

plot_heatmap(
  model_paga_pca,
  expression_source = forebrain_dataset$expression,
  grouping = forebrain_dataset$grouping,
  features_oi = 50
)


plot_heatmap(
  model_paga_cluster_pca,
  expression_source = forebrain_dataset$expression,
  grouping = forebrain_dataset$grouping,
  features_oi = 50
)

## Calculate pseudotime --> uses milestone_begin as root

model_paga_pca <- add_pseudotime(model_paga_pca, pseudotime = calculate_pseudotime(model_paga_pca))

model_paga_cluster_pca <- add_pseudotime(model_paga_cluster_pca, 
                                         pseudotime = calculate_pseudotime(model_paga_cluster_pca))

# Plot the pseudotime 

plot_dimred(
  model_paga_pca, "pseudotime",
  label_milestones = FALSE)

plot_dimred(
  model_paga_cluster_pca, "pseudotime",
  label_milestones = FALSE)

