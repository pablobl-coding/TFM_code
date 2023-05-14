# DJC Add other libraries that Pablo is using
# .libPaths(c('/home/pbeltranl/R/x86_64-pc-linux-gnu-library/4.1/', tail(.libPaths(), 1)))

# Working directory
# setwd("/home/djimenez/NetVolumes/U_BIOINFORMATICA/LABS/U_Bioinformatica/UNIT_MEMBERS/pbeltranl/TFM/Datasets' analysis/6. Ivan_dataset/for_local")




# Here, were going to analyse the neutrophils dataset using trajectory inference methods

#neutrophils_seurat_object <- readRDS("Ivan_Seurat_Obj_no_mutants_PBL.rds")

# Libraries needed
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(SeuratObject)
library(loomR)
library(dyno)
library(tidyverse)
library(cowplot)

#options(dynwrap_backend = c('container'), backend = 'singularity')
# options(dynwrap_backend = c('r_wrapper'))
options(dynwrap_backend = c('container'))
#dynwrap::set_default_config(dynwrap::create_singularity_config(cache_dir = ".singularity/"))






############################################### SEURAT ########################################################


## In this case, we don't have to load the expression matrices, since we are provided
## our own Seurat object with its own data and counts matrices. Here, we also have
## access to the the proper clustering and dimensionality reduction

#neutrophils_seurat_object <- readRDS("data/Ivan_Seurat_Obj_no_mutants_PBL.rds")
neutrophils_seurat_object <- readRDS("Ivan_Seurat_Obj_no_mutants_PBL.rds")


#### DJC - FILTER OUT GENES not expressed in minimum 20 cells
min.value <- 1 #0
min.cells <- 20

num.cells <- Matrix::rowSums(neutrophils_seurat_object@assays[[neutrophils_seurat_object@active.assay]]@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])
#object@data <- object@data[genes.use, ]

print(paste0("Orginal Genes: ", length(num.cells)))
print(paste0("Keep Genes: ", length(genes.use)))

neutrophils_seurat_object <- subset(neutrophils_seurat_object ,features=genes.use)

DefaultAssay(neutrophils_seurat_object) <- "RNA"
neutrophils_seurat_object <- FindVariableFeatures(neutrophils_seurat_object,nfeatures = 2000)
neutrophils_seurat_object <- subset(neutrophils_seurat_object ,features=VariableFeatures(neutrophils_seurat_object))

print(paste0("Final Variable Genes: ", nrow(neutrophils_seurat_object)))

#### DJC 




## This Seurat object has multiple assays, but we are only going to be working with
## one of them, specifically the "RNA" assay.
DefaultAssay(neutrophils_seurat_object) <- "RNA"
DefaultAssay(neutrophils_seurat_object)

# Both matrices are the same, which they shouldn't be
neutrophils_seurat_object@assays$RNA@data[1:10,1:10]
neutrophils_seurat_object@assays$RNA@counts[1:10,1:10]

# Normalize the counts matrix to obtain the expression matrix
neutrophils_seurat_object <- NormalizeData(neutrophils_seurat_object)

# Now, data has the normalized counts matrix (which corresponds to the expression matrix) and counts has the counts matrix
neutrophils_seurat_object@assays$RNA@data[1:10,1:10]
neutrophils_seurat_object@assays$RNA@counts[1:10,1:10]

#neutrophils_seurat_object <- FindVariableFeatures(neutrophils_seurat_object)






############################################### DYNVERSE ########################################################

# Once we have the Seurat object with the counts and expression matrices, we're going to start with dynverse 
# in order to use the other inference methods

# It's important that we remember that Seurat has the genes (features) as rows and the cells as columns, whereas
# dynverse follows the opposite rule (genes=columns; cells=rows), so that's why we have to transpose the matrices

neutrophils_dataset <- wrap_expression(
  counts = Matrix::t(neutrophils_seurat_object@assays$RNA@counts),
  expression = Matrix::t(neutrophils_seurat_object@assays$RNA@data)
)


# Selecting the best methods for our dataset

#guidelines <- guidelines_shiny(neutrophils_dataset)
#methods_selected <- guidelines$methods_selected  
methods_selected <- c("slingshot", "paga_tree", "scorpius", "paga")

# We're also going to be adding our clusters

clusters_v1 <- data.frame(neutrophils_seurat_object$integrated_snn_res.0.2)
clusters_v1$cell_id <- rownames(clusters_v1)
clusters_v1 <- clusters_v1[,c(2,1)]
#clusters <- data.frame(clusters_v1, row.names = c(1:109446))
clusters <- clusters_v1
rownames(clusters) <- NULL
colnames(clusters) <- c("cell_id", "group_id")
neutrophils_dataset <- add_grouping(neutrophils_dataset, clusters)

# Color palette

#neutrophils_color_palette <- c(1= "#d163e6",2= "#008cf9",3= "#ff9287",4= "#878500",
#                                5= "#ebac23",6= "#006e00",7= "#00a76c", 8= "#b80058",
#                                9= "#00c6f8",10= "#5954d6",11= "#b24502",12= "#00bbad")

neutrophils_color_palette <- c("#d163e6","#008cf9","#ff9287","#878500","#ebac23","#006e00","#00a76c","#b80058","#00c6f8","#5954d6","#b24502","#00bbad")
names(neutrophils_color_palette) <- as.character(seq(from=0, to=11))


# We're going to use our own dimensionality reduction 
# In this case, we also have the dimensionality reduction matrix ready

umap_matrix <- neutrophils_seurat_object@reductions[["umap"]]@cell.embeddings

# Since we are going to be using paga_tree (and PAGA), we need to add a start cell. 
# Also add grouping and dimred as prior information. 

neutrophils_dataset <- add_prior_information(
  neutrophils_dataset,
  start_id = "155951-6",
  groups_id = neutrophils_dataset$grouping
)



# REMOVE SEURAT OBJ TO FREE MEM
rm(neutrophils_seurat_object)










############################################### Running the methods ############################################### 

### 3.-Scorpius

# Infer trajectory

model_scorpius <- infer_trajectory(neutrophils_dataset, methods_selected[3])

model_scorpius_root <- infer_trajectory(neutrophils_dataset, methods_selected[3], give_priors = c("start_id"))

model_scorpius_root_cluster <- infer_trajectory(neutrophils_dataset, methods_selected[3], 
                                                give_priors = c("start_id", "groups_id"))


save.image("neutrophils_dataset_further_analysed_scorpius_pre0.RData")

## Plotting the trajectory
# Plot the normal graph

plot_graph(model_scorpius, expression_source = neutrophils_dataset$expression,
           grouping = neutrophils_dataset$grouping, color_cells = "grouping")

plot_graph(model_scorpius_root, expression_source = neutrophils_dataset$expression,
           grouping = neutrophils_dataset$grouping, color_cells = "grouping")

plot_graph(model_scorpius_root_cluster, expression_source = neutrophils_dataset$expression,
           grouping = neutrophils_dataset$grouping, color_cells = "grouping")

# Add the matrix as a dimred

model_scorpius_umap <- add_dimred(model_scorpius, dimred = umap_matrix, 
                                  expression_source = neutrophils_dataset$expression)

model_scorpius_root_umap <- add_dimred(model_scorpius_root, dimred = umap_matrix, 
                                       expression_source = neutrophils_dataset$expression)

model_scorpius_root_cluster_umap <- add_dimred(model_scorpius_root_cluster, dimred = umap_matrix, 
                                               expression_source = neutrophils_dataset$expression)

save.image("neutrophils_dataset_further_analysed_scorpius_pre1.RData")

# Combining a dimensionality reduction, a trajectory model and a cell clustering

png(filename = "neutrophils_scorpius_umap.png", width = 1200, height = 800)
plot_dimred(
  model_scorpius_umap, 
  expression_source = neutrophils_dataset$expression, 
  grouping = neutrophils_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=neutrophils_color_palette)
dev.off()

png(filename = "neutrophils_scorpius_root_umap.png", width = 1200, height = 800)
plot_dimred(
  model_scorpius_root_umap, 
  expression_source = neutrophils_dataset$expression, 
  grouping = neutrophils_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=neutrophils_color_palette)
dev.off()

png(filename = "neutrophils_scorpius_root_cluster_umap.png", width = 1200, height = 800)
plot_dimred(
  model_scorpius_root_cluster_umap, 
  expression_source = neutrophils_dataset$expression, 
  grouping = neutrophils_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=neutrophils_color_palette)
dev.off()

## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_scorpius_umap <- add_root(model_scorpius_umap,
                                root_cell_id = neutrophils_dataset$prior_information$start_id)

model_scorpius_root_umap <- add_root(model_scorpius_root_umap,
                                     root_cell_id = neutrophils_dataset$prior_information$start_id)

model_scorpius_root_cluster_umap <- add_root(model_scorpius_root_cluster_umap,
                                             root_cell_id = neutrophils_dataset$prior_information$start_id)

save.image("neutrophils_dataset_further_analysed_scorpius_pre2.RData")

## Predicting and visualising genes of interest
# A global overview of the most predictive genes

png(filename = "neutrophils_scorpius_umap_heatmap.png", width = 1200, height = 800)
plot_heatmap(
  model_scorpius_umap,
  expression_source = neutrophils_dataset$expression,
  grouping = neutrophils_dataset$grouping,
  features_oi = 50
)
dev.off()


png(filename = "neutrophils_scorpius_root_umap_heatmap.png", width = 1200, height = 800)
plot_heatmap(
  model_scorpius_root_umap,
  expression_source = neutrophils_dataset$expression,
  grouping = neutrophils_dataset$grouping,
  features_oi = 50
)
dev.off()

png(filename = "neutrophils_scorpius_root_cluster_umap_heatmap.png", width = 1200, height = 800)
plot_heatmap(
  model_scorpius_root_cluster_umap,
  expression_source = neutrophils_dataset$expression,
  grouping = neutrophils_dataset$grouping,
  features_oi = 50
)
dev.off()

## Calculate pseudotime --> uses milestone_begin as root

model_scorpius_umap <- add_pseudotime(model_scorpius_umap,pseudotime = calculate_pseudotime(model_scorpius_umap))

model_scorpius_root_umap <- add_pseudotime(model_scorpius_root_umap,
                                           pseudotime = calculate_pseudotime(model_scorpius_root_umap))

model_scorpius_root_cluster_umap <- add_pseudotime(model_scorpius_root_cluster_umap,
                                                   pseudotime = calculate_pseudotime(model_scorpius_root_cluster_umap))


save.image("neutrophils_dataset_further_analysed_scorpius_pre3.RData")

# Plot the pseudotime 

png(filename = "neutrophils_scorpius_umap_pseudotime.png", width = 1200, height = 800)
plot_dimred(
  model_scorpius_umap, "pseudotime",
  label_milestones = FALSE)
dev.off()

png(filename = "neutrophils_scorpius_root_umap_pseudotime.png", width = 1200, height = 800)
plot_dimred(
  model_scorpius_root_umap, "pseudotime",
  label_milestones = FALSE)
dev.off()

png(filename = "neutrophils_scorpius_root_cluster_umap_pseudotime.png", width = 1200, height = 800)
plot_dimred(
  model_scorpius_root_cluster_umap, "pseudotime",
  label_milestones = FALSE)
dev.off()



save.image("neutrophils_dataset_further_analysed_scorpius.RData")














### 1.-Slingshot

# Infer trajectory (we're going to do three different approaches, one where we don't give
# any previous information, one when we add the root cell and one where we add the root cell and the clustering)

model_slingshot <- infer_trajectory(neutrophils_dataset, methods_selected[1])

model_slingshot_root <- infer_trajectory(neutrophils_dataset, methods_selected[1], give_priors = c("start_id"))

model_slingshot_root_cluster <- infer_trajectory(neutrophils_dataset, methods_selected[1], 
                                                 give_priors = c("start_id", "groups_id"))


save.image("neutrophils_dataset_further_analysed_scorpiusslingshot_pre0.RData")

## Plotting the trajectory
# Plot the normal graph

plot_graph(model_slingshot, expression_source = neutrophils_dataset$expression,
           grouping = neutrophils_dataset$grouping, color_cells = "grouping")

plot_graph(model_slingshot_root, expression_source = neutrophils_dataset$expression,
           grouping = neutrophils_dataset$grouping, color_cells = "grouping")

plot_graph(model_slingshot_root_cluster, expression_source = neutrophils_dataset$expression,
           grouping = neutrophils_dataset$grouping, color_cells = "grouping")


# Add the matrix as a dimred

model_slingshot_umap <- add_dimred(model_slingshot, dimred = umap_matrix, 
                                   expression_source = neutrophils_dataset$expression)

model_slingshot_root_umap <- add_dimred(model_slingshot_root, dimred = umap_matrix, 
                                        expression_source = neutrophils_dataset$expression)

model_slingshot_root_cluster_umap <- add_dimred(model_slingshot_root_cluster, dimred = umap_matrix, 
                                                expression_source = neutrophils_dataset$expression)

save.image("neutrophils_dataset_further_analysed_scorpiusslingshot_pre1.RData")

# Combining a dimensionality reduction, a trajectory model and a cell clustering

png(filename = "neutrophils_slingshot_umap.png", width = 1200, height = 800)
plot_dimred(
  model_slingshot_umap, 
  expression_source = neutrophils_dataset$expression, 
  grouping = neutrophils_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=neutrophils_color_palette)
dev.off()

png(filename = "neutrophils_slingshot_root_umap.png", width = 1200, height = 800)
plot_dimred(
  model_slingshot_root_umap, 
  expression_source = neutrophils_dataset$expression, 
  grouping = neutrophils_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=neutrophils_color_palette)
dev.off()

png(filename = "neutrophils_slingshot_root_cluster_umap.png", width = 1200, height = 800)
plot_dimred(
  model_slingshot_root_cluster_umap, 
  expression_source = neutrophils_dataset$expression, 
  grouping = neutrophils_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=neutrophils_color_palette)
dev.off()


## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_slingshot_umap <- add_root(model_slingshot_umap,
                                 root_cell_id = neutrophils_dataset$prior_information$start_id)

model_slingshot_root_umap <- add_root(model_slingshot_root_umap,
                                      root_cell_id = neutrophils_dataset$prior_information$start_id)

model_slingshot_root_cluster_umap <- add_root(model_slingshot_root_cluster_umap,
                                              root_cell_id = neutrophils_dataset$prior_information$start_id)

save.image("neutrophils_dataset_further_analysed_scorpiusslingshot_pre2.RData")

## Predicting and visualising genes of interest
# A global overview of the most predictive genes

png(filename = "neutrophils_slingshot_umap_heatmap.png", width = 1200, height = 800)
plot_heatmap(
  model_slingshot_umap,
  expression_source = neutrophils_dataset$expression,
  grouping = neutrophils_dataset$grouping,
  features_oi = 50)
dev.off()

png(filename = "neutrophils_slingshot_root_umap_heatmap.png", width = 1200, height = 800)
plot_heatmap(
  model_slingshot_root_umap,
  expression_source = neutrophils_dataset$expression,
  grouping = neutrophils_dataset$grouping,
  features_oi = 50) 
dev.off()

png(filename = "neutrophils_slingshot_root_cluster_umap_heatmap.png", width = 1200, height = 800)
plot_heatmap(
  model_slingshot_root_cluster_umap,
  expression_source = neutrophils_dataset$expression,
  grouping = neutrophils_dataset$grouping,
  features_oi = 50) 
dev.off()

## Calculate pseudotime

model_slingshot_umap <- add_pseudotime(model_slingshot_umap,pseudotime = calculate_pseudotime(model_slingshot_umap))

model_slingshot_root_umap <- add_pseudotime(model_slingshot_root_umap,
                                            pseudotime = calculate_pseudotime(model_slingshot_root_umap))

model_slingshot_root_cluster_umap <- add_pseudotime(model_slingshot_root_cluster_umap,
                                                    pseudotime = calculate_pseudotime(model_slingshot_root_cluster_umap))


save.image("neutrophils_dataset_further_analysed_scorpiusslingshot_pre3.RData")

# Plot the pseudotime

png(filename = "neutrophils_slingshot_umap_pseudotime.png", width = 1200, height = 800)
plot_dimred(
  model_slingshot_umap, "pseudotime",
  label_milestones = FALSE)
dev.off()

png(filename = "neutrophils_slingshot_root_umap_pseudotime.png", width = 1200, height = 800)
plot_dimred(
  model_slingshot_root_umap, "pseudotime",
  label_milestones = FALSE)
dev.off()

png(filename = "neutrophils_slingshot_root_cluster_umap_pseudotime.png", width = 1200, height = 800)
plot_dimred(
  model_slingshot_root_cluster_umap, "pseudotime",
  label_milestones = FALSE)
dev.off()




save.image("neutrophils_dataset_further_analysed_scorpiusslingshot.RData")












### 2.-Paga_tree

# Infer trajectory

model_pagatree <- infer_trajectory(neutrophils_dataset, methods_selected[2], give_priors = c("start_id"))

model_pagatree_cluster <- infer_trajectory(neutrophils_dataset, methods_selected[2], 
                                           give_priors = c("start_id", "groups_id"))


save.image("neutrophils_dataset_further_analysed_scorpiusslingshotpagatree_pre0.RData")


## Plotting the trajectory
# Plot the normal graph

plot_graph(model_pagatree, expression_source = neutrophils_dataset$expression,
           grouping = neutrophils_dataset$grouping, color_cells = "grouping")

plot_graph(model_pagatree_cluster, expression_source = neutrophils_dataset$expression,
           grouping = neutrophils_dataset$grouping, color_cells = "grouping")


# Add the matrix as a dimred

model_pagatree_umap <- add_dimred(model_pagatree, dimred = umap_matrix, 
                                  expression_source = neutrophils_dataset$expression)

model_pagatree_cluster_umap <- add_dimred(model_pagatree_cluster, 
                                          dimred = umap_matrix, expression_source = neutrophils_dataset$expression)

save.image("neutrophils_dataset_further_analysed_scorpiusslingshotpagatree_pre1.RData")


# Combining a dimensionality reduction, a trajectory model and a cell clustering

png(filename = "neutrophils_pagatree_umap.png", width = 1200, height = 800)
plot_dimred(
  model_pagatree_umap, 
  expression_source = neutrophils_dataset$expression, 
  grouping = neutrophils_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=neutrophils_color_palette)
dev.off()

png(filename = "neutrophils_pagatree_cluster_umap.png", width = 1200, height = 800)
plot_dimred(
  model_pagatree_cluster_umap, 
  expression_source = neutrophils_dataset$expression, 
  grouping = neutrophils_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=neutrophils_color_palette)
dev.off()


## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_pagatree_umap <- add_root(model_pagatree_umap,
                                root_cell_id = neutrophils_dataset$prior_information$start_id)

model_pagatree_cluster_umap <- add_root(model_pagatree_cluster_umap,
                                        root_cell_id = neutrophils_dataset$prior_information$start_id)

save.image("neutrophils_dataset_further_analysed_scorpiusslingshotpagatree_pre2.RData")


## Predicting and visualising genes of interest
# A global overview of the most predictive genes

png(filename = "neutrophils_pagatree_umap_heatmap.png", width = 1200, height = 800)
plot_heatmap(
  model_pagatree_umap,
  expression_source = neutrophils_dataset$expression,
  grouping = neutrophils_dataset$grouping,
  features_oi = 50
)
dev.off()

png(filename = "neutrophils_pagatree_cluster_umap_heatmap.png", width = 1200, height = 800)
plot_heatmap(
  model_pagatree_cluster_umap,
  expression_source = neutrophils_dataset$expression,
  grouping = neutrophils_dataset$grouping,
  features_oi = 50
)
dev.off()

## Calculate pseudotime 

model_pagatree_umap <- add_pseudotime(model_pagatree_umap,pseudotime = calculate_pseudotime(model_pagatree_umap))

model_pagatree_cluster_umap <- add_pseudotime(model_pagatree_cluster_umap,
                                              pseudotime = calculate_pseudotime(model_pagatree_cluster_umap))


save.image("neutrophils_dataset_further_analysed_scorpiusslingshotpagatree_pre3.RData")


# Plot the pseudotime

png(filename = "neutrophils_pagatree_umap_pseudotime.png", width = 1200, height = 800)
plot_dimred(
  model_pagatree_umap, "pseudotime",
  label_milestones = FALSE)
dev.off()

png(filename = "neutrophils_pagatree_cluster_umap_pseudotime.png", width = 1200, height = 800)
plot_dimred(
  model_pagatree_cluster_umap, "pseudotime",
  label_milestones = FALSE)
dev.off()





save.image("neutrophils_dataset_further_analysed_scorpiusslingshotpagatree.RData")
















### 4.-PAGA

# Infer trajectory

model_paga <- infer_trajectory(neutrophils_dataset, methods_selected[4], give_priors = c("start_id"))

model_paga_cluster <- infer_trajectory(neutrophils_dataset, methods_selected[4], give_priors = c("start_id", "groups_id"))


save.image("neutrophils_dataset_further_analysed_scorpiusslingshotpagatreepaga_pre0.RData")

## Plotting the trajectory
# Plot the normal graph

plot_graph(model_paga, expression_source = neutrophils_dataset$expression,
           grouping = neutrophils_dataset$grouping, color_cells = "grouping")

plot_graph(model_paga_cluster, expression_source = neutrophils_dataset$expression,
           grouping = neutrophils_dataset$grouping, color_cells = "grouping")


# Add the matrix as a dimred

model_paga_umap <- add_dimred(model_paga, dimred = umap_matrix, expression_source = neutrophils_dataset$expression)

model_paga_cluster_umap <- add_dimred(model_paga_cluster, dimred = umap_matrix, 
                                      expression_source = neutrophils_dataset$expression)


save.image("neutrophils_dataset_further_analysed_scorpiusslingshotpagatreepaga_pre1.RData")


# Combining a dimensionality reduction, a trajectory model and a cell clustering

png(filename = "neutrophils_paga_umap.png", width = 1200, height = 800)
plot_dimred(
  model_paga_umap, 
  expression_source = neutrophils_dataset$expression, 
  grouping = neutrophils_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=neutrophils_color_palette)
dev.off()

png(filename = "neutrophils_paga_cluster_umap.png", width = 1200, height = 800)
plot_dimred(
  model_paga_cluster_umap, 
  expression_source = neutrophils_dataset$expression, 
  grouping = neutrophils_dataset$grouping, 
  label_milestones = FALSE)+ scale_fill_manual(values=neutrophils_color_palette)
dev.off()

## For the next steps (heatmap and pseudotime), we're going to add a root cell, which will be the same one
## we've used for paga_tree and paga

model_paga_umap <- add_root(model_paga_umap, root_cell_id = neutrophils_dataset$prior_information$start_id)

model_paga_cluster_umap <- add_root(model_paga_cluster_umap, 
                                    root_cell_id = neutrophils_dataset$prior_information$start_id)


save.image("neutrophils_dataset_further_analysed_scorpiusslingshotpagatreepaga_pre2.RData")



## Predicting and visualising genes of interest
# A global overview of the most predictive genes

png(filename = "neutrophils_paga_umap_heatmap.png", width = 1200, height = 800)
plot_heatmap(
  model_paga_umap,
  expression_source = neutrophils_dataset$expression,
  grouping = neutrophils_dataset$grouping,
  features_oi = 50
)
dev.off()

png(filename = "neutrophils_paga_cluster_umap_heatmap.png", width = 1200, height = 800)
plot_heatmap(
  model_paga_cluster_umap,
  expression_source = neutrophils_dataset$expression,
  grouping = neutrophils_dataset$grouping,
  features_oi = 50
)
dev.off()

## Calculate pseudotime --> uses milestone_begin as root

model_paga_umap <- add_pseudotime(model_paga_umap, pseudotime = calculate_pseudotime(model_paga_umap))

model_paga_cluster_umap <- add_pseudotime(model_paga_cluster_umap, 
                                          pseudotime = calculate_pseudotime(model_paga_cluster_umap))

save.image("neutrophils_dataset_further_analysed_scorpiusslingshotpagatreepaga_pre3.RData")


# Plot the pseudotime 

png(filename = "neutrophils_paga_umap_pseudotime.png", width = 1200, height = 800)
plot_dimred(
  model_paga_umap, "pseudotime",
  label_milestones = FALSE)
dev.off()

png(filename = "neutrophils_paga_cluster_umap_pseudotime.png", width = 1200, height = 800)
plot_dimred(
  model_paga_cluster_umap, "pseudotime",
  label_milestones = FALSE)
dev.off()

save.image("neutrophils_dataset_further_analysed_scorpiusslingshotpagatreepaga.RData")

save.image("neutrophils_dataset_further_analysed.RData")
