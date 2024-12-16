
# MG Trajectory analysis -----------------------------------------------------


# load libraries ----------------------------------------------------------

library(Seurat)
library(tidyverse)
library(cowplot)
# library(Matrix.utils)
# library(magrittr)
# library(Matrix)
# library(purrr)
# library(reshape2)
# library(S4Vectors)
# library(tibble)
library(SingleCellExperiment)
library(pheatmap)
#library(apeglm)
library(png)
#library(DESeq2)
library(RColorBrewer)
library(harmony)
library(monocle3)
# library(plotly)#for interactive plotting
# library(htmlwidgets)
# library(DT)
library(ggpubr)


# load MG data ------------------------------------------------------------




MG<-readRDS("data/MG.rds")

#DimPlot(MG,label=T)
################################################################################

# Convert the MG to Cell_data_set -----------------------------------------

# attempt to make a cell dataset object:
expr_matrix <- GetAssayData(MG, slot='counts')
genes <- data.frame(as.character(rownames(expr_matrix)))
rownames(genes) <- rownames(expr_matrix)
genes <- as.data.frame(cbind(genes,genes))

#to change the column names of genes dataframe
colnames(genes) <- c("GeneSymbol", "gene_short_name")


cds <- new_cell_data_set(
  expr_matrix,
  cell_metadata=MG@meta.data,
  gene_metadata=genes
)



cds <- preprocess_cds(cds,method="PCA",num_dim=50)

plot_pc_variance_explained(cds)

#Remove batch effects with cell alignment.

cds <- align_cds(cds,  preprocess_method="PCA", alignment_group = "sample")
cds <- reduce_dimension(cds, preprocess_method="PCA")
cds <- cluster_cells(cds,reduction_method='UMAP')
colData(cds)<-subset(colData(cds),select=-sample_name)#remove the column with sample variable



###############################################################

#Order cells in pseudotime along a trajectory
## Step 5: Learn a graph


#To obtain a trajectory with loop structure, 
cds <- learn_graph(cds)





##################################################################################
#Automatic root detection in untreated samples
##################################################################################



# WT.1m.Un as the root

get_earliest_principal_node <- function(cds, time_bin="WT.1m.Un"){
  cell_ids <- which(colData(cds)[, "group"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
################################################################################



png("figures/MG_all_pseudotime.png",width=500,height=500)

#plot the graph without trajectory
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           show_trajectory_graph = F)+ggtitle("MG pseudotime")+
  theme(plot.title = element_text(hjust = 0.5))

dev.off()


################################################################################



#extract pseudotime from cds object and added to the metadata of MG seurat object

MG@meta.data$pseudotime<-pseudotime(cds)

#There are inf pseudotimes in the column which has to be removed before binning them
#remove infinite pseudotime

is.na(MG@meta.data)<-sapply(MG@meta.data,is.infinite)
any(is.na(MG@meta.data))

#checking the total number of pseudotime in a column(currently in pseudotime column)
colSums(is.na(MG@meta.data))


#To add a new column named pseudotime_bin with 25 partitions
MG@meta.data$pseudotime_bin<-cut(as.numeric(MG@meta.data$pseudotime),25)



#combining all the proportion values of bins to a single data frame --------
# to get the unique pseudotime_bin
bins <- levels(MG@meta.data$pseudotime_bin)

proportion_df <- data.frame()

for(i in 1:length(bins)){
  bin = bins[i]
  
  # Subset the meta data for the current bin
  cur_meta <- subset(MG@meta.data, pseudotime_bin == bin)
  
  # Check if the current bin is empty, skip if no data
  if(nrow(cur_meta) > 0) {
    # Create a data frame with the proportion of cells in each group
    cur_df <- data.frame(table(cur_meta$group) / nrow(cur_meta))
    
    # Rename columns to proper names
    cur_df <- cur_df %>% dplyr::rename(group = Var1, Cellproportion = Freq)
    
    # Add additional information
    cur_df$bin_num <- i
    cur_df$bin <- bin
    cur_df$Geno <- ifelse(startsWith(as.character(cur_df$group), "WT"), "WT", "Kbk")
    
    # Combine with the overall data frame
    proportion_df <- rbind(proportion_df, cur_df)
  } else {
    # Message if the bin is empty
    message("Skipping bin ", bin, " because it has no data.")
  }
}






#to change the order of levels (proportion_df$group)
proportion_df$group<-factor(proportion_df$group,
                            levels= c("WT.1m.Un","WT.3m.Un","WT.3m.DMSO_PEG","WT.3m.DCF","WT.3m.TCF","Kbk.1m.Un",
                                      "Kbk.3m.Un","Kbk.3m.DMSO_PEG","Kbk.3m.DCF","Kbk.3m.TCF"))





#add a new column(variable) "new_group" to change the name appear as in the graph 
proportion_df$new_group<- factor(proportion_df$group,
                                 labels = c("WT.1m\n(Un)","WT.3m\n(Un)","WT.3m\n(PEG)",
                                            "WT.3m\n(HPBCD)",
                                            "WT.3m\n(TCF)",
                                            "Kbk.1m\n(Un)","Kbk.3m\n(Un)","Kbk.3m\n(PEG)",
                                            "Kbk.3m\n(HPBCD)",
                                            "Kbk.3m\n(TCF)"))


#convert the "Geno" variable  to factor variable
proportion_df$Geno<-factor(proportion_df$Geno,levels=c("WT","Kbk"))




write.csv(proportion_df,file="results/cellProportion_vs_trajectory.csv")

################################################################################
#new theme for ggscatter

new_theme<- theme(
  axis.text.x = element_text(size = 15),#change the x and y axis tick label size
  axis.text.y = element_text(size = 30),
  axis.title.x = element_text(size = 60,hjust = 0),  # Change x-axis label size and position towards the origin
  axis.title.y = element_text(size = 60))+# Adjust the size of y-axis label size
  theme(strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 25, color = "Black"),legend.text = element_text(size = 25),       # Increase legend text size
        #legend.title = element_text(size = 25),      # Increase legend title size
        legend.key.size = unit(2, "lines"),        # Increase size of legend keys
        #legend.background = element_rect(fill = "white", color = "black"),  # Optional: change legend background
        legend.spacing.x = unit(1, 'cm'),# Optional: increase spacing between legend items
        legend.position = "bottom",       # Move legend to the bottom
        legend.direction = "horizontal",  # Arrange legend items in a single row
        legend.title = element_blank())
################################################################################





#change the direcotry to figures_seurat5



png("figures/cellProportion_vs_Trajectory.png",width=2500,height=1250)
p<-ggscatter(
  subset(proportion_df),
  x='bin_num', y='Cellproportion', #col="red",
  color = "Geno",
  
  add='reg.line',
  #facet.by = "new_group",
  add.params=list( fill='lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args=list(method='pearson', label.sep='\n', size=10),
  alpha=0.5, size=3) + xlab('Trajectory') +
  ylab("Cell proportion")+# Add axis labels
  facet_grid(~new_group, scales = "free")+
  theme(panel.spacing = unit(1.5, "lines"))## Increase spacing between panels

# Adjust line colors based on groups
p+new_theme

dev.off()
################################################################################


#pdf format

pdf("figures/cellProportion_vs_Trajectory.pdf",width=30,height=15)
p<-ggscatter(
  subset(proportion_df),
  x='bin_num', y='Cellproportion', #col="red",
  color = "Geno",
  
  add='reg.line',
  #facet.by = "new_group",
  add.params=list( fill='lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args=list(method='pearson', label.sep='\n', size=10),
  alpha=0.5, size=3) + xlab('Trajectory') +
  ylab("Cell proportion")+# Add axis labels
  facet_grid(~new_group, scales = "free")+
  theme(panel.spacing = unit(1.5, "lines"))## Increase spacing between panels

# Adjust line colors based on groups
p+new_theme

dev.off()


# End ---------------------------------------------------------------------


