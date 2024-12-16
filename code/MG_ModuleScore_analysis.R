# module_score vs trajectory(pseudotime graphs)

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
# library(pheatmap)
# library(apeglm)
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
################################################################################



#extract pseudotime from cds object and added to the metadata of MG seurat object

MG@meta.data$pseudotime<-pseudotime(cds)

#There are inf pseudotimes in the column which has to be removed before binning them
#remove infinite pseudotime

is.na(MG@meta.data)<-sapply(MG@meta.data,is.infinite)
any(is.na(MG@meta.data))

#checking the total number of infinite pseudotime in a column(currently in pseudotime column)
colSums(is.na(MG@meta.data))


#To add a new column named pseudotime_bin with 25 partitions
MG@meta.data$pseudotime_bin<-cut(as.numeric(MG@meta.data$pseudotime),25)




# #splitObj is a list with 10 seurat objects ------------------------------


splitObj<-SplitObject(MG,split.by = "group")


# #remove the NA values from the splitObj ---------------------------------



if(any(is.na(splitObj))){
  na.omit(splitObj)
  
}








############################################################
#############################################################################


homeostatic <- c('Hexb', 'Cst', 'Cx3cr1', 'Csf1r', 'Ctss', 'Tmsb4x', 'P2ry12', 'C1qb')
stage1_DAM <- c('Tyrobp', 'Ctsb', 'Apoe', 'B2m', 'Fth1', 'Trem2')
stage2_DAM <- c('Trem2', 'Axl', 'Cst7', 'Ctsl', 'Lpl', 'Cd9', 'Csf1', 'Itgax', 'Clec7a', 'Lilrb4', 'Timp2')





df<-data.frame()
for(i in 1:length(splitObj)){
  SeuObj<-splitObj[[i]]
  
  SeuObj<- AddModuleScore(SeuObj,
                          features=list('stage1_DAM'=stage1_DAM, 'stage2_DAM'=stage2_DAM, 'homeostatic'=homeostatic),
                          pool = rownames(SeuObj), k=F, nbin=24,
                          name=c('stage1_DAM', 'stage2_DAM', 'homeostatic'))
  
  DAM_modules <- select(SeuObj@meta.data, c(pseudotime_bin, stage1_DAM1, stage2_DAM2, homeostatic3,group))
  DAM_modules$pseudotime_bin_num <- as.numeric(DAM_modules$pseudotime_bin)
  DAM_modules<-na.omit(DAM_modules)
  
  df<-bind_rows(df,DAM_modules)
  
}



# Prepare for ggplot ------------------------------------------------------



df<-df %>% pivot_longer( cols = c(2:4), names_to = "features", values_to = "value")

df<-df %>% group_by(group, pseudotime_bin_num,features) %>%
  summarise(value=mean(value))

df$group<- as.character(df$group)
df<-df %>% mutate(Geno=ifelse(startsWith(group,"WT"),"WT","Kbk"))

df$Geno<-factor(df$Geno,levels=c("WT","Kbk"))


#to change the group variable to factor and to change the order of levels (proportion_df$group)
df$group<-factor(df$group, levels= c("WT.1m.Un", "WT.3m.Un","WT.3m.DMSO_PEG","WT.3m.DCF","WT.3m.TCF","Kbk.1m.Un",
                                     "Kbk.3m.Un","Kbk.3m.DMSO_PEG","Kbk.3m.DCF","Kbk.3m.TCF"))


df$new_group<- factor(df$group,labels = c("WT.1m\n(Un)","WT.3m\n(Un)","WT.3m\n(PEG)",
                                          "WT.3m\n(HPBCD)",
                                          "WT.3m\n(TCF)",
                                          "Kbk.1m\n(Un)","Kbk.3m\n(Un)","Kbk.3m\n(PEG)",
                                          "Kbk.3m\n(HPBCD)",
                                          "Kbk.3m\n(TCF)"))


#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/results_seurat5")

write.csv(df,file="results/Trajectory_vs_pseudotime data.csv")


################################################################################

#new theme for ggplot

new_theme <- theme(
  axis.text.x = element_text(size = 15),  # Change x-axis tick label size
  axis.text.y = element_text(size = 30),  # Change y-axis tick label size
  axis.title.x = element_text(size = 60, hjust = 0),  # Adjust x-axis label size and position
  axis.title.y = element_text(size = 60),  # Adjust y-axis label size
  strip.background = element_rect(fill = "white", color = "black"),  # Modify strip background
  strip.text = element_text(size = 25, color = "Black"),  # Modify strip text size
  legend.text = element_text(size = 25),  # Increase legend text size
  legend.key.size = unit(2, "lines"),  # Increase legend key size
  legend.spacing.x = unit(1, 'cm'),  # Increase spacing between legend items
  legend.position = "bottom",  # Move legend to the bottom
  legend.direction = "horizontal",  # Arrange legend items in a row
  legend.title = element_blank()  # Remove legend title
)



################################################################################


# Visualize the modulescore vs psudotime  ---------------------------------


#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/figures_seurat5")
#to separate the groups based on their genotype(Geno)


png("figures/Module score vs pseudotime plot.png", width=2500, height=1250)
p<-ggplot(df, aes(pseudotime_bin_num, value, color=features,fill=Geno)) +
  geom_hline(yintercept=0, linetype='dashed', color='gray', alpha=0.75) +
  geom_point(size=3) +
  geom_smooth(method = "lm", level = 0.95) +
  theme_classic()+
  xlab('Trajectory') + ylab('Module score')+
  facet_grid(~new_group)


p+new_theme#+guides(fill = "none")#to remove the fill legend
dev.off()






# pdf format --------------------------------------------------------------

pdf("figures/Module score vs pseudotime plot.pdf", width=30, height=15)
p<-ggplot(df, aes(pseudotime_bin_num, value, color=features,fill=Geno)) +
  geom_hline(yintercept=0, linetype='dashed', color='gray', alpha=0.75) +
  geom_point(size=3) +
  geom_smooth(method = "lm", level = 0.95) +
  theme_classic()+
  xlab('Trajectory') + ylab('Module score')+
  facet_grid(~new_group)


p+new_theme#+guides(fill = "none")#to remove the fill legend
dev.off()

# end ---------------------------------------------------------------------






# separate the gene modules and plot the graph ----------------------------
#pearson correlation is added


################################################################################


# homeostatic genes module score vs pseudotime ----------------------------


df_homeo<-df %>% filter(features=="homeostatic3")


pdf("figures/homeostatic_genes.pdf",width=30,heigh=15)

p<-ggscatter(
  subset(df_homeo),
  x='pseudotime_bin_num', y='value', #col="red",
  color = "Geno",
  
  add='reg.line',
  add.params=list( fill='lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args=list(method='pearson', label.sep='\n', size=10),
  alpha=0.5, size=3) + xlab('Trajectory') +
  ylab("module score")+# Add axis labels
  ggtitle("Homeostatic_genes") +  # Add title
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold")  # Center the title and increase font size
  )+
  facet_grid(~new_group, scales = "free")+
  theme(panel.spacing = unit(1.5, "lines"))## Increase spacing between panels

# Adjust line colors based on groups
p+new_theme

dev.off()


################################################################################



df_stage1_DAM<-df %>% filter(features=="stage1_DAM1")


pdf("figures/stage1_DAM_genes.pdf",width=30,heigh=15)

p<-ggscatter(
  subset(df_stage1_DAM),
  x='pseudotime_bin_num', y='value', #col="red",
  color = "Geno",
  
  add='reg.line',
  #facet.by = "new_group",
  add.params=list( fill='lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args=list(method='pearson', label.sep='\n', size=10),
  alpha=0.5, size=3) + xlab('Trajectory') +
  ylab("module score")+# Add axis labels
  ggtitle("stage1_DAM_genes") +  # Add title
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold")  # Center the title and increase font size
  )+
  facet_grid(~new_group, scales = "free")+
  theme(panel.spacing = unit(1.5, "lines"))## Increase spacing between panels

# Adjust line colors based on groups
p+new_theme

dev.off()


################################################################################

df_stage2_DAM<-df %>% filter(features=="stage2_DAM2")


pdf("figures/stage2_DAM_genes.pdf",width=30,heigh=15)

p<-ggscatter(
  subset(df_stage2_DAM),
  x='pseudotime_bin_num', y='value', #col="red",
  color = "Geno",
  
  add='reg.line',
  #facet.by = "new_group",
  add.params=list( fill='lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args=list(method='pearson', label.sep='\n', size=10),
  alpha=0.5, size=3) + xlab('Trajectory') +
  ylab("module score")+# Add axis labels
  ggtitle("stage2_DAM_genes") +  # Add title
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold")  # Center the title and increase font size
  )+
  facet_grid(~new_group, scales = "free")+
  theme(panel.spacing = unit(1.5, "lines"))## Increase spacing between panels

# Adjust line colors based on groups
p+new_theme

dev.off()


################################################################################