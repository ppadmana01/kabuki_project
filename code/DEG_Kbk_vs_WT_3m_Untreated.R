
# WT vs Kbk 3m Untreated DEG and heatmaps --------------------------------------


# load libraries ----------------------------------------------------------

######################################################################

library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(harmony)
library(tidyverse)
library(readxl)
library(plotly)#for interactive plotting
library(htmlwidgets)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ComplexHeatmap)

########################################################################################



# Load MG.rds file --------------------------------------------------------




MG<-readRDS("data/MG.rds")


################################################################################


# Subset the DMSO_PEG and DCF cells from Kbk genotype --------------------------
MG_Kbk_WT_3m_Un<-subset(MG,Age=="3m" & Treatment =="Un" & Genotype %in% c("WT","Kbk"))

#reassign to another variable: seurat
seurat<-MG_Kbk_WT_3m_Un

#to check whether I correctly subset the seurat object
table(seurat$unique_name)
length(seurat$unique_name)


################################################################################



# UMAP visualization of MG_Kbk_DCF_TCF seurat_obj ---------------- --------

DimPlot(seurat,label=T)

################################################################################


# Extract raw counts and metadata to create SingleCellExperiment object --------

counts <- seurat[["RNA"]]$counts 

# Extract metadata
metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@meta.data$cell_type)#factor is not essential in this case


# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "Genotype")]

#############################################################################


# Extract celltype: microglia ---------------------------------------------


################################################################################
cluster_names <- levels(colData(sce)$cluster_id)
cluster_names
# Total number of clusters
length(cluster_names)
################################################################################

#  Extract unique names of samples ----------------------------------------



# Extract unique names of samples (= levels of sample_id factor variable)
sample_names <- factor(colData(sce)$unique_name)
sample_names<-levels(sample_names)

#Here, we identify 12 different samples
sample_names

# Total number of samples
length(sample_names)


# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("cluster_id", "unique_name")]
groups$unique_name<-factor(groups$unique_name)
length(groups$cluster_id)

################################################################################


# Generate a new metadata  ------------------------------------------------


####################################################################
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(group,unique_name,sample)

###################################################################


# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]
rownames(metadata)




# The rownames of metadata and colnames of count matrix must be identical-------

# Rename the rownames of metadata
rownames(metadata) <- metadata$unique_name
head(metadata)


#to change the order of rownames so that,
#it will be compatible with that of colnames(counts)
#to make the colnames of aggregated count and the rownames of metadata identical
metadata<-metadata %>% arrange(rownames(metadata))
metadata


t<-as.data.frame(table(colData(sce)$unique_name,
                       colData(sce)$cluster_id))


cell_number<-t$Freq
metadata$cell_number<-cell_number

t
metadata


#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/results_seurat5")
write.csv(metadata,file="results/metadata_Kbk_vs_WT_3m_Un.csv")

################################################################################




# Aggregate the count values of each unique_name(individual sample --------


##################################################
#Create a function to get the rowSums of counts of individual samples
mat_function<-function(x){
  mat<-as.matrix(GetAssayData(x, slot = "counts"))
  #mat<-mat[sig_genes, ]
  #scale the rows
  #mat<-t(scale(t(mat)))
  #get the rowMeans to combine scaled count value of all the cells
  mat<-rowSums(mat)
  #convert back to matrix
  mat<-as.matrix(mat)
  return(mat)
}
#############################################################
#unique_sample names(individual samples)
#All the individual sample containing vector
clusters=factor(unique(seurat@meta.data$unique_name))
clusters
levels(clusters)
#############################################################
#Create an empty dataframe
df<-data.frame()
#for loop to get the individual samples of the seurat_obj to get the combined matrix
for (i in 1 :length(clusters)){
  
  #all the individual sample containing vector
  cur_cluster<-factor(clusters,levels=levels(clusters))
  #extract individual samples from the vector
  cur_seurat_obj<- subset(seurat, unique_name == cur_cluster[i])
  #mat_function extract the matrix,scale the matrix and calculate the rowMeans to make a single vector 
  #get the matrix
  mat<-mat_function(cur_seurat_obj)
  #convert matrix to a data frame
  mat<-as.data.frame(mat)
  #add column names to the matrix
  colnames(mat)<-levels(clusters)[i]
  #transpose the data frame "mat" by t()
  mat<-t(mat)
  df<-rbind(df,mat)#bind the mat data frame to empty df 
  #print(df)
  
}

#dim(df)
#transpose the count matrix
df<-t(df)
colnames(df)
#class(df)
counts<-df
dim(counts)
#counts[1:10,1:12]




#to confirm the colnames(counts) is identical as rownames(metadata)
all(colnames(counts)== rownames(metadata))
################################################################################



# Analysis: The batch effect is removed -------------------

#we add sample(batch information in the following code)
# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(counts, 
                              colData = metadata, 
                              design = ~ sample+group)



##########################################################################
#memory.limit(size = 320000)

# Run DESeq2 differential expression analysis
dds<-DESeq(dds)


saveRDS(dds,"data/dds_obj_Kbk_vs_WT_3m_Un.rds")


################################################################################




# Output results of Wald test for contrast for TCF vs DCF -----------------


WT.3m.Un<-levels(factor(metadata$group))[2]#"WT.3m.Un"
Kbk.3m.Un<-levels(factor(metadata$group))[1]#"Kbk.3m.Un"


contrast <- c("group", Kbk.3m.Un, WT.3m.Un)
library(DESeq2)


res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)


summary(res)

#the number of p_value adjusted genes below 0.05  
sum(res$padj < 0.05, na.rm=TRUE)


################################################################################


# heatmap using pvalue cut off (0.05) -------------------------------------

res_tbl<-data.frame(res@listData)

#gene names stored in the rownames of res object
row.names(res_tbl)<-res@rownames

#create a  new column named "gene" in the res_tbl data frame 
res_tbl$gene<-rownames(res_tbl)




#significant genes based on pvalue cut off 
Kbk_vs_WT_3m_Un_pvalue<-res_tbl %>% filter(pvalue<0.05)


#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/results_seurat5")
write.csv(Kbk_vs_WT_3m_Un_pvalue,file="results/DEG_datasets/DESeq_DEG_Kbk_vs_WT_3m_Un_pvalue_cutoff_0.05.csv")



################################################################################



sig_genes_pvalue<-Kbk_vs_WT_3m_Un_pvalue$gene
normalized_counts <- counts(dds, 
                            normalized = TRUE)

#filter the normalized_counts matrix using significant genes
sig_counts<-normalized_counts[sig_genes_pvalue,]

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")




## Run pheatmap using the metadata data frame for the annotation
sig_counts_s<-t(scale(t(sig_counts)))

sig_counts_s[1:5,1:5]



# change the order of the columns for DMSO_PEG vs DCF order ---------------

new_order<-c(7,8,9,10,11,12,1,2,3,4,5,6)

sig_counts_s<-sig_counts_s[,new_order]

# significant genes based on padj cut off  --------------------------------


Kbk_vs_WT_3m_Un_padj<-res_tbl %>% filter(padj<0.05)
Kbk_vs_WT_3m_Un_padj


genes_of_interest<-rownames(Kbk_vs_WT_3m_Un_padj)
genes_of_interest


################################################################################



# The padj selected genes appear as asterix in the heatmap ---------


row_labels <- ifelse(rownames(sig_counts) %in% genes_of_interest, "*", "")






# Adding 1,2,3 under the column title -------------------------------------


custom_column_names<-rep(1:3,times=4)



################################################################################




# Create the heatmap ------------------------------------------------------

png("figures/Heatmap_Kbk_vs_WT_3m_Un_pvalue_3m_Un.png",width=500,height=500)
ht <- Heatmap(sig_counts_s, 
              name = "Expression",
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 15, fontfamily = "Arial"),  # Font settings for the legend title
                labels_gp = gpar(fontsize = 12, fontfamily = "Arial")  # Font settings for the legend labels
              ),
              #heatmap_legend_param = list(title_gp = gpar(fontsize = 20, fontfamily = "Arial")),  # Set font size and family for the heatmap title
              
              column_split = factor(rep(c("WT.F","WT.M","Kbk.F","Kbk.M"), each = 3), levels = c("WT.F", "WT.M", "Kbk.F", "Kbk.M")),  # Split columns into four groups
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              column_names_side = "top",
              column_title_gp = gpar(fontsize = 15, fontfamily = "Arial"),  # Font settings for column titles
              column_gap = unit(1, "mm"),
              cluster_rows = TRUE,
              show_row_dend = TRUE,
              row_title_gp = gpar(fontsize = 12),
              column_title_rot = 0,  # To make the column title horizontal
              show_heatmap_legend = TRUE,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
              use_raster = TRUE,
              raster_quality = 10,
              row_labels = row_labels,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 30),
              column_labels = rep(1:3, times = 4),  # Custom column names
              column_names_gp = gpar(fontsize = 15, fontfamily = "Arial"),  # Set font size and family for column labels
              column_names_rot = 360)  # Rotation for column labels






# Draw the heatmap
draw(ht, newpage = TRUE)

dev.off()
################################################################################
################################################################################

# heatmap based on padj genes ---------------------------------------------


normalized_counts <- counts(dds, 
                            normalized = TRUE)


# significant genes based on padj cut off  --------------------------------


Kbk_vs_WT_3m_Un_padj<-res_tbl %>% filter(padj<0.05)
Kbk_vs_WT_3m_Un_padj

sig_genes_padj<-Kbk_vs_WT_3m_Un_padj$gene
sig_genes_padj


#filter the normalized_counts matrix using significant genes
sig_counts<-normalized_counts[sig_genes_padj,]

dim(sig_counts)


# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")




## Run pheatmap using the metadata data frame for the annotation
sig_counts_s<-t(scale(t(sig_counts)))

sig_counts_s[1:5,1:5]



# change the order of the columns for DMSO_PEG vs DCF order ---------------

new_order<-c(7,8,9,10,11,12,1,2,3,4,5,6)

new_matrix<-sig_counts_s[,new_order]




genes_of_interest<-rownames(Kbk_vs_WT_3m_Un_padj)
genes_of_interest


################################################################################



# The padj selected genes appear as asterix in the heatmap ---------


row_labels <- rownames(Kbk_vs_WT_3m_Un_padj)






# Adding 1,2,3 under the column title -------------------------------------


custom_column_names<-rep(1:3,times=4)



################################################################################




png("figures/Heatmap_Kbk_vs_WT_3m_Un_padj_3m_Un_padjGenes.png",width=500,height=500)
ht <- Heatmap(new_matrix, 
              name = "Expression",
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 15, fontfamily = "Arial"),  # Font settings for the legend title
                labels_gp = gpar(fontsize = 12, fontfamily = "Arial")  # Font settings for the legend labels
              ),
              #heatmap_legend_param = list(title_gp = gpar(fontsize = 20, fontfamily = "Arial")),  # Set font size and family for the heatmap title
              
              column_split = factor(rep(c("WT.F","WT.M","Kbk.F","Kbk.M"), each = 3), levels = c("WT.F", "WT.M", "Kbk.F", "Kbk.M")),  # Split columns into four groups
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              column_names_side = "top",
              column_title_gp = gpar(fontsize = 15, fontfamily = "Arial"),  # Font settings for column titles
              column_gap = unit(1, "mm"),
              cluster_rows = TRUE,
              show_row_dend = TRUE,
              row_title_gp = gpar(fontsize = 12),
              column_title_rot = 0,  # To make the column title horizontal
              show_heatmap_legend = TRUE,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
              use_raster = TRUE,
              raster_quality = 10,
              row_labels = row_labels,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 10),
              column_labels = rep(1:3, times = 4),  # Custom column names
              column_names_gp = gpar(fontsize = 15, fontfamily = "Arial"),  # Set font size and family for column labels
              column_names_rot = 360)  # Rotation for column labels






# Draw the heatmap
draw(ht, newpage = TRUE)

dev.off()

# end ---------------------------------------------------------------------


