
# Microglia_subset_and_analysis -------------------------------------------



# # Load libraries --------------------------------------------------------
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(magrittr)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(harmony)
library(monocle3)
library(plotly)#for interactive plotting
library(DT)
library(ggpubr)
library(ComplexHeatmap)


# load harmonized data ----------------------------------------------------



kbk_combined.filtered_harmony<-readRDS("data/kbk_combined.filtered_harmony.rds")



# save cell_type in percentage as a csv file ------------------------------

################################################################################

cell_types<-data.frame(table(kbk_combined.filtered_harmony@meta.data$cell_type))
cell_types<-cell_types%>% rename(cell_type=Var1, cell_num=Freq) %>%
  mutate(percent=round((cell_num/sum(cell_num))*100,2))



write.csv(cell_types, "results/cell_types_percentage.csv")


################################################################################



# Subset microglia --------------------------------------------------------

MG<-subset(kbk_combined.filtered_harmony,cell_type=="MG")


# Repeat the preprocessing steps ------------------------------------------


set.seed(123)
MG<-NormalizeData(MG)
MG<-FindVariableFeatures(MG,nfeatures=3000)
MG<-ScaleData(MG)
ElbowPlot(MG)
MG<-RunPCA(MG,dims=1:20)



MG <- RunHarmony(MG, group.by.vars = "orig.ident", dims.use = 1:20, max.iter.harmony = 50)
MG <- RunUMAP(MG, reduction = "harmony", dims = 1:20)
#MG <- RunTSNE(MG, reduction = "harmony", dims = 1:20)
MG <- FindNeighbors(MG, reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.15)


#pdf("Microglia_allSamples.pdf",width=5,height=5)
#cluster the microglia dataset

png("figures/Microglia_Allsamples.png")
DimPlot(MG,reduction="umap",label=TRUE)+
  ggtitle("Microglia")+
  theme(plot.title = element_text(hjust = 0.5))
  
dev.off()






# Add a new column MG_subsets to the metadata -----------------------------

MG@meta.data$MG_subsets<-paste(Idents(MG),"MG",sep="-")#This code join the cluster names(0 to 6) to a string(character) of MG
MG@meta.data$MG_subsets<-as.factor(MG@meta.data$MG_subsets)# this code change a character vector to factor vector
Idents(MG)<-MG@meta.data$MG_subsets
#Idents(MG)
levels(MG)

#visualize the MG cluster
DimPlot(MG,label=T)+NoLegend()



# MG cells per group ------------------------------------------------------



MG_num_df<-data.frame(table(MG@meta.data$group))
names(MG_num_df)<-c("group","number")


MG_num_df$group<-factor(MG_num_df$group,levels=c("WT.1m.Un","WT.3m.Un","WT.3m.DMSO_PEG","WT.3m.DCF","WT.3m.TCF",
                                                 "Kbk.1m.Un","Kbk.3m.Un","Kbk.3m.DMSO_PEG","Kbk.3m.DCF","Kbk.3m.TCF"))
#save as csv file
write.csv(MG_num_df,file="results/Microglia num per group.csv")

write.xlsx(MG_num_df, file="results/Microglia_number_per_group.xlsx", rowNames = TRUE)

#########################################################################


#genotype and SEX
MG_num_per_G_S<-data.frame(table(MG@meta.data$Genotype,MG@meta.data$SEX))


# rename the column names -------------------------------------------------




MG_num_per_G_S<-MG_num_per_G_S %>%
  rename(Genotype=Var1,
         SEX=Var2, 
         cell_num=Freq)

write.xlsx(MG_num_per_G_S,
           file="results/MG_cellnum_split_Gender_and_Sex.xlsx", rowNames=TRUE)


# cluster proportoin of individual subsets --------------------------------




######################################################
#create a function to calculate the cluster proportoin of individual subsets
#######################################################
cluster_proportion<-function(cluster){
  cluster<- subset(MG@meta.data,MG_subsets == cluster)# subset a data frame ("cluster")  from "MG" seuratobj metadata
  cells_percent <-table(cluster$group)/table(MG@meta.data$group)*100#calculate the percentage of each cluster compared to total cells in the group
  percent<-cells_percent/sum(cells_percent)*100# the cell percent devided by sum of cell_percent *100
  percent<-as.data.frame(percent)
  return(percent)
}


clusters <- unique(MG$MG_subsets)
all_proportions<-data.frame()
for (i in 1:length(clusters)){
  cls<-cluster_proportion(clusters[i])
  cls$cluster<-clusters[i]
  cls$group<-as.factor(names(table(MG@meta.data$group)))
  all_proportions<-rbind(all_proportions,cls)
  #print(all_proportions)
}



# Save the cluster proportion to result folder as csv format --------------



write.csv(all_proportions,file="results/microglia_subset_clusterProportion.csv")

################################################################################



#######################################################################
#Changing the levels of group variable
#######################################################################

levels(all_proportions$group)
#Initially the level was as follows

# [1] "Kbk.1m.Un"       "Kbk.3m.DCF"      "Kbk.3m.DMSO_PEG"
# [4] "Kbk.3m.TCF"      "Kbk.3m.Un"       "WT.1m.Un"       
# [7] "WT.3m.DCF"       "WT.3m.DMSO_PEG"  "WT.3m.TCF"      
# [10] "WT.3m.Un" 

#Changing the levels of group variable
vector=c("WT.1m.Un","Kbk.1m.Un","WT.3m.Un","Kbk.3m.Un","WT.3m.DMSO_PEG",
         "Kbk.3m.DMSO_PEG","WT.3m.DCF","Kbk.3m.DCF","WT.3m.TCF","Kbk.3m.TCF")

#add a group variable to all_proportions data frame
all_proportions$group<-factor(all_proportions$group,levels=vector)

#To check whether the level change has been successful
levels(all_proportions$group)
###############################################################################






#visualization of microglia subset cluster proportion
png("figures/MG_subsets_Group_proportions.png",width=639,height=403)
ggplot(all_proportions, aes(y=Freq, x=cluster, fill=group)) +
  geom_bar(stat='identity') +
  scale_fill_brewer(palette = "Paired")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank())+
  RotatedAxis()+
  ggtitle('MG_subsets_Group_proportions_by_cluster')+
  theme(plot.title = element_text(hjust = 0.5))
  
dev.off()


################################################################################


# save harmonized MG seurat data in data folder ---------------------------

#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/data_seurat5")

saveRDS(MG,file="data/MG.rds")


# : Microglia cell number per group and gender ----------------------------






df<-as.data.frame(table(MG@meta.data$Genotype,
                        MG@meta.data$SEX,
                        MG@meta.data$group))

df<-df %>% filter(Freq>0)%>% 
  arrange(Var2)%>%
  rename(Genotype=Var1,Sex=Var2,group=Var3,Cell_num=Freq)


write.csv(df,file="results/Microglia cell number per group and gender.csv")




MG_df<-as.data.frame(table(MG@meta.data$Genotype,
                        MG@meta.data$SEX))

MG_df<-MG_df %>% rename(Genotype=Var1,Sex=Var2,cell_num=Freq)

write.csv(MG_df,file="results/Microglia cell number per Genotype and gender.csv")

# End ---------------------------------------------------------------------





