library(Seurat)
#packageVersion("Seurat")
library(tidyverse)

################################################################################
# Quality control of kbk_combined object -----------------------------


# load the unfiltered seurat object(with joinedlayers) ---------------------------------------

kbk_combined<- readRDS("data/kbk_combined.rds")




#select only the singlets from combined data ----------------------------


kbk_combined<-subset(kbk_combined,HTO_classification.global=="Singlet")
dim(kbk_combined)


#change the working directory to save the figures
#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/figures_seurat5")

pdf("figures/ViolinPlot_feature_healthy_and_unhealthy.pdf",height=6,width=32)
VlnPlot(kbk_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(kbk_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,split.by="sample",pt.size=0)
VlnPlot(kbk_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by="sample",pt.size=0)

VlnPlot(kbk_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,split.by="Genotype",pt.size=0)
VlnPlot(kbk_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by="Genotype",pt.size=0)
dev.off()





 #change the working directory to save the data -----------------------



#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/data_seurat5")

# Filter and keep only the healthy cells ----------------------------------

kbk_combined.filtered<-subset(kbk_combined,nFeature_RNA>500&nFeature_RNA<8000&nCount_RNA<30000&percent.mt<10)





################################################################################


# Normalize, Findvaraiblefeatures, scaeledata and RunPCA------------------------



kbk_combined.filtered<-NormalizeData(kbk_combined.filtered,assay="RNA")
kbk_combined.filtered<-FindVariableFeatures(kbk_combined.filtered,nfeatures=3000)
kbk_combined.filtered<-ScaleData(kbk_combined.filtered)
kbk_combined.filtered<-RunPCA(kbk_combined.filtered,npcs=50)




# save the kbk_combined.filtered as rds file ------------------------------


saveRDS(kbk_combined.filtered,file="data/kbk_combined.filtered.rds")


################################################################################


# Batch correction by harmony ---------------------------------------------



library(harmony)
set.seed(123)
kbk_combined.filtered_harmony <- RunHarmony(kbk_combined.filtered, group.by.vars = "orig.ident", dims.use = 1:20, max.iter.harmony = 50)
kbk_combined.filtered_harmony <- RunUMAP(kbk_combined.filtered_harmony, reduction = "harmony", dims = 1:20)
kbk_combined.filtered_harmony <- RunTSNE(kbk_combined.filtered_harmony, reduction = "harmony", dims = 1:20)
kbk_combined.filtered_harmony <- FindNeighbors(kbk_combined.filtered_harmony, reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.6)





# save the figures in figures_seurat5 folder ------------------------------


#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/figures_seurat5")
pdf("figures/UMAP and TSNE clustering of harmonized sample.pdf")
DimPlot(kbk_combined.filtered_harmony,label=T,reduction="umap",size=10)
DimPlot(kbk_combined.filtered_harmony,label=T,reduction="tsne",size=10)

dev.off()



# Cluster annotation ------------------------------------------------------

#MG<-subset(kbk_combined.filtered_harmony,seurat_clusters %in% c(0,2,7,19))



################################################################################


#Idents(kbk_combined.filtered_harmony)<-kbk_combined.filtered_harmony@meta.data$seurat_clusters

#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/figures_seurat5")

################################################################################
# identify the clusters and assign the names ------------------------------



pdf("figures/Annotated_clusters.pdf",height=10,width=15)
#Identify microglial clusters
Microglia_markers <- c("Tmem119","Cx3cr1","Ptprc","Fcrls","P2ry12","Aif1","Hexb","Slc2a5","Siglech", "Trem2","Itgam")
DotPlot(kbk_combined.filtered_harmony,features = Microglia_markers)+
  ggtitle("Microglia_markers") +
  theme(plot.title = element_text(hjust = 0.5))




#Identify astrocyte clusters
Astrocyte_markers<-c("Gfap","Aqp4","Gldc","Slc1a3","Slc1a2","Aldh1l1","Dab1","Fam107a","Agt")




DotPlot(kbk_combined.filtered_harmony,features = Astrocyte_markers)+
  ggtitle("Astrocyte_markers") +
  theme(plot.title = element_text(hjust = 0.5))




#oligodendrocytes

Oligodendrocyte_markers<-c("Mbp","Sox10","Mog","Car2","Cnp","Plp1","Opalin","Mobp","Olig1","Olig2","Cldn11","Mag","Galc")
DotPlot(kbk_combined.filtered_harmony,features = Oligodendrocyte_markers)+
  ggtitle("Oligodendrocyte_markers") +
  theme(plot.title = element_text(hjust = 0.5))






#dendritic cells 
DCs.sig <- c("Flt3", "Zbtb46")
DotPlot(kbk_combined.filtered_harmony,features=DCs.sig)+
  ggtitle("DCs.sig") +
  theme(plot.title = element_text(hjust = 0.5))





#Oligodendrocyte precursor cells
OPC_markers<-c("Pdgfra","Cspg4","St8sia3","Olig2","Sox10","St8sia1","Cd9")
DotPlot(kbk_combined.filtered_harmony,features = OPC_markers)+
  ggtitle("OPC_markers") +
  theme(plot.title = element_text(hjust = 0.5))






#Neurons
Neuron_markers<- c("Neurod1","Bsn","Dcx","Gap43","Ncam1","Syp","Rbfox3","Nes","Tubb3",
                   "Tbr1","Stmn1",'Slc5a7','Chat','Slc17a8','Slc17a6','Gad1','Gad2',"Slc6a5","Syn1","Syn2","Syn3")

DotPlot(kbk_combined.filtered_harmony,features =Neuron_markers)+
  ggtitle("Neuron_markers") +
  theme(plot.title = element_text(hjust = 0.5))





#pericyte markers

# Pericytes ---------------------------------------------------------------
Pericyte_markers<-c("Pdgfrb","Rgs5","Cspg4","Vtn", "Higd1b", "S1pr3", "Mcam", "Ifitm1",  "Ehd3","Atp13a5")

DotPlot(kbk_combined.filtered_harmony,features =Pericyte_markers)+
  ggtitle("Pericyte_markers") +
  theme(plot.title = element_text(hjust = 0.5))




#Endothelial cell specific genes

endothelial_cell_markers<-c("Cldn5", "Cdh5", "Slc2a1", "Abcb1a", "Vwf", "Mfsd2a","Lef1","Sox17")
DotPlot(kbk_combined.filtered_harmony,features =endothelial_cell_markers)+
  ggtitle("endothelial_cell_markers") +
  theme(plot.title = element_text(hjust = 0.5))





immune_markers<-c("Cd19","Cd79a","Actb","Cd3d","Cd3e","Cd3g","Ptprc","Ms4a1","Ccr7","Ms4a7")
DotPlot(kbk_combined.filtered_harmony,features =immune_markers)+
  ggtitle("immune_markers") +
  theme(plot.title = element_text(hjust = 0.5))





# vascular leptomeningeal cells (VLMC) --------------------------------------------------------------------
VLMC_markers<-c("Col1a2","Col1a1","Mgp","Lum","Dcn","Pdgfra","Serpinf1")
DotPlot(kbk_combined.filtered_harmony,features = VLMC_markers)+
  ggtitle("VLMC_markers") +
  theme(plot.title = element_text(hjust = 0.5))




#choroid plexus markers
choroidPlexus_marker_genes<-c("Aqp1","Folr1","Klk1b5","Lama1","S100b","Otx2","Slc4a5","Igf2")
DotPlot(kbk_combined.filtered_harmony,features = choroidPlexus_marker_genes)+
  ggtitle("choroidPlexus_markers") +
  theme(plot.title = element_text(hjust = 0.5))




dev.off()





# Generate a table with canonical markers used in our study ---------------


library(tidyverse)
library(readxl)



markers_df<-data.frame(marker_genes=Microglia_markers)
markers_df$cell_type<- "MG"
markers_df

astro_markers_df<-data.frame(marker_genes=Astrocyte_markers)
astro_markers_df$cell_type<-"Astro"
markers_df<-rbind(markers_df,astro_markers_df)

OD_markers_df<-data.frame(marker_genes=Oligodendrocyte_markers)
OD_markers_df$cell_type<-"Oligodendrocytes"
markers_df<-rbind(markers_df,OD_markers_df)

DC_df<-data.frame(marker_genes=DCs.sig)
DC_df$cell_type<-"Dendritic_cells"
markers_df<-rbind(markers_df,DC_df)


OPC_df<-data.frame(marker_genes=OPC_markers)
OPC_df$cell_type<-"OPC"  
markers_df<-rbind(markers_df,OPC_df)

neuron_df<-data.frame(marker_genes=Neuron_markers)
neuron_df$cell_type<-"Neurons"
markers_df<-rbind(markers_df,neuron_df)


pericyte_df<-data.frame(marker_genes=Pericyte_markers)
pericyte_df$cell_type<-"Pericytes"
markers_df<-rbind(markers_df,pericyte_df)


endo_df<-data.frame(marker_genes=endothelial_cell_markers)
endo_df$cell_type<-"Endo"
markers_df<-rbind(markers_df,endo_df)


immune_df<-data.frame(marker_genes=immune_markers)
immune_df$cell_type<-"Immune"
markers_df<-rbind(markers_df,immune_df)


VLMC_df<-data.frame(marker_genes=VLMC_markers)
VLMC_df$cell_type<-"VLMC"
markers_df<-rbind(markers_df,VLMC_df)


choroidPlexus_df<-data.frame(marker_genes=choroidPlexus_marker_genes)
choroidPlexus_df$cell_type<-"Choroid_plexus"
markers_df<-rbind(markers_df,choroidPlexus_df)



write_excel_csv(markers_df,"results/tables/canonical_markers_table_used_for_cluster_identification.xls")

################################################################################

#Cluster annotation
################################################################################

Idents(kbk_combined.filtered_harmony)<-kbk_combined.filtered_harmony@meta.data$seurat_clusters

DimPlot(kbk_combined.filtered_harmony,label=T)



#kbk_combined.filtered_harmony<-readRDS("C:\\Users\\kppra\\OneDrive\\Desktop\\Prasad\\DE_analysis_Vivek's code\\kbk_combined.filtered_harmony.rds")


# Change the name of clusters ---------------------------------------------

pdf("figures/Change the cluster names.pdf", width=10, height=10)


DimPlot(kbk_combined.filtered_harmony,reduction="umap",label=T,label.size = 5)+
  ggtitle("Umap_before cluster assignemene")+
  theme(plot.title = element_text(hjust = 0.5))

DimPlot(kbk_combined.filtered_harmony,reduction="tsne",label=T,label.size = 5)+
  ggtitle("tsne_before cluster assignemene")+
  theme(plot.title = element_text(hjust = 0.5))


dev.off()
################################################################################

#change the idents to seurat_clusters
Idents(kbk_combined.filtered_harmony)<-kbk_combined.filtered_harmony@meta.data$seurat_clusters
kbk_combined.filtered_harmony <- RenameIdents(kbk_combined.filtered_harmony, 
                                              "0" = "MG1",
                                              "1" = "Endo1",
                                              "2" = "MG2",
                                              "3" = "Pericytes1",
                                              "4" = "Endo2",
                                              "5" = "Astro1",
                                              "6" = "CP1",
                                              "7" = "MG3",
                                              "8" = "OPC",
                                              "9" = "OD1",
                                              "10" = "Immune_cells1",  # Ms4a7 expressed
                                              "11" = "Pericytes2",
                                              "12" = "Pericytes3",
                                              "13" = "Astro2",
                                              "14" = "Endo3",
                                              "15" = "Astro3",
                                              "16" = "Endo4",
                                              "17" = "CP2",
                                              "18" = "DC",
                                              "19" = "MG4",
                                              "20" = "OD2",
                                              "21" = "Endo5",
                                              "22" = "VLMC",
                                              "23" = "Immune_cells2",
                                              "24" = "Immune_cells3",
                                              "25" = "NPC/Neurons",
                                              "26" = "Astro4",
                                              "27" = "Pericytes4")

#add a new variable cluster_annotation to the metadata
kbk_combined.filtered_harmony@meta.data$cluster_annotation<-Idents(kbk_combined.filtered_harmony)

pdf("figures/Cluster_annotation_varible_named.pdf",height = 10,width=10)

DimPlot(kbk_combined.filtered_harmony,label=T)+
  ggtitle("cluster_annotation")+
  theme(plot.title = element_text(hjust = 0.5))

DimPlot(kbk_combined.filtered_harmony,label=T,reduction="tsne")+
  ggtitle("cluster_annotation")+
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

#save as an rds object
#saveRDS(kbk_combined.filtered_harmony,file="kbk_combined.filtered_harmony.rds")

#add a cell_type in kbk_combined.filtered_harmony

#step 1
#change the idents to seurat_clusters
Idents(kbk_combined.filtered_harmony)<-kbk_combined.filtered_harmony@meta.data$seurat_clusters

#step2 combine similar cell types together
kbk_combined.filtered_harmony <- RenameIdents(kbk_combined.filtered_harmony, 
                                              "0" = "MG",
                                              "1" = "Endo",
                                              "2" = "MG",
                                              "3" = "Pericytes",
                                              "4" = "Endo",
                                              "5" = "Astro",
                                              "6" = "CP",
                                              "7" = "MG",
                                              "8" = "OPC",
                                              "9" = "OD",
                                              "10" = "Immune_cells",  # Ms4a7 expressed
                                              "11" = "Pericytes",
                                              "12" = "Pericytes",
                                              "13" = "Astro",
                                              "14" = "Endo",
                                              "15" = "Astro",
                                              "16" = "Endo",
                                              "17" = "CP",
                                              "18" = "DC",
                                              "19" = "MG",
                                              "20" = "OD",
                                              "21" = "Endo",
                                              "22" = "VLMC",
                                              "23" = "Immune_cells",
                                              "24" = "Immune_cells",
                                              "25" = "NPC/Neuron",
                                              "26" = "Astro",
                                              "27" = "Pericytes")

#add a new variable cell_types to the metadata
kbk_combined.filtered_harmony@meta.data$cell_type<-Idents(kbk_combined.filtered_harmony)

png("figures/clusters_with_cell_type_annotation.png",
    width=2500,height=1250,res = 300)
DimPlot(kbk_combined.filtered_harmony,label=T)+
  ggtitle("cell_types")+
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

pdf("figures/clusters_with_cell_type_annotation.pdf", width=10,height=10)

DimPlot(kbk_combined.filtered_harmony,label=T)+
  ggtitle("cell_types")+
  theme(plot.title = element_text(hjust = 0.5))

DimPlot(kbk_combined.filtered_harmony,label=T,reduction="tsne")+
  ggtitle("cell_types")+
  theme(plot.title = element_text(hjust = 0.5))

dev.off()


# #png file format --------------------------------------------------------



# png("Clustering.png", width=2500,height=1250,res = 300)
# 
# DimPlot(kbk_combined.filtered_harmony,label=T)+
#   ggtitle("scRNAseq")+
#   theme(plot.title = element_text(hjust = 0.5))
# dev.off()





# #batch corrected, cluster annotated #save as rds object -----------------


#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/data_seurat5")

saveRDS(kbk_combined.filtered_harmony,file="data/kbk_combined.filtered_harmony.rds")

################################################################################
# The cell counts of different brain cells  -------------------------------



cell_num_df<-data.frame(table(kbk_combined.filtered_harmony@meta.data$cell_type,kbk_combined.filtered_harmony@meta.data$group))

names(cell_num_df)[2:3]<-c("group","cell_num")

cell_num_df<-pivot_wider(cell_num_df,names_from = "group",values_from = "cell_num")

# convert them to a dataframe ---------------------------------------------

cell_num_df<-as.data.frame(cell_num_df)

#add rownames

rownames(cell_num_df)<-cell_num_df$Var1

#Optionally remove the Var1 column
cell_num_df$Var1<-NULL



# new order for the cell_num_df -------------------------------------------

new_order<-c(6,10,8,7,9,1,5,3,2,4)

cell_type_num_df<-cell_num_df[,new_order]








#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/results_seurat5")

write.csv(cell_type_num_df, file="results/cell_type numbers per group.csv")


# to save the file in excel format ----------------------------------------


library(openxlsx)

write.xlsx(cell_type_num_df, file="results/cell_type numbers per group.xlsx",
           rowNames = TRUE)



######################################################################


# Total cell number per genotype and gender -------------------------------

Gender_df<-as.data.frame(table(kbk_combined.filtered_harmony@meta.data$Genotype,
                               kbk_combined.filtered_harmony@meta.data$SEX))




#rename the columns and add a new variable percent to Gender_df
Gender_df<-Gender_df %>%
  rename(Genotype=Var1,Sex=Var2,cell_num=Freq)%>%
  mutate(percent=(cell_num/sum(cell_num))*100)

write.csv(Gender_df, "results/Total cell number per genotype and gender.csv")
write.xlsx(Gender_df, file="results/Total cell number per genotype and gender.xlsx")



# Microglia cell number per group and gender ------------------------------

df<-as.data.frame(table(kbk_combined.filtered_harmony@meta.data$Genotype,
      kbk_combined.filtered_harmony@meta.data$SEX,
      kbk_combined.filtered_harmony@meta.data$group))

df<-df %>% filter(Freq>0)%>% 
  arrange(Var2)%>%
  rename(Genotype=Var1,Sex=Var2,group=Var3,Cell_num=Freq)


write.csv(df,file="results/Hippocampus cell number per group and gender.csv")
write.xlsx(df,file="results/Hippocampus cell number per group and gender.xlsx")


# End ---------------------------------------------------------------------

################################################################################
# make a table with total cell type numbers
################################################################################

data<-data.frame(table(kbk_combined.filtered_harmony$cell_type)) %>% rename(cell_type=Var1,cell_num=Freq) %>%
  mutate(percent=cell_num/sum(cell_num)*100) 
data

rownames(data)<-data$cell_type
data$cell_type<-NULL

# Calculate totals
cell_num <- sum(data$cell_num)
percent <- sum(data$percent)

# Add totals above the column headings
totals_row <- data.frame(cell_num = cell_num, percent = percent)
totals_row

data_with_totals <- rbind(totals_row, data)
data_with_totals

data_with_totals$percent<-round(data_with_totals$percent,1)

# Update the row names to make it clear
#row.names(data_with_totals) <- c("Total","", rownames(data))



# Insert an empty row after the first row
empty_row <- data.frame(cell_num = "", percent = "")  # Create an empty row
data_with_totals <- rbind(data_with_totals[1, ], empty_row, data_with_totals[-1, ])
#rownames(data_with_totals)[2]<- ""


# Update the row names to make it clear
#row.names(data_with_totals) <- c("Total", as.character(seq_len(nrow(data))))


# Update the row names to make it clear
row.names(data_with_totals) <- c("Total_cell_num","", rownames(data))



# Print the result
print(data_with_totals)

#setwd("C:/Users/kppra/OneDrive/Desktop/scRNAseq_Seurat5/results_seurat5/tables")

write.csv(data_with_totals, "results/tables/cell number per cell types across the dataset.csv")
################################################################################





# end ---------------------------------------------------------------------


