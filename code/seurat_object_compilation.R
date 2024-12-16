
#.libPaths()

# load libraries ----------------------------------------------------------
library(Seurat)
library(tidyverse)


# load data ---------------------------------------------------------------

#Generate seurat objects of all batches(batch1,2,3,4,5,6,7,8)
batch1.data<-Read10X("data/batch1_filtered_feature_bc_matrix")
batch1<-CreateSeuratObject(batch1.data$`Gene Expression`, project="kbk_batch1")
batch1[["HTO"]]<-CreateAssayObject(batch1.data$`Antibody Capture`)
batch1

batch1_HTOs<-as.data.frame(batch1$HTO@counts[1:8,])
batch1[["HTO"]]<-CreateAssayObject(batch1_HTOs)
batch1<-NormalizeData(batch1,assay="HTO",normalization.method="CLR")
batch1<-HTODemux(batch1,positive.quantile=0.99)



#HTO_maxID can be individual sample ID
#Add Genotype column to batch1 metadata and convert Genotype variable to factor
batch1@meta.data$Genotype<-batch1@meta.data$HTO_maxID
batch1@meta.data$Genotype<-recode(batch1@meta.data$Genotype,"M-HTO-1"="Kbk",
                                  "M-HTO-2"="Kbk","M-HTO-3"="Kbk","M-HTO-4"="WT","M-HTO-5"="WT","M-HTO-6"="Kbk",
                                  "M-HTO-7"="Kbk","M-HTO-8"="WT")

batch1@meta.data$Genotype<-factor(batch1@meta.data$Genotype,levels=c("WT","Kbk"))



#Add Age column to the batch1 metadata
batch1@meta.data$Age<-batch1@meta.data$HTO_maxID
batch1@meta.data$Age<-recode(batch1@meta.data$Age,"M-HTO-1"="3m",
                             "M-HTO-2"="3m","M-HTO-3"="3m","M-HTO-4"="3m","M-HTO-5"="3m","M-HTO-6"="3m",
                             "M-HTO-7"="1m","M-HTO-8"="1m")
batch1@meta.data$Age<-factor(batch1@meta.data$Age,levels=c("1m","3m"))

#Add treatment group to the batch1 metadata

batch1@meta.data$Treatment<-batch1@meta.data$HTO_maxID
batch1@meta.data$Treatment<-recode(batch1@meta.data$Treatment,"M-HTO-1"="TCF",
                                   "M-HTO-2"="Un","M-HTO-3"="Un","M-HTO-4"="Un","M-HTO-5"="Un","M-HTO-6"="TCF",
                                   "M-HTO-7"="Un","M-HTO-8"="Un")
batch1@meta.data$Treatment<-factor(batch1@meta.data$Treatment,levels=c("Un","TCF"))



#Add SEX variable to the batch1 metadata
batch1@meta.data$SEX<-batch1@meta.data$HTO_maxID
batch1@meta.data$SEX<-recode(batch1@meta.data$SEX,"M-HTO-1"="M",
                             "M-HTO-2"="F","M-HTO-3"="M","M-HTO-4"="F","M-HTO-5"="M","M-HTO-6"="F",
                             "M-HTO-7"="F","M-HTO-8"="M")
batch1@meta.data$SEX<-factor(batch1@meta.data$SEX,levels=c("M","F"))


#Add sample name(mouse name) to the batch1 metadata
batch1@meta.data$sample_name<-batch1@meta.data$HTO_maxID
batch1@meta.data$sample_name<-recode(batch1@meta.data$sample_name,"M-HTO-1"="Kbk_3m_TCF_M",
                                     "M-HTO-2"="Kbk_3m_Un_F","M-HTO-3"="Kbk_3m_Un_M","M-HTO-4"="WT_3m_Un_F","M-HTO-5"="WT_3m_Un_M","M-HTO-6"="Kbk_3m_TCF_F",
                                     "M-HTO-7"="Kbk_1m_Un_F","M-HTO-8"="WT_1m_Un_M")
batch1@meta.data$unique_name<-batch1@meta.data$HTO_maxID
batch1@meta.data$unique_name<-recode(batch1@meta.data$unique_name,"M-HTO-1"="Kbk_3m_TCF_M_1",
                                     "M-HTO-2"="Kbk_3m_Un_F_1","M-HTO-3"="Kbk_3m_Un_M_1","M-HTO-4"="WT_3m_Un_F_1","M-HTO-5"="WT_3m_Un_M_1","M-HTO-6"="Kbk_3m_TCF_F_1",
                                     "M-HTO-7"="Kbk_1m_Un_F_1","M-HTO-8"="WT_1m_Un_M_1")

#Save the batch1 file with all the new variables added
saveRDS(batch1, file="data/batch1.rds")



# Batch2 to seuratobject --------------------------------------------------

batch2.data<-Read10X("data/batch2_filtered_feature_bc_matrix")
batch2<-CreateSeuratObject(batch2.data$"Gene Expression",project="kbk_batch2")
batch2[["HTO"]]<-CreateAssayObject(batch2.data$`Antibody Capture`)
batch2_HTOs<-as.data.frame(batch2$HTO@counts[1:8,])
batch2[["HTO"]]<-CreateAssayObject(batch2_HTOs)
batch2<-NormalizeData(batch2,assay="HTO",normalization.method="CLR")
batch2<-HTODemux(batch2,positive.quantile=0.99)


batch2@meta.data$Genotype<-batch2@meta.data$HTO_maxID
batch2@meta.data$Genotype<-recode(batch2@meta.data$Genotype,"M-HTO-1"="Kbk",
                                  "M-HTO-2"="Kbk","M-HTO-3"="Kbk","M-HTO-4"="Kbk","M-HTO-5"="Kbk","M-HTO-6"="Kbk",
                                  "M-HTO-7"="Kbk","M-HTO-8"="WT")
batch2@meta.data$Genotype<-factor(batch2@meta.data$Genotype,levels=c("WT","Kbk"))

#Add Age column to the batch2 metadata
batch2@meta.data$Age<-batch2@meta.data$HTO_maxID
batch2@meta.data$Age<-recode(batch2@meta.data$Age,"M-HTO-1"="3m",
                             "M-HTO-2"="3m","M-HTO-3"="3m","M-HTO-4"="3m","M-HTO-5"="3m","M-HTO-6"="3m",
                             "M-HTO-7"="1m","M-HTO-8"="1m")
batch2@meta.data$Age<-factor(batch2@meta.data$Age,levels=c("1m","3m"))

#Add treatment group to the batch2 metadata

batch2@meta.data$Treatment<-batch2@meta.data$HTO_maxID
batch2@meta.data$Treatment<-recode(batch2@meta.data$Treatment,"M-HTO-1"="TCF",
                                   "M-HTO-2"="TCF","M-HTO-3"="TCF","M-HTO-4"="DCF","M-HTO-5"="DCF","M-HTO-6"="DMSO_PEG",
                                   "M-HTO-7"="Un","M-HTO-8"="Un")
batch2@meta.data$Treatment<-factor(batch2@meta.data$Treatment,levels=c("Un","DMSO_PEG" ,"DCF","TCF"))



#Add SEX variable to the batch2 metadata
batch2@meta.data$SEX<-batch2@meta.data$HTO_maxID
batch2@meta.data$SEX<-recode(batch2@meta.data$SEX,"M-HTO-1"="M",
                             "M-HTO-2"="F","M-HTO-3"="F","M-HTO-4"="F","M-HTO-5"="M","M-HTO-6"="M",
                             "M-HTO-7"="M","M-HTO-8"="M")
batch2@meta.data$SEX<-factor(batch2@meta.data$SEX,levels=c("M","F"))


#Add sample name(mouse name) to the batch2 metadata
batch2@meta.data$sample_name<-batch2@meta.data$HTO_maxID
batch2@meta.data$sample_name<-recode(batch2@meta.data$sample_name,"M-HTO-1"="Kbk_3m_TCF_M",
                                     "M-HTO-2"="Kbk_3m_TCF_F","M-HTO-3"="Kbk_3m_TCF_F","M-HTO-4"="Kbk_3m_DCF_F",
                                     "M-HTO-5"="Kbk_3m_DCF_M","M-HTO-6"="Kbk_3m_DMSO_PEG_M",
                                     "M-HTO-7"="Kbk_1m_Un_M","M-HTO-8"="WT_1m_Un_M")

batch2@meta.data$unique_name<-batch2@meta.data$HTO_maxID
batch2@meta.data$unique_name<-recode(batch2@meta.data$unique_name,"M-HTO-1"="Kbk_3m_TCF_M_2",
                                     "M-HTO-2"="Kbk_3m_TCF_F_2","M-HTO-3"="Kbk_3m_TCF_F_3","M-HTO-4"="Kbk_3m_DCF_F_1",
                                     "M-HTO-5"="Kbk_3m_DCF_M_1","M-HTO-6"="Kbk_3m_DMSO_PEG_M_1",
                                     "M-HTO-7"="Kbk_1m_Un_M_1","M-HTO-8"="WT_1m_Un_M_2")


#Save the batch2 file with all the new variables added
saveRDS(batch2, file="data/batch2.rds")



# Batch3 to Seuratobject --------------------------------------------------

#reading batch 3 sequence

batch3.data<-Read10X("data/batch3_filtered_feature_bc_matrix")
batch3<-CreateSeuratObject(batch3.data$"Gene Expression",project="kbk_batch3")
batch3[["HTO"]]<-CreateAssayObject(batch3.data$`Antibody Capture`)
batch3_HTOs<-as.data.frame(batch3$HTO@counts[39:46,])
batch3[["HTO"]]<-CreateAssayObject(batch3_HTOs)
batch3<-NormalizeData(batch3,assay="HTO",normalization.method="CLR")
batch3<-HTODemux(batch3,positive.quantile=0.99)


batch3@meta.data$Genotype<-batch3@meta.data$HTO_maxID
batch3@meta.data$Genotype<-recode(batch3@meta.data$Genotype,"M-HTO-1"="Kbk",
                                  "M-HTO-2"="Kbk","M-HTO-3"="Kbk","M-HTO-4"="Kbk","M-HTO-5"="Kbk","M-HTO-6"="Kbk",
                                  "M-HTO-7"="WT","M-HTO-8"="WT")
batch3@meta.data$Genotype<-factor(batch3@meta.data$Genotype,levels=c("WT","Kbk"))



#Add Age column to the batch3 metadata
#All the mice are of 3 month age
batch3@meta.data$Age<-batch3@meta.data$HTO_maxID
batch3@meta.data$Age<-recode(batch3@meta.data$Age,"M-HTO-1"="3m",
                             "M-HTO-2"="3m","M-HTO-3"="3m","M-HTO-4"="3m","M-HTO-5"="3m","M-HTO-6"="3m",
                             "M-HTO-7"="3m","M-HTO-8"="3m")
batch3@meta.data$Age<-factor(batch3@meta.data$Age,levels="3m")

#Add treatment group to the batch3 metadata

batch3@meta.data$Treatment<-batch3@meta.data$HTO_maxID
batch3@meta.data$Treatment<-recode(batch3@meta.data$Treatment,"M-HTO-1"="TCF",
                                   "M-HTO-2"="Un","M-HTO-3"="Un","M-HTO-4"="DCF","M-HTO-5"="DMSO_PEG","M-HTO-6"="DCF",
                                   "M-HTO-7"="Un","M-HTO-8"="Un")
batch3@meta.data$Treatment<-factor(batch3@meta.data$Treatment,levels=c("Un","DMSO_PEG" ,"DCF","TCF"))



#Add SEX variable to the batch3 metadata
batch3@meta.data$SEX<-batch3@meta.data$HTO_maxID
batch3@meta.data$SEX<-recode(batch3@meta.data$SEX,"M-HTO-1"="M",
                             "M-HTO-2"="F","M-HTO-3"="M","M-HTO-4"="F","M-HTO-5"="F","M-HTO-6"="M",
                             "M-HTO-7"="M","M-HTO-8"="F")
batch3@meta.data$SEX<-factor(batch3@meta.data$SEX,levels=c("M","F"))


#Add sample name(mouse name) to the batch3 metadata
batch3@meta.data$sample_name<-batch3@meta.data$HTO_maxID
batch3@meta.data$sample_name<-recode(batch3@meta.data$sample_name,"M-HTO-1"="Kbk_3m_TCF_M",
                                     "M-HTO-2"="Kbk_3m_Un_F","M-HTO-3"="Kbk_3m_Un_M","M-HTO-4"="Kbk_3m_DCF_F","M-HTO-5"="Kbk_3m_DMSO_PEG_F","M-HTO-6"="Kbk_3m_DCF_M",
                                     "M-HTO-7"="WT_3m_Un_M","M-HTO-8"="WT_3m_Un_F")



#Add unique name(mouse name) to the batch3 metadata
batch3@meta.data$unique_name<-batch3@meta.data$HTO_maxID
batch3@meta.data$unique_name<-recode(batch3@meta.data$unique_name,"M-HTO-1"="Kbk_3m_TCF_M_3",
                                     "M-HTO-2"="Kbk_3m_Un_F_2","M-HTO-3"="Kbk_3m_Un_M_2","M-HTO-4"="Kbk_3m_DCF_F_2","M-HTO-5"="Kbk_3m_DMSO_PEG_F_1","M-HTO-6"="Kbk_3m_DCF_M_2",
                                     "M-HTO-7"="WT_3m_Un_M_2","M-HTO-8"="WT_3m_Un_F_2")


#Save the batch3 file with all the new variables added
saveRDS(batch3, file="data/batch3.rds")





# Load batch4 filtered matrix  --------------------------------------------



batch4.data<-Read10X("data/batch4_filtered_feature_bc_matrix")
batch4<-CreateSeuratObject(batch4.data$`Gene Expression`,project="kbk_batch4")#batch9 is rerun of batch4
batch4[["HTO"]]<-CreateAssayObject(batch4.data$`Antibody Capture`)
batch4_HTOs<-as.data.frame(batch4$HTO@counts[1:8,])
batch4[["HTO"]]<-CreateAssayObject(batch4_HTOs)
batch4<-NormalizeData(batch4,assay="HTO",normalization.method="CLR")
batch4<-HTODemux(batch4,positive.quantile=0.99)

#Add genotype to the batch4
batch4@meta.data$Genotype<-batch4@meta.data$HTO_maxID
batch4@meta.data$Genotype<-recode(batch4@meta.data$Genotype,"M-HTO-1"="WT",
                                  "M-HTO-2"="WT","M-HTO-3"="WT","M-HTO-4"="WT","M-HTO-5"="WT","M-HTO-6"="WT",
                                  "M-HTO-7"="WT","M-HTO-8"="WT")
batch4@meta.data$Genotype<-factor(batch4@meta.data$Genotype,levels="WT")



#Add Age column to the batch4 metadata
#All the mice are of 3 month age
batch4@meta.data$Age<-batch4@meta.data$HTO_maxID
batch4@meta.data$Age<-recode(batch4@meta.data$Age,"M-HTO-1"="3m",
                             "M-HTO-2"="3m","M-HTO-3"="3m","M-HTO-4"="3m","M-HTO-5"="3m","M-HTO-6"="3m",
                             "M-HTO-7"="3m","M-HTO-8"="1m")
batch4@meta.data$Age<-factor(batch4@meta.data$Age,levels=c("1m","3m"))

#Add treatment group to the batch4 metadata

batch4@meta.data$Treatment<-batch4@meta.data$HTO_maxID
batch4@meta.data$Treatment<-recode(batch4@meta.data$Treatment,"M-HTO-1"="Un",
                                   "M-HTO-2"="TCF","M-HTO-3"="TCF","M-HTO-4"="TCF","M-HTO-5"="DCF","M-HTO-6"="DMSO_PEG",
                                   "M-HTO-7"="DMSO_PEG","M-HTO-8"="Un")
batch4@meta.data$Treatment<-factor(batch4@meta.data$Treatment,levels=c("Un","DMSO_PEG","DCF","TCF"))



#Add SEX variable to the batch4 metadata
batch4@meta.data$SEX<-batch4@meta.data$HTO_maxID
batch4@meta.data$SEX<-recode(batch4@meta.data$SEX,"M-HTO-1"="M",
                             "M-HTO-2"="F","M-HTO-3"="F","M-HTO-4"="M","M-HTO-5"="F","M-HTO-6"="F",
                             "M-HTO-7"="M","M-HTO-8"="F")
batch4@meta.data$SEX<-factor(batch4@meta.data$SEX,levels=c("M","F"))


#Add sample name(mouse name) to the batch4 metadata
batch4@meta.data$sample_name<-batch4@meta.data$HTO_maxID
batch4@meta.data$sample_name<-recode(batch4@meta.data$sample_name,"M-HTO-1"="WT_3m_Un_M",
                                     "M-HTO-2"="WT_3m_TCF_F","M-HTO-3"="WT_3m_TCF_F","M-HTO-4"="WT_3m_TCF_M","M-HTO-5"="WT_3m_DCF_F",
                                     "M-HTO-6"="WT_3m_DMSO_PEG_F",
                                     "M-HTO-7"="WT_3m_DMSO_PEG_M","M-HTO-8"="WT_1m_Un_F")



#Add unique name(mouse name) to the batch4 metadata
batch4@meta.data$unique_name<-batch4@meta.data$HTO_maxID
batch4@meta.data$unique_name<-recode(batch4@meta.data$unique_name,"M-HTO-1"="WT_3m_Un_M_3",
                                     "M-HTO-2"="WT_3m_TCF_F_1","M-HTO-3"="WT_3m_TCF_F_2","M-HTO-4"="WT_3m_TCF_M_1","M-HTO-5"="WT_3m_DCF_F_1",
                                     "M-HTO-6"="WT_3m_DMSO_PEG_F_1",
                                     "M-HTO-7"="WT_3m_DMSO_PEG_M_1","M-HTO-8"="WT_1m_Un_F_1")



#Save the batch4 file with all the new variables added
saveRDS(batch4, file="data/batch4.rds")



# Load batch5 filtered matrix  --------------------------------------------




batch5.data<-Read10X("data/batch5_filtered_feature_bc_matrix")
batch5<-CreateSeuratObject(batch5.data$`Gene Expression`,project="kbk_batch5")
batch5[["HTO"]]<-CreateAssayObject(batch5.data$`Antibody Capture`)
batch5_HTOs<-as.data.frame(batch5$HTO@counts[39:46,])
batch5[["HTO"]]<-CreateAssayObject(batch5_HTOs)
batch5<-NormalizeData(batch5,assay="HTO",normalization.method="CLR")
batch5<-HTODemux(batch5,positive.quantile=0.99)

batch5@meta.data$Genotype<-batch5@meta.data$HTO_maxID
batch5@meta.data$Genotype<-recode(batch5@meta.data$Genotype,"M-HTO-1"="Kbk",
                                  "M-HTO-2"="Kbk","M-HTO-3"="WT","M-HTO-4"="WT","M-HTO-5"="WT","M-HTO-6"="Kbk",
                                  "M-HTO-7"="Kbk","M-HTO-8"="WT")
batch5@meta.data$Genotype<-factor(batch5@meta.data$Genotype,levels=c("WT","Kbk"))



#Add Age column to the batch5 metadata
#All the mice are of 3 month age
batch5@meta.data$Age<-batch5@meta.data$HTO_maxID
batch5@meta.data$Age<-recode(batch5@meta.data$Age,"M-HTO-1"="3m",
                             "M-HTO-2"="3m","M-HTO-3"="3m","M-HTO-4"="3m","M-HTO-5"="3m","M-HTO-6"="1m",
                             "M-HTO-7"="1m","M-HTO-8"="1m")
batch5@meta.data$Age<-factor(batch5@meta.data$Age,levels=c("1m","3m"))

#Add treatment group to the batch5 metadata

batch5@meta.data$Treatment<-batch5@meta.data$HTO_maxID
batch5@meta.data$Treatment<-recode(batch5@meta.data$Treatment,"M-HTO-1"="DMSO_PEG",
                                   "M-HTO-2"="DMSO_PEG","M-HTO-3"="DCF","M-HTO-4"="DMSO_PEG","M-HTO-5"="DCF","M-HTO-6"="Un",
                                   "M-HTO-7"="Un","M-HTO-8"="Un")
batch5@meta.data$Treatment<-factor(batch5@meta.data$Treatment,levels=c("Un","DMSO_PEG","DCF"))



#Add SEX variable to the batch5 metadata
batch5@meta.data$SEX<-batch5@meta.data$HTO_maxID
batch5@meta.data$SEX<-recode(batch5@meta.data$SEX,"M-HTO-1"="F",
                             "M-HTO-2"="M","M-HTO-3"="F","M-HTO-4"="F","M-HTO-5"="M","M-HTO-6"="F",
                             "M-HTO-7"="M","M-HTO-8"="F")
batch5@meta.data$SEX<-factor(batch5@meta.data$SEX,levels=c("M","F"))


#Add sample name(mouse name) to the batch5 metadata
batch5@meta.data$sample_name<-batch5@meta.data$HTO_maxID
batch5@meta.data$sample_name<-recode(batch5@meta.data$sample_name,"M-HTO-1"="Kbk_3m_DMSO_PEG_F",
                                     "M-HTO-2"="Kbk_3m_DMSO_PEG_M","M-HTO-3"="WT_3m_DCF_F","M-HTO-4"="WT_3m_DMSO_PEG_F","M-HTO-5"="WT_3m_DCF_M","M-HTO-6"="Kbk_1m_Un_F",
                                     "M-HTO-7"="Kbk_1m_Un_M","M-HTO-8"="WT_1m_Un_F")



#Add unique name(mouse name) to the batch5 metadata
batch5@meta.data$unique_name<-batch5@meta.data$HTO_maxID
batch5@meta.data$unique_name<-recode(batch5@meta.data$unique_name,"M-HTO-1"="Kbk_3m_DMSO_PEG_F_2",
                                     "M-HTO-2"="Kbk_3m_DMSO_PEG_M_2","M-HTO-3"="WT_3m_DCF_F_2","M-HTO-4"="WT_3m_DMSO_PEG_F_2","M-HTO-5"="WT_3m_DCF_M_1","M-HTO-6"="Kbk_1m_Un_F_2",
                                     "M-HTO-7"="Kbk_1m_Un_M_2","M-HTO-8"="WT_1m_Un_F_2")






#Save the batch5 file with all the new variables added
saveRDS(batch5, file="data/batch5.rds")




# load batch6 filtered matrix  --------------------------------------------


batch6.data<-Read10X("data/batch6_filtered_feature_bc_matrix")
batch6<-CreateSeuratObject(batch6.data$`Gene Expression`,project="kbk_batch6")
batch6[["HTO"]]<-CreateAssayObject(batch6.data$`Antibody Capture`)
batch6_HTOs<-as.data.frame(batch6$HTO@counts[1:8,])
batch6[["HTO"]]<-CreateAssayObject(batch6_HTOs)
batch6<-NormalizeData(batch6,assay="HTO",normalization.method="CLR")
batch6<-HTODemux(batch6,positive.quantile=0.99)

batch6@meta.data$Genotype<-batch6@meta.data$HTO_maxID
batch6@meta.data$Genotype<-recode(batch6@meta.data$Genotype,"M-HTO-1"="Kbk",
                                  "M-HTO-2"="Kbk","M-HTO-3"="Kbk","M-HTO-4"="WT","M-HTO-5"="WT","M-HTO-6"="WT",
                                  "M-HTO-7"="WT","M-HTO-8"="WT")
batch6@meta.data$Genotype<-factor(batch6@meta.data$Genotype,levels=c("WT","Kbk"))



#Add Age column to the batch6 metadata
#All the mice are of 3 month age
batch6@meta.data$Age<-batch6@meta.data$HTO_maxID
batch6@meta.data$Age<-recode(batch6@meta.data$Age,"M-HTO-1"="3m",
                             "M-HTO-2"="3m","M-HTO-3"="3m","M-HTO-4"="3m","M-HTO-5"="3m","M-HTO-6"="3m",
                             "M-HTO-7"="3m","M-HTO-8"="1m")
batch6@meta.data$Age<-factor(batch6@meta.data$Age,levels=c("1m","3m"))

#Add treatment group to the batch6 metadata

batch6@meta.data$Treatment<-batch6@meta.data$HTO_maxID
batch6@meta.data$Treatment<-recode(batch6@meta.data$Treatment,"M-HTO-1"="Un",
                                   "M-HTO-2"="DCF","M-HTO-3"="DMSO_PEG","M-HTO-4"="TCF","M-HTO-5"="TCF","M-HTO-6"="DMSO_PEG",
                                   "M-HTO-7"="DCF","M-HTO-8"="Un")
batch6@meta.data$Treatment<-factor(batch6@meta.data$Treatment,levels=c("Un","DMSO_PEG","DCF","TCF"))



#Add SEX variable to the batch6 metadata
batch6@meta.data$SEX<-batch6@meta.data$HTO_maxID
batch6@meta.data$SEX<-recode(batch6@meta.data$SEX,"M-HTO-1"="M",
                             "M-HTO-2"="M","M-HTO-3"="M","M-HTO-4"="F","M-HTO-5"="M","M-HTO-6"="F",
                             "M-HTO-7"="M","M-HTO-8"="F")
batch6@meta.data$SEX<-factor(batch6@meta.data$SEX,levels=c("M","F"))


#Add sample name(mouse name) to the batch6 metadata
batch6@meta.data$sample_name<-batch6@meta.data$HTO_maxID
batch6@meta.data$sample_name<-recode(batch6@meta.data$sample_name,"M-HTO-1"="Kbk_3m_Un_M",
                                     "M-HTO-2"="Kbk_3m_DCF_M","M-HTO-3"="Kbk_3m_DMSO_PEG_M","M-HTO-4"="WT_3m_TCF_F","M-HTO-5"="WT_3m_TCF_M","M-HTO-6"="WT_3m_DMSO_PEG_F",
                                     "M-HTO-7"="WT_3m_DCF_M","M-HTO-8"="WT_1m_Un_F")





#Add unique name(mouse name) to the batch6 metadata
batch6@meta.data$unique_name<-batch6@meta.data$HTO_maxID
batch6@meta.data$unique_name<-recode(batch6@meta.data$unique_name,"M-HTO-1"="Kbk_3m_Un_M_3",
                                     "M-HTO-2"="Kbk_3m_DCF_M_3","M-HTO-3"="Kbk_3m_DMSO_PEG_M_3","M-HTO-4"="WT_3m_TCF_F_3","M-HTO-5"="WT_3m_TCF_M_2","M-HTO-6"="WT_3m_DMSO_PEG_F_3",
                                     "M-HTO-7"="WT_3m_DCF_M_2","M-HTO-8"="WT_1m_Un_F_3")




#batch6@meta.data$sample_name<-factor(batch6@meta.data$sample_name,levels=c("M","F"))


#Save the batch6 file with all the new variables added
saveRDS(batch6, file="data/batch6.rds")




# load batch7 filtered matrix ---------------------------------------------


batch7.data<-Read10X("data/batch7_filtered_feature_bc_matrix")
batch7<-CreateSeuratObject(batch7.data$`Gene Expression`,project="kbk_batch7")
batch7[["HTO"]]<-CreateAssayObject(batch7.data$`Antibody Capture`)
batch7_HTOs<-as.data.frame(batch7$HTO@counts[1:8,])
batch7[["HTO"]]<-CreateAssayObject(batch7_HTOs)
batch7<-NormalizeData(batch7,assay="HTO",normalization.method="CLR")
batch7<-HTODemux(batch7,positive.quantile=0.99)



batch7@meta.data$Genotype<-batch7@meta.data$HTO_maxID
batch7@meta.data$Genotype<-recode(batch7@meta.data$Genotype,"M-HTO-1"="Kbk",
                                  "M-HTO-2"="WT","M-HTO-3"="Kbk","M-HTO-4"="WT","M-HTO-5"="Kbk","M-HTO-6"="WT",
                                  "M-HTO-7"="Kbk","M-HTO-8"="Kbk")
batch7@meta.data$Genotype<-factor(batch7@meta.data$Genotype,levels=c("WT","Kbk"))



#Add Age column to the batch7 metadata
#All the mice are of 3 month age
batch7@meta.data$Age<-batch7@meta.data$HTO_maxID
batch7@meta.data$Age<-recode(batch7@meta.data$Age,"M-HTO-1"="1m",
                             "M-HTO-2"="3m","M-HTO-3"="1m","M-HTO-4"="3m","M-HTO-5"="3m","M-HTO-6"="1m",
                             "M-HTO-7"="3m","M-HTO-8"="3m")
batch7@meta.data$Age<-factor(batch7@meta.data$Age,levels=c("1m","3m"))

#Add treatment group to the batch7 metadata

batch7@meta.data$Treatment<-batch7@meta.data$HTO_maxID
batch7@meta.data$Treatment<-recode(batch7@meta.data$Treatment,"M-HTO-1"="Un",
                                   "M-HTO-2"="TCF","M-HTO-3"="Un","M-HTO-4"="Un","M-HTO-5"="Un","M-HTO-6"="Un",
                                   "M-HTO-7"="DCF","M-HTO-8"="DMSO_PEG")
batch7@meta.data$Treatment<-factor(batch7@meta.data$Treatment,levels=c("Un","DMSO_PEG","DCF","TCF"))



#Add SEX variable to the batch7 metadata
batch7@meta.data$SEX<-batch7@meta.data$HTO_maxID
batch7@meta.data$SEX<-recode(batch7@meta.data$SEX,"M-HTO-1"="F",
                             "M-HTO-2"="M","M-HTO-3"="M","M-HTO-4"="F","M-HTO-5"="F","M-HTO-6"="M",
                             "M-HTO-7"="F","M-HTO-8"="F")
batch7@meta.data$SEX<-factor(batch7@meta.data$SEX,levels=c("M","F"))


#Add sample name(mouse name) to the batch7 metadata
batch7@meta.data$sample_name<-batch7@meta.data$HTO_maxID
batch7@meta.data$sample_name<-recode(batch7@meta.data$sample_name,"M-HTO-1"="Kbk_1m_Un_F",
                                     "M-HTO-2"="WT_3m_TCF_M","M-HTO-3"="Kbk_1m_Un_M","M-HTO-4"="WT_3m_Un_F","M-HTO-5"="Kbk_3m_Un_F","M-HTO-6"="WT_1m_Un_M",
                                     "M-HTO-7"="Kbk_3m_DCF_F","M-HTO-8"="Kbk_3m_DMSO_PEG_F")


#Add unique name(mouse name) to the batch7 metadata
batch7@meta.data$unique_name<-batch7@meta.data$HTO_maxID
batch7@meta.data$unique_name<-recode(batch7@meta.data$unique_name,"M-HTO-1"="Kbk_1m_Un_F_3",
                                     "M-HTO-2"="WT_3m_TCF_M_3","M-HTO-3"="Kbk_1m_Un_M_3","M-HTO-4"="WT_3m_Un_F_3","M-HTO-5"="Kbk_3m_Un_F_3","M-HTO-6"="WT_1m_Un_M_3",
                                     "M-HTO-7" = "Kbk_3m_DCF_F_3","M-HTO-8"="Kbk_3m_DMSO_PEG_F_3")








#Save the batch7 file with all the new variables added
saveRDS(batch7, file="data/batch7.rds")




# load batch8 filtered matrix ---------------------------------------------



#batch8

batch8.data<-Read10X("data/batch8_filtered_feature_bc_matrix")
batch8<-CreateSeuratObject(batch8.data$`Gene Expression`,project="kbk_batch8")#batch9 is rerun of batch8
batch8[["HTO"]]<-CreateAssayObject(batch8.data$`Antibody Capture`)
batch8_HTOs<-as.data.frame(batch8$HTO@counts[1:4,])
batch8[["HTO"]]<-CreateAssayObject(batch8_HTOs)
batch8<-NormalizeData(batch8,assay="HTO",normalization.method="CLR")
batch8<-HTODemux(batch8,positive.quantile=0.99)

#Add genotype to the batch8
batch8@meta.data$Genotype<-batch8@meta.data$HTO_maxID
batch8@meta.data$Genotype<-recode(batch8@meta.data$Genotype,"M-HTO-1"="WT",
                                  "M-HTO-2"="WT","M-HTO-3"="WT","M-HTO-4"="WT")
batch8@meta.data$Genotype<-factor(batch8@meta.data$Genotype,levels="WT")



#1Add Age column to the batch8 metadata
#All the mice are of 3 month age
batch8@meta.data$Age<-batch8@meta.data$HTO_maxID
batch8@meta.data$Age<-recode(batch8@meta.data$Age,"M-HTO-1"="3m",
                             "M-HTO-2"="3m","M-HTO-3"="3m","M-HTO-4"="3m")
batch8@meta.data$Age<-factor(batch8@meta.data$Age,levels=c("3m"))

#Add treatment group to the batch8 metadata

batch8@meta.data$Treatment<-batch8@meta.data$HTO_maxID
batch8@meta.data$Treatment<-recode(batch8@meta.data$Treatment,"M-HTO-1"="DCF",
                                   "M-HTO-2"="DCF","M-HTO-3"="DMSO_PEG","M-HTO-4"="DMSO_PEG")
batch8@meta.data$Treatment<-factor(batch8@meta.data$Treatment,levels=c("DMSO_PEG","DCF"))



#Add SEX variable to the batch8 metadata
batch8@meta.data$SEX<-batch8@meta.data$HTO_maxID
batch8@meta.data$SEX<-recode(batch8@meta.data$SEX,"M-HTO-1"="F",
                             "M-HTO-2"="M","M-HTO-3"="M","M-HTO-4"="M")
batch8@meta.data$SEX<-factor(batch8@meta.data$SEX,levels=c("M","F"))


#Add sample name(mouse name) to the batch8 metadata
batch8@meta.data$sample_name<-batch8@meta.data$HTO_maxID
batch8@meta.data$sample_name<-recode(batch8@meta.data$sample_name,"M-HTO-1"="WT_3m_DCF_F",
                                     "M-HTO-2"="WT_3m_DCF_M","M-HTO-3"="WT_3m_DMSO_PEG_M","M-HTO-4"="WT_3m_DMSO_PEG_M")






#Add unique name(mouse name) to the batch8 metadata
batch8@meta.data$unique_name<-batch8@meta.data$HTO_maxID
batch8@meta.data$unique_name<-recode(batch8@meta.data$unique_name,"M-HTO-1"="WT_3m_DCF_F_3",
                                     "M-HTO-2"="WT_3m_DCF_M_3","M-HTO-3"="WT_3m_DMSO_PEG_M_2","M-HTO-4"="WT_3m_DMSO_PEG_M_3")









#batch8@meta.data$sample_name<-factor(batch8@meta.data$sample_name,levels=c("M","F"))


#Save the batch8 file with all the new variables added
saveRDS(batch8, file="data/batch8.rds")











# #merge all batches ------------------------------------------------------



kbk_combined<-merge(batch1,y=c(batch2,batch3,batch4,batch5,batch6,batch7,batch8),add.cell.id=c("B1","B2","B3","B4","B5","B6","B7","B8"),project ="kbk_all_batches")
kbk_combined



# join the layers ---------------------------------------------------------

kbk_combined<-JoinLayers(kbk_combined)







#Add percentage of mitochondrial genes in kbk_combined ------------------

kbk_combined[["percent.mt"]]<-PercentageFeatureSet(kbk_combined,pattern="^mt")


# Add number of genes per UMI for each cell to metadata -------------------

kbk_combined$log10GenesPerUMI <- log10(kbk_combined$nFeature_RNA) / log10(kbk_combined$nCount_RNA)



# Add additional information to the metadata ------------------------------

#extract meta.data and assign to a varaible called metadata
metadata<-kbk_combined@meta.data


# Add cell IDs to metadata
metadata$cells <- rownames(metadata)


# Create  column named "batches"; it is designated as sample

metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^B1_"))] <- "Batch1"
metadata$sample[which(str_detect(metadata$cells, "^B2_"))] <- "Batch2"
metadata$sample[which(str_detect(metadata$cells, "^B3_"))] <- "Batch3"
metadata$sample[which(str_detect(metadata$cells, "^B4_"))] <- "Batch4"
metadata$sample[which(str_detect(metadata$cells, "^B5_"))] <- "Batch5"
metadata$sample[which(str_detect(metadata$cells, "^B6_"))] <- "Batch6"
metadata$sample[which(str_detect(metadata$cells, "^B7_"))] <- "Batch7"
metadata$sample[which(str_detect(metadata$cells, "^B8_"))] <- "Batch8"



# Add metadata back to Seurat object --------------------------------------

kbk_combined@meta.data<-metadata



# Quality control ---------------------------------------------------------

#Quality control step


kbk_combined@meta.data$group<-paste(kbk_combined@meta.data$Genotype,kbk_combined@meta.data$Age,kbk_combined@meta.data$Treatment,sep=".")

# save the new kbk_combined  with joinlayers rds object ------------------------------------
## save the  object to the data_seurat5 folder -------------------------------

saveRDS(kbk_combined,file="data/kbk_combined.rds")











################################################################################
