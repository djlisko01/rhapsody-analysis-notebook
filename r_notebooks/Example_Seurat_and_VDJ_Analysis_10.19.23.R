#################################################################################################
library(Seurat)
library(dplyr)
library(ggplot2)
library(devtools)

devtools::install_github("satijalab/seurat-data")
devtools::install_github("satijalab/azimuth")

library(SeuratData)
library(Azimuth)

#Read RDS object
cbmc <- readRDS("~/gitrepos/rhapsody-analysis-notebook/data/gex/copied-FASTQs-expected-cells-30K-v2-rerun-MC_Seurat.rds")

#Subset out undetermined / multiplet cells
Idents(cbmc) <- cbmc@meta.data$Sample_Name
cbmc <- subset(cbmc, idents = c("Undetermined", "Multiplet"), invert = T)

#Calculate Percent.MT
cbmc[["percent.mt"]] <- PercentageFeatureSet(cbmc, pattern = "^MT.")

#Visualize QC metrics / set threhsolds
VlnPlot(cbmc, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0, group.by = "orig.ident")

#Subset high quality cells based on QC thresholds
cbmc <- subset(cbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 7.5)

#Put AbSeq and RNA targets into separate assays
adt <- cbmc@assays$RNA@counts[grep("pAbO", rownames(cbmc)),]
rna <- cbmc@assays$RNA@counts[-(grep("pAbO", rownames(cbmc))),]
cbmc[["ADT"]] <- CreateAssayObject(adt)
cbmc[["RNA"]] <- CreateAssayObject(rna)

#Normalize mRNA and AbSeq Separately
DefaultAssay(cbmc) <- "RNA"
cbmc_normalized <- SCTransform(cbmc)

DefaultAssay(cbmc_normalized) <- "ADT"
cbmc_normalized <- NormalizeData(cbmc_normalized, normalization.method = "CLR")

#Set Default Assay to SCTransformed RNA data
DefaultAssay(cbmc_normalized) <- "SCT"
#Identify Variable Features
cbmc_normalized <- FindVariableFeatures(cbmc_normalized)

#Run PCA 
cbmc_normalized <- RunPCA(cbmc_normalized, verbose = FALSE)

#Generate Elbow Plot
ElbowPlot(cbmc_normalized, ndims = 50)

#Nearest Neighbor Graph Construction
cbmc_normalized <- FindNeighbors(cbmc_normalized, dims = 1:12)

#Run UMAP
cbmc_normalized <- RunUMAP(cbmc_normalized, dims = 1:12)

#Find Clustering
cbmc_normalized <- FindClusters(cbmc_normalized, resolution = 0.5, verbose = F)

#Run Azimuth
cbmc_normalized <- RunAzimuth(cbmc_normalized, reference = "pbmcref")

#Visualize Azimuth Results
DimPlot(cbmc_normalized, reduction = "umap", group.by = "predicted.celltype.l2", label = T, repel = T)

#Perform Differential Expression based on cell type
Idents(cbmc_normalized) <- cbmc_normalized@meta.data$predicted.celltype.l2
cd16_monocyte_markers <- FindMarkers(cbmc_normalized, ident.1 = "CD14 Mono")
cd16_vs_cd14_monocyte_markers <- FindMarkers(cbmc_normalized, ident.1 = "CD14 Mono", ident.2 = "CD16 Mono")

#Perform Differential Expression of an individual cell type by sample
cd14_monocytes <- cbmc_normalized[,which(cbmc_normalized@meta.data$predicted.celltype.l2 == "CD14 Mono")]
Idents(cd14_monocytes) <- cd14_monocytes@meta.data$Sample_Name
samp5_vs_samp6_cd14_monocyte_markers <- FindMarkers(cd14_monocytes, ident.1 = "SampleTag05_hs", ident.2 = "SampleTag06_hs")

#Perform Pathway Analysis
BiocManager::install("clusterProfiler")
library(clusterProfiler)
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

original_gene_list <- cd16_monocyte_markers$avg_log2FC
names(original_gene_list) <- rownames(cd16_monocyte_markers)
gene_list <- na.omit(original_gene_list)
gene_list = sort(original_gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH")



require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


#################################################################################################
#Perform VDJ Analysis
#################################################################################################

#NOTE: All file paths will need to be changed to local paths where your data is stored.
devtools::install_github("ncborcherding/scRepertoire")
library(Seurat)
library(scRepertoire)



full_airr <- read.table(file = "~/gitrepos/rhapsody-analysis-notebook/data/vdj/copied-FASTQs-expected-cells-30K-v2-rerun-MC_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("~/gitrepos/rhapsody-analysis-notebook/data/sample_tags/copied-FASTQs-expected-cells-30K-v2-rerun-MC_Sample_Tag_Calls.csv", header = T, sep = ",")
demo_seurat_object <- readRDS("~/gitrepos/rhapsody-analysis-notebook/data/gex/copied-FASTQs-expected-cells-30K-v2-rerun-MC_Seurat.rds")


tcr_paired <- subset(demo_seurat_object, TCR_Paired_Chains == "True")
bcr_paired <- subset(demo_seurat_object, BCR_Paired_Chains == "True")

#NOTE: YOU WILL NEED TO CHANGE the colnames() function to either "tcr_paired" or "bcr_paired" in order to split out T / B cells. 
#Also be sure to change the directory where you are writing these files if you wish to split T/B cells or they will overwrite.
for(i in 1:length(levels(factor(demo_stcalls$Sample_Name)))){
  ##Identify which cell indices are assigned to which sample tags
  x <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == levels(factor(demo_stcalls$Sample_Name))[i])]
  #Subset the full AIRR file into sample-specific files and remove experimental cell type column
  sub_df <- full_airr[which(full_airr$cell_id %in% x & full_airr$cell_id %in% colnames(tcr_paired)),-2]
  #Write separated sample AIRR formatted files to a folder (screpertoire_data_tcr)
  write.table(x = sub_df, file = paste("~/gitrepos/rhapsody-analysis-notebook/data/vdj/", as.character(levels(factor(demo_stcalls$Sample_Name))[i]), 
                                       "_VDJ_Dominant_Contigs_AIRR.tsv", sep = ""), sep = "\t", quote = F, row.names = F)
}

TCR <- loadContigs(dir = "~/gitrepos/rhapsody-analysis-notebook/data/vdj/", format = "BD")
TCR_combined <- combineTCR(TCR, samples = c("Multiplet", "SampleTag05", "SampleTag06", "Undetermined"), removeNA = F, removeMulti = F, filterMulti = F)

subset <- subsetContig(TCR_combined, name = "sample", 
                       variables = c("SampleTag05", "SampleTag06"))

quantContig(subset, cloneCall="gene", scale = T)
abundanceContig(subset, cloneCall = "gene", scale = F)
lengthContig(subset, cloneCall="aa", chain = "both") 

compareClonotypes(subset, 
                  numbers = 10, 
                  #samples = c("SampleTag_05", "SampeTag_06"), 
                  cloneCall="aa", 
                  graph = "alluvial")

vizGenes(subset, gene = "V", 
         chain = "TRB", 
         plot = "heatmap", 
         order = "variance", 
         scale = F)



clonalHomeostasis(subset, cloneCall = "gene", 
                  cloneTypes = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1))


clonalProportion(subset, cloneCall = "gene",
                 split = c(10, 100, 1000, 10000, 30000, 1e+05)) 

clonalOverlap(subset, 
              cloneCall = "gene+nt", 
              method = "morisita")

clonesizeDistribution(subset, 
                      cloneCall = "gene+nt", 
                      method="ward.D2")

clonalDiversity(subset, 
                cloneCall = "gene", 
                group.by = "sample", 
                #x.axis = "sample", 
                n.boots = 100)

scatterClonotype(subset, 
                 cloneCall ="gene", 
                 x.axis = "SampleTag05", 
                 y.axis = "SampleTag06",
                 dot.size = "total",
                 graph = "proportion")
