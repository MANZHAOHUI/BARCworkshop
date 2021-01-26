# R v3.6.0 via UNC OnDemand
# 8 hrs with 128GB RAM, n=1 node

#Set working directory
setwd("/pine/scr/j/m/jmsimon/BARCworkshop_012621")

#Libraries
library(rlang,lib="/proj/jmsimon/Rlibs36")
library(vctrs,lib="/proj/jmsimon/Rlibs36")
library(backports,lib="/proj/jmsimon/Rlibs36")
library(cli,lib="/proj/jmsimon/Rlibs36")
library(Seurat,lib="/proj/jmsimon/Rlibs36")
library(rstudioapi,lib="/proj/jmsimon/Rlibs36")
library(withr,lib="/proj/jmsimon/Rlibs36")
library(tidyverse,lib="/proj/jmsimon/Rlibs36")
library(ggplot2,lib="/proj/jmsimon/Rlibs36")
library(sctransform,lib="/proj/jmsimon/Rlibs36")
library(farver,lib="/proj/jmsimon/Rlibs36")
library(usethis,lib="/proj/jmsimon/Rlibs36")


FAhealthy = Read10X(data.dir="/proj/esmblab/HTSF/201214_UNC41-A00434_0203_AHWWW7DRXX/HWWW7DRXX/counting/1_FA_HEALTHY/outs/filtered_feature_bc_matrix/")
FAinjured = Read10X(data.dir="/proj/esmblab/HTSF/201214_UNC41-A00434_0203_AHWWW7DRXX/HWWW7DRXX/counting/2_FA_INJURED/outs/filtered_feature_bc_matrix/")
OzoneHealthy = Read10X(data.dir="/proj/esmblab/HTSF/201214_UNC41-A00434_0203_AHWWW7DRXX/HWWW7DRXX/counting/3_OZONE_HEALTHY/outs/filtered_feature_bc_matrix/")
OzoneInjured = Read10X(data.dir="/proj/esmblab/HTSF/201214_UNC41-A00434_0203_AHWWW7DRXX/HWWW7DRXX/counting/4_OZONE_INJURED/outs/filtered_feature_bc_matrix/")

FAH.seurat = CreateSeuratObject(FAhealthy,project="Bahnson_10X", min.cells = 100, min.features = 100)
FAH.seurat = RenameCells(FAH.seurat,add.cell.id = "FAH")

FAI.seurat = CreateSeuratObject(FAinjured,project="Bahnson_10X", min.cells = 100, min.features = 100)
FAI.seurat = RenameCells(FAI.seurat,add.cell.id = "FAI")

OH.seurat = CreateSeuratObject(OzoneHealthy,project="Bahnson_10X", min.cells = 100, min.features = 100)
OH.seurat = RenameCells(OH.seurat,add.cell.id = "OH")

OI.seurat = CreateSeuratObject(OzoneInjured,project="Bahnson_10X", min.cells = 100, min.features = 100)
OI.seurat = RenameCells(OI.seurat,add.cell.id = "OI")

#Compute mitochondrial contamination and subset
FAH.seurat <- PercentageFeatureSet(object = FAH.seurat, pattern = "^mm10---mt-", col.name = "percent.mt")
FAH.seurat <- subset(FAH.seurat, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)

FAI.seurat <- PercentageFeatureSet(object = FAI.seurat, pattern = "^mm10---mt-", col.name = "percent.mt")
FAI.seurat <- subset(FAI.seurat, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)

OH.seurat <- PercentageFeatureSet(object = OH.seurat, pattern = "^mm10---mt-", col.name = "percent.mt")
OH.seurat <- subset(OH.seurat, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)

OI.seurat <- PercentageFeatureSet(object = OI.seurat, pattern = "^mm10---mt-", col.name = "percent.mt")
OI.seurat <- subset(OI.seurat, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)

dim(FAH.seurat)
dim(FAI.seurat)
dim(OH.seurat)
dim(OI.seurat)

bahnson.list = as.list(c(FAH.seurat,FAI.seurat,OH.seurat,OI.seurat))
names(bahnson.list) = c("FAhealthy","FAinjured","OzoneHealthy","OzoneInjured")

#Normalize all samples
for (i in 1:length(x = bahnson.list)) {
	bahnson.list[[i]] <- PercentageFeatureSet(object = bahnson.list[[i]], pattern = "^mm10---mt-", col.name = "percent.mt")
	bahnson.list[[i]] <- SCTransform(bahnson.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE,return.only.var.genes=FALSE)
}

#Save image at point of object setup
save.image("Bahnson_10X_cellranger-seurat_scTransform_012621.Rdata")


# Set up integration, then integrate data
### NOTE: I turned verbose to TRUE so I can monitor the process
options(future.globals.maxSize=5242880000)		# This was needed otherwise I got an error about something being too large/above limits. This particular value is 5000*1024^2, based off of this page: https://stackoverflow.com/questions/40536067/how-to-adjust-future-global-maxsize-in-r
bahnson.features <- SelectIntegrationFeatures(object.list = bahnson.list, nfeatures = 5000)
bahnson.list <- PrepSCTIntegration(object.list = bahnson.list, anchor.features = bahnson.features, verbose = TRUE)
bahnson.anchors <- FindIntegrationAnchors(object.list = bahnson.list, normalization.method = "SCT", anchor.features = bahnson.features, verbose = TRUE)
bahnson.integrated <- IntegrateData(anchorset = bahnson.anchors, normalization.method = "SCT", verbose = TRUE)

# Save image at point of object integration
saveRDS(bahnson.integrated,"Bahnson_10X_cellranger-seurat_scTransform_012621_integrated.rds")

# Reset metadata, since it got wiped during integration
condition = gsub("_.+","",colnames(bahnson.integrated))
rep = rep(1,length(bahnson.integrated@meta.data$orig.ident))
bahnson.integrated@meta.data$condition = condition	
bahnson.integrated@meta.data$rep = rep

#Run PCA
bahnson.integrated <- RunPCA(bahnson.integrated, verbose = FALSE, npcs = 100)
bahnson.integrated <- RunUMAP(bahnson.integrated, dims = 1:100)
bahnson.integrated <- FindNeighbors(bahnson.integrated, dims = 1:100, verbose = FALSE)
bahnson.integrated <- FindClusters(bahnson.integrated, verbose = FALSE, resolution = 2, algorithm=2)



#Clustering plots
pdf("Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered.pdf",width=12)
DimPlot(bahnson.integrated,reduction = "umap", label = TRUE)
DimPlot(bahnson.integrated,reduction = "umap", group.by = "condition")
DimPlot(bahnson.integrated,reduction = "umap", group.by = "rep")
dev.off()

#Comparison plots
pdf("Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_comparison.pdf",width=20)
DimPlot(bahnson.integrated,reduction = "umap", label = TRUE, split.by = "condition")
DimPlot(bahnson.integrated,reduction = "umap", label = TRUE, split.by = "rep")
dev.off()


# QC plots
#Make violin plots
features = c("nCount_RNA","nFeature_RNA","percent.mt")
plotwidth = 50
plotheight = 15
pdf("Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_vlnplot_QC.pdf",width=plotwidth,height=plotheight)
for (i in features) {
    print(VlnPlot(bahnson.integrated, features = i, split.by = "condition"))
}
dev.off()



# Compute cell proportions by cluster by replicate

condition <- unique(bahnson.integrated@meta.data$condition)
rep <- unique(bahnson.integrated@meta.data$rep)

comp.prop <- vector()
names.prop <- vector()
for (i in 1:length(condition)) {
        #Subset bahnson.integrated object
        sub <- bahnson.integrated[,bahnson.integrated$condition == condition[i]]
        
        #Generate proportion table by replicate
        sub.table <- prop.table(table(Idents(sub), sub$rep), margin = 2)
        
        #Create table to add to comp (code necessary when not all clusters are in subset)
        sub.final <- array(0, dim=c(length(levels(bahnson.integrated)),length(unique(sub$rep))))
        rownames(sub.final) <- levels(bahnson.integrated)
        colnames(sub.final) <- c(1:length(unique(sub$rep)))
        sub.final[rownames(sub.final) %in% rownames(sub.table),] <- sub.table
        
        #Add to compilation and colnames vector
        comp.prop <- cbind(comp.prop,sub.final)
        names.prop <- c(names.prop, paste0(condition[i],"_",colnames(sub.table)))
}

colnames(comp.prop) <- names.prop
write.table(comp.prop,"Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_proportions_replicates.txt",col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")


# import using tidy's readr function
tbl = as_tibble(rownames_to_column(as.data.frame(comp.prop),var = "Cluster"))


# turn data matrix into tidy tibble
tibble = tbl %>% gather(key="Sample",value="Proportion",-Cluster) %>% separate(Sample,sep="_(?=[:digit:])",into=c("Condition","Replicate"),convert=T)

# boxplot
pdf("Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_proportions_replicates_boxplots.pdf",width=14,height=10)
tibble %>%
    mutate(Cluster = fct_relevel(Cluster, unique(tibble$Cluster))) %>%
    mutate(Genotype = fct_relevel(Condition, "FAH", "FAI","OH","OI")) %>%
    ggplot(aes(fill = Genotype, x= Genotype, y=Proportion)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(pch = 21, position = position_jitterdodge(jitter.width=0.9)) +
    facet_wrap(~Cluster, scales="free") +
    scale_fill_manual(values=c("#f01000","#F8766D","#405900","#7CAE00")) +
    xlab("")
dev.off()

# Save progress
saveRDS(bahnson.integrated,"Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered.rds")

# Find markers of each cluster
#for (i in levels(bahnson.integrated)) {		# This will take many hours to run
for (i in 0:4) {								# Just do first 5 clusters for today's workshop
    markers <- FindConservedMarkers(bahnson.integrated, assay="RNA", slot="counts", ident.1=i,grouping.var="condition",only.pos=T)
    write.table(markers,paste0("Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_ConservedMarkers_RNA_",i,".txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
}

# Extract out top markers for dotplot

markerfiles = Sys.glob("Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_*Markers_RNA_*.txt")
unionmarkers = vector()

for (i in 1:length(markerfiles)){
	markers = read_tsv(markerfiles[i]) %>% 
		select(X1) %>% 
		head(n=5) %>% 
		pull(X1)
	unionmarkers = c(unionmarkers,markers)
}
unionmarkers = unique(unionmarkers)

pdf("Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_dotplot.pdf",width=30,height=12)
DotPlot(bahnson.integrated,features=unique(unionmarkers),assay="SCT",cols = c("blue","red")) + RotatedAxis()
dev.off()


# Save progress for today

saveRDS(bahnson.integrated,"Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_withMarkers.rds")
save.image("Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_withMarkers.Rdata")



# Now plot dotplot with Bahnson genes of interest

bahnson.integrated = readRDS("Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_withMarkers.rds")

goi = unique(c("Acta2","Myh11","Tagln","Cnn1","Cdh5","Pecam1","Icam1","Vcam1","Ptprc","Cd68","Adgre1","Ly6c1","Cd2","Cd3g","S100a4","Ddr2","Thy1","Pdgfrb","Pdgfra","Ly6a","Cd34","Mcam"))
goi2 = paste0("mm10---",goi)

pdf("Bahnson_10X_cellranger-seurat-scTransform_012621_integrated_clustered_dotplot_goi.pdf",width=10,height=10)
DotPlot(bahnson.integrated,features=goi2,assay="SCT",cols = c("blue","red")) + RotatedAxis()
dev.off()

