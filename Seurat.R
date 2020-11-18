# Loading Library

library(dplyr)

library(Seurat)

#Loading Data

DSS.data<-Read10X(data.dir = "~/3x/raw_gene_bc_matrix/mm10/")

DSS <- CreateSeuratObject(raw.data = DSS.data, project = "Mouse", min.cells = 5)

DSS@meta.data$Stimulation <- "DSS"

rm(DSS.data)

# Quality Control

mito.genes <- grep(pattern = "^mt-", x = rownames(x = DSS@raw.data), value = TRUE)

percent.mito <- Matrix::colSums(DSS@raw.data[mito.genes, ])/Matrix::colSums(DSS@raw.data)

DSS <- AddMetaData(object = DSS, metadata = percent.mito, col.name = "percent.mito")

GenePlot(object = DSS, gene1 = "nUMI", gene2 = "percent.mito")

GenePlot(object = DSS, gene1 = "nUMI", gene2 = "nGene")

DSS <- FilterCells(DSS, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(5000, 0.17))

# Variable Feature Identification

DSS <- NormalizeData(DSS)

DSS <- ScaleData(DSS, display.progress = F,vars.to.regress = c("nUMI", "percent.mito"))

DSS <- FindVariableGenes(DSS, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# Dimension Reduction

DSS<- RunPCA(DSS, pc.genes = DSS@var.genes, do.print = FALSE, pcs.compute =  50,reduction.name = "pca",reduction.key = "PC")

DSS <- JackStraw(object = DSS, num.replicate = 100, display.progress = FALSE,num.pc = 50,do.par=TRUE,num.cores = 4)   

JackStrawPlot(object = DSS, PCs = 1:50)

DSS <- RunTSNE(DSS, reduction = "pca", dims.use = c(1:34),  dim.embed = 2,do.fast = T,force.recalc=TRUE)

# Identify Clusters

DSS<- FindClusters(object = DSS, reduction.type = "pca", dims.use = c(1:34), resolution = 0.35, print.output = 0, save.SNN = TRUE)

TSNEPlot(DSS)

# Identify Marker Genes

cluster.markers.s <- FindAllMarkers(object = DSS, only.pos = FALSE, min.pct = 0.1, thresh.use = 0.25, grouping.var = "ident",print.bar = FALSE)

write.csv(cluster.markers.s,file="~/Markers.csv")

# Rename Clusters

new.cluster.ids<-c("Cluster 1","Cluster 2","Cluster 3","Cluster 5","Myofibroblast","MHC II+ Stromal Cell","Mesothelial Cell","Telocyte","Epithelial Cell","Interstitial Cell of Cajal")

current.cluster.ids <- levels(DSS@ident)

DSS@ident <- plyr::mapvalues(x = DSS@ident, from = current.cluster.ids, to = new.cluster.ids)

# Draw ViolinPlot

Data<-SubsetData(DSS,ident.use=c("Cluster 1","Cluster 2","Cluster 3","Cluster 5"))

VlnPlot(object =Data, features = c("Cxcl14","Agt","Figf","Ptn"), use.raw = TRUE, y.log = TRUE,point.size.use = 0,x.lab.rot = TRUE,nCol=4,size.x.use =10,remove.legend = TRUE)

# Draw FeaturePlot

FeaturePlot(object = DSS, features.plot = c("Rspo1"), cols.use = c("grey93", "blue"))

# Draw Heatmap

Data1<-SubsetData(DSS,ident.use=c("Cluster 2","Cluster 5"))

DoHeatmap(object = pbmc, genes.use = "S1pr3", "Tnc",
                        "Hoxb9",
                        "Ces1g",
                        "Sdc1",
                        "Ddit4",
                        "Pappa",
                        "Mab21l",
                        "Cxcl14",
                        "Slc40a1",
                        "Lbh",
                        "Mafb",
                        "Ch25h",
                        "Frzb",
                        "Cxcl9",
                        "Rspo1",
                        "Pcolce2",
                        "Ebf1",
                        "Inmt",
                        "Sfrp4",
                        "Fez1",
                        "Cilp",
                        "Penk",
                        "Lrrn4cl",
                        "Pi16",
                        "Cd81",
                        "Figf",
                        "Il33",
                        "Cd34", slim.col.label = TRUE, remove.key = TRUE)

# Calculate Cluster Gene Average Expression

Data2<-AverageExpression(DSS,use.raw =T)
write.csv(Data2,file="~/Expression.csv")

