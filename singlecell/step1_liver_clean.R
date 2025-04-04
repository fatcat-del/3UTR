options(future.globals.maxSize = 56000 * 1024^2)

library(Seurat)
library(cowplot)
library(harmony)
library(ggplot2)
library(ggpubr)
library(scDblFinder)
library(tidyverse)

colorset <- c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72","#B17BA6","#FF7F00",
              "#FDB462","#E7298A","#E78AC3","#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02",
              "#7570B3","#BEAED4","#666666","#999999","#aa8282","#d4b7b7","#8600bf","#ba5ce3",
              "#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d",
              "#ffff00", "#0077BB", "#33BBEE", "#009988", "#EE7733", "#CC3311", "#EE3377", "#BBBBBB",
              "#027EB6", "#746FB2", "#9651A0", "#C8008F", "#ee64a4",
              "#EE0220", "#D93F00", "#795549", "#6F7C4D", "#008A25",
              "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
              "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

#merge control HFAT single cell 
path.dir <- "/data/fengw/utr/code/data/single_cell_raw"
setwd(path.dir)

ob.list <- list()
samples <- c("control", "HFD")

for (each in samples){
  h5.path <- paste0(path.dir, each, "/filtered_feature_bc_matrix")
  print (h5.path)
  ob <- Read10X(h5.path)
  ob <- CreateSeuratObject(ob, project = each, min.cells = 3, min.features = 200)
  ob <- RenameCells(ob, add.cell.id = each)
  ob$stim <- each
  ob[["percent.mt"]] <- PercentageFeatureSet(ob, pattern = "^mt-")
  ob$log10GenesPerUMI <- log10(ob$nFeature_RNA) / log10(ob$nCount_RNA)
  ob <- subset(ob, subset = nFeature_RNA>=250&percent.mt<50&nCount_RNA>500&log10GenesPerUMI>0.8)
  ob.list[[each]] <- ob
}

dir.create("Integrate")
setwd("./Integrate")

names(ob.list) <- samples
system.time(saveRDS(ob.list, file = "liver.rds"))

sce <- merge(ob.list[[1]], ob.list[[2]])

#remove cell doublets
sce <- JoinLayers(sce)
sce <- NormalizeData(sce)
sce <- as.SingleCellExperiment(sce)
sce <- scDblFinder(sce, samples="stim", dbr.sd=1)
table(sce$scDblFinder.class)
sce <- as.Seurat(sce)
sce <- subset(sce, subset = scDblFinder.class=="singlet")
saveRDS(sce, "liversc.rds")

sce <- SCTransform(sce, vars.to.regress = c("nCount_RNA","percent.mt"))
sce <- RunPCA(sce, npcs=50, verbose=FALSE)
sce <- RunHarmony(sce, group.by.vars="orig.ident", 
                  assay.use="SCT", max.iter.harmony = 20)

sce <- RunTSNE(sce, reduction="harmony", dims=1:50) %>% 
  RunUMAP(reduction="harmony", dims=1:50 ,seed.use = 0)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:50) %>% 
  FindClusters( resolution = 0.6, random.seed = 0)

markergene <- c("Pecam1", "Stab2", "Stab1", #LSECs
                "Alb", "Apoc3", "Hnf4a", #HEPs
                "Cd163", "Clec4f", #KCs
                "Dcn", "Lrat", #HSCs
                "Cd3g", "Cd3e", "Nkg7", "Ncr1", #T cells/NK cells
                "Cd79a", "Ms4a1", #B cells
                "Ly6g", "S100a8", #Neutophils
                "Itgam", "Ccr2", #Monocytes
                "Bst2", #pDCs
                "Sox9", "Krt18", #cholangiocytes
                "Vwf", #macrovascularECs
                "Fbln1", "Msln" #myofibroblast
)


DimPlot(sce, cols = colorset, group.by = "seurat_clusters",label = T,repel = T,label.size = 6)

library(scRNAtoolVis)

jjDotPlot(object = sce,
          gene = markergene,
          id = "seurat_clusters",
          ytree = F)



sce$seurat_clusters <- factor(sce$seurat_clusters, 
                              levels = c(0,3,4,
                                         2,5,6,7,11,12,13,18,
                                         1,10,
                                         9,
                                         8,16,
                                         15,
                                         17,
                                         14,
                                         20,
                                         21,
                                         19,
                                         22))
jjDotPlot(object = sce,
          gene = markergene,
          id = "seurat_clusters",
          ytree = F)


sce$celltype <- ""
sce$celltype <- ifelse(sce$seurat_clusters %in% c(0,3,4), "LSECs",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(2,5,6,7,11,12,13,18), "HEPs",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(1,10), "Kupffer cells",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(9), "HSCs",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(8,16), "T/NK cells",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(15), "B cells",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(17), "Neutrophils",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(14), "Monocytes",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(20), "pDCs",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(21), "Cholangiocytes",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(19), "macrovascularECs",sce$celltype)
sce$celltype <- ifelse(sce$seurat_clusters %in% c(22), "myofibroblast",sce$celltype)
sce$celltype <- factor(sce$celltype, 
                       levels = c("LSECs", "HEPs", "Kupffer cells",
                                  "HSCs", "T/NK cells",
                                  "B cells", "Neutrophils", "Monocytes",
                                  "pDCs","Cholangiocytes",
                                  "macrovascularECs","myofibroblast"))
jjDotPlot(object = sce,
          gene = markergene,
          id = "celltype",
          ytree = F)

saveRDS(sce,"liversc_clean.rds")
