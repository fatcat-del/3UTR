library(Seurat)
library(SCAPE)
library(magrittr)
library(ggplot2)


setwd("/data/fengw/utr/code/data/")
find_df <- function(genes, sceobj, ident1, ident2){
  temp <- c()
  genes <- intersect(rownames(sceobj), genes)
  celltype <- as.data.frame(table(sceobj$celltype, sceobj$orig.ident))
  t1 <- celltype[celltype$Var2 %in% c(ident1,ident2),]
  t1 <- t1[t1$Freq >3, ]
  celltype <- intersect(t1[t1$Var2 %in% c(ident1),]$Var1, t1[t1$Var2 %in% c(ident2),]$Var1 )
  print(celltype)
  for (i in celltype){
    print(i)
    subsce <- sceobj[,sceobj$celltype %in% c(i)]
    #subsce <- PrepSCTFindMarkers(subsce)
    df <- Seurat::FindMarkers(subsce,logfc.threshold = 0,assay= "apa",
                              group.by = "orig.ident", ident.1 = ident1, min.pct = 0,
                              ident.2 = ident2, features = genes)
    df$celltype <- i
    df$gene <- rownames(df)
    temp <- rbind(temp,df)
  }
  return(temp)
}



#Seurat V5 requires modifying the loadData function of SCAPE
loadData <- function (fileList, collapsePa, matrix = FALSE, cores = 1) 
{
  pa_group_info <- data.table::fread(collapsePa)
  pa_group_info$pa <- paste(pa_group_info$chrom, pa_group_info$collapse_pa, 
                            pa_group_info$strand, sep = ":")
  objs <- parallel::mclapply(names(fileList), function(fileLabel) {
    message(fileLabel)
    tmp <- data.table::fread(fileList[[fileLabel]])
    mtx <- as.matrix(tmp[, -1])
    rownames(mtx) <- tmp[["V1"]]
    rm(tmp)
    pa <- data.frame(V1 = rownames(mtx), stringsAsFactors = F)
    pa$V2 <- plyr::mapvalues(from = pa_group_info$pa_site, 
                             to = pa_group_info$pa, x = pa$V1, warn_missing = F)
    pa_dup <- pa[pa$V2 %in% pa$V2[duplicated(pa$V2)], ]
    mtx_dup <- do.call(rbind, lapply(split(pa_dup$V1, pa_dup$V2), 
                                     function(x) {
                                       colSums(mtx[x, ])
                                     }))
    mtx <- rbind(mtx[setdiff(pa$V1, pa_dup$V1), ], mtx_dup)
    rownames(mtx) <- plyr::mapvalues(from = pa$V1, to = pa$V2, 
                                     x = rownames(mtx), warn_missing = F)
    colnames(mtx) <- paste(fileLabel, colnames(mtx), sep = ".")
    dat <- Seurat::CreateSeuratObject(mtx, project = fileLabel, 
                                      names.delim = "[.]")
    return(dat)
  }, mc.cores = cores)
  objs <- Reduce(merge, objs)
  if (isTRUE(matrix)) {
    DefaultAssay(objs) <- "RNA"
    objs <- JoinLayers(objs, overwrite = T)
    return(Seurat::GetAssayData(objs, "RNA", slot = "counts"))
  }
  return(objs)
}


#load scRNA data
sce <- readRDS("./single_cell_raw/Integrate/liversc_clean.rds")

#get APA counts matrix
exp_file <- c("./apa_data/Chow/pasite.csv.gz","./apa_data/HFD/pasite.csv.gz")
names(exp_file) <- c("control", "HFD")
collapse_pa <- c("D:/文章修改/LPGAT1/apa_data/collapse_pa.tsv.gz")

pa_mtx <- loadData(
  fileList = exp_file,
  collapsePa = collapse_pa,
  matrix = TRUE,
  cores = 1
)

#Only these pA sites whcih expressed in more than 50 cell were kept.
binary_filter <- Matrix::rowSums(+pa_mtx)
pa_mtx <- pa_mtx[binary_filter > 50, ]

colnames(pa_mtx) <- gsub("\\.", "_", colnames(pa_mtx))

sce[['apa']] <- CreateAssayObject(pa_mtx[, colnames(sce)])
sce <- NormalizeData(sce, assay = 'apa')
sce <- ScaleData(sce, assay = 'apa')

saveRDS(sce, "./single_cell_raw/Integrate/liver.apa.rds")
