#加载分析所需要的包
library(Seurat)
library(hdf5r)
library(tidyverse)
library(readxl)
library(patchwork) 
library(reshape2)
library(ROGUE)
library(DoubletFinder)
library(cluster)
library(future)
library(parallel)
works_num = 10
args<-commandArgs(T)
dir_dta = paste("/home/ud202180737/stAtlas/2.Result/",args[1],sep = "")
file_list <-list.files(path=dir_dta,pattern = "*_filtered_feature_bc_matrix.h5",full.names=T)
file_list
for (i in c(1:length(file_list))) { 
print(file_list[i])
data = Read10X_h5(file_list[i])
#QC
pbmc <- CreateSeuratObject(counts = data, project = args[1], min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.CD3D"]] <- PercentageFeatureSet(pbmc, pattern = "CD3D" )
pbmc[["percent.CD3G"]] <- PercentageFeatureSet(pbmc, pattern = "CD3G" )
min = unname(quantile(pbmc@meta.data$nFeature_RNA, probs = 0.2,na.rm = TRUE))
max = unname(quantile(pbmc@meta.data$nFeature_RNA, probs = 0.95,na.rm = TRUE))
pbmc <- subset(pbmc, subset = nFeature_RNA > min & nFeature_RNA < max  & percent.mt < 15 & percent.CD3D > 0 & percent.CD3G >0 )
#数据标准化
#默认情况下，使用全局缩放归一化方法“LogNormalize”，用总表达量对每个细胞的基因表达式进行归一化，再乘以一个缩放因子(默认为10,000)，
#然后对结果进行log转换。标准化的数值存储在pbmc[["RNA"]]@data中。
plan("multiprocess", workers = works_num)
pbmc <- NormalizeData(pbmc)
#鉴定高变基因
#接下来，计算数据集中表现出高细胞间变异的特征基因(即，它们在某些细胞中高表达，而在其他细胞中低表达)。这些基因有助于突出单细胞数据集中的生物信号。
#用FindVariableFeatures()函数实现。默认情况下，每个数据集返回2000个features 。这些将用于下游分析，如PCA。
plan("multiprocess", workers = works_num)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#数据缩放
#线性变换(“缩放”)，是在PCA降维之前的一个标准预处理步骤。ScaleData()函数功能:
#转换每个基因的表达值，使每个细胞的平均表达为0
#转换每个基因的表达值，使细胞间的方差为1
#此步骤在下游分析中具有相同的权重，因此高表达的基因不会占主导地位
#结果存储在pbmc[["RNA"]]@scale.data中
all.genes <- rownames(pbmc)
plan("multiprocess", workers = works_num)
pbmc <- ScaleData(pbmc, features = all.genes)
plan("multiprocess", workers = works_num)
suppressMessages(pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)))


# Determine percent of variation associated with each PC
#pct <- pbmc[["pca"]]@stdev / sum( pbmc[["pca"]]@stdev) * 100
pct <- pbmc[["pca"]]@stdev
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
co1 <- which(cumu > 75 & pct < 5)[1]
print(co1)
if(is.na(co1))  co1 = 50 else  print("NA is OK")
if(co1>30)  co1 = 30 else  print("Lager is OK")
if(co1<15)  co1 = 15 else  print("Small is OK") 
pcs = co1
pc.num=1:pcs
plan("multiprocess", workers = works_num)
pbmc <- RunUMAP(pbmc, dims = pc.num)
plan("multiprocess", workers = works_num)
if(dim(pbmc@meta.data)[1]<200){perplexity_x = 10}
if(dim(pbmc@meta.data)[1]>=200){perplexity_x = 30}
pbmc <- RunTSNE(pbmc, dims = pc.num, perplexity = perplexity_x)
plan("multiprocess", workers = works_num)
pbmc <- FindNeighbors(pbmc, dims = pc.num)
plan("multiprocess", workers = works_num)
pca.data <- pbmc@reductions$pca@cell.embeddings
dd <- dist(pca.data)


#使用Rouge包确定最佳resolution
se_cluster_function = function(.x){
        print(.x)
        pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = .x)
        length(unique(as.numeric(pbmc$seurat_clusters)))
      }

silhouette_function = function(.x){
        print(.x)
        pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = .x)
        se_cluster <- as.numeric(pbmc$seurat_clusters)
        if(length(unique(as.numeric(pbmc$seurat_clusters)))==1){
        print(-1)
        }
        else
            {summary(silhouette(as.numeric(pbmc$seurat_clusters), dd))$avg.width}
      }

rogue_function = function(.x){
        print(.x)
        pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = .x)
        se_cluster <- as.numeric(pbmc$seurat_clusters)
        tmp.res <- rogue(pbmc@assays$RNA@counts, labels = se_cluster, samples = rep("a", length(se_cluster)),platform = "UMI")
        mean(as.matrix(tmp.res)[1,])
      }
mc <- getOption("mc.cores", 50)
tibble(
  resolution = seq(0.6,1.2,by=0.1)
) %>%
      dplyr::mutate(
    se_cluster = mclapply(
      resolution,
      se_cluster_function,
      mc.cores = mc
    )
) %>%
      dplyr::mutate(
    silhouette = mclapply(
      resolution,
      silhouette_function,
      mc.cores = mc
    )
)  %>%
    dplyr::mutate(
    rogue = mclapply(
      resolution,
      rogue_function,
      mc.cores = mc
    ))-> res
res$rogue = as.numeric(res$rogue)
res$silhouette = as.numeric(res$silhouette)
res$se_cluster = as.numeric(res$se_cluster)
res[is.na(res)] = 0
res$rogue.scale <- (res$rogue-min(res$rogue))/(max(res$rogue)-min(res$rogue))
res$silhouette.scale <- (res$silhouette-min(res$silhouette))/(max(res$silhouette)-min(res$silhouette))
res$`rogue.scale`[is.nan(res$`rogue.scale`)] = 1
res$avg = 0
for(m in c(1:dim(res)[1])){
if(res$rogue.scale[m] != 0 & res$silhouette.scale[m] != 0)
{
res$avg[m] = (2*(res$rogue.scale[m] * res$silhouette.scale[m]))/(res$rogue.scale[m] + res$silhouette.scale[m])
print(m)
}}
#确定最佳的resolution
best_resolution = res$resolution[which(res$avg == max(res$avg))]
best_resolution = best_resolution[1]
print("best_resolution")
print(best_resolution)
#使用最佳的resolution进行聚类，重置metedata信息
plan("multiprocess", workers = works_num)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = best_resolution)

#去除双细胞
plan("multiprocess", workers = works_num)
sweep.res.list <- paramSweep_v3(pbmc, PCs = pc.num, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
plan("multiprocess", workers = works_num)
suppressMessages(bcmvn <- find.pK(sweep.stats))
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
## Homotypic Doublet Proportion Estimate
plan("multiprocess", workers = works_num)
homotypic.prop <- modelHomotypic(pbmc$seurat_clusters)
#双细胞比率
plan("multiprocess", workers = works_num)
DoubletRate = 0.05
##计算双细胞比例
plan("multiprocess", workers = works_num)
nExp_poi <- round(DoubletRate*ncol(pbmc)) 
##使用同源双细胞比例对计算的双细胞比例进行校正 
plan("multiprocess", workers = works_num)
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 
## Run DoubletFinder with varying classification stringencies
plan("multiprocess", workers = works_num)
pbmc <- doubletFinder_v3(pbmc, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
colnames(pbmc@meta.data)[dim(pbmc@meta.data)[2]] = "DoubletFinder"
pbmc <- subset(pbmc, subset = DoubletFinder  == "Singlet")

# 找出每个cluster的标记与所有剩余的细胞相比较，只报告阳性细胞
plan("multiprocess", workers = works_num)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top500 = pbmc.markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC)


RDS_name = paste("/home/ud202180737/stAtlas/2.Result/",strsplit(file_list[i],'/')[[1]][6],"/",strsplit(strsplit(file_list[i],'/')[[1]][7],"_")[[1]][1],".rds",sep="")
maker_name = paste("/home/ud202180737/stAtlas/2.Result/",strsplit(file_list[i],'/')[[1]][6],"/",strsplit(strsplit(file_list[i],'/')[[1]][7],"_")[[1]][1],".csv",sep="")
saveRDS(pbmc,file = RDS_name)
write_csv(top500,file=maker_name)
print("Done")
}
file_list