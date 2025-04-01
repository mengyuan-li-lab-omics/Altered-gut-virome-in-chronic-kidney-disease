
library(Seurat)                                #加载Seurat
library(dplyr)                                 #加载dplyr包
#library(SingleR)                               #加载SingleR包，细胞注释包
#library(celldex)                               #加载celldex包

#多个样本合并
# 导入Seurat包
library(Seurat)
# 获取数据文件夹下的所有样本文件列表
samples <- list.files("H:/data/")

# 创建一个空的列表来存储Seurat对象
seurat_list <- list()

# 读取每个样本的10x数据并创建Seurat对象
for (sample in samples) {
# 拼接文件路径
  data.path <- paste0("H:/data/", sample)

# 读取10x数据，data.dir参数指定存放文件的路径
  seurat_data <- Read10X(data.dir = data.path)

# 创建Seurat对象，并指定项目名称为样本文件名
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   project = sample,
                                   min.features = 200,
                                   min.cells = 3)

# 将Seurat对象添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}

# 打印所有的Seurat对象列表
seurat_list

# 合并Seurat对象，将所有Seurat对象合并到一个对象中
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)
# 打印合并后的Seurat对象
print(seurat_combined)

pbmc=seurat_combined

head(pbmc@meta.data) 



#计算线粒体基因表达比例
# 注意MT大小写，小写表示另一个物种
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "mt")

#绘制Feature,Count,mt比例图

pdf("H:/绘制Feature,Count,mt比例图.pdf")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)   
dev.off()
#细胞过滤
pbmc <- subset(pbmc,subset = nFeature_RNA <=2000 & nCount_RNA <= 5000 & percent.mt <= 5) 

#绘制过滤后的图片
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = 'orig.ident', ncol = 3)

#数据标准化（Normalizing the data）
# 去除测序深度带来的影响
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#识别高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# plot1 + plot2


#数据标准化 Scale
# PCA降维要求数据为正态分布，即平均值为0，方差为1。

#对所有基因进行标准化
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)

#只对上述2000个高变基因进行标准化，通常仅对高变基因进行标准化和降维
pbmc <- ScaleData(pbmc)                       

#PCA降维
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc) )

DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc)

str(pbmc)

# 细胞聚类
pbmc <- FindNeighbors(pbmc, dims = 1:10)        
pbmc <- FindClusters(pbmc, resolution = 0.5) # 分辨率与最终得到的cluster数量成正比


save(pbmc,file = 'H:/pbmc.Rdata')
load(file = "H:/pbmc.RData")

#UMAP降维
pdf("H:/UMAP.pdf",5,5)
pbmc1 <- RunUMAP(pbmc, dims = 1:10,label=1)
DimPlot(pbmc1, reduction = "umap")
dev.off()
dev.off()

#tSNE降维
pdf("H:/tSNE.pdf",8,8)
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")
dev.off()

save(pbmc1,file = 'H:/pbmc.Rdata')
load(file = "H:/pbmc.RData")

#差异分析
# find all markers of cluster 0（寻找cluster0的基因）
cluster0.markers <- FindMarkers(JoinLayers(pbmc), ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)
#可以根据需要，将每个聚类的差异基因都找一遍，方便后续分析
#head(cluster0.markers, n = 5)，可以将n = 5进行调整，比如换成10，即可显示10个差异基因
###参数解释
#ident.1 设置待分析的细胞类别
#min.pct 在两组细胞中的任何一组中检测到的最小百分
#thresh.test 在两组细胞间以一定数量的差异表达（平均）
#max.cells.per.ident 通过降低每个类的采样值，提高计算速度
 
#寻找cluster1的差异基因
cluster1.markers <- FindMarkers(JoinLayers(pbmc), ident.1 = 1, min.pct = 0.25)
#找到cluster5和cluster0、3之间的markgene
cluster5.markers <- FindMarkers(JoinLayers(pbmc), ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
 
#找到每个cluster相比于其余cluster的markgene，只报道阳性的markgene
pbmc.markers <- FindAllMarkers(JoinLayers(pbmc1), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
## 所有基因先分组，再根据avg_log2FC进行排序
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

#对marker_gene进行筛选p_val<0.05
pbmc.markers %>% subset(p_val<0.05)
## 所有基因先分组，再根据avg_log2FC进行排序，选出每组前十个
list_marker <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
print(list_marker,n = Inf)#打印list_marker所有值
#保存差异分析结果到csv
df_marker=data.frame(p_val = list_marker$p_val,
                     avg_log2FC = list_marker$avg_log2FC,
                     pct.1 = list_marker$pct.1,
                     pct.2 = list_marker$pct.2,
                     p_val_adj = list_marker$p_val_adj,
                     cluster = list_marker$cluster,
                     gene = list_marker$gene)
write.csv(df_marker,"H:/marker-all.csv")
 

pdf("H:/差异基因.pdf",8,8)
DoHeatmap(pbmc1, features = top10$gene) + NoLegend()
dev.off()

# scedata=pbmc1
# table(scedata$orig.ident)#查看各组细胞数
# prop.table(table(Idents(scedata)))
#pbmc1=pbmc
#分样本
#cell maker注释后
library(tidyverse)
library(monocle3)
library(Seurat)
library(pheatmap)
# #添加细胞注释信息,通过CellMarker注释每一个cluster代表的细胞类群


new.cluster.ids <- c("CD8+ T cell","Undefined","Dendritic cell", "Granulocyte-monocyte progenitor","Th cell", 
"Epithelial cells","Epithelial cells","Monocyte", "Neutrophil",
"Undefined" ,"Neutrophil","Macrophage", "Epithelial cells", "B cell")
names(new.cluster.ids) <- levels(pbmc1)
pbmc1 <- RenameIdents(pbmc1, new.cluster.ids)

colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
pdf("H:/注释-总.pdf",7,5)
DimPlot(pbmc1, reduction = 'umap', 
              cols = colour , #设置颜色
              repel = T,  #label不重叠
              label.size = 4, #调整label的大小
              label = TRUE,  #是否展示label
              pt.size = 0.5) 
			  dev.off()


pbmc1$Anno = "NA"
celltype =c("CD8+ T cell","Undefined","Dendritic cell", "Granulocyte-monocyte progenitor","Th cell", 
"Proximal tubular cell","Epithelial cells","Granulocyte-monocyte progenitor",  "Neutrophil",
"Undefined" ,"Neutrophil","Macrophage", "Epithelial cells", "B cell")

#for循环添加
sub_length = length(unique(pbmc1$seurat_clusters)) - 1
for (i in 0:sub_length){ 
  pbmc1$Anno[pbmc1$seurat_clusters==i] = celltype[i+1]
  }

# 绘制注释umap图
DimPlot(pbmc1, reduction = 'umap', group.by='Anno',
label = TRUE, pt.size = 0.5) + NoLegend()

pdf("H:/注释-总-分组排列1.pdf",12,5)
DimPlot(pbmc1,reduction = "umap",label=T,cols = colour , split.by ='orig.ident')
dev.off() 