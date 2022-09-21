# ---------------- macrophage ----------------

mye <- readRDS('wclc/rds/myeloid.tumor.sub.rds')
macro <- subset(mye, cells = rownames(subset(mye@meta.data, clMidwayPr == 'Macro')))

macro <- NormalizeData(object = macro, normalization.method = "LogNormalize", scale.factor = 10000)
macro <- FindVariableFeatures(macro,selection.method="vst",nfeatures=2000)
macro <- ScaleData(macro, features = VariableFeatures(macro))
macro <- RunPCA(object = macro, features = VariableFeatures(macro))
# Run PCA and Determine Dimensions for 90% Variance
PCDeterminators <- function(object){
  stdev <- object@reductions$pca@stdev
  var <- stdev^2
  endVar <- 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    exp.var <- numerator/total
    if(endVar == 0){
      if(exp.var > 0.9){
        endVar <- endVar + 1
        PC.num <- i
      }
    }
  }
  sum(var[1:PC.num])/sum(var)
  return(PC.num)
}
PC.num <- PCDeterminators(macro)
# Find neighbors and clusters with harmony batch correction
macro <- FindNeighbors(macro, dims = 1:PC.num, reduction = "pca")
macro <- FindClusters(macro, resolution = 0.2)
macro <- RunUMAP(macro, dims = 1:PC.num, reduction = "pca")
macro <- RunTSNE(macro,dims = 1:PC.num)

saveRDS(macro, file = 'wclc/rds/macro.rds')

p <- DimPlot(macro,label = T,label.size = 5,split.by = 'TissueSiteSimple',ncol = 1) + 
  NoLegend() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'macro.dr.pdf',p,width = 4,height = 8)


rm(list = ls());gc()
macro <- readRDS('wclc/rds/macro.rds')

mypal <- RColorBrewer::brewer.pal(10,"Paired")
p <- DimPlot(macro,cols = mypal,label = T,label.size = 5) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/macro/seurat.cluster.pdf',p,width = 5,height = 5)

p <- ggplot(macro@meta.data, aes(x = TissueSiteSimple, fill = seurat_clusters)) +
  geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values = mypal) +
  mytheme +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(filename = 'wclc/macro/cell.prop2.pdf',p,width = 8,height = 8)

pathway <- xlsx::read.xlsx("wclc/xlsx/pathway.tmp.xlsx",sheetIndex = 2)
p <- ggplot(pathway, aes(x=Tissue,y=Description,color=Genes,size=-log10(pvalue))) +
  geom_point() + theme_bw() + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 14,colour = "black"),
        legend.position = "left",
        axis.title = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20),limits=c(1,4))
ggsave(filename = "wclc/macro/pathway.lcc.rcc.dotplot.pdf",p,width = 8,height = 6)


###### granulocyte migration #########

lcc.mac <- subset(macro, cells = rownames(subset(macro@meta.data, TissueSiteSimple == 'left')))
rcc.mac <- subset(macro, cells = rownames(subset(macro@meta.data, TissueSiteSimple == 'right')))

# PPBP/CXCL7 CCL23  CCL2 CCL8
genes <- c("PPBP","CCL2","CCL8","CCL23","S100A12")
mypal <- c("grey88","DarkCyan")

for(i in 1:length(genes)){
  p <- FeaturePlot(lcc.macro,features = genes[i],cols = mypal) +
    ggtitle('') +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none',
          panel.border = element_rect(size = 1,colour = "black"))
  ggsave(filename = paste0("wclc/macro/gra.migration.features/LCC.",genes[i],".pdf"),p,width = 4,height = 4)
  p <- FeaturePlot(rcc.macro,features = genes[i],cols = mypal) +
    ggtitle('') +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none',
          panel.border = element_rect(size = 1,colour = "black"))
  ggsave(filename = paste0("wclc/macro/gra.migration.features/RCC.",genes[i],".pdf"),p,width = 4,height = 4)
}


p <- VlnPlot(macro,features = "PPBP",group.by = "TissueSiteSimple",y.max = 7) + NoLegend() + 
  theme(panel.border = element_rect(colour = "black",size = 1.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0.5)) +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/macro/gra.migration.features/PPBP.pdf",width = 4,height = 5)

p <- VlnPlot(macro,features = "CCL23",group.by = "TissueSiteSimple",
             y.max = 8,pt.size = 0,cols = c("#4E9BD2","#DB4740")) + 
  NoLegend() + 
  theme(panel.border = element_rect(colour = "black",size = 1.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0.5)) +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/macro/gra.migration.features/CCL23.pdf",width = 4,height = 5)

p <- VlnPlot(macro,features = "CCL2",group.by = "TissueSiteSimple",y.max = 6.5) + 
  NoLegend() + 
  theme(panel.border = element_rect(colour = "black",size = 1.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0.5)) +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/macro/gra.migration.features/CCL2.pdf",width = 4,height = 5)

p <- VlnPlot(macro,features = "CCL8",group.by = "TissueSiteSimple",y.max = 6) + 
  NoLegend() + 
  theme(panel.border = element_rect(colour = "black",size = 1.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0.5)) +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/macro/gra.migration.features/CCL8.pdf",width = 4,height = 5)

p <- VlnPlot(macro,features = "S100A12",group.by = "TissueSiteSimple",y.max = 6,
             pt.size = 0, cols = c("#4E9BD2","#DB4740")) + 
  NoLegend() + 
  theme(panel.border = element_rect(colour = "black",size = 1.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0.5)) +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/macro/gra.migration.features/S100A12.pdf",width = 4,height = 5)

Idents(macro) <- 'TissueSiteSimple'
rlt <- AverageExpression(macro,features = genes)
expr <- rlt$RNA
expr <- log2(expr)
pheatmap(expr,border_color = 'black',filename = 'wclc/macro/gra.migra.pdf',
         cellwidth = 40,cellheight = 40)

# GENE PROPORTION
# PPBP
tag1 <- lcc.mac@assays$RNA@counts['PPBP',]
tag2 <- rcc.mac@assays$RNA@counts['PPBP',]
lcc.mac$GroupI <- ifelse(tag1 > 0, "PPBP+",'PPBP-')
rcc.mac$GroupI <- ifelse(tag2 > 0, "PPBP+",'PPBP-')
table(lcc.mac$GroupI)
sum(unname(table(lcc.mac$GroupI)))
unname(table(lcc.mac$GroupI))[2]/sum(unname(table(lcc.mac$GroupI)))
table(rcc.mac$GroupI)
sum(unname(table(rcc.mac$GroupI)))
unname(table(rcc.mac$GroupI))[2]/sum(unname(table(rcc.mac$GroupI)))
# CCL23
tag1 <- lcc.mac@assays$RNA@counts['CCL23',]
tag2 <- rcc.mac@assays$RNA@counts['CCL23',]
lcc.mac$GroupI <- ifelse(tag1 > 0, "CCL23+",'CCL23-')
rcc.mac$GroupI <- ifelse(tag2 > 0, "CCL23+",'CCL23-')
table(lcc.mac$GroupI)
sum(unname(table(lcc.mac$GroupI)))
unname(table(lcc.mac$GroupI))[2]/sum(unname(table(lcc.mac$GroupI)))
table(rcc.mac$GroupI)
sum(unname(table(rcc.mac$GroupI)))
unname(table(rcc.mac$GroupI))[2]/sum(unname(table(rcc.mac$GroupI)))
# CCL2
tag1 <- lcc.mac@assays$RNA@counts['CCL2',]
tag2 <- rcc.mac@assays$RNA@counts['CCL2',]
lcc.mac$GroupI <- ifelse(tag1 > 0, "CCL2+",'CCL2-')
rcc.mac$GroupI <- ifelse(tag2 > 0, "CCL2+",'CCL2-')
table(lcc.mac$GroupI)
sum(unname(table(lcc.mac$GroupI)))
unname(table(lcc.mac$GroupI))[2]/sum(unname(table(lcc.mac$GroupI)))
table(rcc.mac$GroupI)
sum(unname(table(rcc.mac$GroupI)))
unname(table(rcc.mac$GroupI))[2]/sum(unname(table(rcc.mac$GroupI)))
# CCL8
tag1 <- lcc.mac@assays$RNA@counts['CCL8',]
tag2 <- rcc.mac@assays$RNA@counts['CCL8',]
lcc.mac$GroupI <- ifelse(tag1 > 0, "CCL8+",'CCL8-')
rcc.mac$GroupI <- ifelse(tag2 > 0, "CCL8+",'CCL8-')
table(lcc.mac$GroupI)
sum(unname(table(lcc.mac$GroupI)))
unname(table(lcc.mac$GroupI))[2]/sum(unname(table(lcc.mac$GroupI)))
table(rcc.mac$GroupI)
sum(unname(table(rcc.mac$GroupI)))
unname(table(rcc.mac$GroupI))[2]/sum(unname(table(rcc.mac$GroupI)))

# all.markers <- FindAllMarkers(macro, only.pos = T)
# 
# load('data/signature_collection.rda')
# setdiff(signature_collection$M1_Azizi_et_al,rownames(macro))
# # [1] "iNOS"          "IL12"          "CD64 (FcγR1A)" "CD64 (FcγR1B)" "CD64 (FcγR1C)"
# # [6] "CD80 (B7-1)"   "CXCR10"        "IL23"          "CXCL11"        "CD86 (B7-2)"  
# # [11] "TNFa"          "MHCII"         "CD40"          "KYNU"
# diff.gene <- c('FCGR1A','FCGR1B','CD80','CXCL10','IL23A','CD86','TNF')
# m1.sig <- c(diff.gene,intersect(signature_collection$M1_Azizi_et_al,rownames(macro)))
# setdiff(m1.sig,rownames(macro))
# 
# setdiff(signature_collection$M2_Azizi_et_al,rownames(macro))
# # [1] "ARG1/2 (Arginase)"  "CD32"               "CD23 (FCER2)"       "PD-L2 (PDCD1LG2)"  
# # [5] "PD-L1 (CD274)"      "CSF1R"              "CD206 (MRC1)"       "Il1RA (IL1RN)"     
# # [9] "Il1R2"              "VEGFC"              "VEGFD"              "EGF"               
# # [13] "Cathepsin A (CTSA)" "CSTC"               "MMP19"              "WNT7b"             
# # [17] "FASL"               "TNFSF12"            "CD276 (B7-H3)"      "VTCN1 (BH-H4)"     
# # [21] "MSR1 (CD204)"
# diff.gene <- c('ARG1','ARG2','FCGR2A','FCER2','PDCD1LG2','CD274','IL1RN','IL1R2',
#                'CTSA','CTSC','FASLG','MSR1')
# m2.sig <- c(diff.gene,intersect(signature_collection$M2_Azizi_et_al,rownames(macro)))
# setdiff(m2.sig,rownames(macro))
# gene.match <- c(m1.sig,m2.sig)
# tmp <- macro[gene.match,]
# tmp <- ScaleData(tmp,features = rownames(tmp))
# mat <- GetAssayData(tmp, assay = "RNA", slot = "scale.data")
# 
# mat[mat > 5] = 5
# mat[mat < -2] = -2
# cls.info <- sort(tmp$seurat_clusters)
# 
# mat <- as.matrix(mat[gene.match, names(cls.info)])
# gene <- unique(gene.match)
# gene.pos <- which(rownames(mat) %in% gene)
# row.anno <- ComplexHeatmap::rowAnnotation(gene = ComplexHeatmap::anno_mark(at = gene.pos, labels = rownames(mat)[which(rownames(mat) %in% gene)]))
# cols <- RColorBrewer::brewer.pal(10,"Paired")
# names(cols) <- levels(cls.info)
# top.anno <- ComplexHeatmap::HeatmapAnnotation(
#   cluster = ComplexHeatmap::anno_block(gp = grid::gpar(fill=cols), labels = levels(cls.info), labels_gp = grid::gpar(cex=0.5,col='white'))
# )
# col.fun <- circlize::colorRamp2(seq(min(mat),max(mat),length=3), c('#377EB8','white','#E41A1C'))
# pdf("wclc/diff.cls.heatmap.pdf",width = 10,height = 7)
# ComplexHeatmap::Heatmap(matrix = mat, cluster_rows = F, cluster_columns = F,
#                         show_column_names = F, show_row_names = F,
#                         column_split = cls.info, right_annotation = row.anno,
#                         column_title = NULL, top_annotation = top.anno,
#                         heatmap_legend_param = list(title = 'Expression', title_position = 'leftcenter-rot'),col = col.fun)
# dev.off()


rm(list = ls());gc()
macro <- readRDS('wclc/rds/macro.rds')
mypal <- RColorBrewer::brewer.pal(10,"Paired")
p1 <- VlnPlot(macro,features = 'C1QC',cols = mypal,pt.size = 0) + NoLegend() +
  coord_flip()  + xlab('Clusters') +
  theme(axis.text.x = element_text(angle = 0,hjust = .5),
        panel.border = element_rect(size = 1,colour = 'black'))
p2 <- VlnPlot(macro,features = 'C1QA',pt.size = 0,cols = mypal) + NoLegend() +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = .5),
        panel.border = element_rect(size = 1,colour = 'black'))
p3 <- VlnPlot(macro,features = 'APOE',pt.size = 0,cols = mypal) + NoLegend() +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = .5),
        panel.border = element_rect(size = 1,colour = 'black'))
pc <- p1 | p2 | p3
ggsave("wclc/macro/C1Q.pdf",pc,width = 5,height = 6)

macro$celltype <- "C1QC+TAM"
saveRDS(macro,file = 'wclc/rds/macro.rds')

p <- DimPlot(macro,group.by = 'celltype',pt.size = 0.8,cols = '#FF7F0E') +
  NoLegend() + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank())
ggsave("wclc/macro/celltype.pdf",p,width = 5,height = 5)


macro.query <- readRDS('wclc/rds/macro.rds')
macro.ref <- readRDS('validation/macro.rds')

anchors <- FindTransferAnchors(reference = macro.ref, query = macro.query,dims = 1:30,k.filter = 500)
predictions <- TransferData(anchorset = anchors, refdata = macro.ref$Cell_subtype, dims = 1:30)
macro.query <- AddMetaData(macro.query, metadata = predictions)
macro.query$prediction.match <- macro.query$predicted.id == macro.ref$Cell_subtype
saveRDS(macro.query, file = 'wclc/rds/macro.celltype.predicted.rds')

macro <- readRDS('wclc/rds/macro.celltype.predicted.rds')

p <- ggplot(macro@meta.data, aes(x = TissueSiteSimple, fill = predicted.id)) +
  geom_bar(position = "fill") + theme_bw() + 
  labs(y = "Cell Fraction") + coord_flip()