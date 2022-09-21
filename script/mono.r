# ---------------- monocytes ----------------

rm(list = ls());gc()
mye.T <- readRDS('wclc/rds/myeloid.tumor.sub.rds')

# mono
mono <- subset(mye.T, cells = rownames(subset(mye.T@meta.data, clMidwayPr == 'Mono')))
mono <- NormalizeData(object = mono, normalization.method = "LogNormalize", scale.factor = 10000)
mono <- FindVariableFeatures(mono,selection.method="vst",nfeatures=2000)
mono <- ScaleData(mono, features = VariableFeatures(mono))
mono <- RunPCA(object = mono, features = VariableFeatures(mono))
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
PC.num <- PCDeterminators(mono)
# Find neighbors and clusters with harmony batch correction
mono <- FindNeighbors(mono, dims = 1:PC.num, reduction = "pca")
mono <- FindClusters(mono, resolution = 0.2)
mono <- RunUMAP(mono, dims = 1:PC.num, reduction = "pca")
mono <- RunTSNE(mono,dims = 1:PC.num)
p <- DimPlot(mono,label = T,label.size = 5) + NoLegend() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/mono/mono.dr.pdf',p,width = 4,height = 4)

p <- DimPlot(mono,label = T,label.size = 5,split.by = 'TissueSiteSimple',ncol = 1) + NoLegend() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/mono/mono.splitBy.rcc.lcc.dr.pdf',p,width = 4,height = 8)

saveRDS(mono, file = 'wclc/rds/mono.rds')

mypal <- RColorBrewer::brewer.pal(7,"Paired")

p <- DimPlot(mono,cols = mypal,label = T,label.size = 5) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/mono/seurat.cluster.pdf',p,width = 5,height = 5)

p <- ggplot(mono@meta.data, aes(x = TissueSiteSimple, fill = seurat_clusters)) +
  geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values = mypal) +
  mytheme +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(filename = 'wclc/mono/cell.prop.pdf',p,width = 5,height = 10)
p <- ggplot(mono@meta.data, aes(x = TissueSiteSimple, fill = seurat_clusters)) +
  geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values = mypal) +
  mytheme +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(filename = 'wclc/mono/cell.prop2.pdf',p,width = 8,height = 8)

mypal <- RColorBrewer::brewer.pal(7,"Paired")
p <- DimPlot(mono,cols = mypal,label = T,label.size = 5,reduction = "tsne") +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_blank())
ggsave(filename = 'wclc/mono/seurat.cluster.tsne.pdf',p,width = 5,height = 5)

pathway <- xlsx::read.xlsx("wclc/xlsx/pathway.tmp.xlsx",sheetIndex = 1)
p <- ggplot(pathway, aes(x=Tissue,y=Description,color=Genes,size=-log10(pvalue))) +
  geom_point() + theme_bw() + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 14,colour = "black"),
        legend.position = "left",
        axis.title = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20),limits=c(1,4))
ggsave(filename = "wclc/mono/pathway.lcc.rcc.dotplot.pdf",p,width = 8,height = 6)


library(Seurat)
library(ggplot2)

mono <- readRDS('wclc/rds/mono.rds')

# Mono_CD14: FCN1, S100A8, S100A9
# Mono_CD16: FCGR3A, LST1, LILRB2
mypal <- RColorBrewer::brewer.pal(7,"Paired")
p1 <- VlnPlot(mono,features = 'CD14',pt.size = 0,cols = mypal) + NoLegend() +
  coord_flip() + xlab('Clusters') +
  theme(axis.text.x = element_text(angle = 0,hjust = .5),
        panel.border = element_rect(size = 1,colour = 'black'))
p2 <- VlnPlot(mono,features = 'FCGR3B',pt.size = 0,cols = mypal) + NoLegend() +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = .5),
        panel.border = element_rect(size = 1,colour = 'black'))
pc <- p1|p2
ggsave("wclc/mono/CD14.CD16.pdf",pc,width = 4,height = 6)

# CD14 <- mono@assays$RNA@counts['CD14',]
# CD16 <- mono@assays$RNA@counts['FCGR3B',]
# mono$CD14 <- ifelse(CD14 > 0, "CD14+",'CD14-')
# mono$CD16 <- ifelse(CD16 > 0, 'CD16+','CD14+')
# mono$celltype <- ifelse(CD16>0 & CD14>0, "CD14+CD16+",ifelse(CD14>0,"CD14+CD16-",ifelse(CD16>0,"CD14-CD16+",'CD14-CD16-')))

# markers <- xlsx::read.xlsx("data/NIHMS910854-supplement-Supplementary_Tables_1-16.xlsx",sheetIndex = 5,startRow = 3)
# markers$Cluster.ID <- gsub("Mono1","Mono1(CD14+ monocyte)",markers$Cluster.ID)
# markers$Cluster.ID <- gsub("Mono2","Mono1(CD16+ monocyte)",markers$Cluster.ID)
# markers <- markers[,c("Gene.ID","Cluster.ID")]
# colnames(markers) <- c("markers","celltype")

old.ident <- c(0:6)
new.ident <- c('CD14+ monocytes',"CD14+CD16+ monocytes", "CD16+ monocytes", "CD14+CD16+ monocytes",
               "CD14+ monocytes",'CD14+ monocytes','CD14+ monocytes')
mono$celltype <- plyr::mapvalues(mono$seurat_clusters, from = old.ident, to = new.ident)
saveRDS(mono,file = 'wclc/rds/mono.rds')

# c('#1F77B4', '#39D486', '#43315C', '#00C8FF', '#DE4343', '#FF7F0E', '#E377C2')
mypal <- c('#00C8FF','#DE4343','#FF7F0E')
p <- DimPlot(mono,group.by = 'celltype',pt.size = 0.8,reduction = 'tsne',cols = mypal) +
  NoLegend() + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank())
ggsave("wclc/mono/celltype.pdf",p,width = 5,height = 5)

cell.prop <- as.matrix(table(mono$TissueSiteSimple,mono$celltype))
cell.prop <- reshape2::melt(cell.prop)

lcc.cell.prop <- cell.prop[cell.prop$Var1=="left",]
rcc.cell.prop <- cell.prop[cell.prop$Var1=="right",]
label <- paste0(round(lcc.cell.prop$value/sum(lcc.cell.prop$value)*100,1),'%')
label <- paste0(round(rcc.cell.prop$value/sum(rcc.cell.prop$value)*100,1),'%')

library(ggpubr)
p1 <- ggdonutchart(lcc.cell.prop,'value',label=label,fill='Var2',
             palette = c('#FF7F0E','#00C8FF','#DE4343'),lab.pos = 'in') +
  theme(legend.position = 'none')
p2 <- ggdonutchart(rcc.cell.prop,'value',label=label,fill='Var2',
                   palette = c('#FF7F0E','#00C8FF','#DE4343'),lab.pos = 'in') +
  theme(legend.position = 'none')
pc <- p1 | p2
ggsave('wclc/mono/cell.prop.circle.pdf',pc,width = 6,height = 4)

p <- ggplot(mono@meta.data, aes(x = TissueSiteSimple, fill = celltype)) +
  geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values = mypal) +
  labs(y = "Cell Fraction") +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14,colour = 'black'),
        axis.title.y = element_text(size = 16,colour = 'black'),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/mono/cell.prop.2.pdf',p,width = 3,height = 5)

library(ROGUE)
mono <- readRDS('wclc/rds/mono.rds')
mono$celltype <- as.character(mono$celltype)
old.ident <- names(table(mono$celltype))
new.ident <- c('CD14+CD16- monocytes','CD14+CD16+ monocytes','CD14-CD16+ monocytes')
mono$celltype <- plyr::mapvalues(mono$celltype,from = old.ident,to = new.ident)
saveRDS(mono,file = 'wclc/rds/mono.rds')

Idents(mono) <- 'celltype'
# all.markers <- FindAllMarkers(mono,only.pos = T)
# CD14.markers <- all.markers[all.markers$cluster=="CD14+CD16- monocytes",]
# CD16.markers <- all.markers[all.markers$cluster=="CD14-CD16+ monocytes",]
# CD14.CD16.markers <- all.markers[all.markers$cluster=="CD14+CD16+ monocytes",]
# CD14.markers <- CD14.markers[order(CD14.markers$avg_logFC,decreasing = T),]

CD14.markers <- FindMarkers(mono,group.by = 'TissueSiteSimple',
                            ident.1 = 'right',ident.2 = 'left',
                            subset.ident = 'CD14+CD16- monocytes')
CD14.markers <- CD14.markers[!grepl("^AC[0-9]",rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl("^RP",rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl("^MT",rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^RN[0-9]',rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^RNU[0-9]',rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^RNA[0-9]',rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^CTD',rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^CTC',rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^Y-RNA',rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^GS1',rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^AL928768',rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^XXbac',rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^AP00',rownames(CD14.markers)),]
CD14.markers <- CD14.markers[!grepl('^CTA-',rownames(CD14.markers)),]

CD14.markers <- CD14.markers[order(CD14.markers$avg_logFC,decreasing = T),]
CD14.markers$celltype <- "CD14+CD16- monocytes"
CD14.markers$group <- ifelse(CD14.markers$avg_logFC>0,"RCC","LCC")
CD14.markers$gene <- rownames(CD14.markers)

CD16.markers <- FindMarkers(mono,group.by = 'TissueSiteSimple',
                            ident.1 = 'right',ident.2 = 'left',
                            subset.ident = 'CD14-CD16+ monocytes')
CD16.markers <- CD16.markers[!grepl("^AC[0-9]",rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl("^RP",rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl("^MT",rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^RN[0-9]',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^RNU[0-9]',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^RNA[0-9]',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^CTD',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^CTC',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^Y-RNA',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^CMB9',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^CTA',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^CTB',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^AL[0-9]',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^AP[0-9]',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^XX',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^CH17',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^BX47',rownames(CD16.markers)),]
CD16.markers <- CD16.markers[!grepl('^KB-',rownames(CD16.markers)),]

CD16.markers <- CD16.markers[order(CD16.markers$avg_logFC,decreasing = T),]
CD16.markers$celltype <- "CD14-CD16+ monocytes"
CD16.markers$group <- ifelse(CD16.markers$avg_logFC>0,"RCC","LCC")
CD16.markers$gene <- rownames(CD16.markers)

CD14.CD16.markers <- FindMarkers(mono,group.by = 'TissueSiteSimple',
                                 ident.1 = 'right',ident.2 = 'left',
                                 subset.ident = 'CD14+CD16+ monocytes')
CD14.CD16.markers <- CD14.CD16.markers[!grepl("^AC[0-9]",rownames(CD14.CD16.markers)),]
CD14.CD16.markers <- CD14.CD16.markers[!grepl("^RP",rownames(CD14.CD16.markers)),]
CD14.CD16.markers <- CD14.CD16.markers[!grepl("^MT",rownames(CD14.CD16.markers)),]
CD14.CD16.markers <- CD14.CD16.markers[!grepl('^RN[0-9]',rownames(CD14.CD16.markers)),]
CD14.CD16.markers <- CD14.CD16.markers[!grepl('^RNU[0-9]',rownames(CD14.CD16.markers)),]
CD14.CD16.markers <- CD14.CD16.markers[!grepl('^RNA[0-9]',rownames(CD14.CD16.markers)),]
CD14.CD16.markers <- CD14.CD16.markers[!grepl('^CTD',rownames(CD14.CD16.markers)),]
CD14.CD16.markers <- CD14.CD16.markers[!grepl('^CTC',rownames(CD14.CD16.markers)),]
CD14.CD16.markers <- CD14.CD16.markers[!grepl('^Y-RNA',rownames(CD14.CD16.markers)),]
CD14.CD16.markers <- CD14.CD16.markers[!grepl('^GS1',rownames(CD14.CD16.markers)),]

CD14.CD16.markers <- CD14.CD16.markers[order(CD14.CD16.markers$avg_logFC,decreasing = T),]
CD14.CD16.markers$celltype <- "CD14+CD16+ monocytes"
CD14.CD16.markers$group <- ifelse(CD14.CD16.markers$avg_logFC>0,"RCC","LCC")
CD14.CD16.markers$gene <- rownames(CD14.CD16.markers)

all.markers <- rbind(CD14.markers,CD14.CD16.markers,CD16.markers)


library(dplyr)
library(ggrepel)
library(ggsci)
top5pos <- all.markers %>% group_by(celltype) %>% top_n(n = 5, wt = avg_logFC)
top5neg <- all.markers %>% group_by(celltype) %>% top_n(n = -5, wt = avg_logFC)
top10 <- rbind(top5pos, top5neg)

all.markers$celltype <- factor(all.markers$celltype,
                               levels = c('CD14+CD16- monocytes',
                                          'CD14+CD16+ monocytes',
                                          'CD14-CD16+ monocytes'))
openxlsx::write.xlsx(all.markers,file = 'wclc/mono/markers.xlsx')
p <- ggplot(all.markers, aes(x = pct.2 - pct.1, y = avg_logFC)) +
  geom_point(color='grey80') +
  geom_hline(yintercept = c(-0.25,0.25),lty = 'dashed', size = 1, color = 'grey50') +
  geom_text_repel(data = top5pos,aes(x = pct.2-pct.1, y=avg_logFC, label=gene,color=celltype),
                  show.legend = F, direction = 'y', hjust = 1, nudge_y = 0.25, force = 5,
                  nudge_x = 0.8 - (top5pos$pct.2 - top5pos$pct.1)) +
  geom_text_repel(data = top5neg, aes(x = pct.2-pct.1, y=avg_logFC, label=gene,color=celltype),
                  show.legend = F, direction = 'y', hjust = 0, nudge_y = 0, force = 2.5, 
                  nudge_x = -0.8 - (top5neg$pct.2 - top5neg$pct.1)) +
  geom_point(data = top10, show.legend = F,aes(x=pct.2-pct.1, y=avg_logFC, color=celltype)) +
  scale_color_manual(values = c("#FF7F0E","#00C8FF","#DE4343")) +
  scale_y_continuous(limits = c(-2,2),breaks = seq(-2,2,1)) +
  scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_rect(colour = NA, fill = 'grey90')) +
  xlab(expression(Delta~'Percentage Difference')) +
  ylab('Log2-Fold Change') +
  facet_wrap(~celltype, nrow = 1, scales = 'fixed')
ggsave('wclc/mono/DEgenes.pdf',p,width = 7,height = 4)


### Enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
EnrichmentAnalysis <- function(DEgenes, enrich.method = c("GO","KEGG")){
  enrich.method <- match.arg(enrich.method)
  if(enrich.method == "GO") {
    ego_BP <- enrichGO(gene = DEgenes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_BP@result$Description <- substring(ego_BP@result$Description, 1, 70)
    return(ego_BP)
  } else if(enrich.method == "KEGG") {
    genelist <- bitr(DEgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    genelist <- pull(genelist, ENTREZID)
    ekegg <- enrichKEGG(gene = genelist, organism = "hsa")
    return(ekegg)
  }
}

CD14.markers <- all.markers[all.markers$celltype=="CD14+CD16- monocytes",]
RCC.CD14.gene <- CD14.markers[CD14.markers$group=="RCC",]$gene
LCC.CD14.gene <- CD14.markers[CD14.markers$group=="LCC",]$gene
CD14.RCC.GO <- EnrichmentAnalysis(DEgenes = RCC.CD14.gene, enrich.method = 'GO')
CD14.RCC.GO <- CD14.RCC.GO@result
CD14.LCC.GO <- EnrichmentAnalysis(DEgenes = LCC.CD14.gene, enrich.method = 'GO')
CD14.LCC.GO <- CD14.LCC.GO@result
CD14.RCC.KEGG <- EnrichmentAnalysis(DEgenes = RCC.CD14.gene, enrich.method = 'KEGG')
CD14.RCC.KEGG <- CD14.RCC.KEGG@result
CD14.LCC.KEGG <- EnrichmentAnalysis(DEgenes = LCC.CD14.gene, enrich.method = 'KEGG')
CD14.LCC.KEGG <- CD14.LCC.KEGG@result

openxlsx::write.xlsx(CD14.LCC.GO,file = 'wclc/mono/EnrichmentResults.xlsx')
xlsx::write.xlsx(CD14.RCC.GO,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD14.RCC.GO',row.names = F)
xlsx::write.xlsx(CD14.LCC.KEGG,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD14.LCC.KEGG',row.names = F)
xlsx::write.xlsx(CD14.RCC.KEGG,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD14.RCC.KEGG',row.names = F)

CD16.markers <- all.markers[all.markers$celltype=="CD14-CD16+ monocytes",]
RCC.CD16.gene <- CD16.markers[CD16.markers$group=="RCC",]$gene
LCC.CD16.gene <- CD16.markers[CD16.markers$group=="LCC",]$gene
CD16.RCC.GO <- EnrichmentAnalysis(DEgenes = RCC.CD16.gene, enrich.method = 'GO')
CD16.RCC.GO <- CD16.RCC.GO@result
CD16.LCC.GO <- EnrichmentAnalysis(DEgenes = LCC.CD16.gene, enrich.method = 'GO')
CD16.LCC.GO <- CD16.LCC.GO@result
CD16.RCC.KEGG <- EnrichmentAnalysis(DEgenes = RCC.CD16.gene, enrich.method = 'KEGG')
CD16.RCC.KEGG <- CD16.RCC.KEGG@result
CD16.LCC.KEGG <- EnrichmentAnalysis(DEgenes = LCC.CD16.gene, enrich.method = 'KEGG')
CD16.LCC.KEGG <- CD16.LCC.KEGG@result

xlsx::write.xlsx(CD16.LCC.GO,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD16.LCC.GO',row.names = F)
xlsx::write.xlsx(CD16.RCC.GO,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD16.RCC.GO',row.names = F)
xlsx::write.xlsx(CD16.LCC.KEGG,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD16.LCC.KEGG',row.names = F)
xlsx::write.xlsx(CD16.RCC.KEGG,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD16.RCC.KEGG',row.names = F)

CD14.CD16.markers <- all.markers[all.markers$celltype=="CD14+CD16+ monocytes",]
RCC.CD14.CD16.gene <- CD14.CD16.markers[CD14.CD16.markers$group=="RCC",]$gene
LCC.CD14.CD16.gene <- CD14.CD16.markers[CD14.CD16.markers$group=="LCC",]$gene
CD14.CD16.RCC.GO <- EnrichmentAnalysis(DEgenes = RCC.CD14.CD16.gene, enrich.method = 'GO')
CD14.CD16.RCC.GO <- CD14.CD16.RCC.GO@result
CD14.CD16.LCC.GO <- EnrichmentAnalysis(DEgenes = LCC.CD14.CD16.gene, enrich.method = 'GO')
CD14.CD16.LCC.GO <- CD14.CD16.LCC.GO@result
CD14.CD16.RCC.KEGG <- EnrichmentAnalysis(DEgenes = RCC.CD14.CD16.gene, enrich.method = 'KEGG')
CD14.CD16.RCC.KEGG <- CD14.CD16.RCC.KEGG@result
CD14.CD16.LCC.KEGG <- EnrichmentAnalysis(DEgenes = LCC.CD14.CD16.gene, enrich.method = 'KEGG')
CD14.CD16.LCC.KEGG <- CD14.CD16.LCC.KEGG@result

xlsx::write.xlsx(CD14.CD16.LCC.GO,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD14.CD16.LCC.GO',row.names = F)
xlsx::write.xlsx(CD14.CD16.RCC.GO,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD14.CD16.RCC.GO',row.names = F)
xlsx::write.xlsx(CD14.CD16.LCC.KEGG,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD14.CD16.LCC.KEGG',row.names = F)
xlsx::write.xlsx(CD14.CD16.RCC.KEGG,file = 'wclc/mono/EnrichmentResults.xlsx',append = T,sheetName = 'CD14.CD16.RCC.KEGG',row.names = F)

rm(list = ls());gc()

go.rst <- xlsx::read.xlsx('wclc/mono/EnrichmentResults.InShort.xlsx',sheetIndex = 1)
go.rst$`-Log10(pvalue)` <- -log10(go.rst$pvalue)

p <- ggdotplot(go.rst,x = 'Description',y = '-Log10(pvalue)',color = 'black',
          fill = 'group',size = 3,palette = RColorBrewer::brewer.pal(6, "Paired")) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        panel.border = element_rect(size = 1,colour = 'black',fill = NA),
        axis.title.x = element_blank())
ggsave("wclc/mono/enrich.go.pdf",p,width = 8,height = 6)


### correlation analysis 
mono <- readRDS('wclc/rds/mono.rds')
Idents(mono) <- 'celltype'

# avg.exp <- AverageExpression(mono,assays = 'RNA')$RNA
# cg <- names(tail(sort(apply(avg.exp,1,sd)),50))
# cor.T <- cor(avg.exp[cg,],method = 'spearman')
# mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(7,"RdBu"))(11))
# pheatmap::pheatmap(cor.T,color = mypal,border_color = NA,width = 5.5,height = 3.5,
#                    show_colnames = F,display_numbers = T,fontsize_number = 12,
#                    filename = 'wclc/mono/celltype.cor.pdf')
# 
# top.exp <- mono@assays$RNA@counts[top5.pct.var.genes,]
# top.exp <- t(as.matrix(top.exp))
# cor.T2 <- cor(top.exp,method = 'spearman')
# a <- pheatmap(cor.T2,color = mypal,border_color = NA,show_colnames = F)
# cor.T2 <- cor.T2[a$tree_row$order, a$tree_col$order]
# pheatmap::pheatmap(cor.T2,color = mypal,border_color = NA,width = 5,height = 3.5,
#                    show_colnames = F,
#                    filename = 'wclc/mono/top.var.cor.pdf')
# pheatmap::pheatmap(cor.T2,color = mypal,border_color = NA,width = 30,height = 29,
#                    show_colnames = F,
#                    filename = 'wclc/mono/top.var.cor.2.pdf')

lcc.mono <- subset(mono,cells=rownames(subset(mono@meta.data,TissueSiteSimple=='left')))
rcc.mono <- subset(mono,cells=rownames(subset(mono@meta.data,TissueSiteSimple=='right')))
lcc.mono <- FindVariableFeatures(lcc.mono)
rcc.mono <- FindVariableFeatures(rcc.mono)

mypal.1 <- rev(colorRampPalette(RColorBrewer::brewer.pal(7,"RdBu"))(7))
mypal.2 <- rev(colorRampPalette(RColorBrewer::brewer.pal(7,"RdYlGn"))(7))

lcc.avg.exp <- AverageExpression(lcc.mono,assays = 'RNA')$RNA
cg <- names(tail(sort(apply(lcc.avg.exp,1,sd)),100))
lcc.cor.T <- cor(lcc.avg.exp[cg,],method = 'spearman')
pheatmap(lcc.cor.T,color = mypal.1,border_color = NA,width = 5.5,height = 3.5,
                   show_colnames = F,display_numbers = T,fontsize_number = 12,
                   filename = 'wclc/mono/lcc.celltype.cor.pdf')

rcc.avg.exp <- AverageExpression(rcc.mono,assays = 'RNA')$RNA
cg <- names(tail(sort(apply(rcc.avg.exp,1,sd)),100))
rcc.cor.T <- cor(rcc.avg.exp[cg,],method = 'spearman')
pheatmap(rcc.cor.T,color = mypal.2,border_color = NA,width = 5.5,height = 3.5,
                   show_colnames = F,display_numbers = T,fontsize_number = 12,
                   filename = 'wclc/mono/rcc.celltype.cor.pdf')

lcc.top5.var.genes <- head(VariableFeatures(lcc.mono),100)
rcc.top5.var.genes <- head(VariableFeatures(rcc.mono),100)

lcc.top.exp <- lcc.mono@assays$RNA@counts[lcc.top5.var.genes,]
lcc.top.exp <- t(as.matrix(lcc.top.exp))
cor.T1 <- cor(lcc.top.exp,method = 'spearman')
a <- pheatmap(cor.T1,color = mypal.1,border_color = NA,show_colnames = F)
cor.T1 <- cor.T1[a$tree_row$order, a$tree_col$order]
pheatmap(cor.T1,color = mypal.1,border_color = NA,cellwidth = 2,cellheight = 2,
         show_colnames = F,filename = 'wclc/mono/lcc.mono.cor.var.exp.1.pdf')
pheatmap(cor.T1,color = mypal.1,border_color = NA,width = 30,height = 29,
         show_colnames = F,filename = 'wclc/mono/lcc.mono.cor.var.exp.2.pdf')

rcc.top.exp <- rcc.mono@assays$RNA@counts[rcc.top5.var.genes,]
rcc.top.exp <- t(as.matrix(rcc.top.exp))
cor.T2 <- cor(rcc.top.exp,method = 'spearman')
a <- pheatmap(cor.T2,color = mypal.2,border_color = NA,show_colnames = F)
cor.T2 <- cor.T2[a$tree_row$order, a$tree_col$order]
pheatmap(cor.T2,color = mypal.2,border_color = NA,cellheight = 2,
         cellwidth = 2,show_colnames = F,filename = 'wclc/mono/rcc.mono.cor.var.exp.1.pdf')
pheatmap(cor.T2,color = mypal.2,border_color = NA,width = 30,height = 29,
         show_colnames = F,filename = 'wclc/mono/rcc.mono.cor.var.exp.2.pdf')

save(cor.T1,cor.T2,file = 'wclc/mono/lcc.rcc.cor.top100.var.rda')



# ---------------------------------------------------------------------------------------------

mono <- readRDS('wclc/rds/mono.rds')

genes <- c('PPBP','PF4V1','CCL3','CXCL3')
cols <- c('#00C8FF','#DE4343')
p.list <- list()
for(i in 1:length(genes)){
  p.list[[i]] <- VlnPlot(mono,features = genes[i],group.by = 'TissueSiteSimple',
                         cols = cols,pt.size = 1e-10,y.max = 7.5) + NoLegend() +
    theme(axis.line = element_blank(), axis.title.x = element_blank(),
          panel.border = element_rect(size = 1,colour = 'black',fill = NA),
          axis.text.x = element_text(angle = 0,hjust = .5)) +
    ggpubr::stat_compare_means(comparisons = list(c('left','right')))
  ggsave(paste0('wclc/mono/',genes[i],'.pdf'),p.list[[i]],width = 3,height = 4)
}
p <- cowplot::plot_grid(plotlist = p.list,nrow = 1)
ggsave('wclc/mono/migration.features.pdf',p,width = 12,height = 4)

mypal <- c('#FF7F0E','#00C8FF','#DE4343')
p <- DimPlot(mono,split.by = 'TissueSiteSimple',reduction = 'tsne',
        group.by = 'celltype',cols = mypal) +
  theme(legend.position = 'top',plot.title = element_blank())
ggsave('wclc/mono/celltype.split.by.tissue.pdf',p,width = 8,height = 5)

p.list <- list()
for(i in 1:length(genes)){
  p.list[[i]] <- FeaturePlot(mono,features = genes[i],reduction = 'tsne',
                             cols = c('grey90','#DE4343')) +
    theme(axis.line = element_blank(),
          panel.border = element_rect(size = 1,colour = 'black',fill = NA))
  ggsave(paste0('wclc/mono/',genes[i],'.tsne.pdf'),p.list[[i]],width = 4.2,height = 4)
}
p <- cowplot::plot_grid(plotlist = p.list,nrow = 1)
ggsave('wclc/mono/migration.features.tsne.pdf',p,width = 14,height = 3.2)

