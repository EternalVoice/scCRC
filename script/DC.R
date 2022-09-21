# ----------------  DC ----------------
rm(list = ls());gc()
library(Seurat)

mye <- readRDS('wclc/rds/myeloid.tumor.sub.rds')
DC <- subset(mye, cells = rownames(subset(mye@meta.data, clMidwayPr == "DC")))
DC <- NormalizeData(object = DC, normalization.method = "LogNormalize", scale.factor = 10000)
DC <- FindVariableFeatures(DC,selection.method="vst",nfeatures=2000)
DC <- ScaleData(DC, features = VariableFeatures(DC))
DC <- RunPCA(object = DC, features = VariableFeatures(DC))
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
PC.num <- PCDeterminators(DC)
# Find neighbors and clusters with harmony batch correction
DC <- FindNeighbors(DC, dims = 1:PC.num, reduction = "pca")
DC <- FindClusters(DC, resolution = 0.2)
DC <- RunUMAP(DC, dims = 1:PC.num, reduction = "pca")
DC <- RunTSNE(DC,dims = 1:PC.num)
p <- DimPlot(DC,label = T,label.size = 5) + NoLegend() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/DC/DC.dr.pdf',p,width = 4,height = 4)

p <- DimPlot(DC,label = T,label.size = 5,reduction = "tsne") + NoLegend() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/DC/DC.dr.tsne.pdf',p,width = 6,height = 6)

p <- DimPlot(DC,label = T,label.size = 5,split.by = 'TissueSiteSimple',ncol = 1) + NoLegend() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/DC/DC.splitBy.rcc.lcc.dr.pdf',p,width = 4,height = 8)

mypal <- RColorBrewer::brewer.pal(7,"Paired")
p <- DimPlot(DC,cols = mypal,label = T,group.by = "celltype",reduction = "tsne") +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_blank())
ggsave(filename = "wclc/DC/celltype.tsne.pdf",p,width = 5,height = 5)

pathway <- xlsx::read.xlsx("wclc/xlsx/pathway.tmp.xlsx",sheetIndex = 4)
p <- ggplot(pathway, aes(x=Tissue,y=Description,color=Genes,size=-log10(pvalue))) +
  geom_point() + theme_bw() + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 14,colour = "black"),
        legend.position = "left",
        axis.title = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20),limits=c(1,4))
ggsave(filename = "wclc/DC/pathway.lcc.rcc.dotplot.pdf",p,width = 8,height = 6)


# ACKR1
p1 <- FeaturePlot(DC, features = 'ACKR1',cols = c("grey90","red")) + NoLegend() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/DC/ACKR1.pdf',p1,width = 4,height = 4)

DC$ACKR1 <- DC@assays$RNA@data["ACKR1",]
DC$GroupI <- ifelse(DC$ACKR1 > 0, "ACKR1+","ACKR1-")

saveRDS(DC, file = 'wclc/rds/DC.rds')

# mycols <- c('#1f77b4', '#39d486', '#43315c', '#00c8ff', '#de4343', '#ff7f0e', '#e377c2')
# scales::show_col(mycols)
p <- DimPlot(DC, group.by = 'GroupI', cols = c('#00c8ff','#de4343')) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(size = 1,colour = 'black')) +
  ggtitle('')
ggsave(filename = 'wclc/DC/ACKR1.pos.neg.pdf',p,width = 4,height = 4.2)

p <- DimPlot(DC, group.by = 'GroupI', cols = c('#00c8ff','#de4343'),reduction = "tsne") + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(size = 1,colour = 'black')) +
  ggtitle('')
ggsave(filename = 'wclc/DC/ACKR1.pos.neg.tsne.pdf',p,width = 6,height = 6.2)

all.markers <- FindAllMarkers(DC)

save(all.markers, file = 'wclc/DC/all.markers.c0-7.rda')

c1.markers <- all.markers[all.markers$cluster == 1,]

for(i in 1:100){
  p1 <- FeaturePlot(DC, features = c1.markers$gene[i],cols = c("grey90","red")) + NoLegend() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(size = 1,colour = 'black'))
  ggsave(paste0("wclc/DC/markers/c1/",i,'.',c1.markers$gene[i],".pdf"),plot = p1, width = 4,height = 4)
}

mycol <- RColorBrewer::brewer.pal(8, "Paired")
p <- VlnPlot(DC, features = c('ACKR1','SLC16A11','CD1A'),ncol = 1,
             pt.size = 0,group.by = "celltype",cols = mycol) + NoLegend()
ggsave(filename = 'wclc/DC/ACKR1.SLC16A11.CD1A.ct.pdf',p,width = 6,height = 8)

# # DC
# load('wclc/DC/all.markers.c0-7.rda')
# DC <- readRDS('wclc/rds/DC.rds')
# 
# for(i in 0:7){
#   if(!dir.exists(paste0("wclc/DC/markers/VlnPlot/c",i))) dir.create(paste0("wclc/DC/markers/VlnPlot/c",i),recursive = T)
# }
# all.markers <- all.markers[order(all.markers[,2],decreasing = T),]
# c1.markers <- all.markers[all.markers$cluster == 1,]
# 
# for(i in 1:50){
#   p1 <- FeaturePlot(DC, features = c7.markers$gene[i],cols = c("grey90","red")) + NoLegend() +
#     theme(axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           axis.line = element_blank(),
#           axis.title = element_blank(),
#           panel.border = element_rect(size = 1,colour = 'black'))
#   ggsave(paste0("wclc/DC/markers/c7/",i,'.',c7.markers$gene[i],".pdf"),plot = p1, width = 4,height = 4)
#   p2 <- VlnPlot(DC, features = c7.markers$gene[i],pt.size = 1e-5) + NoLegend() +
#     theme(axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           axis.line = element_blank(),
#           axis.title = element_blank(),
#           panel.border = element_rect(size = 1,colour = 'black'))
#   ggsave(paste0("wclc/DC/markers/VlnPlot/c7/",i,'.',c7.markers$gene[i],".pdf"),plot = p2, width = 6,height = 4)
# }

### ACKR1
rm(list = ls());gc()
DC <- readRDS('wclc/rds/DC.rds')

Idents(DC) <- DC$GroupI
markers <- FindMarkers(DC, group.by = 'GroupI', ident.1 = 'ACKR1+', ident.2 = 'ACKR1-')
# remove ribosomal and mitochondrial genes
# markers <- markers[!grepl("^RP[SL]",rownames(markers)),]
# markers <- markers[!grepl("^MT-",rownames(markers)),]
markers$cluster <- ifelse(markers$avg_logFC > 0, "ACKR1+","ACKR1-")
save(markers, file = 'wclc/DC/DC.ACKR1.pos_vs_neg.DEgenes.rda')
nrow(markers[markers$avg_logFC > 0,]) # 277
nrow(markers[markers$avg_logFC < 0,]) # 209

markers$Significance <- ifelse(markers$p_val < 0.01, TRUE, FALSE)
markers <- markers[order(markers[,2],decreasing = T),]
top10 <- rbind(head(markers,10),tail(markers,10))

mytheme <- theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 16),
                 axis.text = element_text(size = 12),
                 legend.title = element_text(size = 12),
                 legend.text = element_text(size = 12)) 

p.all.markers <- ggplot(markers, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')
ggsave(filename = 'wclc/DC/ACKR1.pos_vs_neg.DEgenes.pdf',p.all.markers,width = 5,height = 4)

EnrichmentAnalysis <- function(DEgenes, enrich.method = c("GO","KEGG")) {
  
  pacman::p_load(clusterProfiler,org.Hs.eg.db,dplyr,patchwork,enrichplot,ggplot2)
  
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
    p <- dotplot(ekegg, showCategory = 20)
    ggsave("KEGG.pdf",plot = p, device = "pdf", width = 12, height = 10)
    return(ekegg)
  }
}

# 
up <- rownames(markers[markers$avg_logFC > 0,])
dn <- rownames(markers[markers$avg_logFC < 0,])
ACKR1pos.up <- EnrichmentAnalysis(DEgenes = up, enrich.method = 'GO')
openxlsx::write.xlsx(ACKR1pos.up$ego_BP@result, file = 'wclc/DC/ACKR1pos.up.GO.xlsx')
ACKR1pos.dn <- EnrichmentAnalysis(DEgenes = dn, enrich.method = 'GO')
openxlsx::write.xlsx(ACKR1pos.dn$ego_BP@result, file = 'wclc/DC/ACKR1pos.dn.GO.xlsx')

## RCC
rcc.DC <- subset(DC, cells = rownames(subset(DC@meta.data, TissueSiteSimple == "right")))
p <- DimPlot(rcc.DC, group.by = 'GroupI', cols = c('#00c8ff','#de4343')) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(size = 1,colour = 'black')) +
  ggtitle('')
ggsave(filename = 'wclc/DC/rcc.ACKR1.pos.neg.pdf',p,width = 4,height = 4.2)

p <- DimPlot(rcc.DC, group.by = 'GroupI', cols = c('#00c8ff','#de4343'),reduction = "tsne") +
  ggtitle('')  +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/DC/rcc.ACKR1.pos.neg.tsne.pdf',p,width = 4,height = 4.2)

markers <- FindMarkers(rcc.DC, group.by = 'GroupI', ident.1 = 'ACKR1+', ident.2 = 'ACKR1-')
# remove ribosomal and mitochondrial genes
markers <- markers[!grepl("^RP[SL]",rownames(markers)),]
markers <- markers[!grepl("^MT-",rownames(markers)),]

nrow(markers[markers$avg_logFC > 0,]) # 269
nrow(markers[markers$avg_logFC < 0,]) # 208
markers$cluster <- ifelse(markers$avg_logFC > 0, "ACKR1+","ACKR1-")
save(markers, file = 'wclc/DC/RCC.DC.ACKR1.pos_vs_neg.DEgenes.rda')
markers$Significance <- ifelse(markers$p_val < 0.01, TRUE, FALSE)
markers <- markers[order(markers[,2],decreasing = T),]
top10 <- rbind(head(markers,10),tail(markers,10))
p.all.markers <- ggplot(markers, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')
ggsave(filename = 'wclc/DC/rcc.ACKR1.pos_vs_neg.DEgenes.pdf',p.all.markers,width = 5,height = 5)

up <- rownames(markers[markers$avg_logFC > 0,])
dn <- rownames(markers[markers$avg_logFC < 0,])
ACKR1pos.up <- EnrichmentAnalysis(DEgenes = up, enrich.method = 'GO')
openxlsx::write.xlsx(ACKR1pos.up$ego_BP@result, file = 'wclc/DC/rcc.ACKR1pos.up.GO.xlsx')
ACKR1pos.dn <- EnrichmentAnalysis(DEgenes = dn, enrich.method = 'GO')
openxlsx::write.xlsx(ACKR1pos.dn$ego_BP@result, file = 'wclc/DC/rcc.ACKR1pos.dn.GO.xlsx')

# LCC
lcc.DC <- subset(DC, cells = rownames(subset(DC@meta.data, TissueSiteSimple == "left")))
p <- DimPlot(lcc.DC, group.by = 'GroupI', cols = c('#00c8ff','#de4343')) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(size = 1,colour = 'black')) +
  ggtitle('')
ggsave(filename = 'wclc/DC/lcc.ACKR1.pos.neg.pdf',p,width = 4,height = 4.2)

p <- DimPlot(lcc.DC, group.by = 'GroupI', cols = c('#00c8ff','#de4343'),reduction = "tsne") + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(size = 1,colour = 'black')) +
  ggtitle('')
ggsave(filename = 'wclc/DC/lcc.ACKR1.pos.neg.tsne.pdf',p,width = 4,height = 4.2)

markers <- FindMarkers(lcc.DC, group.by = 'GroupI', ident.1 = 'ACKR1+', ident.2 = 'ACKR1-')
# remove ribosomal and mitochondrial genes
markers <- markers[!grepl("^RP[SL]",rownames(markers)),]
markers <- markers[!grepl("^MT-",rownames(markers)),]

nrow(markers[markers$avg_logFC > 0,]) # 322
nrow(markers[markers$avg_logFC < 0,]) # 193
markers$cluster <- ifelse(markers$avg_logFC > 0, "ACKR1+","ACKR1-")
save(markers, file = 'wclc/DC/LCC.DC.ACKR1.pos_vs_neg.DEgenes.rda')
markers$Significance <- ifelse(markers$p_val < 0.01, TRUE, FALSE)
markers <- markers[order(markers[,2],decreasing = T),]
top10 <- rbind(head(markers,10),tail(markers,10))
p.all.markers <- ggplot(markers, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')
ggsave(filename = 'wclc/DC/lcc.ACKR1.pos_vs_neg.DEgenes.pdf',p.all.markers,width = 5,height = 5)

up <- rownames(markers[markers$avg_logFC > 0,])
dn <- rownames(markers[markers$avg_logFC < 0,])
ACKR1pos.up <- EnrichmentAnalysis(DEgenes = up, enrich.method = 'GO')
openxlsx::write.xlsx(ACKR1pos.up$ego_BP@result, file = 'wclc/DC/lcc.ACKR1pos.up.GO.xlsx')
ACKR1pos.dn <- EnrichmentAnalysis(DEgenes = dn, enrich.method = 'GO')
openxlsx::write.xlsx(ACKR1pos.dn$ego_BP@result, file = 'wclc/DC/lcc.ACKR1pos.dn.GO.xlsx')


pathway <- xlsx::read.xlsx("wclc/DC/GO.summary.xlsx",sheetIndex = 1)
p <- ggplot(pathway, aes(x=Group,y=Description,color=Genes,size=-log10(pvalue))) +
  geom_point() + theme_bw() + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 14,colour = "black"),
        legend.position = "left",
        axis.title = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20),limits=c(1,4))
ggsave(filename = "wclc/DC/GO.summary.pdf",p,width = 8,height = 6)

pathway <- xlsx::read.xlsx("wclc/DC/GO.summary.xlsx",sheetIndex = 2)
p <- ggplot(pathway, aes(x=Group,y=Description,color=Genes,size=-log10(pvalue))) +
  geom_point() + theme_bw() + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 14,colour = "black"),
        legend.position = "left",
        axis.title = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20),limits=c(1,4))
ggsave(filename = "wclc/DC/LCC.ACKR1.GO.pdf",p,width = 8,height = 6)

pathway <- xlsx::read.xlsx("wclc/DC/GO.summary.xlsx",sheetIndex = 3)
p <- ggplot(pathway, aes(x=Group,y=Description,color=Genes,size=-log10(pvalue))) +
  geom_point() + theme_bw() + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 14,colour = "black"),
        legend.position = "left",
        axis.title = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20),limits=c(1,4))
ggsave(filename = "wclc/DC/RCC.ACKR1.GO.pdf",p,width = 8,height = 6)


load('wclc/DC/ACKR1/DC.ACKR1.pos_vs_neg.DEgenes.rda')
fea1 <- c('ACKR1','SLC16A11','C1QA','CD1A','C1QC','EPHA8')
fea2 <- c('PLLP','C9orf152','GZMH','LACTB')
fea <- c(fea1,fea2)
mypal <- rev(RColorBrewer::brewer.pal(11,"RdYlGn"))
p <- DotPlot(DC,features = fea, group.by = "GroupI") + 
  coord_flip() + 
  scale_color_gradientn(colours = mypal) +
  theme(panel.border = element_rect(size = 1,colour = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid.major = element_line(linetype = 3,colour = "grey50"))
ggsave("wclc/DC/ACKR1/markers.pdf",p,width = 4,height = 4)

# DC
# pDC: IL3RA, LILRA4, CLEC4C, PLD4, PHEX, PTCRA, IRF8, IRF7, GZMB, CXCR3
# pDC: PTPRC

# cDC1s: XCR1, CD103(ITGAE), CD8
# cDC2s: CD172, CD11b, CD4
# 
# iDC: CD1B, VASH1, F13A1, CD1E, MMP12M, FABP4, CLEC10A, SYT17, MS4A6A, CTNS
# 
# GUCA1A, CARD9, ABCG2, CD1A, PPARG, RAP1GAP, SLC7A8, GSTT1, FZD2, CSF1R
# HS3ST2, CH25H, LMAN2L, SLC26A6, BLVRB, NUDT9, PREP, TM7SF4, TACSTD2, CD1C

# aDC: CCL1, EBI3, INDO, LAMP3, OAS3


# pDC: CCR7, CD45RA(PTPRC), CD209, CLEC4C, LILRA4, NRP1, B220(PTPRC), SiglecH(SIGLEC7), IFN-α(IFNA1), IRF4, IRF7, IRF8
# cDC1: BTLA, CADM1, CD8A, CLEC9A, ITGAE, ITGAX, LY75, THBD, XCR1, BATF3, ID2, IRF8, ZBTB46
# cDC2: CD14, CD163, CLEC10A, NOTCH2, ITGAM, SIRPA, CX3CR1, CD1C, CD2, ID2, IRF4, KLF4, ZBTB46
# moDC: CD14, CD1A, CD1C, CD209, FCER1, ITGAM, MRC1, SIRPA, IRF4, KLF4, ZBTB46
# LC: CD1A, CD207, ID2
# markers <- c('IL3RA','CCR7','PTPRC','CD209','CLEC4C','LILRA4','NRP1', 'PTPRC', 'SIGLEC7', 'IFNA1', 'IRF4', 'IRF7', 'IRF8',
#              'BTLA', 'CADM1', 'CD8A', 'CLEC9A', 'ITGAE', 'ITGAX', 'LY75', 'THBD', 'XCR1', 'BATF3', 'ID2', 'IRF8', 'ZBTB46',
#              'CD14', 'CD163', 'CLEC10A', 'NOTCH2', 'ITGAM', 'SIRPA', 'CX3CR1', 'CD1C', 'CD2', 'ID2', 'IRF4', 'KLF4', 'ZBTB46',
#              'CD14', 'CD1A', 'CD1C', 'CD209', 'FCER1G', 'ITGAM', 'MRC1', 'SIRPA', 'IRF4', 'KLF4', 'ZBTB46',
#              'CD1A', 'CD207', 'ID2')


rm(list = ls());gc()
DC <- readRDS('wclc/rds/DC.rds')

# FeaturePlot(DC, features = 'IL3RA')
# FeaturePlot(DC, features = 'LILRA4')
# 
# FeaturePlot(DC, features = 'C1QA')
# FeaturePlot(DC, features = 'C1QC')
# FeaturePlot(DC, features = 'C1QL3')

DC$cl295v11SubFull <- as.character(DC$cl295v11SubFull)
old.ident <- names(table(DC$cl295v11SubFull))
new.ident <- c('DC1','DC2','C1Q+DC2','IL22RA2+DC','pDC','AS-DC','mregDC')
DC$celltype <- plyr::mapvalues(DC$cl295v11SubFull, from = old.ident, to = new.ident)

mycol <- c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3")

old.ident <- names(table(DC$seurat_clusters))
new.ident <- c('C0','C1','C2','C3','C4','C5','C6','C7')
DC$tmp <- plyr::mapvalues(DC$seurat_clusters, from = old.ident, to = new.ident)

p <- DimPlot(DC, group.by = 'tmp', cols = mycol,label = T,label.size = 5) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/DC/seurat_clusters.pdf',p,width = 5,height = 5)

mycol <- RColorBrewer::brewer.pal(8, "Paired")

p <- DimPlot(DC, group.by = 'celltype', cols = mycol) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/DC/celltype.pdf',p,width = 5,height = 5.2)

pDC <- subset(DC, cells = rownames(subset(DC@meta.data, celltype == 'pDC')))

lcc.pDC <- subset(pDC, cells = rownames(subset(DC@meta.data, TissueSiteSimple == 'left')))
rcc.pDC <- subset(pDC, cells = rownames(subset(DC@meta.data, TissueSiteSimple == 'right')))
# CDHR5
tag1 <- lcc.pDC@assays$RNA@counts['CDHR5',]
tag2 <- rcc.pDC@assays$RNA@counts['CDHR5',]
lcc.pDC$GroupI <- ifelse(tag1 > 0, "CDHR5+",'CDHR5-')
rcc.pDC$GroupI <- ifelse(tag2 > 0, "CDHR5+",'CDHR5-')
table(lcc.pDC$GroupI)
sum(unname(table(lcc.pDC$GroupI)))
unname(table(lcc.pDC$GroupI))[2]/sum(unname(table(lcc.pDC$GroupI)))
table(rcc.pDC$GroupI)
sum(unname(table(rcc.pDC$GroupI)))
unname(table(rcc.pDC$GroupI))[2]/sum(unname(table(rcc.pDC$GroupI)))
# ACKR1
tag1 <- lcc.pDC@assays$RNA@counts['ACKR1',]
tag2 <- rcc.pDC@assays$RNA@counts['ACKR1',]
lcc.pDC$GroupI <- ifelse(tag1 > 0, "ACKR1+",'ACKR1-')
rcc.pDC$GroupI <- ifelse(tag2 > 0, "ACKR1+",'ACKR1-')
table(lcc.pDC$GroupI)
sum(unname(table(lcc.pDC$GroupI)))
unname(table(lcc.pDC$GroupI))[2]/sum(unname(table(lcc.pDC$GroupI)))
table(rcc.pDC$GroupI)
sum(unname(table(rcc.pDC$GroupI)))
unname(table(rcc.pDC$GroupI))[2]/sum(unname(table(rcc.pDC$GroupI)))

# DC cell proportion
mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "black"),
                 panel.grid = element_blank(),
                 axis.text = element_text(size = 12,face = "bold",colour = "black"),
                 axis.title.x = element_blank(),
                 axis.title = element_text(size = 14,face = "bold",colour = "black"),
                 legend.position = 'none')

p <- ggplot(DC@meta.data, aes(x = TissueSiteSimple, fill = celltype)) +
  geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values = mycol) +
  mytheme +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(filename = 'wclc/DC/cell.prop.pdf',p,width = 5,height = 10)

saveRDS(DC, file = 'wclc/rds/DC.rds')

# cell correlation
DC <- readRDS('wclc/rds/DC.rds')
Idents(DC) <- "celltype"
avg <- AverageExpression(LCC.DC)$RNA
cg <- names(tail(sort(apply(avg,1,sd)),1000))
cols <- rev(colorRampPalette(RColorBrewer::brewer.pal(7,"RdYlBu"))(256))
pheatmap::pheatmap(cor(avg[cg,],method = "spearman"),
                   color = cols,border_color = NA,
                   filename = "wclc/DC/DC.subset.corr.pdf",
                   width = 5.4,height = 5)

Idents(RCC.DC) <- "celltype"
rcc.avg <- AverageExpression(RCC.DC)$RNA
cg <- names(tail(sort(apply(rcc.avg,1,sd)),1000))
pheatmap::pheatmap(cor(rcc.avg[cg,],method = "spearman"),
                   color = cols,border_color = NA)


### Pseudotime

library(Seurat)
library(monocle)
table(lcc.DC$celltype)

lcc.dat <- as(as.matrix(lcc.DC@assays$RNA@counts), "sparseMatrix")
pd <- new('AnnotatedDataFrame', data = lcc.DC@meta.data)
fData <- data.frame(gene_short_name = rownames(lcc.dat), row.names = rownames(lcc.dat))
fd <- new('AnnotatedDataFrame', data = fData)
lcc.mycds <- newCellDataSet(lcc.dat, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
lcc.mycds <- estimateSizeFactors(lcc.mycds)
lcc.mycds <- estimateDispersions(lcc.mycds,relative_expr = T)
# 使用monocle选择的高变基因
disp.tbl <- dispersionTable(lcc.mycds)
disp.genes <- subset(disp.tbl, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
lcc.mycds <- setOrderingFilter(lcc.mycds, disp.genes)
p1 <- plot_ordering_genes(lcc.mycds)
# 使用disp.genes展开后续分析
lcc.mycds <- reduceDimension(lcc.mycds, max_components = 2, method = "DDRTree")
lcc.mycds <- orderCells(lcc.mycds,root_state = 2)
saveRDS(lcc.mycds, file = "wclc/rds/LCC.DC.monocle2.rds")

mycol <- RColorBrewer::brewer.pal(8, "Paired")
source("wclc/pseudotime_heatmap.R")
p1 <- plot_cell_trajectory(lcc.mycds, color_by = "celltype") + 
  scale_color_manual(values = mycol) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("wclc/DC/lcc.DC.trajectory.pdf",p1,width = 6,height = 4)

p2 <- plot_cell_trajectory(lcc.mycds, color_by = "Pseudotime") +
  scale_color_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256))) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("wclc/DC/lcc.DC.pseudo.pdf",p2,width = 6,height = 4)

# RCC
rcc.dat <- as(as.matrix(rcc.DC@assays$RNA@counts), "sparseMatrix")
pd <- new('AnnotatedDataFrame', data = rcc.DC@meta.data)
fData <- data.frame(gene_short_name = rownames(rcc.dat), row.names = rownames(rcc.dat))
fd <- new('AnnotatedDataFrame', data = fData)
rcc.mycds <- newCellDataSet(rcc.dat, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
rcc.mycds <- estimateSizeFactors(rcc.mycds)
rcc.mycds <- estimateDispersions(rcc.mycds,relative_expr = T)
# 使用monocle选择的高变基因
disp.tbl <- dispersionTable(rcc.mycds)
disp.genes <- subset(disp.tbl, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
rcc.mycds <- setOrderingFilter(rcc.mycds, disp.genes)
p1 <- plot_ordering_genes(rcc.mycds)
# 使用disp.genes展开后续分析
rcc.mycds <- reduceDimension(rcc.mycds, max_components = 2, method = "DDRTree")
rcc.mycds <- orderCells(rcc.mycds,root_state = 2)
saveRDS(rcc.mycds, file = "wclc/rds/RCC.DC.monocle2.rds")

p1 <- plot_cell_trajectory(rcc.mycds, color_by = "celltype") + 
  scale_color_manual(values = mycol) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("wclc/DC/rcc.DC.trajectory.pdf",p1,width = 6,height = 4)

p2 <- plot_cell_trajectory(rcc.mycds, color_by = "Pseudotime") +
  scale_color_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256))) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("wclc/DC/rcc.DC.pseudo.pdf",p2,width = 6,height = 4)

## DE analysis
Idents(lcc.DC) <- "celltype"
lcc.all.markers <- FindAllMarkers(lcc.DC,only.pos = T)

# mregDC
markers.mregDC <- lcc.all.markers[lcc.all.markers$cluster=="mregDC",]
lcc.mregDC.go <- EnrichmentAnalysis(DEgenes = markers.mregDC$gene,enrich.method = "GO")
# Highlight: HDAC1,TNIP1,CX3CL1,TLR8,BIRC3,LITAF,PRDX4,BIRC2
xlsx::write.xlsx(lcc.mregDC.go@result, file = "wclc/DC/enrichment/LCC.GO.xlsx",sheetName = "mregDC",append = F,row.names = F)

# DC2
markers.DC2 <- lcc.all.markers[lcc.all.markers$cluster=="DC2",]
lcc.DC2.go <- EnrichmentAnalysis(DEgenes = markers.DC2$gene,enrich.method = "GO")
# Highlight: CD1A/HLA-DPA1/TAP1/HLA-DOA/HLA-DQA1/HLA-DRB1/HLA-DRA
xlsx::write.xlsx(lcc.DC2.go@result, file = "wclc/DC/enrichment/LCC.GO.xlsx",sheetName = "DC2",append = T,row.names = F)

# C1Q+DC2
markers.C1Q_DC2 <- lcc.all.markers[lcc.all.markers$cluster=="C1Q+DC2",]
lcc.C1Q_DC2.go <- EnrichmentAnalysis(DEgenes = markers.C1Q_DC2$gene,enrich.method = "GO")
# Highlight: C1QC/C1QA/NCR3/CD1D/CD1A/CD1B/PRDX1/CRP/HSPD1/HFE/IGHG1/IGHD
xlsx::write.xlsx(lcc.C1Q_DC2.go@result, file = "wclc/DC/enrichment/LCC.GO.xlsx",sheetName = "C1Q+DC2",append = T,row.names = F)

# IL22RA2+DC
markers.IL22RA2_DC2 <- lcc.all.markers[lcc.all.markers$cluster=="IL22RA2+DC",]
lcc.IL22RA2_DC2.go <- EnrichmentAnalysis(DEgenes = markers.IL22RA2_DC2$gene,enrich.method = "GO")
# Highlight: TNF/CD74/HSPA1B/CRP/CD14
xlsx::write.xlsx(lcc.IL22RA2_DC2.go@result, file = "wclc/DC/enrichment/LCC.GO.xlsx",sheetName = "IL22RA2+DC",append = T,row.names = F)

# DC1
markers.DC1 <- lcc.all.markers[lcc.all.markers$cluster=="DC1",]
lcc.DC1.go <- EnrichmentAnalysis(DEgenes = markers.DC1$gene,enrich.method = "GO")
# Highlight: HLA-DQB2/HLA-DPA1/HLA-DOA/PSMB9/HLA-DQA1/KIF2A/HLA-C/PSMD5/TAP2
xlsx::write.xlsx(lcc.DC1.go@result, file = "wclc/DC/enrichment/LCC.GO.xlsx",sheetName = "DC1",append = T,row.names = F)

# pDC
markers.pDC <- lcc.all.markers[lcc.all.markers$cluster=="pDC",]
lcc.pDC.go <- EnrichmentAnalysis(DEgenes = markers.pDC$gene,enrich.method = "GO")
# Highlight: GPSM3/GPR18/VEGFB/AKIRIN1/RAC2/EDN2/IL23A/SPICE1/CEP63/PLK2/CHMP5/BRCA2
xlsx::write.xlsx(lcc.pDC.go@result, file = "wclc/DC/enrichment/LCC.GO.xlsx",sheetName = "pDC",append = T,row.names = F)

# AS-DC
markers.AS.DC <- lcc.all.markers[lcc.all.markers$cluster=="AS-DC",]
lcc.AS.DC.go <- EnrichmentAnalysis(DEgenes = markers.AS.DC$gene,enrich.method = "GO")
# Highlight: TERF2/MAPK3/NBN
xlsx::write.xlsx(lcc.AS.DC.go@result, file = "wclc/DC/enrichment/LCC.GO.xlsx",sheetName = "AS-DC",append = T,row.names = F)

# heatmap
gene <- c(
  'HDAC1','TNIP1','CX3CL1','TLR8','BIRC3','LITAF','PRDX4','BIRC2',  # mregDC
  'TNF','CD74','HSPA1B','CRP','CD14',  # IL22RA2+DC
  'HLA-DQB2','HLA-DPA1','HLA-DOA','PSMB9','HLA-DQA1','KIF2A','HLA-C','PSMD5','TAP2',  # DC1
  'CD1A','HLA-DPA1','TAP1','HLA-DOA','HLA-DQA1','HLA-DRB1','HLA-DRA',  # DC2
  'C1QC','C1QA','NCR3','CD1D','CD1A','CD1B','PRDX1','CRP','HSPD1','HFE','IGHG1','IGHD',  # C1Q+DC2
  'TERF2','MAPK3','NBN',  # AS-DC
  'GPSM3','GPR18','VEGFB','AKIRIN1','RAC2','EDN2','IL23A','SPICE1','CEP63','PLK2','CHMP5','BRCA2'  # pDC
)

Idents(lcc.DC) <- "celltype"
cols <- rev(colorRampPalette(RColorBrewer::brewer.pal(10,"RdYlBu"))(256))
mycol <- RColorBrewer::brewer.pal(8, "Paired")
p <- DoHeatmap(lcc.DC,features = lcc.all.markers$gene,group.by = "celltype",group.colors = mycol) +
  scale_fill_gradientn(colours = cols)
ggsave(filename = "wclc/DC/LCC.DEgenes.heatmap.pdf",p,width = 10,height = 10)

# RCC
Idents(rcc.DC) <- "celltype"
rcc.all.markers <- FindAllMarkers(rcc.DC,only.pos = T)

DEgenes <- xlsx::read.xlsx("wclc/DC/DEgene.count.xlsx",sheetIndex = 1)

p <- ggpubr::ggbarplot(DEgenes,x = "celltype",y = "Degenes",fill = "Tissue",combine = T,
                       palette = c('#00c8ff','#de4343')) + coord_flip()
ggsave(filename = "wclc/DC/DEgene.count.pdf",p,width = 4,height = 6)

save(lcc.all.markers,rcc.all.markers, file = "wclc/DC/lcc.rcc.all.markers.rda")

# mregDC
markers.mregDC <- rcc.all.markers[rcc.all.markers$cluster=="mregDC",]
rcc.mregDC.go <- EnrichmentAnalysis(DEgenes = markers.mregDC$gene,enrich.method = "GO")
# Highlight: HDAC1/TNIP1/CX3CL1/LITAF/TLR8/BIRC3/PRDX4/BIRC2
# I-kappaB kinase/NF-kappaB signaling
xlsx::write.xlsx(rcc.mregDC.go@result, file = "wclc/DC/enrichment/RCC.GO.xlsx",sheetName = "mregDC",append = F,row.names = F)

# DC2
markers.DC2 <- rcc.all.markers[rcc.all.markers$cluster=="DC2",]
rcc.DC2.go <- EnrichmentAnalysis(DEgenes = markers.DC2$gene,enrich.method = "GO")
# Highlight: CD1A/HLA-DOA/HLA-DQA1/HLA-DPA1/HLA-DRB1/TAP1/CD1B/HLA-DRA
# antigen processing and presentation
xlsx::write.xlsx(rcc.DC2.go@result, file = "wclc/DC/enrichment/RCC.GO.xlsx",sheetName = "DC2",append = T,row.names = F)

# C1Q+DC2
markers.C1Q_DC2 <- rcc.all.markers[rcc.all.markers$cluster=="C1Q+DC2",]
rcc.C1Q_DC2.go <- EnrichmentAnalysis(DEgenes = markers.C1Q_DC2$gene,enrich.method = "GO")
# Highlight: CD1B/CD1A/CD1D/TAP1/HLA-DOA/HLA-DQA1/HLA-DRB1/HFE/HLA-DRA
# antigen processing and presentation
xlsx::write.xlsx(rcc.C1Q_DC2.go@result, file = "wclc/DC/enrichment/RCC.GO.xlsx",sheetName = "C1Q+DC2",append = T,row.names = F)

# IL22RA2+DC
markers.IL22RA2_DC2 <- lcc.all.markers[rcc.all.markers$cluster=="IL22RA2+DC",]
rcc.IL22RA2_DC2.go <- EnrichmentAnalysis(DEgenes = markers.IL22RA2_DC2$gene,enrich.method = "GO")
# Highlight: IL1RL1/TNF
# regulation of chemokine secretion
xlsx::write.xlsx(rcc.IL22RA2_DC2.go@result, file = "wclc/DC/enrichment/RCC.GO.xlsx",sheetName = "IL22RA2+DC",append = T,row.names = F)

# DC1
markers.DC1 <- rcc.all.markers[rcc.all.markers$cluster=="DC1",]
rcc.DC1.go <- EnrichmentAnalysis(DEgenes = markers.DC1$gene,enrich.method = "GO")
# Highlight: HLA-DQB2/HLA-DPA1/HLA-DOA/HLA-DQA1/PSMB9/HLA-DRB1/HLA-DRA/HLA-C/TAP1
# antigen processing and presentation
xlsx::write.xlsx(rcc.DC1.go@result, file = "wclc/DC/enrichment/RCC.GO.xlsx",sheetName = "DC1",append = T,row.names = F)

# pDC
markers.pDC <- rcc.all.markers[rcc.all.markers$cluster=="pDC",]
rcc.pDC.go <- EnrichmentAnalysis(DEgenes = markers.pDC$gene,enrich.method = "GO")
# Highlight: GPR18/GPSM3/VEGFB/RAC2/IL23A/CCR1/BST1/CCL3/SPICE1/CEP63/PLK2/BRCA2/CHMP5
# regulation of leukocyte chemotaxis
# centrosome duplication
xlsx::write.xlsx(rcc.pDC.go@result, file = "wclc/DC/enrichment/RCC.GO.xlsx",sheetName = "pDC",append = T,row.names = F)

# AS-DC
markers.AS.DC <- rcc.all.markers[rcc.all.markers$cluster=="AS-DC",]
rcc.AS.DC.go <- EnrichmentAnalysis(DEgenes = markers.AS.DC$gene,enrich.method = "GO")
# Highlight: MAPK3/PNKP/TERF2/ERCC1
# telomere capping
xlsx::write.xlsx(rcc.AS.DC.go@result, file = "wclc/DC/enrichment/RCC.GO.xlsx",sheetName = "AS-DC",append = T,row.names = F)

# heatmap
gene <- c(
  'HDAC1','TNIP1','CX3CL1','LITAF','TLR8','BIRC3','PRDX4','BIRC2',  # mregDC
  'IL1RL1','TNF',  # IL22RA2+DC
  'HLA-DQB2','HLA-DPA1','HLA-DOA','HLA-DQA1','PSMB9','HLA-DRB1','HLA-DRA','HLA-C','TAP1',  # DC1
  'CD1A','HLA-DOA','HLA-DQA1','HLA-DPA1','HLA-DRB1','TAP1','CD1B','HLA-DRA',  # DC2
  'CD1B','CD1A','CD1D','TAP1','HLA-DOA','HLA-DQA1','HLA-DRB1','HFE','HLA-DRA',  # C1Q+DC2
  'MAPK3','PNKP','TERF2','ERCC1',  # AS-DC
  'GPR18','GPSM3','VEGFB','RAC2','IL23A','CCR1','BST1','CCL3','SPICE1','CEP63','PLK2','BRCA2','CHMP5'  # pDC
)

DC <- readRDS("wclc/rds/DC.rds")
rcc.DC <- subset(DC, cells = rownames(subset(DC@meta.data, TissueSiteSimple == "right")))
Idents(rcc.DC) <- "celltype"
cols <- rev(colorRampPalette(RColorBrewer::brewer.pal(10,"RdYlBu"))(256))
mycol <- RColorBrewer::brewer.pal(8, "Paired")
p <- DoHeatmap(rcc.DC,features = rcc.all.markers$gene,group.by = "celltype",group.colors = mycol) +
  scale_fill_gradientn(colours = cols)
ggsave(filename = "wclc/DC/RCC.DEgenes.heatmap.pdf",p,width = 10,height = 10)

# pDC
pDC <- subset(DC,cells = rownames(subset(DC@meta.data, celltype=="pDC")))

# CCL3-CCR1
p <- VlnPlot(pDC, features = "CCL3",group.by = "TissueSiteSimple") + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0),
        panel.border = element_rect(size = 1,colour = "black"))
ggsave(filename = "wclc/DC/CCL3.pdf",p,width = 3,height = 4)

p <- VlnPlot(pDC, features = "CCR1",group.by = "TissueSiteSimple") + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0),
        panel.border = element_rect(size = 1,colour = "black"))
ggsave(filename = "wclc/DC/CCR1.pdf",p,width = 3,height = 4)

## Tissue preference
rm(list = ls());gc()
DC <- readRDS("wclc/rds/DC.rds")

chiq.L <- chisq.test(table(DC$celltype,DC$TissueSiteSimple)[,1])
chiq.R <- chisq.test(table(DC$celltype,DC$TissueSiteSimple)[,2])
Left <- chiq.L$observed/chiq.L$expected
Right <- chiq.R$observed/chiq.R$expected
result <- rbind(Left,Right)
result <- t(result)
wilcox.test(result)$p.value

dat <- reshape2::melt(result)
colnames(dat) <- c("CellTypes","Group","Roe")
p <- ggplot(dat, aes(x=CellTypes,y = Roe,fill=Group)) + 
  geom_bar(stat = "identity",position = "dodge",color="black",size=1.3) +
  geom_text(aes(label=signif(Roe,3),size=6), position = position_dodge2(width = 0.9,preserve = 'single'),
            vjust = -0.20, hjust = 0.5) +
  scale_fill_manual(values = c('#4E9BD2','#DE4343')) +
  theme_bw() + coord_flip() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14,angle =0,hjust = .5,face = "bold",colour = 'black'),
        axis.title.x = element_text(size = 16,colour = 'black'),
        axis.text.y = element_text(size = 14, face = "bold",colour = 'black'),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black'),
        legend.title = element_text(size = 14,colour = 'black'),
        legend.text = element_text(size = 13,colour = 'black')) +
  annotate("text",x=1.3,y=1.5,label=paste0("p = ",wilcox.test(result)$p.value),size=5)
ggsave("wclc/DC/Roe.pdf",p,width = 6,height = 8)


## rogue
# library(ROGUE)
# library(tibble)
# DC <- readRDS("wclc/rds/DC.rds")
# expr <- as.matrix(DC@assays$RNA@counts)
# meta <- DC@meta.data
# filt.expr <- matr.filter(expr,min.cells = 10,min.genes = 10)
# rogue.res <- rogue(filt.expr,labels = meta$celltype, 
#                    samples = meta$orig.ident,
#                    platform = 'UMI', span = 0.6)
# rogue.boxplot(rogue.res)

tmp <- table(DC$orig.ident, DC$celltype)
tmp.2 <- table(DC$orig.ident, DC$TissueSiteSimple)
tmp.3 <- table(DC$orig.ident)
dat <- cbind(tmp.3,tmp.2,tmp)
dat <- as.data.frame(dat)
colnames(dat)[1] <- "Total"
dat$left <- ifelse(dat$left == "0", "RCC","LCC")
dat$right <- NULL
colnames(dat)[2] <- "TissueSite"

STAT <- c()
for(j in 1:nrow(dat)){
  stat <- NULL
  for(i in 3:ncol(dat)){
    stat[i] <- dat[j,i]/dat[j,1]*100
  }
  STAT <- rbind(STAT, stat)
}
colnames(STAT) <- colnames(dat)
STAT <- as.data.frame(STAT)
STAT$TissueSite <- dat$TissueSite
STAT$Total <- NULL

library(ggpubr)

for(i in 2:8){
  p <- ggpaired(STAT,x = 'TissueSite',y = colnames(STAT)[i],
                color = 'TissueSite',add = 'jitter',
                palette = c("#4E9BD3","#C72719"),size = 1.2,
                ylab = "Fraction of cells (%)",
                title = colnames(STAT)[i],
                line.color = 'gray',
                line.size = 0.4) +
    theme(legend.position = 'top',axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14,colour = 'black'),
          axis.text = element_text(size = 12,colour = 'black'),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0("wclc/DC/subset/",colnames(STAT)[i],'.2.pdf'),p,width = 3,height = 3.5)
}

# for(i in 2:8){
#   p <- ggboxplot(STAT,x = 'TissueSite',y = colnames(STAT)[i],color = 'TissueSite',
#                  add = 'jitter',
#                  palette = c("#4E9BD3","#C72719"),size = 1.2,
#                  ylab = "Fraction of cells (%)",
#                  title = colnames(STAT)[i]) +
#     theme(legend.position = 'none',axis.title.x = element_blank(),
#           axis.title.y = element_text(size = 14,colour = 'black'),
#           axis.text = element_text(size = 12,colour = 'black'),
#           plot.title = element_text(hjust = 0.5)) +
#     stat_compare_means(comparisons = list(c("NR","R")))
#   ggsave(paste0("wclc/DC/subset/",colnames(STAT)[i],'.pdf'),p,width = 3,height = 3.5)
# }

p <- VlnPlot(DC,features = 'TREM2',group.by = 'celltype') + NoLegend() +
  theme(axis.title.x = element_blank())
ggsave('wclc/DC/TREM2.pdf',p,width = 4,height = 3)
# interferon-stimulated gene
isg <- c('MX2','ISG15','IRF7','BST2','IFITM2','IFI27')
DC <- AddModuleScore(DC,features = list(isg),name = 'ISG')
mypal <- RColorBrewer::brewer.pal(7,"Paired")
mycomparison <- list(
  c('pDC','AS-DC'),
  c('pDC','C1Q+DC2'),
  c('pDC','DC1'),
  c('pDC','DC2'),
  c('pDC','IL22RA2+DC'),
  c('pDC','mregDC')
)
p <- VlnPlot(DC,features = 'ISG1',group.by = 'celltype',pt.size = 0,cols = mypal,y.max = 1.7) + 
  NoLegend() + ylab('ISG Score') +
  geom_boxplot(width = 0.1, fill = 'white') +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank()) +
  stat_compare_means(comparisons = mycomparison,label = 'p.signif')
ggsave("wclc/DC/ISG.pdf",p,width = 6,height = 4)

saveRDS(DC,file = 'wclc/rds/DC.rds')


# ------------------------ACKR1+ DC------------------------------------------------

DC <- readRDS('wclc/rds/DC.rds')
mypal <- RColorBrewer::brewer.pal(7,"Paired")

customVlnPlot <- function(obj,features){
  VlnPlot(obj,features = features,pt.size = 0,cols = mypal,group.by = 'celltype') + 
    NoLegend() + coord_flip() +
    theme(axis.text.x = element_text(angle = 0),
          axis.title = element_blank(),
          panel.border = element_rect(size = 1,colour = 'black',fill = NA),
          axis.line = element_blank())
}
p1  <- customVlnPlot(DC,'ACKR1')
p2 <- customVlnPlot(DC,'SLC16A11')
p3 <- customVlnPlot(DC,'CD1A')
p <- p1 | p2 | p3
ggsave('wclc/DC/ACKR1/ACKR1.SLC16A11.CD1A.pdf',p,width = 8,height = 6)

p <- DimPlot(DC,cells.highlight = rownames(DC@meta.data[DC@meta.data$GroupI=='ACKR1+',]),reduction = 'tsne') + 
  NoLegend() + ggtitle('ACKR1+ DC cells') +
  theme(axis.line = element_blank(),
        plot.title = element_text(hjust = .5),
        panel.border = element_rect(size = 1,color = 'black',fill = NA))
ggsave('wclc/DC/ACKR1/ACKR1.tsne.pdf',p,width = 4,height = 4.1)

mypal <- RColorBrewer::brewer.pal0(7,"Paired")
p <- DimPlot(DC,cols = mypal,label = T,group.by = "celltype",
             reduction = "tsne",label.size = 4) +
  ggtitle('') +
  theme(axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1,color = 'black',fill = NA))
ggsave(filename = "wclc/DC/celltype.tsne2.pdf",p,width = 4,height = 4)

cells1 <- rownames(subset(DC@meta.data,celltype=='DC2'))
cells2 <- rownames(subset(DC@meta.data,celltype=='C1Q+DC2'))
DC2 <- subset(DC,cells = c(cells1,cells2))

DC2 <- NormalizeData(object = DC2, normalization.method = "LogNormalize", scale.factor = 10000)
DC2 <- FindVariableFeatures(DC2,selection.method="vst",nfeatures=2000)
DC2 <- ScaleData(DC2, features = VariableFeatures(DC2))
DC2 <- RunPCA(object = DC2, features = VariableFeatures(DC2))
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
PC.num <- PCDeterminators(DC2)
# Find neighbors and clusters with harmony batch correction
DC2 <- FindNeighbors(DC2, dims = 1:20, reduction = "pca")
DC2 <- FindClusters(DC2, resolution = 0.2)
DC2 <- RunUMAP(DC2, dims = 1:20, reduction = "pca")
DC2 <- RunTSNE(DC2,dims = 1:20)
saveRDS(DC2,file = 'wclc/rds/DC2.rds')

p <- DimPlot(DC2,cells.highlight = rownames(DC2@meta.data[DC2@meta.data$GroupI=='ACKR1+',]),
             reduction = 'tsne') + 
  NoLegend() + ggtitle('ACKR1+ DC cells') +
  theme(axis.line = element_blank(),
        plot.title = element_text(hjust = .5),
        panel.border = element_rect(size = 1,color = 'black',fill = NA))
ggsave('wclc/DC/ACKR1/ACKR1.DC2.pdf',p,width = 4,height = 4.1)

mypal <- RColorBrewer::brewer.pal(2,"Paired")
p <- DimPlot(macro,cols = mypal,label = T,label.size = 5) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1,colour = 'black'))

# LCC.DC <- subset(DC,cells = rownames(subset(DC@meta.data, TissueSiteSimple == 'left')))
# RCC.DC <- subset(DC,cells = rownames(subset(DC@meta.data, TissueSiteSimple == 'right')))

proprtion <- openxlsx::read.xlsx('wclc/DC/ACKR1/proportion.xlsx',sheet = 1)
p <- ggbarplot(proprtion,x = 'celltype',y = 'proportion',
          color = 'tissue',fill = 'tissue', palette = c("#4E9BD3","#C72719"),
          merge = T) + coord_flip() +
  ylab('ACKR1+ cells (%)') +
  theme(axis.line = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black',fill = NA))
ggsave('wclc/DC/ACKR1/proportion.pdf',p,width = 4,height = 6)


# ------------- Append --------------------
library(Seurat)
library(ggplot2)
DC <- readRDS('wclc/rds/DC.rds')
ACKRposDC <- subset(DC, GroupI == 'ACKR1+')

tmp <- cbind(table(DC$orig.ident,DC$TissueSiteSimple),table(DC$orig.ident,DC$GroupI))
tmp <- as.data.frame(tmp)
tmp$Sample <- rownames(tmp)
tmp$group <- ifelse(tmp$left > 0, 'LCC','RCC')
tmp$cells <- tmp$left + tmp$right
tmp$prop <- tmp$`ACKR1+` / tmp$cells * 100
p <- ggpaired(tmp,x = 'group',y = 'prop',color = 'group',palette = c('#4E9BD2','#DB4740'),
         line.color = 'gray',line.size = 0.4,short.panel.labs = F) +
  stat_compare_means(method = 'anova') + ylab('ACKR1+ DC proportion (%)') +
  theme(legend.position = 'none', axis.title.x = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black',fill = NA))
ggsave('wclc/DC/ACKR1/ACKR1.prop.pdf',p,width = 3,height = 3.6)
saveRDS(tmp, file = 'wclc/DC/ACKR1/ACKR1.prop.rds')

genes <- c('CD8A','ACKR1','SLC16A11','C1QA','CD1A','C1QC','EPHA8')
Idents(DC) <- 'orig.ident'
expr <- AverageExpression(DC,features = genes)$RNA
expr <- as.data.frame(t(expr))
p.list <- list()
for(i in 1:(length(genes)-1)){
  p.list[[i]] <- ggplot(expr, aes_string(x = colnames(expr)[1], y = colnames(expr)[i+1])) +
    geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25,'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
}

p <- cowplot::plot_grid(plotlist = p.list,nrow = 1,ncol = 6)
ggsave('wclc/DC/ACKR1/sigs_vs_CD8A.pdf',p,width = 18,height = 3)

library(openxlsx)
library(survival)
library(survminer)

coad.expr <- read.table('data/TCGA/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt', sep = "\t", header = T)
coad.expr <- coad.expr[!duplicated(coad.expr$Hugo_Symbol),]
rownames(coad.expr) <- coad.expr$Hugo_Symbol
coad.expr <- coad.expr[, -c(1,2)]
colnames(coad.expr) <- gsub("\\.", "-", colnames(coad.expr))
coad.expr <- data.frame(t(coad.expr))

clin <- read.table('data/TCGA/data_clinical_patient.txt',skip = 4, header = T,sep = '\t')
small_clin <- data.frame(clin$PATIENT_ID, clin$OS_MONTHS, clin$OS_STATUS, clin$PRIMARY_SITE_PATIENT)
small_clin <- small_clin[small_clin$clin.OS_STATUS != "[Not Available]",]
small_clin <- small_clin[small_clin$clin.PRIMARY_SITE_PATIENT != "[Not Available]",]
small_clin <- small_clin[small_clin$clin.PRIMARY_SITE_PATIENT != "[Discrepancy]",]
small_clin$clin.PRIMARY_SITE_PATIENT <- as.character(small_clin$clin.PRIMARY_SITE_PATIENT)

dat0 <- coad.expr[,genes[-1]]
index <- match(gsub('-01','',rownames(dat0)),small_clin$clin.PATIENT_ID)
df <- cbind.data.frame(small_clin[index[!is.na(index)],],dat0[-which(is.na(index)),])
rownames(df) <- df$clin.PATIENT_ID
df <- df[,-c(1,4)]
df$clin.OS_STATUS <- stringr::str_split(df$clin.OS_STATUS,':',simplify = T)[,1]
names(df)[1:2] <- c("OS",'Status')
df$Status <- as.numeric(df$Status)
df$OS <- as.numeric(df$OS)

b <- c()
for(i in 1:nrow(df)){
  a <- sum(df[i,3:8])/6
  b <- c(b,a)
}
df$ACKR1pos <- b
df.new <- df[,c(1,2,9)]
cox_model <- coxph(Surv(OS,Status)~., data = df.new)
risk_score <- predict(cox_model, type = 'risk',newdata = df.new)
risk_tab <- cbind.data.frame(df.new, risk_score)
res.cut <- surv_cutpoint(risk_tab, time = 'OS', event = 'Status', variables = 'ACKR1pos')
res.cat <- surv_categorize(res.cut)
df.new$group <- res.cat$ACKR1pos
fit <- survfit(Surv(OS,Status)~group, data = df.new)
pdf('wclc/DC/ACKR1/ACKR1pos.surv.pdf',width = 4,height = 4)
ggsurvplot(fit,data = df.new, conf.int = FALSE,
           risk.table = F,
           palette = c('#E63E40','#367DB6'),
           surv.median.line = "hv",pval = T,
           ggtheme = theme_classic())
dev.off()

# TISMO
ACKR1_DC <- c('Ackr1','Slc16a11','C1qa','Cd1d1','C1qc','Epha8')
tismo <- readRDS('data/TISMO_expressionvivo_profiles.RDS')
lst <- list(ACKR1posDC = ACKR1_DC)
gsva_matrix <- GSVA::gsva(as.matrix(tismo),lst,method = 'ssgsea',kcdf = 'Gaussian', abs.ranking = T)
df <- reshape2::melt(gsva_matrix)
colnames(df) <- c('Celltype','Sample','Score')
tismo.meta <- readRDS('data/TISMO_vivosample_annotations.rds')
tismo.meta$Cancer_type <- as.character(tismo.meta$Cancer_type)
tismo.meta$ICBtreatment <- as.character(tismo.meta$ICBtreatment)
tismo.meta$Responder <- as.character(tismo.meta$Responder)
tismo.meta$Tumor <- as.character(tismo.meta$Tumor)
df$Cancer <- plyr::mapvalues(df$Sample,from = tismo.meta$SampleName,to = tismo.meta$Cancer_type)
df$Treatment <- plyr::mapvalues(df$Sample,from = tismo.meta$SampleName,to = tismo.meta$ICBtreatment)
df$Responder <- plyr::mapvalues(df$Sample,from = tismo.meta$SampleName,to = tismo.meta$Responder)
df$Tumor <- plyr::mapvalues(df$Sample,from = tismo.meta$SampleName,to = tismo.meta$Tumor)

df.2 <- df[which(df$Responder == 'NR' | df$Responder == 'R'),]
df.2$Responder <- as.character(df.2$Responder)
saveRDS(df.2, file = 'wclc/DC/ACKR1/tismo/tismo.dat.rds')


customPairplot <- function(dat, cancer){
  ggpaired(dat, x = 'Responder', y = 'Score', color = 'Responder', palette = 'jco',
           line.size = 0.4, line.color = 'gray') +
    ggtitle(paste0(cancer,' antiPD1 treatment biopsy')) + 
    ylab('ACKR1+ DC Score') +
    stat_compare_means() +
    theme(legend.position = 'none',axis.title.x = element_blank(),
          axis.text = element_text(size = 14,colour = 'black'),
          axis.title.y = element_text(size = 16,colour = 'black'),
          plot.title = element_text(hjust = .5))
}
# crc
crc <- df.2[which(df.2$Cancer == 'Colorectal carcinoma'),]
crc <- crc[crc$Score > 0,]
p <- customPairplot(crc,'Clon Cancer')
ggsave('wclc/DC/ACKR1/tismo/crc.pdf',p,width = 3.7,height = 3)


customBoxplot <- function(dat,cancer,treatment){
  ggboxplot(dat, x = 'Responder', y = 'Score', color = 'Responder',
            palette = 'jco', add = 'jitter', size = 1) +
    stat_compare_means() + ylab('ACKR1+ DC Score') +
    ggtitle(paste(cancer,treatment,'treatment biopsy')) +
    theme(legend.position = 'none',axis.title.x = element_blank(),
          axis.text = element_text(size = 14,colour = 'black'),
          axis.title.y = element_text(size = 16,colour = 'black'),
          plot.title = element_text(hjust = .5))
}

# Lung cancer
lcc <- df.2[which(df.2$Cancer == 'Lung carcinoma'),]
lcc <- lcc[lcc$Score > 0,]
lcc$Responder <- factor(lcc$Responder, levels = c('R','NR'))
p1 <- customBoxplot(lcc,'Lung cancer','antiPD1')
# melanoma
mela <- df.2[which(df.2$Cancer == 'Melanoma'),]
mela <- mela[mela$Score > 0,]
mela$Responder <- factor(mela$Responder, levels = c('R','NR'))
p2 <- customBoxplot(mela,'Melanoma','antiPD1')
# pdac
# pdac <- df.2[which(df.2$Cancer == 'Pancreatic ductal adenocarcinoma'),]
# pdac <- pdac[pdac$Score > 0,]
# customBoxplot(pdac,'PDAC','antiPD1')
# Glioblastoma multiforme
# Mammary carcinoma
# mc <- df.2[which(df.2$Cancer == 'Mammary carcinoma'),]
# mc <- mc[mc$Score > 0,]
# mc$Responder <- factor(mc$Responder, levels = c('R','NR'))
# customBoxplot(mc,'Mammary carcinoma','antiPD1')
# Mammary cancer, NOS
# mcNOS <- df.2[which(df.2$Cancer == 'Mammary cancer, NOS'),]
# mcNOS <- mcNOS[mcNOS$Score > 0,]
# mcNOS$Responder <- factor(mcNOS$Responder, levels = c('R','NR'))
# customBoxplot(mcNOS,'Mammary cancer, NOS','antiPD1')
# Sarcoma
Sarcoma <- df.2[which(df.2$Cancer == 'Sarcoma'),]
Sarcoma <- Sarcoma[Sarcoma$Score > 0,]
Sarcoma$Responder <- factor(Sarcoma$Responder, levels = c('R','NR'))
p3 <- customBoxplot(Sarcoma,'Sarcoma','antiPD1')
# Mammary adenocarcinoma
# ma <- df.2[which(df.2$Cancer == 'Mammary adenocarcinoma'),]
# ma <- ma[ma$Score > 0,]
# customBoxplot(ma,'Mammary adenocarcinoma','antiPD1')
# Hepatocellular carcinoma
hcc <- df.2[which(df.2$Cancer == 'Hepatocellular carcinoma'),]
hcc <- hcc[hcc$Score > 0,]
hcc$Responder <- factor(hcc$Responder, levels = c('R','NR'))
p4 <- customBoxplot(hcc,'Liver cancer','antiPD1')
# Mesothelioma
# mesothelioma <- df.2[which(df.2$Cancer == 'Mesothelioma'),]
# mesothelioma <- mesothelioma[mesothelioma$Score > 0,]
# mesothelioma$Responder <- factor(mesothelioma$Responder, levels = c('R','NR'))
# customBoxplot(mesothelioma,'Mesothelioma','antiPD1')
# Renal adenocarcinoma
ra <- df.2[which(df.2$Cancer == 'Renal adenocarcinoma'),]
ra <- ra[ra$Score > 0,]
ra$Responder <- factor(ra$Responder, levels = c('R','NR'))
p5 <- customBoxplot(ra,'Renal adenocarcinoma','antiPD1')
# Head and neck squamous cell carcinoma
# hnscc <- df.2[which(df.2$Cancer == 'Head and neck squamous cell carcinoma'),]
# hnscc <- hnscc[hnscc$Score > 0,]
# customBoxplot(hnscc,'HNSCC','antiPD1')
# Gastric adenocarcinoma
gc <- df.2[which(df.2$Cancer == 'Gastric adenocarcinoma'),]
gc <- gc[gc$Score > 0,]
gc$Responder <- factor(gc$Responder, levels = c('R','NR'))
p6 <- customBoxplot(gc,'Gastric cancer','antiPD1')

p.list <- list(p1,p2,p3,p4,p5,p6)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 3)
ggsave('wclc/DC/ACKR1/tismo/pancancer2.pdf',p,width = 7,height = 6)
