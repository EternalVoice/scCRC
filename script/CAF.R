# ----------------- CAF ---------------------
rm(list = ls());gc()

caf <- readRDS('wclc/rds/caf.tumor.rds')

caf <- NormalizeData(object = caf, normalization.method = "LogNormalize", scale.factor = 10000)
caf <- FindVariableFeatures(caf,selection.method="vst",nfeatures=2000)
caf <- ScaleData(caf, features = rownames(caf))
caf <- RunPCA(object = caf, features = VariableFeatures(caf))
# Find neighbors and clusters with harmony batch correction
caf <- FindNeighbors(caf, dims = 1:30, reduction = "pca")
caf <- FindClusters(caf, resolution = 0.2)
# myeloid <- myeloid %>% RunHarmony("SampleName",plot_convergence = F)
caf <- RunUMAP(caf, dims = 1:30, reduction = "pca")
caf <- RunTSNE(caf,dims = 1:30)

old.ident <- names(table(caf$cl295v11SubFull))
new.ident <- c('CXCL14+CAF','GREM1+CAF','MMP3+CAF','CCL8+CAF')
caf$celltype <- plyr::mapvalues(caf$cl295v11SubFull, from = old.ident, to = new.ident)

saveRDS(caf, file = 'wclc/rds/caf.tumor.rds')

# vis
mypal <- c('#1F78B4','#E31A1C','#FF7F00','#FB9A99')
p <- DimPlot(caf, group.by = 'celltype', cols = mypal,pt.size = 1) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/CAF/celltype.pdf',p,width = 6,height = 6.2)

lcc.caf <- subset(caf, cells = rownames(subset(caf@meta.data, TissueSiteSimple == 'left')))
p <- DimPlot(lcc.caf, group.by = 'celltype', cols = mypal,pt.size = 1) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/CAF/lcc.celltype.pdf',p,width = 6,height = 6)

rcc.caf <- subset(caf, cells = rownames(subset(caf@meta.data, TissueSiteSimple == 'right')))
p <- DimPlot(rcc.caf, group.by = 'celltype', cols = mypal,pt.size = 1) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/CAF/rcc.celltype.pdf',p,width = 6,height = 6)

mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "black"),
                 panel.grid = element_blank(),
                 axis.text = element_text(size = 12,face = "bold",colour = "black"),
                 axis.title.y = element_blank(),
                 axis.title = element_text(size = 14,face = "bold",colour = "black"),
                 legend.position = 'bottom')

p <- ggplot(caf@meta.data, aes(x = TissueSiteSimple, fill = celltype)) +
  geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values = mypal) +
  coord_flip() +
  mytheme +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(filename = 'wclc/CAF/cell.prop.pdf',p,width = 8,height = 4)

# -- DE analysis
library(ggrepel)
caf.markers <- FindMarkers(caf, group.by = 'TissueSiteSimple',ident.1 = 'right',ident.2 = 'left')
# remove ribosomal and mitochondrial genes
caf.markers <- caf.markers[!grepl("^RP[SL]",rownames(caf.markers)),]
caf.markers <- caf.markers[!grepl("^MT-",rownames(caf.markers)),]
caf.markers$Significance <- ifelse(caf.markers$p_val < 0.01, TRUE, FALSE)
caf.markers <- caf.markers[order(caf.markers[,2],decreasing = T),]
top10 <- rbind(head(caf.markers,10),tail(caf.markers,10))

mytheme <- theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 16),
                 axis.text = element_text(size = 12),
                 legend.title = element_text(size = 12),
                 legend.text = element_text(size = 12)) 

p.all.markers <- ggplot(caf.markers, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')
ggsave(filename = 'wclc/CAF/DEgenes.rcc.lcc.volc.pdf',p.all.markers,width = 5,height = 5)

up.genes <- rownames(caf.markers[caf.markers$avg_logFC > 0,])
dn.genes <- rownames(caf.markers[caf.markers$avg_logFC < 0,])

save(caf.markers, up.genes, dn.genes, file = 'wclc/CAF/all.caf.DEgenes.lcc.rcc.rda')

caf.up <- EnrichmentAnalysis(DEgenes = up.genes, enrich.method = 'GO')
openxlsx::write.xlsx(caf.up$ego_BP@result, file = 'wclc/CAF/RCC.up.GO.xlsx')
caf.dn <- EnrichmentAnalysis(DEgenes = dn.genes, enrich.method = 'GO')
openxlsx::write.xlsx(caf.dn$ego_BP@result, file = 'wclc/CAF/LCC.up.GO.xlsx')

## Visualization of DE pathway over LCC and RCC
DEpathway <- xlsx::read.xlsx("wclc/CAF/DEpathway.sel.rcc.vs.lcc.xlsx",sheetIndex = 1)
DEpathway$value <- -log10(DEpathway$pvalue)
DEpathway$value[DEpathway$Group=="LCC.CAF"] <- -DEpathway$value[DEpathway$Group=="LCC.CAF"]
DEpathway$hjust <- ifelse(DEpathway$Group == "RCC.CAF",1,0)
p <- ggplot(data = DEpathway,aes(factor(Description,levels=Description),value,color=Group,fill=Group)) +
  geom_bar(stat = "identity", width = 0.5) + theme_bw() + coord_flip() +
  scale_color_manual(values = c("RCC.CAF"="#DB4740","LCC.CAF"="#4E9BD2")) +
  scale_fill_manual(values = c("RCC.CAF"="#DB4740","LCC.CAF"="#4E9BD2")) +
  geom_text(data = DEpathway, mapping = aes(label=Description, y=0), 
            hjust=DEpathway$hjust, size = 5, color = "black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(DEpathway$value))),ceiling(max(abs(DEpathway$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(filename = "wclc/CAF/DEpathway.LCC.vs.RCC.pdf",p,width = 6,height = 6)

## violin plot of interested genes

# granulocyte migration: PPBP, PF4V1, CCL7, CXCL8, CXCL3
# immune response: IGFBP2, IGFBP1, LCN2, ALDH1A1, CASP1
mycol <- c("#4E9BD2","#DB4740")
mytheme <- theme(axis.title.x = element_blank(),
                 axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0),
                 panel.border = element_rect(size = 1,colour = "black"))

p <- VlnPlot(CAF,features = "PPBP",group.by = "TissueSiteSimple",
             y.max = 5.5,cols = mycol,pt.size = 1e-5) + 
  NoLegend() + mytheme +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/CAF/features/PPBP.pdf",p,width = 3,height = 4)
p <- VlnPlot(CAF,features = "PF4V1",group.by = "TissueSiteSimple",y.max = 7,
             pt.size = 1e-5,cols = c("#4E9BD2","#DB4740")) + 
  NoLegend() + mytheme +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/CAF/features/PF4V1.pdf",p,width = 3,height = 4)
p <- VlnPlot(CAF,features = "CCL7",group.by = "TissueSiteSimple",y.max = 6,
             cols = c("#4E9BD2","#DB4740"),pt.size = 1e-5) + 
  NoLegend() + mytheme +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/CAF/features/CCL7.pdf",p,width = 3,height = 4)
p <- VlnPlot(CAF,features = "CXCL8",group.by = "TissueSiteSimple",y.max = 5,
             cols = c("#4E9BD2","#DB4740"),pt.size = 1e-5) + 
  NoLegend() + mytheme +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/CAF/features/CXCL8.pdf",p,width = 3,height = 4)
p <- VlnPlot(CAF,features = "CXCL3",group.by = "TissueSiteSimple",y.max = 5,
             cols = c("#4E9BD2","#DB4740"),pt.size = 1e-5) + 
  NoLegend() + mytheme +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/CAF/features/CXCL3.pdf",p,width = 3,height = 4)

# IGFBP2, IGFBP1, LCN2, ALDH1A1, CASP1
p <- VlnPlot(CAF,features = "IGFBP2",group.by = "TissueSiteSimple",y.max = 7,
             cols = c("#4E9BD2","#DB4740"),pt.size = 1e-5) + 
  NoLegend() + mytheme +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/CAF/features/IGFBP2.pdf",p,width = 3,height = 4)
p <- VlnPlot(CAF,features = "IGFBP1",group.by = "TissueSiteSimple",y.max = 6,
             cols = c("#4E9BD2","#DB4740"),pt.size = 1e-5) + 
  NoLegend() + mytheme +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/CAF/features/IGFBP1.pdf",p,width = 3,height = 4)
p <- VlnPlot(CAF,features = "LCN2",group.by = "TissueSiteSimple",y.max = 4,
             cols = c("#4E9BD2","#DB4740"),pt.size = 1e-5) + 
  NoLegend() + mytheme +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/CAF/features/LCN2.pdf",p,width = 3,height = 4)
p <- VlnPlot(CAF,features = "ALDH1A1",group.by = "TissueSiteSimple",y.max = 4.5,
             cols = c("#4E9BD2","#DB4740"),pt.size = 1e-5) + 
  NoLegend() + mytheme +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/CAF/features/ALDH1A1.pdf",p,width = 3,height = 4)

p <- VlnPlot(CAF,features = "CASP1",group.by = "TissueSiteSimple",
        y.max = 3.5,cols = c("#4E9BD2","#DB4740"),pt.size = 1e-5) + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0),
        panel.border = element_rect(size = 1,colour = "black")) +
  ggpubr::stat_compare_means(comparisons = list(c("left","right")))
ggsave(filename = "wclc/CAF/features/CASP1.pdf",p,width = 3,height = 4)

Idents(caf) <- caf$celltype
# CCL8+CAF
ccl8.caf.markers <- FindMarkers(caf, group.by = 'TissueSiteSimple',ident.1 = 'right',ident.2 = 'left',subset.ident = 'CCL8+CAF')
# remove ribosomal and mitochondrial genes
ccl8.caf.markers <- ccl8.caf.markers[!grepl("^RP[SL]",rownames(ccl8.caf.markers)),]
ccl8.caf.markers <- ccl8.caf.markers[!grepl("^MT-",rownames(ccl8.caf.markers)),]
ccl8.caf.markers$Significance <- ifelse(ccl8.caf.markers$p_val < 0.01, TRUE, FALSE)
ccl8.caf.markers <- ccl8.caf.markers[order(ccl8.caf.markers[,2],decreasing = T),]
top10 <- rbind(head(ccl8.caf.markers,10),tail(ccl8.caf.markers,10))
p.all.markers <- ggplot(ccl8.caf.markers, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')
ggsave(filename = 'wclc/CAF/ccl8.caf.rcc.lcc.volc.pdf',p.all.markers,width = 5,height = 5)

up.genes <- rownames(ccl8.caf.markers[ccl8.caf.markers$avg_logFC > 0,])
dn.genes <- rownames(ccl8.caf.markers[ccl8.caf.markers$avg_logFC < 0,])

save(ccl8.caf.markers, up.genes, dn.genes, file = 'wclc/CAF/ccl8.caf.DEgenes.lcc.rcc.rda')

caf.up <- EnrichmentAnalysis(DEgenes = up.genes, enrich.method = 'GO')
openxlsx::write.xlsx(caf.up$ego_BP@result, file = 'wclc/CAF/CCL8.CAF.RCC.up.GO.xlsx')
caf.dn <- EnrichmentAnalysis(DEgenes = dn.genes, enrich.method = 'GO')
openxlsx::write.xlsx(caf.dn$ego_BP@result, file = 'wclc/CAF/CCL8.CAF.LCC.up.GO.xlsx')

# p <- barplot(caf.up$ego_BP, showCategory = 20, title = "TOP 20 enriched GO terms") +
#   scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")))
# ggsave(filename = 'wclc/CAF/CCL8.CAF.RCC.up.GO.pdf',p,width = 8,height = 10)

# CXCL14+CAF
cxcl14.caf.markers <- FindMarkers(caf, group.by = 'TissueSiteSimple',ident.1 = 'right',ident.2 = 'left',subset.ident = 'CXCL14+CAF')
# remove ribosomal and mitochondrial genes
cxcl14.caf.markers <- cxcl14.caf.markers[!grepl("^RP[SL]",rownames(cxcl14.caf.markers)),]
cxcl14.caf.markers <- cxcl14.caf.markers[!grepl("^MT-",rownames(cxcl14.caf.markers)),]
cxcl14.caf.markers$Significance <- ifelse(cxcl14.caf.markers$p_val < 0.01, TRUE, FALSE)
cxcl14.caf.markers <- cxcl14.caf.markers[order(cxcl14.caf.markers[,2],decreasing = T),]
top10 <- rbind(head(cxcl14.caf.markers,10),tail(cxcl14.caf.markers,10))
p.all.markers <- ggplot(cxcl14.caf.markers, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')
ggsave(filename = 'wclc/CAF/cxcl14.caf.rcc.lcc.volc.pdf',p.all.markers,width = 5,height = 5)

up.genes <- rownames(cxcl14.caf.markers[cxcl14.caf.markers$avg_logFC > 0,])
dn.genes <- rownames(cxcl14.caf.markers[cxcl14.caf.markers$avg_logFC < 0,])

save(cxcl14.caf.markers, up.genes, dn.genes, file = 'wclc/CAF/cxcl14.caf.DEgenes.lcc.rcc.rda')

caf.up <- EnrichmentAnalysis(DEgenes = up.genes, enrich.method = 'GO')
openxlsx::write.xlsx(caf.up$ego_BP@result, file = 'wclc/CAF/CXCL14.CAF.RCC.up.GO.xlsx')
caf.dn <- EnrichmentAnalysis(DEgenes = dn.genes, enrich.method = 'GO')
openxlsx::write.xlsx(caf.dn$ego_BP@result, file = 'wclc/CAF/CXCL14.CAF.LCC.up.GO.xlsx')

# GREM1+CAF
grem1.caf.markers <- FindMarkers(caf, group.by = 'TissueSiteSimple',ident.1 = 'right',ident.2 = 'left',subset.ident = 'GREM1+CAF')
# remove ribosomal and mitochondrial genes
grem1.caf.markers <- grem1.caf.markers[!grepl("^RP[SL]",rownames(grem1.caf.markers)),]
grem1.caf.markers <- grem1.caf.markers[!grepl("^MT-",rownames(grem1.caf.markers)),]
grem1.caf.markers$Significance <- ifelse(grem1.caf.markers$p_val < 0.01, TRUE, FALSE)
grem1.caf.markers <- grem1.caf.markers[order(grem1.caf.markers[,2],decreasing = T),]
top10 <- rbind(head(grem1.caf.markers,10),tail(grem1.caf.markers,10))
p.all.markers <- ggplot(grem1.caf.markers, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')
ggsave(filename = 'wclc/CAF/grem1.caf.rcc.lcc.volc.pdf',p.all.markers,width = 5,height = 5)

up.genes <- rownames(grem1.caf.markers[grem1.caf.markers$avg_logFC > 0,])
dn.genes <- rownames(grem1.caf.markers[grem1.caf.markers$avg_logFC < 0,])

save(grem1.caf.markers, up.genes, dn.genes, file = 'wclc/CAF/grem1.caf.DEgenes.lcc.rcc.rda')

caf.up <- EnrichmentAnalysis(DEgenes = up.genes, enrich.method = 'GO')
openxlsx::write.xlsx(caf.up$ego_BP@result, file = 'wclc/CAF/GREM1.CAF.RCC.up.GO.xlsx')
caf.dn <- EnrichmentAnalysis(DEgenes = dn.genes, enrich.method = 'GO')
openxlsx::write.xlsx(caf.dn$ego_BP@result, file = 'wclc/CAF/GREM1.CAF.LCC.up.GO.xlsx')

# MMP3+CAF
mmp3.caf.markers <- FindMarkers(caf, group.by = 'TissueSiteSimple',ident.1 = 'right',ident.2 = 'left',subset.ident = 'MMP3+CAF')
# remove ribosomal and mitochondrial genes
mmp3.caf.markers <- mmp3.caf.markers[!grepl("^RP[SL]",rownames(mmp3.caf.markers)),]
mmp3.caf.markers <- mmp3.caf.markers[!grepl("^MT-",rownames(mmp3.caf.markers)),]
mmp3.caf.markers$Significance <- ifelse(mmp3.caf.markers$p_val < 0.01, TRUE, FALSE)
mmp3.caf.markers <- mmp3.caf.markers[order(mmp3.caf.markers[,2],decreasing = T),]
top10 <- rbind(head(mmp3.caf.markers,10),tail(mmp3.caf.markers,10))
p.all.markers <- ggplot(mmp3.caf.markers, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')
ggsave(filename = 'wclc/CAF/mmp3.caf.rcc.lcc.volc.pdf',p.all.markers,width = 5,height = 5)

up.genes <- rownames(mmp3.caf.markers[mmp3.caf.markers$avg_logFC > 0,])
dn.genes <- rownames(mmp3.caf.markers[mmp3.caf.markers$avg_logFC < 0,])

save(mmp3.caf.markers, up.genes, dn.genes, file = 'wclc/CAF/mmp3.caf.DEgenes.lcc.rcc.rda')

caf.up <- EnrichmentAnalysis(DEgenes = up.genes, enrich.method = 'GO')
openxlsx::write.xlsx(caf.up$ego_BP@result, file = 'wclc/CAF/MMP3.CAF.RCC.up.GO.xlsx')
caf.dn <- EnrichmentAnalysis(DEgenes = dn.genes, enrich.method = 'GO')
openxlsx::write.xlsx(caf.dn$ego_BP@result, file = 'wclc/CAF/MMP3.CAF.LCC.up.GO.xlsx')


load("wclc/CAF/ccl8.caf.DEgenes.lcc.rcc.rda")
load("wclc/CAF/cxcl14.caf.DEgenes.lcc.rcc.rda")
load("wclc/CAF/grem1.caf.DEgenes.lcc.rcc.rda")
load("wclc/CAF/mmp3.caf.DEgenes.lcc.rcc.rda")

ccl8.caf.markers$celltype <- "CCL8+CAF"
ccl8.caf.markers$gene <- c(rownames(head(ccl8.caf.markers,10)),
                           rep("",nrow(ccl8.caf.markers)-20),
                           rownames(tail(ccl8.caf.markers,10)))
cxcl14.caf.markers$celltype <- "CXCL14+CAF"
cxcl14.caf.markers$gene <- c(rownames(head(cxcl14.caf.markers,10)),
                             rep("",nrow(cxcl14.caf.markers)-20),
                             rownames(tail(cxcl14.caf.markers,10)))
grem1.caf.markers$celltype <- "GREM1"
grem1.caf.markers$gene <- c(rownames(head(grem1.caf.markers,10)),
                            rep("",nrow(grem1.caf.markers)-20),
                            rownames(tail(grem1.caf.markers,10)))
mmp3.caf.markers$celltype <- "MMP3+CAF"
mmp3.caf.markers$gene <- c(rownames(head(mmp3.caf.markers,10)),
                            rep("",nrow(mmp3.caf.markers)-20),
                            rownames(tail(mmp3.caf.markers,10)))

all.markers <- rbind(ccl8.caf.markers,cxcl14.caf.markers,grem1.caf.markers,mmp3.caf.markers)

all.markers$bar.1 <- 0
all.markers$bar.2 <- 0
for(i in 1:length(unique(all.markers$celltype))){
  tmp <- which(all.markers$celltype==unique(all.markers$celltype)[i])[1]
  all.markers$bar.1[tmp] <- 0.2
  all.markers$bar.2[tmp] <- -0.2
}
all.markers$celltype <- factor(all.markers$celltype,levels = c("CCL8+CAF","CXCL14+CAF","GREM1","MMP3+CAF"))
p <- ggplot(all.markers, aes(celltype,avg_logFC, color=Significance)) +
  geom_point(aes(col = Significance)) +
  geom_jitter() + scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(aes(celltype, avg_logFC, label = gene),color="black") +
  scale_y_continuous(limits = c(-2,2)) +
  theme_classic() +
  ylab("Average log2FC") +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9,0.9),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size=5)))

p2 <- p + geom_bar(aes(celltype,bar.1,fill=celltype),stat = 'identity',color=NA) +
  geom_bar(aes(celltype,bar.2,fill=celltype),stat = 'identity',color=NA) +
  scale_fill_manual(values = c('#1F78B4','#E31A1C','#FF7F00','#FB9A99')) +
  theme(legend.position = "none")

ggsave(filename = "wclc/CAF/all.DEgenes.pdf",p2,width = 8,height = 6)

pathway <- xlsx::read.xlsx("wclc/CAF/DEpathway.top5.xlsx",sheetIndex = 1)

p <- ggplot(pathway, aes(x=Group,y=Description)) + 
  geom_point(aes(size=GeneRatio,color = -log10(pvalue))) +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_blank()) +
  scale_color_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(10,"RdBu"))(10)))
ggsave(filename = "wclc/CAF/top5.pathway.pdf",p,width = 20,height = 6)

### Evolution trajectory of CAFs
rm(list = ls());gc()
library(Seurat)
library(monocle)

CAF <- readRDS("wclc/rds/caf.tumor.rds")
lcc.CAF <- subset(CAF, cells = rownames(subset(CAF@meta.data, TissueSiteSimple == "left")))
rcc.CAF <- subset(CAF, cells = rownames(subset(CAF@meta.data, TissueSiteSimple == "right")))

mypal <- c('#1F78B4','#E31A1C','#FF7F00','#FB9A99')

table(lcc.CAF$celltype)
lcc.dat <- as(as.matrix(lcc.CAF@assays$RNA@counts), "sparseMatrix")
pd <- new('AnnotatedDataFrame', data = lcc.CAF@meta.data)
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
lcc.mycds <- orderCells(lcc.mycds,root_state = 8)
saveRDS(lcc.mycds, file = "wclc/rds/LCC.CAF.monocle2.rds")

source("script/pseudotime_heatmap.R")
p1 <- plot_cell_trajectory(lcc.mycds, color_by = "celltype") + 
  scale_color_manual(values = mypal) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("wclc/CAF/trajectory/lcc.CAF.trajectory.pdf",p1,width = 5,height = 4)

p2 <- plot_cell_trajectory(lcc.mycds, color_by = "Pseudotime") +
  scale_color_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256))) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("wclc/CAF/trajectory/lcc.CAF.pseudo.pdf",p2,width = 5,height = 4)

# RCC
table(rcc.CAF$celltype)
rcc.dat <- as(as.matrix(rcc.CAF@assays$RNA@counts), "sparseMatrix")
pd <- new('AnnotatedDataFrame', data = rcc.CAF@meta.data)
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
rcc.mycds <- orderCells(rcc.mycds)
saveRDS(rcc.mycds, file = "wclc/rds/RCC.CAF.monocle2.rds")

p1 <- plot_cell_trajectory(rcc.mycds, color_by = "celltype") + 
  scale_color_manual(values = mypal) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("wclc/CAF/trajectory/rcc.CAF.trajectory.pdf",p1,width = 5,height = 4)

p2 <- plot_cell_trajectory(rcc.mycds, color_by = "Pseudotime") +
  scale_color_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256))) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey45", fill = NA),
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
ggsave("wclc/CAF/trajectory/rcc.CAF.pseudo.pdf",p2,width = 5,height = 4)


## Tissue preference
rm(list = ls());gc()
CAF <- readRDS("wclc/rds/caf.tumor.rds")

chiq.L <- chisq.test(table(CAF$celltype,CAF$TissueSiteSimple)[,1])
chiq.R <- chisq.test(table(CAF$celltype,CAF$TissueSiteSimple)[,2])
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
ggsave("wclc/CAF/Roe.pdf",p,width = 6,height = 8)

## Pan-cancer and survival

# SIGNATURES
# CCL8+CAF: ADAMDEC1, APOE, CCL11, CCL13, TCF21
# CXCL14+CAF:CXCL14, AGT, NSG1, MEST, EMID1, CST1, BMP4, WNT4, INHBA
# GREM1+CAF: COL10A1, GAS1, RSPO3, COL11A1, FAP, INHBA
# MMP3+CAF: MMP10, CCL20, IL1B, CSF2, STC1, INHBA
rm(list = ls());gc()
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)

CCL8 <- c('ADAMDEC1','APOE','CCL11','CCL13','TCF21')
CXCL14 <- c('CXCL14','AGT','NSG1','MEST','EMID1','CST1','BMP4','WNT4','INHBA')
CXCL14 <- gsub("NSG1","D4S234E",CXCL14)
GREM1 <- c('COL10A1','GAS1','RSPO3','COL11A1','FAP','INHBA')
MMP3 <- c('MMP10','CCL20','IL1B','CSF2','STC1','INHBA')
## Pancacer
load("data/TCGA/TCGA_36cancers.rnaseq.RData")
sig <- list(CCL8,CXCL14,GREM1,MMP3)
cancers <- gsub(".rnaseq",'',names(TCGA_36cancers.rnaseq))

df <- as.data.frame(matrix(nrow = 4, ncol = 0))
for(i in 1:length(TCGA_36cancers.rnaseq)){
  expr <- TCGA_36cancers.rnaseq[[i]]
  rownames(expr) <- expr$bcr_patient_barcode
  expr$bcr_patient_barcode <- NULL
  expr <- log2(t(expr) + 1)
  expr <- expr[, order(colnames(expr))]
  gsva_matrix <- gsva(as.matrix(expr), sig, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = T)
  
  x1 <- sum(gsva_matrix[1,])/ncol(gsva_matrix)
  x2 <- sum(gsva_matrix[2,])/ncol(gsva_matrix)
  x3 <- sum(gsva_matrix[3,])/ncol(gsva_matrix)
  x4 <- sum(gsva_matrix[4,])/ncol(gsva_matrix)
  x <- c(x1,x2,x3,x4)
  names(x) <- rownames(gsva_matrix)
  x <- data.frame(x)
  colnames(x) <- cancers[i]
  df <- cbind(df, x)
}

df$celltype <- c("CCL8+CAF","CXCL14+CAF","GREM1+CAF","MMP3+CAF")
df.2 <- reshape2::melt(df)
colnames(df.2) <- c("celltype","cancers","score")
p <- ggplot(df.2,aes(x=cancers,y=celltype,fill=score,color=score)) +
  geom_point(size=5) + theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 12,colour = "black")) +
  scale_color_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(7, "Reds"))(100))
ggsave("wclc/CAF/CAF.sub.score.in.pancancer.pdf",p,width = 12,height = 2.5)
openxlsx::write.xlsx(df.2,file = "wclc/CAF/CAF.sub.score.xlsx")

## LCC vs RCC in COAD
tcga.exp.df <- read.table("data/TCGA/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",header = T)
tcga.exp.df <- tcga.exp.df[!duplicated(tcga.exp.df$Hugo_Symbol),]
rownames(tcga.exp.df) <- tcga.exp.df$Hugo_Symbol
tcga.exp.df <- tcga.exp.df[,-c(1,2)]
colnames(tcga.exp.df) <- gsub("\\.","-",colnames(tcga.exp.df))

tcga.exp.df <- data.frame(t(tcga.exp.df))
colnames(tcga.exp.df) <- gsub("\\.","-",colnames(tcga.exp.df))
Type <- ifelse(substr(rownames(tcga.exp.df),14,15) == "11", "Normal", "Tumor")
tcga.exp.df <- data.frame(tibble::add_column(tcga.exp.df, Type, .before = 2))
tcga.exp.df <- tcga.exp.df[tcga.exp.df$Type=="Tumor",]
clin <- read.table("data/TCGA/data_clinical_patient.txt",header = T,sep = "\t",skip = 4)
clin <- data.frame(clin$PATIENT_ID, clin$OS_MONTHS, clin$OS_STATUS, clin$PRIMARY_SITE_PATIENT)
clin <- clin[clin$clin.OS_STATUS != "[Not Available]",]
clin <- clin[clin$clin.PRIMARY_SITE_PATIENT != "[Not Available]",]
clin <- clin[clin$clin.PRIMARY_SITE_PATIENT != "[Discrepancy]",]
clin$clin.PRIMARY_SITE_PATIENT <- as.character(clin$clin.PRIMARY_SITE_PATIENT)

colon <- data.frame(TissueSite = names(table(clin$clin.PRIMARY_SITE_PATIENT)),
                    TissueSiteSimple = c('right','right','left','right','left','left',
                                         'left','left','left'),
                    stringsAsFactors = F)
clin$TissueSiteSimple <- plyr::mapvalues(
  x = clin$clin.PRIMARY_SITE_PATIENT, 
  from = colon$TissueSite, 
  to = colon$TissueSiteSimple
)
tcga.exp.df <- t(tcga.exp.df)
tcga.exp.df <- tcga.exp.df[, order(colnames(tcga.exp.df))]
gsva_matrix <- gsva(as.matrix(tcga.exp.df), sig, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = T)
gsva_matrix <- data.frame(t(gsva_matrix))
gsva_matrix$PATIENT_ID <- rownames(gsva_matrix)
gsva_matrix$PATIENT_ID <- gsub("-01",'',gsva_matrix$PATIENT_ID)
colnames(gsva_matrix)[1:4] <- c("CCL8+CAF","CXCL14+CAF","GREM1+CAF","MMP3+CAF")
merged.data <- merge(clin,gsva_matrix,by.x="clin.PATIENT_ID",by.y="PATIENT_ID")
merged.data <- merged.data[order(merged.data$TissueSiteSimple,decreasing = F),]
xlsx::write.xlsx(merged.data,file = "wclc/CAF/CAF.sub.score.xlsx",sheetName = "COAD",append = T,row.names = F)

library(ggpubr)
p.list <- list()
for(i in 6:ncol(merged.data)){
  p.list[[i]] <- ggboxplot(merged.data,x = colnames(merged.data)[5], 
                           y = colnames(merged.data)[i],
                           color = colnames(merged.data)[5],
                           palette = c('#4E9BD2','#DE4343'),
                           size = 1.5, add = 'jitter') + 
    ggtitle(colnames(merged.data)[i]) +
    ylab("Enrichment Score") + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 14,colour = "black"),
          axis.text = element_text(size = 13,colour = "black")) +
    stat_compare_means(comparisons = list(c("left","right")))
  ggsave(paste0("wclc/CAF/",gsub("\\+","_",colnames(merged.data)[i]),".left.vs.right.pdf"),
         p.list[[i]],width = 3.5,height = 3.6)
}

pc <- cowplot::plot_grid(plotlist = p.list,nrow = 2)
ggsave(filename = "wclc/CAF/CAF.sub.COAD.score.pdf",pc,width = 6,height = 6)


## Survival
clin <- data.frame(clin$PATIENT_ID, clin$OS_MONTHS, clin$OS_STATUS)
clin <- clin[clin$clin.OS_STATUS != "[Not Available]",]
Dead <- rep(NA, dim(clin)[1])
for(i in 1:length(clin$clin.OS_STATUS)){
  if(stringr::str_split(clin$clin.OS_STATUS,":",simplify = T)[,2][i]=="DECEASED"){
    Dead[i] <- 1
  }else{
    Dead[i] <- 0
  }
}
clin$Dead <- Dead
clin$clin.OS_MONTHS <- as.integer(clin$clin.OS_MONTHS)
## CCL8
new.tcga.exp <- data.frame(tcga.exp.df[,sig[[1]]])
removeColAllNA <- function(x){x[,apply(x,2,function(y) any(!is.na(y)))]}
new.tcga.exp <- removeColAllNA(new.tcga.exp)
new.tcga.exp <- data.frame(scale(new.tcga.exp))
new.tcga.exp$PATIENT_ID <- gsub("-01","",rownames(new.tcga.exp))
avg.exp <- aggregate(.~PATIENT_ID,FUN = mean, data = new.tcga.exp)
merged.data <- merge(clin,avg.exp,by.x="clin.PATIENT_ID",by.y="PATIENT_ID")
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS,Dead)~.,
                          data = merged.data[,c(2,5:ncol(merged.data))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix(merged.data[,6:ncol(merged.data)] - colMeans(merged.data[,6:ncol(merged.data)])) %*% coefs
group <- rep(NA,dim(merged.data)[1])
for(i in 1:dim(merged.data)[1]){
  if(risk.score[i] > 0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data$group <- group

mytheme <- theme(legend.title = element_text(size = 14), 
                 legend.text = element_text(size = 14),
                 axis.text.x = element_text(size = 14),
                 axis.text.y = element_text(size = 14),
                 axis.line.x = element_line(size = 1),
                 axis.line.y = element_line(size = 1),
                 axis.title.x = element_text(size = 16),
                 axis.title.y = element_text(size = 16))

mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead)~group,data=merged.data)
p.Surv <- survminer::ggsurvplot(mysurv,conf.int = T,risk.table = T)
p1 <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 40, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 358"), x = 20, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "CCL8+CAF\nsignatures")) +
  mytheme + 
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave(filename = "wclc/CAF/CCL8_CAF.surv.pdf",p1,width = 6,height = 6)
openxlsx::write.xlsx(merged.data,file = "wclc/CAF/CAF.sub.surv.xlsx")

## CXCL14
sig[[2]] <- gsub("D4S234E","NSG1",sig[[2]])
new.tcga.exp <- data.frame(tcga.exp.df[,sig[[2]]])
removeColAllNA <- function(x){x[,apply(x,2,function(y) any(!is.na(y)))]}
new.tcga.exp <- removeColAllNA(new.tcga.exp)
new.tcga.exp <- data.frame(scale(new.tcga.exp))
new.tcga.exp$PATIENT_ID <- gsub("-01","",rownames(new.tcga.exp))
avg.exp <- aggregate(.~PATIENT_ID,FUN = mean, data = new.tcga.exp)
merged.data <- merge(clin,avg.exp,by.x="clin.PATIENT_ID",by.y="PATIENT_ID")
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS,Dead)~.,
                          data = merged.data[,c(2,5:ncol(merged.data))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix(merged.data[,6:ncol(merged.data)] - colMeans(merged.data[,6:ncol(merged.data)])) %*% coefs
group <- rep(NA,dim(merged.data)[1])
for(i in 1:dim(merged.data)[1]){
  if(risk.score[i] > 0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data$group <- group

mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead)~group,data=merged.data)
p.Surv <- survminer::ggsurvplot(mysurv,conf.int = T,risk.table = T)
p2 <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 40, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 358"), x = 20, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "CXCL14+CAF\nsignatures")) +
  mytheme + 
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave(filename = "wclc/CAF/CXCL14_CAF.surv.pdf",p2,width = 6,height = 6)
xlsx::write.xlsx(merged.data,file = "wclc/CAF/CAF.sub.surv.xlsx",sheetName = "CXCL14+CAF",append = T)

## GREM1
new.tcga.exp <- data.frame(tcga.exp.df[,sig[[3]]])
removeColAllNA <- function(x){x[,apply(x,2,function(y) any(!is.na(y)))]}
new.tcga.exp <- removeColAllNA(new.tcga.exp)
new.tcga.exp <- data.frame(scale(new.tcga.exp))
new.tcga.exp$PATIENT_ID <- gsub("-01","",rownames(new.tcga.exp))
avg.exp <- aggregate(.~PATIENT_ID,FUN = mean, data = new.tcga.exp)
merged.data <- merge(clin,avg.exp,by.x="clin.PATIENT_ID",by.y="PATIENT_ID")
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS,Dead)~.,
                          data = merged.data[,c(2,5:ncol(merged.data))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix(merged.data[,6:ncol(merged.data)] - colMeans(merged.data[,6:ncol(merged.data)])) %*% coefs
group <- rep(NA,dim(merged.data)[1])
for(i in 1:dim(merged.data)[1]){
  if(risk.score[i] > 0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data$group <- group
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead)~group,data=merged.data)
p.Surv <- survminer::ggsurvplot(mysurv,conf.int = T,risk.table = T)
p3 <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 40, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 358"), x = 20, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "GREM1+CAF\nsignatures")) +
  mytheme + 
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave(filename = "wclc/CAF/GREM1_CAF.surv.pdf",p3,width = 6,height = 6)
xlsx::write.xlsx(merged.data,file = "wclc/CAF/CAF.sub.surv.xlsx",sheetName = "GREM1+CAF",append = T)
## MMP3
new.tcga.exp <- data.frame(tcga.exp.df[,sig[[4]]])
removeColAllNA <- function(x){x[,apply(x,2,function(y) any(!is.na(y)))]}
new.tcga.exp <- removeColAllNA(new.tcga.exp)
new.tcga.exp <- data.frame(scale(new.tcga.exp))
new.tcga.exp$PATIENT_ID <- gsub("-01","",rownames(new.tcga.exp))
avg.exp <- aggregate(.~PATIENT_ID,FUN = mean, data = new.tcga.exp)
merged.data <- merge(clin,avg.exp,by.x="clin.PATIENT_ID",by.y="PATIENT_ID")
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS,Dead)~.,
                          data = merged.data[,c(2,5:ncol(merged.data))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix(merged.data[,6:ncol(merged.data)] - colMeans(merged.data[,6:ncol(merged.data)])) %*% coefs
group <- rep(NA,dim(merged.data)[1])
for(i in 1:dim(merged.data)[1]){
  if(risk.score[i] > 0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data$group <- group
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead)~group,data=merged.data)
p.Surv <- survminer::ggsurvplot(mysurv,conf.int = T,risk.table = T)
p4 <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 40, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 358"), x = 20, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "MMP3+CAF\nsignatures")) +
  mytheme + 
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave(filename = "wclc/CAF/MMP3_CAF.surv.pdf",p4,width = 6,height = 6)
xlsx::write.xlsx(merged.data,file = "wclc/CAF/CAF.sub.surv.xlsx",sheetName = "MMP3+CAF",append = T)

l1 <- ggpubr::get_legend(p1)
l2 <- ggpubr::get_legend(p2)
l3 <- ggpubr::get_legend(p3)
l4 <- ggpubr::get_legend(p4)

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")

p.list <- list(p1,p2,p3,p4)
pc <- cowplot::plot_grid(plotlist = p.list,nrow = 2)
ggsave(filename = "wclc/CAF/CAF.sub.surv.pdf",pc,width = 8,height = 8)
p.legend <- list(l1,l2,l3,l4)
pc.l <- cowplot::plot_grid(plotlist = p.legend,nrow = 4)
ggsave(filename = "wclc/CAF/CAF.sub.surv.legend.pdf",pc.l,width = 8,height = 4)


### CAF assignment and Immune Score Calculation

rm(list = ls());gc()
CAF <- readRDS('wclc/rds/caf.tumor.rds')

## Assign CAF subset to myCAF & iCAF
iCAF.markers <- c('IL6','CFD','LMNA','AGTR1','HAS1','CXCL1','CXCL2','CXCL8')
myCAF.markers <- c('ACTA2','TAGLN','MMP11','MYL9','POSTN','TPM1','TPM2')

library(garnett)
library(org.Hs.eg.db)
dat <- as(as.matrix(CAF@assays$RNA@counts),'sparseMatrix')
pd <- new('AnnotatedDataFrame',data = CAF@meta.data)
fData <- data.frame(gene_short_name = rownames(dat), row.names = rownames(dat))
fd <- new('AnnotatedDataFrame',data = fData)
mycds <- newCellDataSet(dat,phenoData = pd, featureData = fd)
mycds <- estimateSizeFactors(mycds)
marker_file_path <- "data/CAF.markers.txt"
marker_check <- check_markers(
  mycds, marker_file_path, db = org.Hs.eg.db,
  cds_gene_id_type = 'SYMBOL',
  marker_file_gene_id_type = "SYMBOL"
)

plot_markers(marker_check)

mycds_classifier <- train_cell_classifier(
  cds = mycds, marker_file = marker_file_path,
  db = org.Hs.eg.db, cds_gene_id_type = 'SYMBOL',
  marker_file_gene_id_type = 'SYMBOL',
  num_unknown = 50
)

mycds <- classify_cells(
  mycds, mycds_classifier,
  db = org.Hs.eg.db, cluster_extend = T,
  cds_gene_id_type = 'SYMBOL'
)

table(pData(mycds)$cluster_ext_type, pData(mycds)$celltype)

mypal <- c('#1F78B4','#E31A1C','#FF7F00','#FB9A99')
p1 <- VlnPlot(CAF,features = 'CXCL8',combine = T,pt.size = 1e-5,
              group.by = 'celltype',cols = mypal) + 
  NoLegend() + ggtitle("IL8(iCAF)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p2 <- VlnPlot(CAF,features = 'ACTA2',combine = T,pt.size = 1e-5,
              group.by = 'celltype',cols = mypal) + 
  ggtitle("α-SMA(myCAF)") +
  NoLegend() + theme(axis.title.x = element_blank())
pc <- p1/p2
ggsave("wclc/CAF/iCAF.myCAF.pdf",pc,width = 5,height = 5)

CAF$CAF <- ifelse(CAF$celltype=="CXCL14+CAF","myCAF",'iCAF')
saveRDS(CAF,file = "wclc/rds/caf.tumor.rds")

## Immune Score Calculation
rm(list = ls()[-grep("CAF",ls())]);gc()
load('data/signature_collection.rda')
Tac <- signature_collection$T_cell_accumulation_Peng_et_al
Tac <- intersect(Tac, rownames(CAF))
Tex <- signature_collection$T_cell_exhaustion_Peng_et_al
Tex <- intersect(Tex, rownames(CAF))
Treg <- signature_collection$T_cell_regulatory_Peng_et_al
Treg <- intersect(Treg, rownames(CAF))
ICB <- signature_collection$ICB_resistance_Peng_et_al
ICB <- intersect(ICB, rownames(CAF))
# Roh Immune Score --- Roh et al., 2017
Roh_IS <- c(
  'GZMA','GZMB','PRF1','GNLY',   # CYT, cytolytic activity
  'HLA-A','HLA-B','HLA-C','HLA-E','HLA-F','HLA-G','HLA-H',
  'HLA-DMA','HLA-DMB','HLA-DOA','HLA-DOB','HLA-DPA1',
  'HLA-DPB1','HLA-DQA1','HLA-DQA2','HLA-DQB1','HLA-DRA',
  'HLA-DRB1',   # HLA molecules
  'IFNG','IFNGR1','IFNGR2','IRF1','STAT1','PSMB9',  # IFN-γ pathway genes
  'CCR5','CCL3','CCL4','CCL5','CXCL9','CXCL10','CXCL11',  # chemokines
  'ICAM1','ICAM2','ICAM3','ICAM4','ICAM5','VCAM1'  # adhesion molecules
)
Roh_IS <- intersect(Roh_IS, rownames(CAF))

library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
expr <- CAF@assays$RNA@counts
expr <- log2(expr + 1)
expr <- expr[, order(colnames(expr))]
lst <- list(Tac=Tac,Tex=Tex, Treg=Treg, ICB=ICB,Roh_IS=Roh_IS)
gsva_matrix <- gsva(as.matrix(expr), lst, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = T)

tmp <- as.data.frame(t(gsva_matrix))
tmp$celltype <- plyr::mapvalues(rownames(tmp),from = colnames(CAF),to = CAF$celltype)
tmp <- tmp[order(tmp$celltype,decreasing = F),]
sum.T <- NULL
for(i in 1:(ncol(tmp)-1)){
  sum.T[i] <- data.frame(tapply(tmp[,i],tmp$celltype, sum))
}
names(sum.T) <- colnames(tmp)[-ncol(tmp)]
sum.T <- data.frame(sum.T)
rownames(sum.T) <- names(tapply(tmp[,1],tmp$celltype,sum))
sum.T$celltype <- unname(table(tmp$celltype))
for(i in 1:(ncol(sum.T)-1)){
  for(j in 1:nrow(sum.T)){
    sum.T[j,i] <- sum.T[j,i]/sum.T$celltype[j]
  }
}
sum.T$celltype <- NULL
tmp <- apply(sum.T,2,scale)
rownames(tmp) <- rownames(sum.T)
mypal <- colorRampPalette(RColorBrewer::brewer.pal(9,"GnBu"))(10)
mypal <- mypal[-c(9,10)]
pheatmap::pheatmap(tmp,border_color = NA, color = mypal,
                   filename = "wclc/CAF/immScore.pdf",
                   width = 4,height = 2.5)

xlsx::write.xlsx(sum.T,file = 'wclc/CAF/immScore.xlsx')


###----- Transcription regulation -----------

# ---------- pySCENIC pipe ---------------
rm(list = ls());gc()
#-------- SCENIC preparation -----------

library(SCENIC)
library(SCopeLoomR)
CAF <- readRDS("wclc/rds/caf.tumor.rds")

if(!file.exists("scenicOptions.rds")){
  mydbdir <- "/home/tongqiang/reference/cisTarget/mc9nr"
  mydbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
             "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
  names(mydbs) <- c("500bp","10kb")
  scenicOptions <- initializeScenic(org = "hgnc",nCores = 9, dbDir = mydbdir,
                                    dbs = mydbs, datasetTitle = "CRC")
  saveRDS(scenicOptions,file = "wclc/CAF/GRN/scenicOptions.rds")
}else{
  scenicOptions <- readRDS("scenicOptions.rds")
}

if(!file.exists("expr.mat.filt.rds")){
  CAF <- readRDS("caf.tumor.rds")
  expr.mat <- as.matrix(CAF@assays$RNA@counts)
  runFiltering <- function(expr.mat,stringent = F,scenicOptions){
    if(stringent){
      genesKept <- geneFiltering(expr.mat, scenicOptions,
                                 minCountsPerGene = 3*0.01*ncol(expr.mat),
                                 minSamples = ncol(expr.mat)*0.01)
    }else{
      genesKept <- geneFiltering(expr.mat,scenicOptions = scenicOptions,
                                 minCountsPerGene = 1, minSamples = 20)
    }
    return(expr.mat[genesKept,])
  }
  expr.mat.filt <- runFiltering(expr.mat,scenicOptions = scenicOptions)
  saveRDS(expr.mat.filt, file = "wclc/CAF/GRN/expr.mat.filt.rds")
}else{
  expr.mat.filt <- readRDS("wclc/CAF/GRN/expr.mat.filt.rds")
}


# --------------------- 1.FULL CAF ---------------------

#-------- pySCENIC preparation ---------------

cellInfo <- data.frame(CAF@meta.data)
colnames(cellInfo)[which(colnames(cellInfo) == "orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo) == "seurat_clusters")] <- "cluster"
cellInfo <- cellInfo[, c("sample","TissueSiteSimple","cluster","celltype")]
cellInfo$CellState <- paste(cellInfo$celltype, cellInfo$cluster, sep = "_")
saveRDS(cellInfo, file = "wclc/CAF/GRN/cellInfo.rds")

## step 1

AvgN <- function(data, cell.info, N = 20, seed = 1024, field = "CellState"){
  cell.types <- table(cell.info[[field]])
  cell.types <- cell.types[cell.types >= N]
  cell.types <- sort(cell.types, decreasing = T)
  new.data <- lapply(names(cell.types), function(x){
    new.cell.info <- subset(cell.info, get(field) == x)
    set.seed(seed)
    new.cell.info <- new.cell.info[permute::shuffle(rownames(new.cell.info)),]
    new.cell.info$PoolID <- floor(0:(nrow(new.cell.info)-1)/N)
    pool.ids <- table(new.cell.info$PoolID)
    pool.ids <- pool.ids[pool.ids >= N]
    pool.ids <- sort(names(pool.ids))
    avgN.data <- lapply(pool.ids, function(y){
      select.cells <- rownames(subset(new.cell.info, PoolID == y))
      new.data <- data[, select.cells]
      return(rowMeans(new.data))
    })
    avgN.data <- do.call(cbind, avgN.data)
    colnames(avgN.data) <- paste(x, pool.ids, sep = ".")
    return(avgN.data)
  })
  new.data <- do.call(cbind, new.data)
  return(new.data)
}
# average 20
avg20.rep1 <- AvgN(expr.mat.filt, cellInfo, seed = 1)
avg20.rep2 <- AvgN(expr.mat.filt, cellInfo, seed = 2)
avg20.rep3 <- AvgN(expr.mat.filt, cellInfo, seed = 3)
# average 2
avg2.rep1 <- AvgN(expr.mat.filt, cellInfo, N = 2,seed = 1)
avg2.rep2 <- AvgN(expr.mat.filt, cellInfo, N = 2,seed = 2)
avg2.rep3 <- AvgN(expr.mat.filt, cellInfo, N = 2,seed = 3)

addCellInfo <- function(loom, cellAnno){
  cellAnno <- data.frame(cellAnno)
  if(any(c("nGene","nUMI") %in% colnames(cellAnno))){
    warning("'Columns' 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnno <- cellAnno[, colnames(cellAnno) != "nGene", drop = F]
    cellAnno <- cellAnno[, colnames(cellAnno) != "nUMI", drop = F]
  }
  if(ncol(cellAnno) <= 0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnno))){
    stop("Cell IDs are missing in the annotation")
  }
  cellAnno <- cellAnno[get_cell_ids(loom),,drop=F]
  for(cn in colnames(cellAnno)){
    add_col_attr(loom = loom, key = cn, value = cellAnno[,cn])
  }
  invisible(loom)
}

saveLoom <- function(expr.mat, output){
  cellInfo <- data.frame(
    row.names = colnames(expr.mat),
    celltype <- sub("\\_[0-9]+.*","",colnames(expr.mat))
  )
  loom <- build_loom(output, dgem = expr.mat)
  loom <- addCellInfo(loom,cellInfo)
  close_loom(loom)
}
# average 20
saveLoom(avg20.rep1, "wclc/CAF/GRN/pyscenic/s1_avg20_rep1.loom")
saveLoom(avg20.rep2, "wclc/CAF/GRN/pyscenic/s1_avg20_rep2.loom")
saveLoom(avg20.rep3, "wclc/CAF/GRN/pyscenic/s1_avg20_rep3.loom")
# average 2
saveLoom(avg2.rep1, "wclc/CAF/GRN/pyscenic/s1_avg2_rep1.loom")
saveLoom(avg2.rep2, "wclc/CAF/GRN/pyscenic/s1_avg2_rep2.loom")
saveLoom(avg2.rep3, "wclc/CAF/GRN/pyscenic/s1_avg2_rep3.loom")

loom <- build_loom("wclc/CAF/GRN/s1_expr.mat.loom",dgem = expr.mat.filt)
loom <- addCellInfo(loom,cellInfo)
close_loom(loom)
## step 2
# s2_runPySCENIC.sh ===>> s2_cmd.sh
## step 3
# s3_postSCENIC.py ===>> s3_postSCENIC.sh
## step 4
options(stringsAsFactors = F)
suppressMessages(pacman::p_load(ggplot2,plyr,ggsci))

cellInfo <- readRDS("wclc/CAF/GRN/cellInfo.rds")
CAF <- readRDS("wclc/rds/caf.tumor.rds")
tsne <- CAF@reductions$tsne@cell.embeddings
umap <- CAF@reductions$umap@cell.embeddings
cellInfo <- do.call(cbind, list(cellInfo, tsne,umap))
saveRDS(cellInfo, file = "wclc/CAF/GRN/pyscenic/step4/s4_cellInfo.rds")

# calculating the position of cluster labels
get_label_pos <- function(data, emb="UMAP",group.by="celltype"){
  new.data <- data[, c(paste(emb,1:2,sep = "_"),group.by)]
  colnames(new.data) <- c("x","y","celltype")
  celltypes <- names(table(new.data$celltype))
  new.pos <- lapply(celltypes,function(i){
    tmp.data <- subset(new.data, celltype == i)
    data.frame(x = median(tmp.data$x), y = median(tmp.data$y), label = i)
  })
  do.call(rbind,new.pos)
}

mypal <- c('#1F78B4','#E31A1C','#FF7F00','#FB9A99')
p <- ggplot(cellInfo, aes(UMAP_1, UMAP_2, color = as.character(celltype))) +
  geom_point(size = 0.8) +
  geom_text(inherit.aes = F, data = get_label_pos(cellInfo, emb = "UMAP"),
            aes(x, y, label = label), size = 6) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_discrete_manual(values = mypal,aesthetics = c("color"))
ggsave(filename = "wclc/CAF/GRN/pyscenic/step4/s4_celltype.umap.pdf",p,width = 5,height = 5)

UMAPPlot.celltype <- function(data, celltype.use, topn,color){
  clusters.use <- table(subset(data, celltype == celltype.use)$cluster)
  if(topn > length(clusters.use)) topn <- length(clusters.use)
  clusters.use <- names(head(sort(clusters.use, decreasing = T), n = topn))
  data$new.cluster <- ifelse(data$celltype == celltype.use & data$cluster %in% 
                               clusters.use, data$celltype, "others")
  data$new.cluster <- factor(data$new.cluster, levels = c(setdiff(names(table(data$new.cluster)),"Others"),"Others"))
  data$pt.size <- ifelse(data$new.cluster == "Others", 0.05, 0.2)
  
  ggplot(data, aes(UMAP_1, UMAP_2, color = new.cluster)) +
    geom_point(size = data$pt.size) +
    scale_color_manual(values = c(color,"grey")) +
    guides(colour = guide_legend(keyheight = .7, keywidth = .1, override.aes = list(size=3))) +
    theme_bw(base_size = 12) +
    ggtitle(celltype.use) +
    theme(legend.justification = c(0,0),
          legend.position = c(0,0),
          legend.title = element_blank(),
          legend.key = element_rect(fill = alpha("white", 0)),
          legend.background = element_rect(fill = alpha("white", 0)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = .5, face = "bold"))
}

p1 <- UMAPPlot.celltype(cellInfo,"CCL8+CAF",4,color = "#1F78B4")
p2 <- UMAPPlot.celltype(cellInfo,"CXCL14+CAF",4,color = "#E31A1C")
p3 <- UMAPPlot.celltype(cellInfo,"GREM1+CAF",4,color = "#FF7F00")
p4 <- UMAPPlot.celltype(cellInfo,"MMP3+CAF",4,color = "#FB9A99")

pc <- (p1 + p2)/(p3+p4)
ggsave(filename = "wclc/CAF/GRN/pyscenic/step4/s4_ct.split.umap.pdf",pc,width = 8,height = 8)


## step 5
# philentropy to calculate JS Divergence: http://blog.fens.me/r-entropy/
suppressMessages(pacman::p_load(data.table,pbapply,philentropy,ggrepel,latex2exp,patchwork))

if(!file.exists("wclc/CAF/GRN/pyscenic/step5")){
  dir.create("wclc/CAF/GRN/pyscenic/step5",recursive = T)
}

samples <- c("avg20_rep1","avg20_rep2","avg20_rep3","avg2_rep1","avg2_rep2","avg2_rep3")

for(i in 1:length(samples)){
  # Load RAS matrix
  # rasMat: rows = cellID, cols = regulon, values = Regulon Activity Score
  rasMat <- fread(paste0("wclc/CAF/GRN/pyscenic/step3/s3_",samples[i],".AUCell.txt"),
                  sep = "\t",header = T,data.table = F)
  rownames(rasMat) <- rasMat$V1
  colnames(rasMat) <- sub("(+)","",colnames(rasMat),fixed = T)
  rasMat <- rasMat[,-1]
  saveRDS(rasMat,file = paste0("wclc/CAF/GRN/pyscenic/step5/s5_",samples[i],".rasMat.rds"))
  
  # load cell info
  cellInfo <- readRDS("wclc/CAF/GRN/pyscenic/step4/s4_cellInfo.rds")
  cell.types <- names(table(cellInfo$celltype))
  ctMat <- lapply(cell.types, function(x){as.numeric(cellInfo$celltype == x)})
  ctMat <- do.call(cbind,ctMat)
  colnames(ctMat) <- cell.types
  rownames(ctMat) <- rownames(cellInfo)
  
  # Calculate RSS matrix (Regulon Specificity Score)
  rssMat <- pblapply(colnames(rasMat),function(i){
    sapply(colnames(ctMat),function(j){
      1 - JSD(rbind(rasMat[,i], ctMat[,j]),unit = 'log2',est.prob = 'empirical')
    })
  })
  rssMat <- do.call(rbind,rssMat)
  rownames(rssMat) <- colnames(rasMat)
  colnames(rssMat) <- colnames(ctMat)
  rownames(rssMat) <- sub("(+)","",rownames(rssMat),fixed = T)
  saveRDS(rssMat, file = paste0("wclc/CAF/GRN/pyscenic/step5/s5_",samples[i],".rssMat.rds"))
  
  binMat <- read.table(paste0("wclc/CAF/GRN/pyscenic/step3/s3_",samples[i],".binary_mtx.txt"),
                       sep = "\t",header = T,row.names = 1,check.names = F)
  colnames(binMat) <- sub("(+)","",colnames(binMat),fixed = T)
  saveRDS(binMat, file = paste0("wclc/CAF/GRN/pyscenic/step5/s5_",samples[i],".binMat.rds"))
  
}

PlotRegulonRank <- function(rssMat, cell.type, topn=5){
  data <- data.frame(
    Regulons = 1:nrow(rssMat),
    RSS = sort(rssMat[, cell.type], decreasing = T),
    label = names(sort(rssMat[,cell.type],decreasing = T))
  )
  
  data$pt.col <- ifelse(data$Regulons <= topn, "#007D9B","#BECEE3")
  data <- head(data,200)
  data.label <- head(data,n = topn)
  ggplot(data,aes(Regulons,RSS)) +
    geom_point(size = 3, color = data$pt.col) +
    ggrepel::geom_text_repel(inherit.aes = F,data = data.label,
                             aes(Regulons, RSS, label=label),size=5) +
    ggtitle(cell.type) + ylab("Specificity score") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black",size = 14),
          plot.title = element_text(hjust = .5),
          panel.background = element_rect(colour = "black",size = 1.5),
          axis.title = element_text(size = 15))
}

for(i in 1:length(samples)){
  rssMat <- readRDS(paste0("wclc/CAF/GRN/pyscenic/step5/s5_",samples[i],".rssMat.rds"))
  p1 <- PlotRegulonRank(rssMat,"CCL8+CAF")
  p2 <- PlotRegulonRank(rssMat,"CXCL14+CAF")
  p3 <- PlotRegulonRank(rssMat,"GREM1+CAF")
  p4 <- PlotRegulonRank(rssMat,"MMP3+CAF")
  pc <- (p1+p2)/(p3+p4)
  ggsave(filename = paste0("wclc/CAF/GRN/pyscenic/step5/s5_",samples[i],".rank.pdf"),
         pc,width = 8,height = 10)
}

# library(SCENIC)
# scenicLoomPath <- "wclc/CAF/GRN/pyscenic/step2/s2_avg2_rep1.pyscenic.loom"
# loom <- open_loom(scenicLoomPath)
# # Read information from loom file
# regulons_incidMat <- get_regulons(loom,column.attr.name = "Regulons")
# regulons <- regulonsToGeneLists(regulons_incidMat)
# regulonAUC <- get_regulons_AUC(loom,column.attr.name = "RegulonsAUC")
# CAF <- readRDS("wclc/rds/caf.tumor.rds")
# 
# test <- c("CREB3L1(+)","NEF2L2(+)","MEF2C(+)","TBX2(+)","MAFK(+)","MEF2D(+)")


# ------------ 2.CAF(LCC.vs.RCC) --------------

### Step 1

CAF <- readRDS("wclc/rds/caf.tumor.rds")
LCC.CAF <- subset(CAF, cells = rownames(subset(CAF@meta.data, TissueSiteSimple=="left")))
RCC.CAF <- subset(CAF, cells = rownames(subset(CAF@meta.data, TissueSiteSimple=="right")))
rm(CAF);gc()
saveRDS(LCC.CAF, file = "wclc/CAF/GRN/LCC.CAF.rds")
saveRDS(RCC.CAF, file = "wclc/CAF/GRN/RCC.CAF.rds")

lcc.expr.mat <- as.matrix(LCC.CAF@assays$RNA@counts)
rcc.expr.mat <- as.matrix(RCC.CAF@assays$RNA@counts)

runFiltering <- function(expr.mat,stringent = F,scenicOptions){
  if(stringent){
    genesKept <- geneFiltering(expr.mat, scenicOptions,
                               minCountsPerGene = 3*0.01*ncol(expr.mat),
                               minSamples = ncol(expr.mat)*0.01)
  }else{
    genesKept <- geneFiltering(expr.mat,scenicOptions = scenicOptions,
                               minCountsPerGene = 1, minSamples = 20)
  }
  return(expr.mat[genesKept,])
}

scenicOptions <- readRDS("wclc/CAF/GRN/scenicOptions.rds")

lcc.expr.mat.filt <- runFiltering(lcc.expr.mat,scenicOptions = scenicOptions)
rcc.expr.mat.filt <- runFiltering(rcc.expr.mat,scenicOptions = scenicOptions)
saveRDS(lcc.expr.mat.filt, file = "wclc/CAF/GRN/lcc.expr.mat.filt.rds")
saveRDS(rcc.expr.mat.filt, file = "wclc/CAF/GRN/rcc.expr.mat.filt.rds")

LCC.cellInfo <- data.frame(LCC.CAF@meta.data)
colnames(LCC.cellInfo)[which(colnames(LCC.cellInfo) == "orig.ident")] <- "sample"
colnames(LCC.cellInfo)[which(colnames(LCC.cellInfo) == "seurat_clusters")] <- "cluster"
LCC.cellInfo <- LCC.cellInfo[, c("sample","TissueSiteSimple","cluster","celltype")]
LCC.cellInfo$CellState <- paste(LCC.cellInfo$celltype, LCC.cellInfo$cluster, sep = "_")
saveRDS(LCC.cellInfo, file = "wclc/CAF/GRN/LCC.cellInfo.rds")

RCC.cellInfo <- data.frame(RCC.CAF@meta.data)
colnames(RCC.cellInfo)[which(colnames(RCC.cellInfo) == "orig.ident")] <- "sample"
colnames(RCC.cellInfo)[which(colnames(RCC.cellInfo) == "seurat_clusters")] <- "cluster"
RCC.cellInfo <- RCC.cellInfo[, c("sample","TissueSiteSimple","cluster","celltype")]
RCC.cellInfo$CellState <- paste(RCC.cellInfo$celltype, RCC.cellInfo$cluster, sep = "_")
saveRDS(RCC.cellInfo, file = "wclc/CAF/GRN/RCC.cellInfo.rds")

AvgN <- function(data, cell.info, N = 20, seed = 1024, field = "CellState"){
  cell.types <- table(cell.info[[field]])
  cell.types <- cell.types[cell.types >= N]
  cell.types <- sort(cell.types, decreasing = T)
  new.data <- lapply(names(cell.types), function(x){
    new.cell.info <- subset(cell.info, get(field) == x)
    set.seed(seed)
    new.cell.info <- new.cell.info[permute::shuffle(rownames(new.cell.info)),]
    new.cell.info$PoolID <- floor(0:(nrow(new.cell.info)-1)/N)
    pool.ids <- table(new.cell.info$PoolID)
    pool.ids <- pool.ids[pool.ids >= N]
    pool.ids <- sort(names(pool.ids))
    avgN.data <- lapply(pool.ids, function(y){
      select.cells <- rownames(subset(new.cell.info, PoolID == y))
      new.data <- data[, select.cells]
      return(rowMeans(new.data))
    })
    avgN.data <- do.call(cbind, avgN.data)
    colnames(avgN.data) <- paste(x, pool.ids, sep = ".")
    return(avgN.data)
  })
  new.data <- do.call(cbind, new.data)
  return(new.data)
}
# lcc
lcc.avg2.rep1 <- AvgN(lcc.expr.mat.filt, LCC.cellInfo, seed = 1,N = 2)
lcc.avg2.rep2 <- AvgN(lcc.expr.mat.filt, LCC.cellInfo, seed = 2,N = 2)
lcc.avg2.rep3 <- AvgN(lcc.expr.mat.filt, LCC.cellInfo, seed = 3,N = 2)
# rcc
rcc.avg2.rep1 <- AvgN(rcc.expr.mat.filt, RCC.cellInfo, seed = 1,N = 2)
rcc.avg2.rep2 <- AvgN(rcc.expr.mat.filt, RCC.cellInfo, seed = 2,N = 2)
rcc.avg2.rep3 <- AvgN(rcc.expr.mat.filt, RCC.cellInfo, seed = 3,N = 2)

addCellInfo <- function(loom, cellAnno){
  cellAnno <- data.frame(cellAnno)
  if(any(c("nGene","nUMI") %in% colnames(cellAnno))){
    warning("'Columns' 'nGene' and 'nUMI' wiill not be added as annotations to the loom file.")
    cellAnno <- cellAnno[, colnames(cellAnno) != "nGene", drop = F]
    cellAnno <- cellAnno[, colnames(cellAnno) != "nUMI", drop = F]
  }
  if(ncol(cellAnno) <= 0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnno))){
    stop("Cell IDs are missing in the annotation")
  }
  cellAnno <- cellAnno[get_cell_ids(loom),,drop=F]
  for(cn in colnames(cellAnno)){
    add_col_attr(loom = loom, key = cn, value = cellAnno[,cn])
  }
  invisible(loom)
}

saveLoom <- function(expr.mat, output){
  cellInfo <- data.frame(
    row.names = colnames(expr.mat),
    celltype <- sub("\\_[0-9]+.*","",colnames(expr.mat))
  )
  loom <- build_loom(output, dgem = expr.mat)
  loom <- addCellInfo(loom,cellInfo)
  close_loom(loom)
}

# average 2
saveLoom(lcc.avg2.rep1, "wclc/CAF/GRN/pyscenic/LCC.RCC/step1/s1_lcc_avg2_rep1.loom")
saveLoom(lcc.avg2.rep2, "wclc/CAF/GRN/pyscenic/LCC.RCC/step1/s1_lcc_avg2_rep2.loom")
saveLoom(lcc.avg2.rep3, "wclc/CAF/GRN/pyscenic/LCC.RCC/step1/s1_lcc_avg2_rep3.loom")

saveLoom(rcc.avg2.rep1, "wclc/CAF/GRN/pyscenic/LCC.RCC/step1/s1_rcc_avg2_rep1.loom")
saveLoom(rcc.avg2.rep2, "wclc/CAF/GRN/pyscenic/LCC.RCC/step1/s1_rcc_avg2_rep2.loom")
saveLoom(rcc.avg2.rep3, "wclc/CAF/GRN/pyscenic/LCC.RCC/step1/s1_rcc_avg2_rep3.loom")

loom <- build_loom("wclc/CAF/GRN/pyscenic/LCC.RCC/step1/s1_lcc.expr.mat.loom",dgem = lcc.expr.mat.filt)
loom <- addCellInfo(loom,LCC.cellInfo)
close_loom(loom)

loom <- build_loom("wclc/CAF/GRN/pyscenic/LCC.RCC/step1/s1_rcc.expr.mat.loom",dgem = rcc.expr.mat.filt)
loom <- addCellInfo(loom,RCC.cellInfo)
close_loom(loom)

## step 2
# s2_runPySCENIC.LCC.sh ===>> s2_cmd.sh
# s2_runPySCENIC.RCC.sh ===>> s2_cmd.sh
## step 3
# s3_postSCENIC.py ===>> s3_postSCENIC.sh

## step 4
options(stringsAsFactors = F)
suppressMessages(pacman::p_load(ggplot2,plyr,ggsci))

cellInfo <- readRDS("wclc/CAF/GRN/cellInfo.rds")
LCC.CAF <- readRDS("wclc/CAF/GRN/LCC.CAF.rds")
RCC.CAF <- readRDS("wclc/CAF/GRN/RCC.CAF.rds")

LCC.cellInfo <- readRDS("wclc/CAF/GRN/LCC.cellInfo.rds")
RCC.cellInfo <- readRDS("wclc/CAF/GRN/RCC.cellInfo.rds")

lcc.tsne <- LCC.CAF@reductions$tsne@cell.embeddings
lcc.umap <- LCC.CAF@reductions$umap@cell.embeddings
LCC.cellInfo <- do.call(cbind, list(LCC.cellInfo, lcc.tsne,lcc.umap))
saveRDS(LCC.cellInfo, file = "wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_lcc.cellInfo.rds")
rcc.tsne <- RCC.CAF@reductions$tsne@cell.embeddings
rcc.umap <- RCC.CAF@reductions$umap@cell.embeddings
RCC.cellInfo <- do.call(cbind, list(RCC.cellInfo, rcc.tsne,rcc.umap))
saveRDS(RCC.cellInfo, file = "wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_rcc.cellInfo.rds")


# calculating the position of cluster labels
get_label_pos <- function(data, emb="UMAP",group.by="celltype"){
  new.data <- data[, c(paste(emb,1:2,sep = "_"),group.by)]
  colnames(new.data) <- c("x","y","celltype")
  celltypes <- names(table(new.data$celltype))
  new.pos <- lapply(celltypes,function(i){
    tmp.data <- subset(new.data, celltype == i)
    data.frame(x = median(tmp.data$x), y = median(tmp.data$y), label = i)
  })
  do.call(rbind,new.pos)
}

mypal <- c('#1F78B4','#E31A1C','#FF7F00','#FB9A99')
p <- ggplot(LCC.cellInfo, aes(UMAP_1, UMAP_2, color = as.character(celltype))) +
  geom_point(size = 0.8) +
  geom_text(inherit.aes = F, data = get_label_pos(cellInfo, emb = "UMAP"),
            aes(x, y, label = label), size = 6) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_discrete_manual(values = mypal,aesthetics = c("color"))
ggsave(filename = "wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_lcc.celltype.umap.pdf",p,width = 5,height = 5)
p <- ggplot(RCC.cellInfo, aes(UMAP_1, UMAP_2, color = as.character(celltype))) +
  geom_point(size = 0.8) +
  geom_text(inherit.aes = F, data = get_label_pos(cellInfo, emb = "UMAP"),
            aes(x, y, label = label), size = 6) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_discrete_manual(values = mypal,aesthetics = c("color"))
ggsave(filename = "wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_rcc.celltype.umap.pdf",p,width = 5,height = 5)

UMAPPlot.celltype <- function(data, celltype.use, topn,color){
  clusters.use <- table(subset(data, celltype == celltype.use)$cluster)
  if(topn > length(clusters.use)) topn <- length(clusters.use)
  clusters.use <- names(head(sort(clusters.use, decreasing = T), n = topn))
  data$new.cluster <- ifelse(data$celltype == celltype.use & data$cluster %in% 
                               clusters.use, data$celltype, "others")
  data$new.cluster <- factor(data$new.cluster, levels = c(setdiff(names(table(data$new.cluster)),"Others"),"Others"))
  data$pt.size <- ifelse(data$new.cluster == "Others", 0.05, 0.2)
  
  ggplot(data, aes(UMAP_1, UMAP_2, color = new.cluster)) +
    geom_point(size = data$pt.size) +
    scale_color_manual(values = c(color,"grey")) +
    guides(colour = guide_legend(keyheight = .7, keywidth = .1, override.aes = list(size=3))) +
    theme_bw(base_size = 12) +
    ggtitle(celltype.use) +
    theme(legend.justification = c(0,0),
          legend.position = c(0,0),
          legend.title = element_blank(),
          legend.key = element_rect(fill = alpha("white", 0)),
          legend.background = element_rect(fill = alpha("white", 0)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = .5, face = "bold"))
}
p1 <- UMAPPlot.celltype(LCC.cellInfo,"CCL8+CAF",4,color = "#1F78B4")
p2 <- UMAPPlot.celltype(LCC.cellInfo,"CXCL14+CAF",4,color = "#E31A1C")
p3 <- UMAPPlot.celltype(LCC.cellInfo,"GREM1+CAF",4,color = "#FF7F00")
p4 <- UMAPPlot.celltype(LCC.cellInfo,"MMP3+CAF",4,color = "#FB9A99")

pc <- (p1 + p2)/(p3+p4)
ggsave(filename = "wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_lcc.celltype.split.umap.pdf",
       pc,width = 8,height = 8)

p1 <- UMAPPlot.celltype(RCC.cellInfo,"CCL8+CAF",4,color = "#1F78B4")
p2 <- UMAPPlot.celltype(RCC.cellInfo,"CXCL14+CAF",4,color = "#E31A1C")
p3 <- UMAPPlot.celltype(RCC.cellInfo,"GREM1+CAF",4,color = "#FF7F00")
p4 <- UMAPPlot.celltype(RCC.cellInfo,"MMP3+CAF",4,color = "#FB9A99")

pc <- (p1 + p2)/(p3+p4)
ggsave(filename = "wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_rcc.celltype.split.umap.pdf",
       pc,width = 8,height = 8)

## step 5
# philentropy to calculate JS Divergence: http://blog.fens.me/r-entropy/
suppressMessages(pacman::p_load(data.table,pbapply,philentropy,ggrepel,latex2exp,patchwork))

if(!file.exists("wclc/CAF/GRN/pyscenic/step5")){
  dir.create("wclc/CAF/GRN/pyscenic/LCC.RCC/step5",recursive = T)
}

## LCC
lcc.samples <- c("lcc_avg2_rep1","lcc_avg2_rep2","lcc_avg2_rep3")
for(i in 1:length(lcc.samples)){
  # Load RAS matrix
  # rasMat: rows = cellID, cols = regulon, values = Regulon Activity Score
  rasMat <- fread(paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step3/s3_",lcc.samples[i],".AUCell.txt"),
                  sep = "\t",header = T,data.table = F)
  rownames(rasMat) <- rasMat$V1
  colnames(rasMat) <- sub("(+)","",colnames(rasMat),fixed = T)
  rasMat <- rasMat[,-1]
  saveRDS(rasMat,file = paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_",lcc.samples[i],".rasMat.rds"))
  
  # load cell info
  cellInfo <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_lcc.cellInfo.rds")
  cell.types <- names(table(cellInfo$celltype))
  ctMat <- lapply(cell.types, function(x){as.numeric(cellInfo$celltype == x)})
  ctMat <- do.call(cbind,ctMat)
  colnames(ctMat) <- cell.types
  rownames(ctMat) <- rownames(cellInfo)
  
  # Calculate RSS matrix (Regulon Specificity Score)
  rssMat <- pblapply(colnames(rasMat),function(i){
    sapply(colnames(ctMat),function(j){
      1 - JSD(rbind(rasMat[,i], ctMat[,j]),unit = 'log2',est.prob = 'empirical')
    })
  })
  rssMat <- do.call(rbind,rssMat)
  rownames(rssMat) <- colnames(rasMat)
  colnames(rssMat) <- colnames(ctMat)
  rownames(rssMat) <- sub("(+)","",rownames(rssMat),fixed = T)
  saveRDS(rssMat, file = paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_",lcc.samples[i],".rssMat.rds"))
  
  binMat <- read.table(paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step3/s3_",lcc.samples[i],".binary_mtx.txt"),
                       sep = "\t",header = T,row.names = 1,check.names = F)
  colnames(binMat) <- sub("(+)","",colnames(binMat),fixed = T)
  saveRDS(binMat, file = paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_",lcc.samples[i],".binMat.rds"))
}

## RCC
rcc.samples <- c("rcc_avg2_rep1","rcc_avg2_rep2","rcc_avg2_rep3")
for(i in 1:length(rcc.samples)){
  # Load RAS matrix
  # rasMat: rows = cellID, cols = regulon, values = Regulon Activity Score
  rasMat <- fread(paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step3/s3_",rcc.samples[i],".AUCell.txt"),
                  sep = "\t",header = T,data.table = F)
  rownames(rasMat) <- rasMat$V1
  colnames(rasMat) <- sub("(+)","",colnames(rasMat),fixed = T)
  rasMat <- rasMat[,-1]
  saveRDS(rasMat,file = paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_",rcc.samples[i],".rasMat.rds"))
  
  # load cell info
  cellInfo <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_rcc.cellInfo.rds")
  cell.types <- names(table(cellInfo$celltype))
  ctMat <- lapply(cell.types, function(x){as.numeric(cellInfo$celltype == x)})
  ctMat <- do.call(cbind,ctMat)
  colnames(ctMat) <- cell.types
  rownames(ctMat) <- rownames(cellInfo)
  
  # Calculate RSS matrix (Regulon Specificity Score)
  rssMat <- pblapply(colnames(rasMat),function(i){
    sapply(colnames(ctMat),function(j){
      1 - JSD(rbind(rasMat[,i], ctMat[,j]),unit = 'log2',est.prob = 'empirical')
    })
  })
  rssMat <- do.call(rbind,rssMat)
  rownames(rssMat) <- colnames(rasMat)
  colnames(rssMat) <- colnames(ctMat)
  rownames(rssMat) <- sub("(+)","",rownames(rssMat),fixed = T)
  saveRDS(rssMat, file = paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_",rcc.samples[i],".rssMat.rds"))
  
  binMat <- read.table(paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step3/s3_",rcc.samples[i],".binary_mtx.txt"),
                       sep = "\t",header = T,row.names = 1,check.names = F)
  colnames(binMat) <- sub("(+)","",colnames(binMat),fixed = T)
  saveRDS(binMat, file = paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_",rcc.samples[i],".binMat.rds"))
}

PlotRegulonRank <- function(rssMat, cell.type, topn=5){
  data <- data.frame(
    Regulons = 1:nrow(rssMat),
    RSS = sort(rssMat[, cell.type], decreasing = T),
    label = names(sort(rssMat[,cell.type],decreasing = T))
  )
  
  data$pt.col <- ifelse(data$Regulons <= topn, "#007D9B","#BECEE3")
  data <- head(data,200)
  data.label <- head(data,n = topn)
  ggplot(data,aes(Regulons,RSS)) +
    geom_point(size = 3, color = data$pt.col) +
    ggrepel::geom_text_repel(inherit.aes = F,data = data.label,
                             aes(Regulons, RSS, label=label),size=5) +
    ggtitle(cell.type) + ylab("Specificity score") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black",size = 14),
          plot.title = element_text(hjust = .5),
          panel.background = element_rect(colour = "black",size = 1.5),
          axis.title = element_text(size = 15))
}

samples <- c(lcc.samples,rcc.samples)

for(i in 1:length(samples)){
  rssMat <- readRDS(paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_",samples[i],".rssMat.rds"))
  p1 <- PlotRegulonRank(rssMat,"CCL8+CAF")
  p2 <- PlotRegulonRank(rssMat,"CXCL14+CAF")
  p3 <- PlotRegulonRank(rssMat,"GREM1+CAF")
  p4 <- PlotRegulonRank(rssMat,"MMP3+CAF")
  pc <- (p1+p2)/(p3+p4)
  ggsave(filename = paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_",samples[i],".rank.pdf"),
         pc,width = 8,height = 8)
}


# LCC.CAF <- readRDS("wclc/CAF/GRN/LCC.CAF.rds")
# RCC.CAF <- readRDS("wclc/CAF/GRN/RCC.CAF.rds")
LCC.cellInfo <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_lcc.cellInfo.rds")
RCC.cellInfo <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_rcc.cellInfo.rds")

LCC.binMat <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_lcc_avg2_rep1.binMat.rds")
RCC.binMat <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_rcc_avg2_rep1.binMat.rds")

LCC.rssMat <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_lcc_avg2_rep1.rssMat.rds")
RCC.rssMat <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_rcc_avg2_rep1.rssMat.rds")

LCC.cellInfo <- cbind(LCC.cellInfo, LCC.binMat[rownames(LCC.cellInfo),])
RCC.cellInfo <- cbind(RCC.cellInfo, RCC.binMat[rownames(RCC.cellInfo),])

Dimplot <- function(cell.info,dim.1="UMAP_1",dim.2="UMAP_2",
                    cell.type=NULL, regulon=NULL, color=NULL){
  if(!is.null(cell.type)){
    data <- cell.info[, c(dim.1, dim.2, "celltype")]
    data$pt.col <- ifelse(data$celltype == cell.type, color, "#DFDFDF")
    data$pt.size <- ifelse(data$celltype == cell.type, 0.4, 0.2)
    title <- paste(cell.type)
    col.title <- color
  }else{
    data <- cell.info[, c(dim.1, dim.2, regulon)]
    data$pt.col <- ifelse(data[,regulon],"#006464", "#DFDFDF")
    data$pt.size <- ifelse(data[,regulon],0.4,0.2)
    title <- paste0("Regulon: ", regulon)
    col.title <- "#006464"
  }
  ggplot(data, aes(get(dim.1), get(dim.2))) +
    geom_point(size = data$pt.size, color = data$pt.col) +
    theme_bw(base_size = 12) + ggtitle("") +
    xlab(dim.1) + ylab(dim.2) +
    annotate("text",x=Inf,y=Inf,hjust=1.1,vjust=1.5,label=title,color=col.title,size=6) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(size = 1.5, colour = "black"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 15))
}
mypal <- c('#1F78B4','#E31A1C','#FF7F00','#FB9A99')
## LCC
p1 <- Dimplot(cell.info = LCC.cellInfo, cell.type = "CCL8+CAF",color = "#1F78B4")
p2 <- PlotRegulonRank(LCC.rssMat,"CCL8+CAF")
p3.1 <- Dimplot(cell.info = LCC.cellInfo, regulon = "HOXD9")
p3.2 <- Dimplot(cell.info = LCC.cellInfo, regulon = "HES6")
p3.3 <- Dimplot(cell.info = LCC.cellInfo, regulon = "MYBL1")
p3.4 <- Dimplot(cell.info = LCC.cellInfo, regulon = "HSF2")
p3.5 <- Dimplot(cell.info = LCC.cellInfo, regulon = "ZNF585B")
p.list <- list(p1,p2,p3.1,p3.2,p3.3,p3.4,p3.5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 7,rel_widths = c(4,3,4,4,4,4,4))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/Regulon/LCC.CCL8_CAF.Regulon.pdf",p,width = 26,height = 4)

p1 <- Dimplot(cell.info = LCC.cellInfo,cell.type = "CXCL14+CAF",color = "#E31A1C")
p2 <- PlotRegulonRank(LCC.rssMat,"CXCL14+CAF")
p3.1 <- Dimplot(cell.info = LCC.cellInfo, regulon = "SIN3A")
p3.2 <- Dimplot(cell.info = LCC.cellInfo, regulon = "JUNB")
p3.3 <- Dimplot(cell.info = LCC.cellInfo, regulon = "TCF7")
p3.4 <- Dimplot(cell.info = LCC.cellInfo, regulon = "ZNF274")
p3.5 <- Dimplot(cell.info = LCC.cellInfo, regulon = "RBBP5")
p.list <- list(p1,p2,p3.1,p3.2,p3.3,p3.4,p3.5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 7,rel_widths = c(4,3,4,4,4,4,4))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/Regulon/LCC.CXCL14_CAF.Regulon.pdf",p,width = 26,height = 4)

p1 <- Dimplot(cell.info = LCC.cellInfo, cell.type = "GREM1+CAF",color = "#FF7F00")
p2 <- PlotRegulonRank(LCC.rssMat,"GREM1+CAF")
p3.1 <- Dimplot(cell.info = LCC.cellInfo, regulon = "ATF6")
p3.2 <- Dimplot(cell.info = LCC.cellInfo, regulon = "TFAP4")
p3.3 <- Dimplot(cell.info = LCC.cellInfo, regulon = "NR2C1")
p3.4 <- Dimplot(cell.info = LCC.cellInfo, regulon = "ERF")
p3.5 <- Dimplot(cell.info = LCC.cellInfo, regulon = "CREB3L1")
p.list <- list(p1,p2,p3.1,p3.2,p3.3,p3.4,p3.5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 7,rel_widths = c(4,3,4,4,4,4,4))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/Regulon/LCC.GREM1_CAF.Regulon.pdf",p,width = 26,height = 4)

p1 <- Dimplot(cell.info = LCC.cellInfo, cell.type = "MMP3+CAF",color = "#FB9A99")
p2 <- PlotRegulonRank(LCC.rssMat,"MMP3+CAF")
p3.1 <- Dimplot(cell.info = LCC.cellInfo, regulon = "CREB1")
p3.2 <- Dimplot(cell.info = LCC.cellInfo, regulon = "IRF3")
p3.3 <- Dimplot(cell.info = LCC.cellInfo, regulon = "TCFL5")
p3.4 <- Dimplot(cell.info = LCC.cellInfo, regulon = "IRF6")
p3.5 <- Dimplot(cell.info = LCC.cellInfo, regulon = "ERF")
p.list <- list(p1,p2,p3.1,p3.2,p3.3,p3.4,p3.5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 7,rel_widths = c(4,3,4,4,4,4,4))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/Regulon/LCC.MMP3_CAF.Regulon.pdf",p,width = 26,height = 4)

## RCC
p1 <- Dimplot(cell.info = RCC.cellInfo, cell.type = "CCL8+CAF",color = "#1F78B4")
p2 <- PlotRegulonRank(RCC.rssMat,"CCL8+CAF")
p3.1 <- Dimplot(cell.info = RCC.cellInfo, regulon = "CREB3")
p3.2 <- Dimplot(cell.info = RCC.cellInfo, regulon = "IRF7")
p3.3 <- Dimplot(cell.info = RCC.cellInfo, regulon = "CREB5")
p3.4 <- Dimplot(cell.info = RCC.cellInfo, regulon = "RARA")
p3.5 <- Dimplot(cell.info = RCC.cellInfo, regulon = "ZNF274")
p.list <- list(p1,p2,p3.1,p3.2,p3.3,p3.4,p3.5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 7,rel_widths = c(4,3,4,4,4,4,4))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/Regulon/RCC.CCL8_CAF.Regulon.pdf",p,width = 26,height = 4)

p1 <- Dimplot(cell.info = RCC.cellInfo,cell.type = "CXCL14+CAF",color = "#E31A1C")
p2 <- PlotRegulonRank(RCC.rssMat,"CXCL14+CAF")
p3.1 <- Dimplot(cell.info = RCC.cellInfo, regulon = "FOXC2")
p3.2 <- Dimplot(cell.info = RCC.cellInfo, regulon = "RARA")
p3.3 <- Dimplot(cell.info = RCC.cellInfo, regulon = "TBX2")
p3.4 <- Dimplot(cell.info = RCC.cellInfo, regulon = "JUNB")
p3.5 <- Dimplot(cell.info = RCC.cellInfo, regulon = "ZNF250")
p.list <- list(p1,p2,p3.1,p3.2,p3.3,p3.4,p3.5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 7,rel_widths = c(4,3,4,4,4,4,4))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/Regulon/RCC.CXCL14_CAF.Regulon.pdf",p,width = 26,height = 4)

p1 <- Dimplot(cell.info = RCC.cellInfo, cell.type = "GREM1+CAF",color = "#FF7F00")
p2 <- PlotRegulonRank(RCC.rssMat,"GREM1+CAF")
p3.1 <- Dimplot(cell.info = RCC.cellInfo, regulon = "FOXN1")
p3.2 <- Dimplot(cell.info = RCC.cellInfo, regulon = "TFEC")
p3.3 <- Dimplot(cell.info = RCC.cellInfo, regulon = "HOXC11")
p3.4 <- Dimplot(cell.info = RCC.cellInfo, regulon = "NR2F1")
p3.5 <- Dimplot(cell.info = RCC.cellInfo, regulon = "RUNX2")
p.list <- list(p1,p2,p3.1,p3.2,p3.3,p3.4,p3.5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 7,rel_widths = c(4,3,4,4,4,4,4))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/Regulon/RCC.GREM1_CAF.Regulon.pdf",p,width = 26,height = 4)

p1 <- Dimplot(cell.info = RCC.cellInfo, cell.type = "MMP3+CAF",color = "#FB9A99")
p2 <- PlotRegulonRank(RCC.rssMat,"MMP3+CAF")
p3.1 <- Dimplot(cell.info = RCC.cellInfo, regulon = "IRF6")
p3.2 <- Dimplot(cell.info = RCC.cellInfo, regulon = "IRF9")
p3.3 <- Dimplot(cell.info = RCC.cellInfo, regulon = "FOSB")
p3.4 <- Dimplot(cell.info = RCC.cellInfo, regulon = "HES5")
p3.5 <- Dimplot(cell.info = RCC.cellInfo, regulon = "MLXIPL")
p.list <- list(p1,p2,p3.1,p3.2,p3.3,p3.4,p3.5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 7,rel_widths = c(4,3,4,4,4,4,4))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/Regulon/RCC.MMP3_CAF.Regulon.pdf",p,width = 26,height = 4)

###### Use SEEK to verify regulon
SeekPlot <- function(seek.res, key.words, regulon){
  seek.res$Related <- grepl(key.words, seek.res$Description, ignore.case = T)
  data.use <- data.frame(
    Dataset = seek.res$Rank,
    P = -log10(seek.res$Coexpression.PValue),
    Related = seek.res$Coexpression.PValue < 0.01
  )
  m <- sum(seek.res$Related); n <- sum(data.use$Related)
  M <- nrow(seek.res); N <- sum(seek.res$Coexpression.PValue < 0.01)
  fisher.res <- fisher.test(matrix(c(n,m-n,N-n,M-N-m+n),ncol = 2),alternative = "greater")
  max.p <- max(data.use$P[is.infinite(data.use$P)])
  data.use$P <- ifelse(is.finite(data.use$P), data.use$P, max.p)
  
  ggplot(data.use, aes(Dataset, P)) +
    geom_point(color = "#4590CE", size = 3) +
    geom_point(inherit.aes = F, data = subset(data.use, Related),
               aes(Dataset, P), color = "#E2AE2D", size = 3) +
    geom_hline(yintercept = 2, color = "grey", size = 1) +
    annotate("text",x=Inf,y=Inf,hjust=1.3,vjust=8.5,label=paste0("(",n," out of ",m,")")) +
    annotate("text",x=Inf,y=Inf,hjust=1.4,vjust=10,label=paste0("p=",signif(fisher.res$p.value,3))) +
    annotate("text",x=Inf,y=Inf,hjust=1.2,vjust=2,label=paste0(regulon,"(+)"),size=5,color="black") +
    ylab(TeX("-log_{10}(p-value)")) + ggtitle("") +
    theme_bw(base_size = 15) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black",size = 14),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(colour = "black",size = 1.5))
}

SeekPlot2 <- function(seek.res, key.words, regulon){
  seek.res$Related <- grepl(key.words, seek.res$Description, ignore.case = T)
  data.use <- data.frame(
    Dataset = seek.res$Rank,
    P = -log10(seek.res$Coexpression.PValue),
    Related = seek.res$Coexpression.PValue < 0.05 & seek.res$Related
  )
  m <- sum(seek.res$Related); n <- sum(data.use$Related)
  M <- nrow(seek.res); N <- sum(seek.res$Coexpression.PValue < 0.01)
  fisher.res <- fisher.test(matrix(c(n,m-n,N-n,M-N-m+n),ncol = 2),alternative = "greater")
  max.p <- max(data.use$P[is.infinite(data.use$P)])
  data.use$P <- ifelse(is.finite(data.use$P), data.use$P, max.p)
  
  ggplot(data.use, aes(Dataset, P)) +
    geom_point(color = "#4590CE", size = 3) +
    geom_point(inherit.aes = F, data = subset(data.use, Related),
               aes(Dataset, P), color = "#E2AE2D", size = 3) +
    geom_hline(yintercept = 2, color = "grey", size = 1) +
    annotate("text",x=Inf,y=Inf,hjust=1.3,vjust=8.5,label=paste0("(",n," out of ",m,")")) +
    annotate("text",x=Inf,y=Inf,hjust=1.4,vjust=10,label=paste0("p=",signif(fisher.res$p.value,3))) +
    annotate("text",x=Inf,y=Inf,hjust=1.2,vjust=2,label=paste0(regulon,"(+)"),size=5,color="black") +
    ylab(TeX("-log_{10}(p-value)")) + ggtitle("") +
    theme_bw(base_size = 15) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black",size = 14),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(colour = "black",size = 1.5))
}

### LCC.CCL8+CAF
# HOXD9
HOXD9 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CCL8_CAF_HOXD9",header = T,sep = "\t")
HOXD9 <- na.omit(HOXD9)
p1 <- SeekPlot(HOXD9,key.words = "Fibroblast", regulon = "HOXD9")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.CCL8_CAF.HOXD9.seek.pdf",p1,width = 4,height = 4.2)
# HES6
HES6 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CCL8_CAF_HES6",header = T,sep = "\t")
HES6 <- na.omit(HES6)
p2 <- SeekPlot(HES6,key.words = "Fibroblast", regulon = "HES6")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.CCL8_CAF.HES6.seek.pdf",p2,width = 4,height = 4.2)
# MYBL1
MYBL1 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CCL8_CAF_MYBL1",header = T,sep = "\t")
MYBL1 <- na.omit(MYBL1)
p3 <- SeekPlot(MYBL1,key.words = "Fibroblast", regulon = "MYBL1")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.CCL8_CAF.MYBL1.seek.pdf",p3,width = 4,height = 4.2)
# HSF2
HSF2 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CCL8_CAF_HSF2",header = T,sep = "\t")
HSF2 <- na.omit(HSF2)
p4 <- SeekPlot(HSF2,key.words = "Fibroblast", regulon = "HSF2")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.CCL8_CAF.HSF2.seek.pdf",p4,width = 4,height = 4.2)
# ZNF585B
ZNF585B <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC_CCL8_CAF_ZNF585B",header = T,sep = "\t")
ZNF585B <- na.omit(ZNF585B)
p5 <- SeekPlot(ZNF585B,key.words = "Fibroblast", regulon = "ZNF585B")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.CCL8_CAF.ZNF585B.seek.pdf",p5,width = 4,height = 4.2)
p.list <- list(p1,p2,p3,p4,p5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 5,rel_widths = c(3,3,3,3,3))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CCL8_CAF.seek.verified.pdf",p,width = 13,height = 3)

### LCC.CXCL14+CAF
rm(p1,p2,p3,p4,p5)
# SIN3A
SIN3A <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CXCL14_CAF_SIN3A",header = T,sep = "\t")
SIN3A <- na.omit(SIN3A)
p1 <- SeekPlot(SIN3A,key.words = "Fibroblast", regulon = "SIN3A")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.CXCL14_CAF.SIN3A.seek.pdf",p1,width = 4,height = 4.2)
# JUNB
JUNB <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CXCL14_CAF_JUNB",header = T,sep = "\t")
JUNB <- na.omit(JUNB)
p2 <- SeekPlot2(JUNB,key.words = "Fibroblast", regulon = "JUNB")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.CXCL14_CAF.JUNB.seek.pdf",p2,width = 4,height = 4.2)
# TCF7
TCF7 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CXCL14_CAF_TCF7",header = T,sep = "\t")
TCF7 <- na.omit(TCF7)
p3 <- SeekPlot(TCF7,key.words = "Fibroblast", regulon = "TCF7")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.CXCL14_CAF.TCF7.seek.pdf",p3,width = 4,height = 4.2)
# ZNF274
ZNF274 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CXCL14_CAF_ZNF274",header = T,sep = "\t")
ZNF274 <- na.omit(ZNF274)
p4 <- SeekPlot(ZNF274,key.words = "Fibroblast", regulon = "ZNF274")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.CXCL14_CAF.ZNF274.seek.pdf",p4,width = 4,height = 4.2)
# RBBP5
RBBP5 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CXCL14_CAF_RBBP5",header = T,sep = "\t")
RBBP5 <- na.omit(RBBP5)
p5 <- SeekPlot2(RBBP5,key.words = "Fibroblast", regulon = "RBBP5")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.CXCL14_CAF.RBBP5.seek.pdf",p5,width = 4,height = 4.2)
p.list <- list(p1,p2,p3,p4,p5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 5,rel_widths = c(3,3,3,3,3))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.CXCL14_CAF.seek.verified.pdf",p,width = 13,height = 3)

### LCC.GREM1+CAF
rm(p1,p2,p3,p4,p5)
# ATF6
ATF6 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.GREM1_CAF_ATF6",header = T,sep = "\t")
ATF6 <- na.omit(ATF6)
p1 <- SeekPlot(ATF6,key.words = "Fibroblast", regulon = "ATF6")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.GREM1_CAF.ATF6.seek.pdf",p1,width = 4,height = 4.2)
# TFAP4
TFAP4 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.GREM1_CAF_TFAP4",header = T,sep = "\t")
TFAP4 <- na.omit(TFAP4)
p2 <- SeekPlot(TFAP4,key.words = "Fibroblast", regulon = "TFAP4")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.GREM1_CAF.TFAP4.seek.pdf",p2,width = 4,height = 4.2)
# NR2C1
NR2C1 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.GREM1_CAF_NR2C1",header = T,sep = "\t")
NR2C1 <- na.omit(NR2C1)
p3 <- SeekPlot(NR2C1,key.words = "Fibroblast", regulon = "NR2C1")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.GREM1_CAF.NR2C1.seek.pdf",p3,width = 4,height = 4.2)
# ERF
ERF <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.GREM1_CAF_ERF",header = T,sep = "\t")
ERF <- na.omit(ERF)
p4 <- SeekPlot2(ERF,key.words = "Fibroblast", regulon = "ERF")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.GREM1_CAF.ERF.seek.pdf",p4,width = 4,height = 4.2)
# CREB3L1
CREB3L1 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.GREM1_CAF_CREB3L1",header = T,sep = "\t")
CREB3L1 <- na.omit(CREB3L1)
p5 <- SeekPlot(CREB3L1,key.words = "Fibroblast", regulon = "CREB3L1")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.GREM1_CAF.CREB3L1.seek.pdf",p5,width = 4,height = 4.2)
p.list <- list(p1,p2,p3,p4,p5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 5,rel_widths = c(3,3,3,3,3))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.GREM1_CAF.seek.verified.pdf",p,width = 13,height = 3)

### LCC.MMP3+CAF
rm(p1,p2,p3,p4,p5)
# CREB1
CREB1 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.MMP3_CAF_CREB1",header = T,sep = "\t")
CREB1 <- na.omit(CREB1)
p1 <- SeekPlot2(CREB1,key.words = "Fibroblast", regulon = "CREB1")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.MMP3_CAF.CREB1.seek.pdf",p1,width = 4,height = 4.2)
# IRF3
IRF3 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.MMP3_CAF_IRF3",header = T,sep = "\t")
IRF3 <- na.omit(IRF3)
p2 <- SeekPlot2(IRF3,key.words = "Fibroblast", regulon = "IRF3")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.MMP3_CAF.IRF3.seek.pdf",p2,width = 4,height = 4.2)
# TCFL5
TCFL5 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.MMP3_CAF_TCFL5",header = T,sep = "\t")
TCFL5 <- na.omit(TCFL5)
p3 <- SeekPlot(TCFL5,key.words = "Fibroblast", regulon = "TCFL5")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.MMP3_CAF.TCFL5.seek.pdf",p3,width = 4,height = 4.2)
# IRF6
IRF6 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.MMP3_CAF_IRF6",header = T,sep = "\t")
IRF6 <- na.omit(IRF6)
p4 <- SeekPlot(IRF6,key.words = "Fibroblast", regulon = "IRF6")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.MMP3_CAF.IRF6.seek.pdf",p4,width = 4,height = 4.2)
# ERF
ERF <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.MMP3_CAF_ERF",header = T,sep = "\t")
ERF <- na.omit(ERF)
p5 <- SeekPlot2(ERF,key.words = "Fibroblast", regulon = "ERF")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/LCC.MMP3_CAF.ERF.seek.pdf",p5,width = 4,height = 4.2)
p.list <- list(p1,p2,p3,p4,p5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 5,rel_widths = c(3,3,3,3,3))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/LCC.MMP3_CAF.seek.verified.pdf",p,width = 13,height = 3)

### RCC.CCL8+CAF
rm(p1,p2,p3,p4,p5)
# CREB3
CREB3 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CCL8_CAF_CREB3",header = T,sep = "\t")
CREB3 <- na.omit(CREB3)
p1 <- SeekPlot(CREB3,key.words = "Fibroblast", regulon = "CREB3")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.CCL8_CAF.CREB3.seek.pdf",p1,width = 4,height = 4.2)
# IRF7
IRF7 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CCL8_CAF_IRF7",header = T,sep = "\t")
IRF7 <- na.omit(IRF7)
p2 <- SeekPlot2(IRF7,key.words = "Fibroblast", regulon = "IRF7")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.CCL8_CAF.IRF7.seek.pdf",p2,width = 4,height = 4.2)
# CREB5
CREB5 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CCL8_CAF_CREB5",header = T,sep = "\t")
CREB5 <- na.omit(CREB5)
p3 <- SeekPlot(CREB5,key.words = "Fibroblast", regulon = "CREB5")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.CCL8_CAF.CREB5.seek.pdf",p3,width = 4,height = 4.2)
# RARA
RARA <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CCL8_CAF_RARA",header = T,sep = "\t")
RARA <- na.omit(RARA)
p4 <- SeekPlot(RARA,key.words = "Fibroblast", regulon = "RARA")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.CCL8_CAF.RARA.seek.pdf",p4,width = 4,height = 4.2)
# ZNF274
ZNF274 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CCL8_CAF_ZNF274",header = T,sep = "\t")
ZNF274 <- na.omit(ZNF274)
p5 <- SeekPlot(ZNF274,key.words = "Fibroblast", regulon = "ZNF274")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.CCL8_CAF.ZNF274.seek.pdf",p5,width = 4,height = 4.2)
p.list <- list(p1,p2,p3,p4,p5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 5,rel_widths = c(3,3,3,3,3))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CCL8_CAF.seek.verified.pdf",p,width = 13,height = 3)

### RCC.CXCL14+CAF
rm(p1,p2,p3,p4,p5)
# FOXC2
FOXC2 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CXCL14_CAF_FOXC2",header = T,sep = "\t")
FOXC2 <- na.omit(FOXC2)
p1 <- SeekPlot(FOXC2,key.words = "Fibroblast", regulon = "FOXC2")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.CXCL14_CAF.FOXC2.seek.pdf",p1,width = 4,height = 4.2)
# RARA
RARA <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CXCL14_CAF_RARA",header = T,sep = "\t")
RARA <- na.omit(RARA)
p2 <- SeekPlot(RARA,key.words = "Fibroblast", regulon = "RARA")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.CXCL14_CAF.RARA.seek.pdf",p2,width = 4,height = 4.2)
# TBX2
TBX2 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CXCL14_CAF_TBX2",header = T,sep = "\t")
TBX2 <- na.omit(TBX2)
p3 <- SeekPlot(TBX2,key.words = "Fibroblast", regulon = "TBX2")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.CXCL14_CAF.TBX2.seek.pdf",p3,width = 4,height = 4.2)
# JUNB
JUNB <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CXCL14_CAF_JUNB",header = T,sep = "\t")
JUNB <- na.omit(JUNB)
p4 <- SeekPlot2(JUNB,key.words = "Fibroblast", regulon = "JUNB")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.CXCL14_CAF.JUNB.seek.pdf",p4,width = 4,height = 4.2)
# ZNF250
ZNF250 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CXCL14_CAF_ZNF250",header = T,sep = "\t")
ZNF250 <- na.omit(ZNF250)
p5 <- SeekPlot(ZNF250,key.words = "Fibroblast", regulon = "ZNF250")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.CXCL14_CAF.ZNF250.seek.pdf",p5,width = 4,height = 4.2)
p.list <- list(p1,p2,p3,p4,p5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 5,rel_widths = c(3,3,3,3,3))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.CXCL14_CAF.seek.verified.pdf",p,width = 13,height = 3)

### RCC.GREM1+CAF
rm(p1,p2,p3,p4,p5)
# FOXN1
FOXN1 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.GREM1_CAF_FOXN1",header = T,sep = "\t")
FOXN1 <- na.omit(FOXN1)
p1 <- SeekPlot2(FOXN1,key.words = "Fibroblast", regulon = "FOXN1")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.GREM1_CAF.FOXN1.seek.pdf",p1,width = 4,height = 4.2)
# TFEC
TFEC <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.GREM1_CAF_TFEC",header = T,sep = "\t")
TFEC <- na.omit(TFEC)
p2 <- SeekPlot2(TFEC,key.words = "Fibroblast", regulon = "TFEC")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.GREM1_CAF.TFEC.seek.pdf",p2,width = 4,height = 4.2)
# HOXC11
HOXC11 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.GREM1_CAF_HOXC11",header = T,sep = "\t")
HOXC11 <- na.omit(HOXC11)
p3 <- SeekPlot2(HOXC11,key.words = "Fibroblast", regulon = "HOXC11")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.GREM1_CAF.HOXC11.seek.pdf",p3,width = 4,height = 4.2)
# NR2F1
NR2F1 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.GREM1_CAF_NR2F1",header = T,sep = "\t")
NR2F1 <- na.omit(NR2F1)
p4 <- SeekPlot(NR2F1,key.words = "Fibroblast", regulon = "NR2F1")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.GREM1_CAF.NR2F1.seek.pdf",p4,width = 4,height = 4.2)
# RUNX2
RUNX2 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.GREM1_CAF_RUNX2",header = T,sep = "\t")
RUNX2 <- na.omit(RUNX2)
p5 <- SeekPlot(RUNX2,key.words = "Fibroblast", regulon = "RUNX2")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.GREM1_CAF.RUNX2.seek.pdf",p5,width = 4,height = 4.2)
p.list <- list(p1,p2,p3,p4,p5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 5,rel_widths = c(3,3,3,3,3))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.GREM1_CAF.seek.verified.pdf",p,width = 13,height = 3)

### RCC.MMP3+CAF
rm(p1,p2,p3,p4,p5)
# IRF6
IRF6 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.MMP3_CAF_IRF6",header = T,sep = "\t")
IRF6 <- na.omit(IRF6)
p1 <- SeekPlot(IRF6,key.words = "Fibroblast", regulon = "IRF6")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.MMP3_CAF.IRF6.seek.pdf",p1,width = 4,height = 4.2)
# IRF9
IRF9 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.MMP3_CAF_IRF9",header = T,sep = "\t")
IRF9 <- na.omit(IRF9)
p2 <- SeekPlot2(IRF9,key.words = "Fibroblast", regulon = "IRF9")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.MMP3_CAF.IRF9.seek.pdf",p2,width = 4,height = 4.2)
# FOSB
FOSB <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.MMP3_CAF_FOSB",header = T,sep = "\t")
FOSB <- na.omit(FOSB)
p3 <- SeekPlot2(FOSB,key.words = "Fibroblast", regulon = "FOSB")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.MMP3_CAF.FOSB.seek.pdf",p3,width = 4,height = 4.2)
# HES5
HES5 <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.MMP3_CAF_HES5",header = T,sep = "\t")
HES5 <- na.omit(HES5)
p4 <- SeekPlot(HES5,key.words = "Fibroblast", regulon = "HES5")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.MMP3_CAF.HES5.seek.pdf",p4,width = 4,height = 4.2)
# MLXIPL
MLXIPL <- read.csv("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.MMP3_CAF_MLXIPL",header = T,sep = "\t")
MLXIPL <- na.omit(MLXIPL)
p5 <- SeekPlot2(MLXIPL,key.words = "Fibroblast", regulon = "MLXIPL")
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/single/RCC.MMP3_CAF.MLXIPL.seek.pdf",p5,width = 4,height = 4.2)
p.list <- list(p1,p2,p3,p4,p5)
p <- cowplot::plot_grid(plotlist = p.list,ncol = 5,rel_widths = c(3,3,3,3,3))
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/seek/RCC.MMP3_CAF.seek.verified.pdf",p,width = 13,height = 3)

## Correlation of AUCell score
LCC.AUCell <- read.table("wclc/CAF/GRN/pyscenic/LCC.RCC/step3/s3_lcc_avg2_rep1.AUCell.txt",
                         header = T,sep = "\t",row.names = 1)
cor.regulon <- cor(LCC.AUCell,method = "spearman")
mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(7,"RdYlBu"))(11))
pheatmap::pheatmap(cor.regulon,color = mypal,
                   border_color = NA,
                   width = 6.5,height = 6,
                   show_rownames = F, show_colnames = F,
                   filename = "wclc/CAF/GRN/pyscenic/LCC.RCC/step3/LCC.126regulon.corr.pdf")

RCC.AUCell <- read.table("wclc/CAF/GRN/pyscenic/LCC.RCC/step3/s3_rcc_avg2_rep1.AUCell.txt",
                         header = T,sep = "\t",row.names = 1)
cor.regulon <- cor(RCC.AUCell,method = "spearman")
pheatmap::pheatmap(cor.regulon,color = mypal,
                   border_color = NA,
                   width = 6.5,height = 6,
                   show_rownames = F, show_colnames = F,
                   filename = "wclc/CAF/GRN/pyscenic/LCC.RCC/step3/RCC.157regulon.corr.pdf")

##### correlation heatmap
rm(list = ls());gc()
suppressMessages(pacman::p_load(grid,pbapply,circlize,ggsci,plyr,ggplot2,ggrepel,dendextend,ComplexHeatmap))

### LCC
lcc.rasMat <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_lcc_avg2_rep1.rasMat.rds")
lcc.regulon <- read.table("wclc/CAF/GRN/pyscenic/LCC.RCC/step3/s3_lcc_avg2_rep1.regulons.txt",sep = "\t")
lcc.regulon.names <- sapply(strsplit(as.character(lcc.regulon$V1),split = "\\("),function(x) x[1])
lcc.regulon.sizes <- sapply(strsplit(as.character(lcc.regulon$V3),split = ","),function(x) length(x))
lcc.regulon.names <- lcc.regulon.names[lcc.regulon.sizes >= 10]

lcc.rasMat <- lcc.rasMat[,lcc.regulon.names]
lcc.pccMat <- cor(lcc.rasMat)

CSI <- function(r1,r2,pccMat){
  delta <- pccMat[r1,r2]
  r.others <- setdiff(colnames(pccMat),c(r1,r2))
  N <- sum(pccMat[r1,r.others] < delta) + sum(pccMat[r2,r.others] < delta)
  M <- length(r.others) * 2
  return(N/M)
}

lcc.csiMat <- pblapply(rownames(lcc.pccMat),function(i){
  sapply(colnames(lcc.pccMat), function(j) CSI(i,j,lcc.pccMat))
})

lcc.csiMat <- do.call(rbind, lcc.csiMat)
rownames(lcc.csiMat) <- rownames(lcc.pccMat)

lcc.csiMat.binary <- matrix(as.numeric(lcc.csiMat >= 0.7), nrow = nrow(lcc.csiMat))
colnames(lcc.csiMat.binary) <- colnames(lcc.csiMat)
rownames(lcc.csiMat.binary) <- rownames(lcc.csiMat)
saveRDS(lcc.csiMat, file = "wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_lcc_avg2_rep1.csiMat.rds")
write.table(lcc.csiMat.binary,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_lcc_avg2_rep1.csiMat.binary.txt",
            sep = "\t",quote = F)

# dendrogram
h = 5
row_dend <- as.dendrogram(hclust(dist(lcc.csiMat),method = "complete"))
clusters <- dendextend::cutree(row_dend, h = h)
row_dend <- dendextend::color_branches(row_dend, h = h, col = ggsci::pal_d3("category10")(6))
plot(row_dend)

col_range <- c(0.7, 1)
col_fun <- circlize::colorRamp2(col_range, c("#FCF8DE","#253177"))

pdf("wclc/CAF/GRN/pyscenic/LCC.RCC/step6/LCC.regulon.corr.heatmap.pdf",width = 6,height = 6.2)

ht <- Heatmap(matrix = lcc.csiMat, col = col_fun, name = "ht1", cluster_rows = T,
              cluster_columns = T, show_column_names = F,
              show_row_names = F, show_heatmap_legend = F)
lgd <- Legend(col_fun = col_fun, title = '', at = col_range,
              labels = c('low','high'), direction = 'horizontal',
              legend_width = unit(1,'in'), border = F)
draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c('bottom'))
cols <- c("#786140","#8C0603","#B2712B","#B48521","#C6EEC4","#DCD3B1")
decorate_heatmap_body("ht1",{
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(i) which(i)[1]
  last_index = function(i) {x = which(i); x[length(x)]}
  clusters = names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind==x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  x1 = x1/length(ind)
  x2 = x2/length(ind)
  grid.rect(x = x1, width = (x2-x1), y = 1-x1, height = (x1-x2),
            hjust = 0, vjust = 0, default.units = "npc",
            gp = gpar(fill = NA, col = cols, lwd = 3))
  grid.text(label = paste0("M", clusters),
            x = x2 - length(clusters)/length(ind),
            y = 1 - x1 - (x2 - x1)/2,
            default.units = "npc", hjust = 1, vjust = 0.5,
            gp = gpar(fontsize = 18, fontface = "bold", col = cols))
})

decorate_column_dend("ht1",{
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(i) which(i)[1]
  last_index = function(i) {x = which(i); x[length(x)]}
  clusters = names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  grid.rect(x = x1/length(ind), width = (x2-x1)/length(ind),
            just = 'left', default.units = 'npc',
            gp = gpar(fill=pal_d3("category10")(6),alpha=.5,col=NA))
})

dev.off()

tree <- column_dend(ht)
ind <- cutree(as.hclust(tree),h = h)[order.dendrogram(tree)]
clusters <- names(table(ind))
regulon.clusters <- data.frame(regulon = names(ind), module = paste0("M",ind))
write.table(regulon.clusters,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_lcc_avg2_rep1.regulon_modules.txt",
            sep = "\t",quote = F,row.names = F)

### UMAP plot of different regulon modules

k <- length(clusters)
LCC.cellInfo <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_lcc.cellInfo.rds")
moduleRasMat <- lapply(paste0("M",1:k),function(x){
  regulon.use <- subset(regulon.clusters, module == x)$regulon
  rowMeans(lcc.rasMat[,regulon.use])
})
names(moduleRasMat) <- paste0("M",1:k)
moduleRasMat <- do.call(cbind, moduleRasMat)
LCC.cellInfo <- cbind(LCC.cellInfo, moduleRasMat[rownames(LCC.cellInfo),])
saveRDS(LCC.cellInfo, file = "wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_lcc.cellInfo.rds")

p.list <- lapply(paste0("M",1:k),function(module){
  data.use <- LCC.cellInfo
  expression.color <- rev(RColorBrewer::brewer.pal(11,"RdYlGn"))
  max.val <- quantile(data.use[,module], 0.99)
  min.val <- quantile(data.use[,module], 0.1)
  data.use[,module] <- ifelse(data.use[,module]>max.val, max.val, data.use[,module])
  ggplot(data.use, aes(UMAP_1, UMAP_2, color = get(module))) +
    geom_point(size = 0.6) +
    theme_bw(base_size = 15) +
    ggtitle(module) +
    scale_color_gradientn(name=NULL, colors = expression.color) +
    theme(legend.position = 'right', legend.title = element_blank(),
          plot.title = element_text(hjust = .5, face = 'bold', size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = 'black'),
          panel.background = element_rect(colour = 'black',size = 1.5))
})
for(i in 1:6){
  ggsave(paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step6/module/s6_lcc_M",i,".pdf"),p.list[[i]],width = 4.5,height = 4)
}
pc <- cowplot::plot_grid(plotlist = p.list, ncol = 3)
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_lcc_module_umap.pdf",pc,width = 12,height = 6)

# cell type rank
clusters <- names(table(LCC.cellInfo$cluster))
moduleScoreInCelltype <- lapply(clusters,function(x){
  cells.use <- rownames(subset(LCC.cellInfo), cluster == x)
  colMeans(LCC.cellInfo[cells.use, paste0('M',1:6)])
})
names(moduleScoreInCelltype) <- clusters
moduleScoreInCelltype <- do.call(rbind, moduleScoreInCelltype)
moduleScoreInCelltype <- as.data.frame(moduleScoreInCelltype)
moduleScoreInCelltype$CellType <- mapvalues(
  x = rownames(moduleScoreInCelltype),
  from = LCC.cellInfo$cluster,
  to = LCC.cellInfo$celltype
)
plotCelltypeRank <- function(data, module, topn = 5){
  data.use <- data
  data.use <- data.use[order(data.use[,module],decreasing = T),]
  data.use$Rank <- 1:nrow(data.use)
  data.use$pt.col <- ifelse(data.use$Rank <= topn, "#007D9B", "#BECEE3")
  data.label <- head(data.use, n = topn)
  data.label$delta <- c(Inf, abs(diff(data.label[,module])))
  
  ggplot(data.use, aes(Rank,get(module))) +
    geom_point(size = 3, color = data.use$pt.col) +
    geom_text_repel(inherit.aes = F, data = data.label,
                    aes(Rank,get(module),label=CellType),
                    size = 4, max.iter = 2e4) +
    ggtitle(module) +
    ylab("Regulon activity score") +
    xlab("Cell Type") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          axis.text = element_text(colour = 'black'),
          plot.title = element_text(hjust = .5, face = 'bold'))
}
p.list <- lapply(paste0('M',1:k),function(x){
  plotCelltypeRank(moduleScoreInCelltype,module = x, topn = 10)
})
pc <- cowplot::plot_grid(plotlist = p.list, ncol = 3)
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_lcc_celltype_score.pdf",pc,width = 12,height = 8)


## RCC
rm(list = ls());gc()
rcc.rasMat <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step5/s5_rcc_avg2_rep1.rasMat.rds")
rcc.regulon <- read.table("wclc/CAF/GRN/pyscenic/LCC.RCC/step3/s3_rcc_avg2_rep1.regulons.txt",sep = "\t")
rcc.regulon.names <- sapply(strsplit(as.character(rcc.regulon$V1),split = "\\("),function(x) x[1])
rcc.regulon.sizes <- sapply(strsplit(as.character(rcc.regulon$V3),split = ","),function(x) length(x))
rcc.regulon.names <- rcc.regulon.names[rcc.regulon.sizes >= 10]

rcc.rasMat <- rcc.rasMat[,rcc.regulon.names]
rcc.pccMat <- cor(rcc.rasMat)

CSI <- function(r1,r2,pccMat){
  delta <- pccMat[r1,r2]
  r.others <- setdiff(colnames(pccMat),c(r1,r2))
  N <- sum(pccMat[r1,r.others] < delta) + sum(pccMat[r2,r.others] < delta)
  M <- length(r.others) * 2
  return(N/M)
}

rcc.csiMat <- pblapply(rownames(rcc.pccMat),function(i){
  sapply(colnames(rcc.pccMat), function(j) CSI(i,j,rcc.pccMat))
})

rcc.csiMat <- do.call(rbind, rcc.csiMat)
rownames(rcc.csiMat) <- rownames(rcc.pccMat)

rcc.csiMat.binary <- matrix(as.numeric(rcc.csiMat >= 0.7), nrow = nrow(rcc.csiMat))
colnames(rcc.csiMat.binary) <- colnames(rcc.csiMat)
rownames(rcc.csiMat.binary) <- rownames(rcc.csiMat)
saveRDS(rcc.csiMat, file = "wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_rcc_avg2_rep1.csiMat.rds")
write.table(rcc.csiMat.binary,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_rcc_avg2_rep1.csiMat.binary.txt",
            sep = "\t",quote = F)

# dendrogram
h = 6
row_dend <- as.dendrogram(hclust(dist(rcc.csiMat),method = "complete"))
clusters <- dendextend::cutree(row_dend, h = h)
row_dend <- dendextend::color_branches(row_dend, h = h, col = ggsci::pal_d3("category10")(5))
plot(row_dend)

col_range <- c(0.7, 1)
col_fun <- circlize::colorRamp2(col_range, c("#FCF8DE","#253177"))

pdf("wclc/CAF/GRN/pyscenic/LCC.RCC/step6/RCC.regulon.corr.heatmap2.pdf",width = 6,height = 6.2)

ht <- Heatmap(matrix = rcc.csiMat, col = col_fun, name = "ht1", cluster_rows = T,
              cluster_columns = T, show_column_names = F,
              show_row_names = F, show_heatmap_legend = F)
lgd <- Legend(col_fun = col_fun, title = '', at = col_range,
              labels = c('low','high'), direction = 'horizontal',
              legend_width = unit(1,'in'), border = F)
draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c('bottom'))
cols <- c("#786140","#8C0603","#B2712B","#B48521","#C6EEC4")
decorate_heatmap_body("ht1",{
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(i) which(i)[1]
  last_index = function(i) {x = which(i); x[length(x)]}
  clusters = names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind==x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  x1 = x1/length(ind)
  x2 = x2/length(ind)
  grid.rect(x = x1, width = (x2-x1), y = 1-x1, height = (x1-x2),
            hjust = 0, vjust = 0, default.units = "npc",
            gp = gpar(fill = NA, col = cols, lwd = 3))
  grid.text(label = paste0("M", clusters),
            x = x2 - length(clusters)/length(ind),
            y = 1 - x1 - (x2 - x1)/2,
            default.units = "npc", hjust = 1, vjust = 0.5,
            gp = gpar(fontsize = 18, fontface = "bold", col = cols))
})

decorate_column_dend("ht1",{
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(i) which(i)[1]
  last_index = function(i) {x = which(i); x[length(x)]}
  clusters = names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  grid.rect(x = x1/length(ind), width = (x2-x1)/length(ind),
            just = 'left', default.units = 'npc',
            gp = gpar(fill=pal_d3("category10")(6),alpha=.5,col=NA))
})

dev.off()

tree <- column_dend(ht)
ind <- cutree(as.hclust(tree),h = h)[order.dendrogram(tree)]
clusters <- names(table(ind))
regulon.clusters <- data.frame(regulon = names(ind), module = paste0("M",ind))
write.table(regulon.clusters,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_rcc_avg2_rep1.regulon_modules.txt",
            sep = "\t",quote = F,row.names = F)

## umap plot of regulon modules
k <- length(clusters)
RCC.cellInfo <- readRDS("wclc/CAF/GRN/pyscenic/LCC.RCC/step4/s4_rcc.cellInfo.rds")
moduleRasMat <- lapply(paste0('M',1:k),function(x){
  regulon.use <- subset(regulon.clusters, module == x)$regulon
  rowMeans(rcc.rasMat[,regulon.use])
})
names(moduleRasMat) <- paste0('M',1:k)
moduleRasMat <- do.call(cbind, moduleRasMat)
RCC.cellInfo <- cbind(RCC.cellInfo, moduleRasMat[rownames(RCC.cellInfo),])
saveRDS(RCC.cellInfo,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_rcc.cellInfo.rds")

p.list <- lapply(paste0('M',1:k),function(module){
  data.use <- RCC.cellInfo
  expression.color <- rev(RColorBrewer::brewer.pal(11,"RdYlGn"))
  max.val <- quantile(data.use[,module],0.99)
  min.val <- quantile(data.use[,module],0.1)
  data.use[,module] <- ifelse(data.use[,module]>max.val, max.val, data.use[,module])
  ggplot(data.use, aes(UMAP_1, UMAP_2, color=get(module))) +
    geom_point(size=0.6) + theme_bw(base_size = 15) + ggtitle(module) +
    scale_color_gradientn(name = NULL, colours = expression.color) +
    theme(legend.position = 'right', legend.title = element_blank(),
          plot.title = element_text(hjust = .5, face = 'bold', size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = 'black'),
          panel.background = element_rect(colour = 'black',size = 1.5))
})

for(i in 1:5){
  ggsave(paste0("wclc/CAF/GRN/pyscenic/LCC.RCC/step6/module/s6_rcc_M",i,".pdf"),p.list[[i]],width = 4.5,height = 4)
}
pc <- cowplot::plot_grid(plotlist = p.list, ncol = 3)
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_rcc_module_umap.pdf",pc,width = 12,height = 6)

## cell type rank
clusters <- names(table(RCC.cellInfo$cluster))
moduleScoreInCelltype <- lapply(clusters,function(x){
  cells.use <- rownames(subset(RCC.cellInfo), cluster == x)
  colMeans(RCC.cellInfo[cells.use, paste0('M',1:k)])
})
names(moduleScoreInCelltype) <- clusters
moduleScoreInCelltype <- do.call(rbind, moduleScoreInCelltype)
moduleScoreInCelltype <- as.data.frame(moduleScoreInCelltype)
moduleScoreInCelltype$CellType <- mapvalues(
  x = rownames(moduleScoreInCelltype),
  from = RCC.cellInfo$cluster,
  to = RCC.cellInfo$celltype
)

p.list <- lapply(paste0('M',1:k),function(x){
  plotCelltypeRank(moduleScoreInCelltype,module = x, topn = 10)
})

pc <- cowplot::plot_grid(plotlist = p.list, ncol = 3)
ggsave("wclc/CAF/GRN/pyscenic/LCC.RCC/step6/s6_rcc_celltype_score.pdf",pc,width = 12,height = 8)

## enrichment analysis of Regulons
EnrichmentAnalysis <- function(DEgenes, enrich.method = c("GO","KEGG")) {
  pacman::p_load(clusterProfiler,org.Hs.eg.db,dplyr,patchwork,enrichplot,ggplot2)
  enrich.method <- match.arg(enrich.method)
  if(enrich.method == "GO") {
    ego_BP <- enrichGO(gene = DEgenes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_BP@result$Description <- substring(ego_BP@result$Description, 1, 70)
    return(ego_BP)
    return(ego_results)
  } else if(enrich.method == "KEGG") {
    genelist <- bitr(DEgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    genelist <- pull(genelist, ENTREZID)
    ekegg <- enrichKEGG(gene = genelist, organism = "hsa")
    return(ekegg)
  }
}

# LCC
lcc.uniq <- read.table("wclc/CAF/GRN/pyscenic/LCC.RCC/OL.LCC.uniq.37.txt",
                       header = F,stringsAsFactors = F)$V1
lcc.uniq.go <- EnrichmentAnalysis(DEgenes = lcc.uniq, enrich.method = 'GO')
lcc.uniq.kegg <- EnrichmentAnalysis(DEgenes = lcc.uniq, enrich.method = 'KEGG')
openxlsx::write.xlsx(lcc.uniq.go@result,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/LCC.uniq.enrich.xlsx")
xlsx::write.xlsx(lcc.uniq.kegg@result,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/LCC.uniq.enrich.xlsx",
                 sheetName = "KEGG",append = T,row.names = F)

# RCC
rcc.uniq <- read.table("wclc/CAF/GRN/pyscenic/LCC.RCC/OL.RCC.uniq.68.txt",
                       header = F,stringsAsFactors = F)$V1
rcc.uniq.go <- EnrichmentAnalysis(DEgenes = rcc.uniq, enrich.method = 'GO')
rcc.uniq.kegg <- EnrichmentAnalysis(DEgenes = rcc.uniq, enrich.method = 'KEGG')
openxlsx::write.xlsx(rcc.uniq.go@result,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/RCC.uniq.enrich.xlsx")
xlsx::write.xlsx(rcc.uniq.kegg@result,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/RCC.uniq.enrich.xlsx",
                 sheetName = "KEGG",append = T,row.names = F)

# common
comm <- read.table("wclc/CAF/GRN/pyscenic/LCC.RCC/OL.LCC.and.RCC.common.89.txt",
                   header = F,stringsAsFactors = F)$V1
comm.go <- EnrichmentAnalysis(DEgenes = comm, enrich.method = 'GO')
comm.kegg <- EnrichmentAnalysis(DEgenes = comm, enrich.method = 'KEGG')
openxlsx::write.xlsx(comm.go@result,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/common.enrich.xlsx")
xlsx::write.xlsx(comm.kegg@result,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/common.enrich.xlsx",
                 sheetName = "KEGG",append = T,row.names = F)

kegg <- xlsx::read.xlsx("wclc/CAF/GRN/pyscenic/LCC.RCC/kegg.forVis.xlsx",sheetIndex = 1)
kegg[,ncol(kegg)] <- NULL
kegg$`-log10(pvalue)` <- -log10(kegg$pvalue)
openxlsx::write.xlsx(kegg,file = "wclc/CAF/GRN/pyscenic/LCC.RCC/kegg.forVis.xlsx")

kegg <- openxlsx::read.xlsx('wclc/CAF/GRN/pyscenic/LCC.RCC/kegg.forVis.xlsx',sheet = 2)

colnames(kegg) <- gsub("\\."," ",colnames(kegg))
colnames(kegg)[9] <- gsub("-beta","β",colnames(kegg)[9])
colnames(kegg)[11] <- gsub("-kappa B","κB",colnames(kegg)[11])
rownames(kegg) <- kegg[,1]
kegg[,1] <- NULL
kegg <- rbind(rep(4,11),rep(0,6),kegg)

colors_border <- c(rgb(0.8549,0.64706,0.12549, 0.9),
                   rgb(1,0.49804,0.31373, 0.9),
                   rgb(0.81569, 0.12549, 0.9))
colors_in <- c(rgb(0.8549,0.64706,0.12549,0.5),
               rgb(1,0.49804, 0.31373, 0.5),
               rgb(0.81569, 0.12549, 0.56471, 0.4))

colors_border <- c(rgb(0.06275,0.41569,0.69804,0.9),
                   rgb(1,0.5255,0.13333,0.9),
                   rgb(0.70588,0.11373,0.13725,0.9))
colors_in <- c(rgb(0.06275,0.41569,0.69804,0.4),
                   rgb(1,0.5255,0.13333,0.5),
                   rgb(0.70588,0.11373,0.13725,0.4))
pdf("wclc/CAF/GRN/pyscenic/LCC.RCC/kegg.pdf",width = 5,height = 5)
radarchart(kegg, axistype = 0, pcol = colors_border,
           pfcol = colors_in, plwd = 1.3, plty = 1, pty = 32,
           cglcol = 'black', cglty = 3, cglwd = 0.6)
legend(x=1.5,y=1,legend=rownames(kegg[-c(1,2),]),bty='n',pch=20,
       col = colors_border, text.col = 'black',cex = 1,pt.cex = 2)
dev.off()




# library(pheatmap)
# binMat <- read.table("wclc/CAF/GRN/pyscenic/step3/s3_avg2_rep1.binary_mtx.txt",
#                      header = T,sep = "\t",row.names = 1)
# 
# colnames(binMat) <- gsub("\\...","(+)",colnames(binMat))
# CAF <- readRDS("wclc/rds/caf.tumor.rds")
# 
# meta <- CAF@meta.data[,c("TissueSiteSimple","celltype")]
# anno_cols <- list(
#   TissueSiteSimple = c(left="#4E9BD2", right="#DB4740"),
#   celltype = c(`CCL8+CAF`="#1F78B4",`CXCL14+CAF`="#E31A1C",`GREM1+CAF`="#FF7F00",`MMP3+CAF`="#FB9A99")
# )
# 
# binMat$TISS <- plyr::mapvalues(rownames(binMat),from = rownames(CAF@meta.data),to = as.character(CAF$TissueSiteSimple))
# binMat$barcode <- rownames(binMat)
# binMat <- binMat[order(binMat[,ncol(binMat)-1]),]
# binMat$barcode <- NULL;binMat$TISS <- NULL
# pheatmap(binMat,show_colnames = F,annotation_row = meta,
#          color = c("#FCF8DE","#253177"),
#          annotation_colors = anno_cols,
#          width = 16,height = 8,show_rownames = F,
#          cluster_rows = F)


# ---------- SCENIC pipe ---------------

# suppressPackageStartupMessages(pacman::p_load(future, future.apply))
# plan("multiprocess", workers = 10)
# options(future.globals.maxSize = 50000 * 1024^2)
# expr.mat.filt <- readRDS("expr.mat.filt.rds")
# library(SCENIC)
# scenicOptions <- readRDS("scenicOptions.rds")
# # compute correlation matrix
# runCorrelation(expr.mat.filt, scenicOptions)
# # TF-targets correlation regression analysis
# expr.mat.filt.log <- log2(expr.mat.filt + 1)
# runGenie3(expr.mat.filt.log, scenicOptions, nParts = 10)
# # Infer co-expression modules
# runSCENIC_1_coexNetwork2modules(scenicOptions)
# # Infer GRNs (regulon)
# runSCENIC_2_createRegulons(scenicOptions,coexMethods = c("top5perTarget"))
# source("script/runSCENIC_3_scoreCells_mod.R")
# runSCENIC_3_scoreCells_mod(ScenicOptions,exp.mat.norm="wclc/CAF/GRN/expr.mat.filt.norm.rds")
# runSCENIC_4_aucell_binarize(scenicOptions, exprMat = expr.mat.filt.log)
# 
# ### Visualization via Seurat
# rm(list = ls());gc()
# CAF <- readRDS("wclc/rds/caf.tumor.rds")
# 
# if(!file.exists("wclc/CAF/GRN/scenic/scenic_seurat")){
#   dir.create("wclc/CAF/GRN/scenic/scenic_seurat",recursive = T)
# }
# 
# # import raw regulonAUC matrix
# AUCmatrix <- readRDS("wclc/CAF/GRN/scenic/int/3.4_regulonAUC.Rds")
# AUCmatrix <- AUCmatrix@assays@data@listData$AUC
# AUCmatrix <- data.frame(t(AUCmatrix), check.names = F)
# RegulonName_AUC <- colnames(AUCmatrix)
# RegulonName_AUC <- gsub(" \\(", "_", RegulonName_AUC)
# RegulonName_AUC <- gsub("\\)", "", RegulonName_AUC)
# colnames(AUCmatrix) <- RegulonName_AUC
# CAF <- AddMetaData(CAF, AUCmatrix)
# saveRDS(CAF, file = "wclc/rds/CAF.scenic.auc.rds")
# 
# rm(CAF);gc()
# 
# CAF <- readRDS("wclc/rds/caf.tumor.rds")
# # import binary regulon matrix
# BINmatrix <- readRDS("wclc/CAF/GRN/scenic/int/4.1_binaryRegulonActivity.Rds")
# BINmatrix <- data.frame(t(BINmatrix), check.names = F)
# RegulonName_BIN <- colnames(BINmatrix)
# RegulonName_BIN <- gsub(" \\(","_",RegulonName_BIN)
# RegulonName_BIN <- gsub("\\)","",RegulonName_BIN)
# colnames(BINmatrix) <- RegulonName_BIN
# CAF <- AddMetaData(CAF, BINmatrix)
# saveRDS(CAF, file = "wclc/rds/CAF.scenic.bin.rds")
# 
# 
# meta <- CAF@meta.data[,c("TissueSiteSimple","celltype")]
# anno_cols <- list(
#   TissueSiteSimple = c(left="#4E9BD2", right="#DB4740"),
#   celltype = c(`CCL8+CAF`="#1F78B4",`CXCL14+CAF`="#E31A1C",`GREM1+CAF`="#FF7F00",`MMP3+CAF`="#FB9A99")
# )
# pheatmap(t(BINmatrix),show_colnames = F,annotation_col = meta,
#          color = c("#FCF8DE","#253177"),
#          annotation_colors = anno_cols,
#          width = 16,height = 8,
#          filename = "wclc/CAF/GRN/scenic/scenic_seurat/Binary.TF.activity.Heatmap.pdf")
# 
# CAF <- readRDS("wclc/rds/CAF.scenic.auc.rds")
# 
# mytheme <- theme(axis.title.x = element_blank(),
#                  axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0),
#                  panel.border = element_rect(size = 1,colour = "black"))
# mypal <- c('#1F78B4','#E31A1C','#FF7F00','#FB9A99')
# 
# for(i in 1:ncol(AUCmatrix)){
#   p <- VlnPlot(CAF,features = colnames(AUCmatrix)[i],group.by = "TissueSiteSimple",
#                cols = c("#4E9BD2","#DB4740"),pt.size = 1e-5) + 
#     NoLegend() + mytheme + ylab("TF Activity") +
#     ggpubr::stat_compare_means(comparisons = list(c("left","right")))
#   ggsave(filename = paste("wclc/CAF/GRN/scenic/scenic_seurat/lcc.vs.rcc_",colnames(AUCmatrix)[i],".pdf"),
#          p,width = 4,height = 4)
#   p <- VlnPlot(CAF,features = colnames(AUCmatrix)[i],group.by = "celltype",cols = mypal,pt.size = 1e-5) + 
#     NoLegend() + mytheme + ylab("TF Activity") +
#     theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))
#   ggsave(filename = paste("wclc/CAF/GRN/scenic/scenic_seurat/celltype_",colnames(AUCmatrix)[i],".pdf"),
#          p,width = 5,height = 4)
# }


### 比较CXCL14+CAF和其他CAF亚群的差异
if(F){
  rm(list = ls());gc()
  caf <- readRDS('wclc/rds/caf.tumor.rds')
  caf$celltypeII <- ifelse(caf$celltype == 'CXCL14+CAF', 'CXCL14+CAF', 'otherCAFs')
  
  p <- DimPlot(caf, group.by = 'celltypeII',cols = c("#C72719","#BEBBBD"),label = T,label.size = 5) +
    NoLegend() + ggtitle('') +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none',
          panel.border = element_rect(size = 1,colour = 'black'))
  ggsave(filename = 'wclc/CAF/CXCL14+CAF.pdf',p,width = 5,height = 5)
  
  CXCL14.markers <- FindMarkers(caf, group.by = 'celltypeII', ident.1 = 'CXCL14+CAF', ident.2 = 'otherCAFs')
  # remove ribosomal and mitochondrial genes
  CXCL14.markers <- CXCL14.markers[!grepl("^RP[SL]",rownames(CXCL14.markers)),]
  CXCL14.markers <- CXCL14.markers[!grepl("^MT-",rownames(CXCL14.markers)),]
  CXCL14.markers$Significance <- ifelse(CXCL14.markers$p_val < 0.01, TRUE, FALSE)
  CXCL14.markers <- CXCL14.markers[order(CXCL14.markers[,2],decreasing = T),]
  top10 <- rbind(head(CXCL14.markers,10),tail(CXCL14.markers,10))
  p.all.markers <- ggplot(CXCL14.markers, aes(avg_logFC,-log10(p_val))) +
    geom_point(aes(col = Significance)) +
    scale_color_manual(values = c("#BEBBBD","#C72719")) +
    geom_text_repel(data = top10, aes(label = rownames(top10)),
                    box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                    segment.color = "black", show.legend = F) +
    theme_classic() + ylab("-log10(P value)") + mytheme +
    theme(legend.position = 'none')
  ggsave(filename = 'wclc/CAF/CXCL14caf_vs_otherCAFs.volc.pdf',p.all.markers,width = 5,height = 5)
  
  up.genes <- rownames(CXCL14.markers[CXCL14.markers$avg_logFC > 0,])
  dn.genes <- rownames(CXCL14.markers[CXCL14.markers$avg_logFC < 0,])
  
  CXCL14pos.up <- EnrichmentAnalysis(DEgenes = up.genes, enrich.method = 'GO')
  CXCL14pos.dn <- EnrichmentAnalysis(DEgenes = dn.genes, enrich.method = 'GO')
  xlsx::write.xlsx(CXCL14pos.up$ego_BP@result, file = 'wclc/CAF/CXCL14pos_vs_otherCAFs.GO.xlsx',sheetName = 'up',append = F)
  xlsx::write.xlsx(CXCL14pos.dn$ego_BP@result, file = 'wclc/CAF/CXCL14pos_vs_otherCAFs.GO.xlsx',sheetName = 'down',append = T)
  
  lcc.caf <- subset(caf, cells = rownames(subset(caf@meta.data, TissueSiteSimple == 'left')))
  rcc.caf <- subset(caf, cells = rownames(subset(caf@meta.data, TissueSiteSimple == 'right')))
  
  p <- DimPlot(lcc.caf, group.by = 'celltypeII',cols = c("#C72719","#BEBBBD"),label = T,label.size = 5) +
    NoLegend() + ggtitle('') +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none',
          panel.border = element_rect(size = 1,colour = 'black'))
  ggsave(filename = 'wclc/CAF/lcc.CXCL14+CAF.pdf',p,width = 5,height = 5)
  p <- DimPlot(rcc.caf, group.by = 'celltypeII',cols = c("#C72719","#BEBBBD"),label = T,label.size = 5) +
    NoLegend() + ggtitle('') +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none',
          panel.border = element_rect(size = 1,colour = 'black'))
  ggsave(filename = 'wclc/CAF/rcc.CXCL14+CAF.pdf',p,width = 5,height = 5)
  
  # LCC
  CXCL14.markers <- FindMarkers(lcc.caf, group.by = 'celltypeII', ident.1 = 'CXCL14+CAF', ident.2 = 'otherCAFs')
  # remove ribosomal and mitochondrial genes
  CXCL14.markers <- CXCL14.markers[!grepl("^RP[SL]",rownames(CXCL14.markers)),]
  CXCL14.markers <- CXCL14.markers[!grepl("^MT-",rownames(CXCL14.markers)),]
  CXCL14.markers$Significance <- ifelse(CXCL14.markers$p_val < 0.01, TRUE, FALSE)
  CXCL14.markers <- CXCL14.markers[order(CXCL14.markers[,2],decreasing = T),]
  top10 <- rbind(head(CXCL14.markers,10),tail(CXCL14.markers,10))
  p.all.markers <- ggplot(CXCL14.markers, aes(avg_logFC,-log10(p_val))) +
    geom_point(aes(col = Significance)) +
    scale_color_manual(values = c("#BEBBBD","#C72719")) +
    geom_text_repel(data = top10, aes(label = rownames(top10)),
                    box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                    segment.color = "black", show.legend = F) +
    theme_classic() + ylab("-log10(P value)") + mytheme +
    theme(legend.position = 'none')
  ggsave(filename = 'wclc/CAF/LCC.CXCL14caf_vs_otherCAFs.volc.pdf',p.all.markers,width = 5,height = 5)
  
  up.genes <- rownames(CXCL14.markers[CXCL14.markers$avg_logFC > 0,])
  dn.genes <- rownames(CXCL14.markers[CXCL14.markers$avg_logFC < 0,])
  
  CXCL14pos.up <- EnrichmentAnalysis(DEgenes = up.genes, enrich.method = 'GO')
  CXCL14pos.dn <- EnrichmentAnalysis(DEgenes = dn.genes, enrich.method = 'GO')
  xlsx::write.xlsx(CXCL14pos.up$ego_BP@result, file = 'wclc/CAF/LCC.CXCL14pos_vs_otherCAFs.GO.xlsx',sheetName = 'up',append = F)
  xlsx::write.xlsx(CXCL14pos.dn$ego_BP@result, file = 'wclc/CAF/LCC.CXCL14pos_vs_otherCAFs.GO.xlsx',sheetName = 'down',append = T)
  
  # RCC
  CXCL14.markers <- FindMarkers(rcc.caf, group.by = 'celltypeII', ident.1 = 'CXCL14+CAF', ident.2 = 'otherCAFs')
  # remove ribosomal and mitochondrial genes
  CXCL14.markers <- CXCL14.markers[!grepl("^RP[SL]",rownames(CXCL14.markers)),]
  CXCL14.markers <- CXCL14.markers[!grepl("^MT-",rownames(CXCL14.markers)),]
  CXCL14.markers$Significance <- ifelse(CXCL14.markers$p_val < 0.01, TRUE, FALSE)
  CXCL14.markers <- CXCL14.markers[order(CXCL14.markers[,2],decreasing = T),]
  top10 <- rbind(head(CXCL14.markers,10),tail(CXCL14.markers,10))
  p.all.markers <- ggplot(CXCL14.markers, aes(avg_logFC,-log10(p_val))) +
    geom_point(aes(col = Significance)) +
    scale_color_manual(values = c("#BEBBBD","#C72719")) +
    geom_text_repel(data = top10, aes(label = rownames(top10)),
                    box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                    segment.color = "black", show.legend = F) +
    theme_classic() + ylab("-log10(P value)") + mytheme +
    theme(legend.position = 'none')
  ggsave(filename = 'wclc/CAF/RCC.CXCL14caf_vs_otherCAFs.volc.pdf',p.all.markers,width = 5,height = 5)
  
  up.genes <- rownames(CXCL14.markers[CXCL14.markers$avg_logFC > 0,])
  dn.genes <- rownames(CXCL14.markers[CXCL14.markers$avg_logFC < 0,])
  
  CXCL14pos.up <- EnrichmentAnalysis(DEgenes = up.genes, enrich.method = 'GO')
  CXCL14pos.dn <- EnrichmentAnalysis(DEgenes = dn.genes, enrich.method = 'GO')
  xlsx::write.xlsx(CXCL14pos.up$ego_BP@result, file = 'wclc/CAF/RCC.CXCL14pos_vs_otherCAFs.GO.xlsx',sheetName = 'up',append = F)
  xlsx::write.xlsx(CXCL14pos.dn$ego_BP@result, file = 'wclc/CAF/RCC.CXCL14pos_vs_otherCAFs.GO.xlsx',sheetName = 'down',append = T)
}

## output results
rm(list = ls());gc()
load('wclc/CAF/all.caf.DEgenes.lcc.rcc.rda')
caf.markers$genes <- rownames(caf.markers)
openxlsx::write.xlsx(caf.markers, file = 'wclc/CAF/DEgenes.xlsx')
load('wclc/CAF/ccl8.caf.DEgenes.lcc.rcc.rda')
ccl8.caf.markers$genes <- rownames(ccl8.caf.markers)
xlsx::write.xlsx(ccl8.caf.markers,file = 'wclc/CAF/DEgenes.xlsx',sheetName = 'CCL8+CAF(RCC_vs_LCC)',append = T)
load('wclc/CAF/cxcl14.caf.DEgenes.lcc.rcc.rda')
cxcl14.caf.markers$genes <- rownames(cxcl14.caf.markers)
xlsx::write.xlsx(cxcl14.caf.markers,file = 'wclc/CAF/DEgenes.xlsx',sheetName = 'CXCL14+CAF(RCC_vs_LCC)',append = T)

load('wclc/CAF/grem1.caf.DEgenes.lcc.rcc.rda')
grem1.caf.markers$genes <- rownames(grem1.caf.markers)
xlsx::write.xlsx(grem1.caf.markers,file = 'wclc/CAF/DEgenes.xlsx',sheetName = 'GREM1+CAF(RCC_vs_LCC)',append = T)

load('wclc/CAF/mmp3.caf.DEgenes.lcc.rcc.rda')
mmp3.caf.markers$genes <- rownames(mmp3.caf.markers)
xlsx::write.xlsx(grem1.caf.markers,file = 'wclc/CAF/DEgenes.xlsx',sheetName = 'MMP3+CAF(RCC_vs_LCC)',append = T)





