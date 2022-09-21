
# ---------------- 1. DimRed ----------------

library(Seurat)
library(ggplot2)
library(ggrepel)
library(ggpubr)

myeloid <- readRDS('data/GSE178341/rds/myeloid.rds')
myeloid$clMidwayPr <- as.character(myeloid$clMidwayPr)

myeloid <- NormalizeData(object = myeloid, normalization.method = "LogNormalize", scale.factor = 10000)
myeloid <- FindVariableFeatures(myeloid,selection.method="vst",nfeatures=2000)
myeloid <- ScaleData(myeloid, features = VariableFeatures(myeloid))
myeloid <- RunPCA(object = myeloid, features = VariableFeatures(myeloid))
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
PC.num <- PCDeterminators(myeloid)
# Find neighbors and clusters with harmony batch correction
myeloid <- FindNeighbors(myeloid, dims = 1:PC.num, reduction = "pca")
myeloid <- FindClusters(myeloid, resolution = 0.2)
# myeloid <- myeloid %>% RunHarmony("SampleName",plot_convergence = F)
myeloid <- RunUMAP(myeloid, dims = 1:PC.num, reduction = "pca")
myeloid <- RunTSNE(myeloid,dims = 1:PC.num)

# vis
p1 <- DimPlot(myeloid, group.by = "SPECIMEN_TYPE",reduction = 'tsne',cols = c("#4E9BD3","#DB4840")) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.title = element_blank(),
        legend.text = element_text(size = 14))
ggsave(filename = 'wclc/mye.N_vs_T.pdf',p1,width = 5,height = 4.8)

p2 <- DimPlot(myeloid, group.by = "clMidwayPr",reduction = 'tsne',cols = c("#224444","#e68a00","#33adff","#439373")) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom",
        plot.title = element_blank(),
        legend.text = element_text(size = 14))
ggsave(filename = 'wclc/mye.celltype.pdf',p2,width = 5,height = 5)

p3 <- DimPlot(myeloid, group.by = "TissueSiteSimple",reduction = 'tsne',cols = c("#4E9BD3","#DB4840")) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.title = element_blank(),
        legend.text = element_text(size = 14))
ggsave(filename = 'wclc/mye.left_right.pdf',p3,width = 5,height = 4.8)

saveRDS(myeloid, file = 'wclc/myeloid.total.rds')

# extract tumor samples
mye.T <- subset(myeloid, cells = rownames(subset(myeloid@meta.data, SPECIMEN_TYPE == 'T')))
mye.T$SPECIMEN_TYPE <- as.character(mye.T$SPECIMEN_TYPE)

p.mye <- DimPlot(mye.T, reduction = 'tsne',split.by = 'TissueSiteSimple',group.by = 'clMidwayPr',cols = c("#224444","#e68a00","#33adff","#439373")) +
  theme(panel.border = element_rect(size = 1,colour = 'black'),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = 'bottom')
ggsave(filename = 'wclc/p.mye.pdf',p.mye,width = 8,height = 5)

saveRDS(mye.T, file = 'wclc/myeloid.tumor.sub.rds')

# visualization
p1 <- DimPlot(mye.T, group.by = "TissueSiteSimple",reduction = 'tsne',cols = c("#4E9BD3","#C72719")) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.title = element_blank(),
        legend.text = element_text(size = 14))
ggsave(filename = 'wclc/mye.T.LCC_RCC.pdf',p1,width = 5,height = 4.6)

p2 <- DimPlot(mye.T, group.by = "clMidwayPr",reduction = 'tsne',cols = c("#224444","#e68a00","#33adff","#439373")) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom",
        plot.title = element_blank(),
        legend.text = element_text(size = 14))
ggsave(filename = 'wclc/mye.T.ct.pdf',p2,width = 5,height = 5)


mytheme.1 <- theme(
  panel.border = element_rect(size = 1.5, colour = "grey75"),
  panel.grid = element_blank(),
  axis.text = element_text(size = 12,face = "bold",colour = "grey25"),
  axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
  legend.position = "top", legend.direction = "horizontal",
  legend.text = element_text(size = 14),legend.title = element_text(size = 15)
)

mycols <- c("#224444","#e68a00","#33adff","#439373")

p3 <- ggplot(mye.T@meta.data, aes(x = TissueSiteSimple, fill = clMidwayPr)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mycols) +
  mytheme.1 +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(filename = 'wclc/mye.T.cell.prop.left_right.pdf',p3,width = 8,height = 3)

dat <- mye.T@meta.data[,c('PatientTypeID','clMidwayPr')]
dat$PatientTypeID <- factor(dat$PatientTypeID, 
                            levels = c('C103_T','C124_T','C125_T','C129_T','C132_T','C134_T',
                                       'C136_T','C140_T','C150_T','C153_T','C157_T','C166_T',
                                       'C173_T','C104_T','C105_T','C106_T','C107_T','C109_T',
                                       'C110_T','C111_T','C112_T','C113_T','C114_T','C115_T',
                                       'C116_T','C118_T','C119_T','C122_T','C123_T','C126_T',
                                       'C130_TA','C130_TB','C133_T','C135_T','C137_T','C138_T',
                                       'C139_T','C142_T','C143_T','C144_T','C145_T','C146_T',
                                       'C147_T','C149_T','C151_T','C152_T','C154_T','C155_T',
                                       'C156_T','C158_T','C159_T','C160_T','C161_T','C162_T',
                                       'C163_T','C164_T','C165_T','C167_T','C168_T','C169_T',
                                       'C170_T','C171_TA','C171_TB','C172_T'),
                            ordered = T)

p4 <- ggplot(dat, aes(x = PatientTypeID, fill = clMidwayPr)) +
  geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values = mycols) +
  mytheme.1 +
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1,size = 10)) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction",x = NULL)

ggsave(filename = 'wclc/mye.T.cell.prop.pdf',p4,width = 20,height = 5)

# tbl <- data.frame(site = c(rep('left',4),rep('right',4)), 
#            celltype = rep(c('DC','Granulo','Macro','Mono'),2), 
#            value = c(0.131,0.012,0.539,0.318,0.119,0.064,0.469,0.349))
# ggplot(tbl, aes(celltype,value,fill=site)) +
#   geom_bar(stat = 'identity',position = 'dodge') +
#   ggsci::scale_fill_jco() +
#   ggpubr::stat_compare_means()

# ---------------- 2. DEAnalysis ----------------

Idents(mye.T) <- mye.T$TissueSiteSimple

all.markers <- FindMarkers(mye.T,ident.1 = 'right',group.by = 'TissueSiteSimple',assay = 'RNA')
# remove ribosomal and mitochondrial genes
all.markers <- all.markers[!grepl("^RP[SL]",rownames(all.markers)),]
all.markers <- all.markers[!grepl("^MT-",rownames(all.markers)),]
all.markers$Significance <- ifelse(all.markers$p_val < 0.01, TRUE, FALSE)
all.markers <- all.markers[order(all.markers[,2],decreasing = T),]
top10 <- rbind(head(all.markers,10),tail(all.markers,10))

mytheme <- theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 16),
                   axis.text = element_text(size = 12),
                   legend.title = element_text(size = 12),
                   legend.text = element_text(size = 12)) 

p.all.markers <- ggplot(all.markers, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')

ggsave(filename = 'wclc/all.markers.volcano.pdf',p.all.markers,width = 5,height = 5)

up.gene <- rownames(all.markers[all.markers$avg_logFC > 0,])
dn.gene <- rownames(all.markers[all.markers$avg_logFC < 0,])

write(up.gene,file = 'wclc/all.RCC.up.genes.txt')
write(dn.gene,file = 'wclc/all.LCC.up.genes.txt')

save(all.markers, file = 'wclc/all.markers.rda')

Idents(mye.T) <- mye.T$clMidwayPr
# mono
markers.mono <- FindMarkers(mye.T,ident.1 = 'right',group.by = 'TissueSiteSimple',assay = 'RNA',subset.ident = 'Mono')
# remove ribosomal and mitochondrial genes
markers.mono <- markers.mono[!grepl("^RP[SL]",rownames(markers.mono)),]
markers.mono <- markers.mono[!grepl("^MT-",rownames(markers.mono)),]
markers.mono$Significance <- ifelse(markers.mono$p_val < 0.01, TRUE, FALSE)
markers.mono <- markers.mono[order(markers.mono[,2],decreasing = T),]
top10 <- rbind(head(markers.mono,10),tail(markers.mono,10))
p.all.markers <- ggplot(markers.mono, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')

ggsave(filename = 'wclc/markers.mono.volcano.pdf',p.all.markers,width = 5,height = 5)

save(markers.mono, file = 'wclc/markers.mono.rda')

# Macro
markers.macro <- FindMarkers(mye.T,ident.1 = 'right',group.by = 'TissueSiteSimple',assay = 'RNA',subset.ident = 'Macro')
# remove ribosomal and mitochondrial genes
markers.macro <- markers.macro[!grepl("^RP[SL]",rownames(markers.macro)),]
markers.macro <- markers.macro[!grepl("^MT-",rownames(markers.macro)),]
markers.macro$Significance <- ifelse(markers.macro$p_val < 0.01, TRUE, FALSE)
markers.macro <- markers.macro[order(markers.macro[,2],decreasing = T),]
top10 <- rbind(head(markers.macro,10),tail(markers.macro,10))
p.all.markers <- ggplot(markers.macro, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#C72719","#BEBBBD")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')

ggsave(filename = 'wclc/markers.macro.volcano.pdf',p.all.markers,width = 5,height = 5)

save(markers.macro, file = 'wclc/markers.macro.rda')

# Granulo
markers.gra <- FindMarkers(mye.T,ident.1 = 'right',group.by = 'TissueSiteSimple',assay = 'RNA',subset.ident = 'Granulo')
# remove ribosomal and mitochondrial genes
markers.gra <- markers.gra[!grepl("^RP[SL]",rownames(markers.gra)),]
markers.gra <- markers.gra[!grepl("^MT-",rownames(markers.gra)),]
markers.gra$Significance <- ifelse(markers.gra$p_val < 0.01, TRUE, FALSE)
markers.gra <- markers.gra[order(markers.gra[,2],decreasing = T),]
top10 <- rbind(head(markers.gra,10),tail(markers.gra,10))
p.all.markers <- ggplot(markers.gra, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')

ggsave(filename = 'wclc/markers.gra.volcano.pdf',p.all.markers,width = 5,height = 5)
save(markers.gra, file = 'wclc/markers.gra.rda')

# DC
markers.DC <- FindMarkers(mye.T,ident.1 = 'right',group.by = 'TissueSiteSimple',assay = 'RNA',subset.ident = 'DC')
# remove ribosomal and mitochondrial genes
markers.DC <- markers.DC[!grepl("^RP[SL]",rownames(markers.DC)),]
markers.DC <- markers.DC[!grepl("^MT-",rownames(markers.DC)),]
markers.DC$Significance <- ifelse(markers.DC$p_val < 0.01, TRUE, FALSE)
markers.DC <- markers.DC[order(markers.DC[,2],decreasing = T),]
top10 <- rbind(head(markers.DC,10),tail(markers.DC,10))
p.all.markers <- ggplot(markers.DC, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme +
  theme(legend.position = 'none')

ggsave(filename = 'wclc/markers.DC.volcano.pdf',p.all.markers,width = 5,height = 5)
save(markers.DC, file = 'wclc/markers.DC.rda')

# heatmap
load('wclc/all.markers.rda')
mat <- GetAssayData(mye.T, assay = "RNA", slot = "scale.data")
mat[mat > 2] = 2
mat[mat < -2] = -2
cls.info <- sort(mye.T$TissueSiteSimple)
gene.match <- intersect(rownames(all.markers),rownames(mat))
mat <- as.matrix(mat[gene.match, names(cls.info)])
gene <- c('PPBP','CCL23','MEPE','IL1A','ZMIZ1','PF4V1','CCL3','SPATA5L1','FRMD8','CLEC7A', 
          'KIAA1683', 'ZNF335', 'DAPK1', 'MMP3', 'TBC1D3D', 'SNRNP27', 'IFITM5',
          'TMEM121','IGLJ3','BTK','CST9','CCDC152','CYP2S1','LYZ','ENAM', 'FBL', 'FLT3LG')
gene <- unique(gene)
gene.pos <- which(rownames(mat) %in% gene)
row.anno <- ComplexHeatmap::rowAnnotation(gene = ComplexHeatmap::anno_mark(at = gene.pos, labels = rownames(mat)[which(rownames(mat) %in% gene)]))
cols <- c("#4E9BD3","#DB4840","#DF78C1")
names(cols) <- levels(cls.info)
top.anno <- ComplexHeatmap::HeatmapAnnotation(
  cluster = ComplexHeatmap::anno_block(gp = grid::gpar(fill=cols), labels = levels(cls.info), labels_gp = grid::gpar(cex=0.5,col='white'))
)
col.fun <- circlize::colorRamp2(seq(min(mat),max(mat),length=3), c('#377EB8','white','#E41A1C'))
pdf("wclc/diff.cls.heatmap.pdf",width = 10,height = 7)
ComplexHeatmap::Heatmap(matrix = mat, cluster_rows = F, cluster_columns = F,
                        show_column_names = F, show_row_names = F,
                        column_split = cls.info, right_annotation = row.anno,
                        column_title = NULL, top_annotation = top.anno,
                        heatmap_legend_param = list(title = 'Expression', title_position = 'leftcenter-rot'),col = col.fun)
dev.off()

rm(list = ls());gc()

# DEPTH <- function(exp, match){
#   #Two files need to be input into this function.
#   exp<-as.data.frame(exp)
#   match <- as.data.frame(match)
#   
#   pos_nor <- which(colnames(exp) %in% match[which(match[, 2] == "N"), 1])
#   pos_dis <- which(colnames(exp) %in% match[which(match[, 2] == "T"), 1])
#   
#   exp_nor <- exp[, pos_nor]
#   exp_dis <- exp[, pos_dis]
#   
#   #Caculate the average values of each gene in normal/disease sample.
#   mean_nor <- as.data.frame(rowMeans(exp_nor,na.rm=T))
#   mean_dis <- as.data.frame(rowMeans(exp_dis,na.rm=T))
#   
#   score <- as.data.frame(matrix(0, nrow <- dim(exp_dis)[1], ncol <- dim(exp_dis)[2]))
#   #Caculate the heterogeneity score of each gene.
#   if(length(pos_nor) > 0){
#     for(u in 1 : dim(exp_dis)[2]){
#       score[,u] <- (exp_dis[,u] - mean_nor[,1])^2
#     }
#   }else if(length(pos_nor) == 0){
#     for(u in 1 : dim(exp_dis)[2]){
#       score[,u] <- (exp_dis[,u] - mean_dis[,1])^2
#     }
#   }
#   
#   colnames(score) <- colnames(exp_dis)
#   rownames(score) <- rownames(exp_dis)
#   
#   heterogeneity_score <- c()
#   
#   for(z in 1:length(colnames(exp_dis))){
#     heterogeneity_score[z] <- sd(score[, z],na.rm=T)
#   }
#   heterogeneity_score <- cbind(colnames(exp_dis), heterogeneity_score)
#   #Caculate the heterogeneity score of each sample.
#   colnames(heterogeneity_score) <- c("Sample", "ITH score")
#   #DEPTH function will output the heterogeneity score of each tumor sample.
#   return(heterogeneity_score)
# }
# mye <- readRDS('wclc/rds/myeloid.total.rds')
# expr <- mye@assays$RNA@data
# identifier <- mye@meta.data[,c('orig.ident','SPECIMEN_TYPE')]
# identifier$orig.ident <- rownames(identifier)
# ITHscore <- DEPTH(expr,identifier)
# load('wclc/rds/ITHscore.rda')
# ITHscore <- data.frame(ITHscore,stringsAsFactors = F)
# ITHscore$ITH.score <- as.numeric(ITHscore$ITH.score)
# identifier <- mye@meta.data[,c('TissueSiteSimple','SPECIMEN_TYPE','orig.ident')]
# identifier <- identifier[identifier$SPECIMEN_TYPE == 'T',]
# identifier$TissueSiteSimple <- as.character(identifier$TissueSiteSimple)
# identifier$orig.ident <- as.character(identifier$orig.ident)
# identifier$SPECIMEN_TYPE <- as.character(identifier$SPECIMEN_TYPE)
# ITHscore$TissueSite <- plyr::mapvalues(x = ITHscore$Sample, 
#                                        from = row.names(identifier), 
#                                        to = identifier$TissueSiteSimple)
# ITHscore$Patient <- plyr::mapvalues(x = ITHscore$Sample, 
#                                        from = row.names(identifier), 
#                                        to = identifier$orig.ident)
# ITHscore <- ITHscore[order(ITHscore$Patient,decreasing = F),]
# sum.T <- data.frame(tapply(ITHscore$ITH.score, ITHscore$Patient, median))
# colnames(sum.T) <- 'ITH score'
# sum.T$TissueSiteSimple <- plyr::mapvalues(rownames(sum.T),
#                                           from = identifier$orig.ident,
#                                           to = identifier$TissueSiteSimple)
# ggpubr::ggboxplot(sum.T,x = 'TissueSiteSimple',y = 'ITH score',
#                   add = 'jitter',color = 'TissueSiteSimple') +
#   ggpubr::stat_compare_means(comparisons = list(c('left','right')))


# ---------------- 3. Enrichment ----------------

EnrichmentAnalysis <- function(DEgenes, enrich.method = c("GO","KEGG")) {
  
  pacman::p_load(clusterProfiler,org.Hs.eg.db,dplyr,patchwork,enrichplot,ggplot2)
  
  enrich.method <- match.arg(enrich.method)
  
  if(enrich.method == "GO") {
    # ego_ALL <- enrichGO(gene = rownames(DEgenes), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL",
    #                     pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    # ego_CC <- enrichGO(gene = DEgenes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC",
    #                    pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    # ego_MF <- enrichGO(gene = DEgenes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF",
    #                    pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_BP <- enrichGO(gene = DEgenes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    # ego_CC@result$Description <- substring(ego_CC@result$Description, 1, 70)
    ego_BP@result$Description <- substring(ego_BP@result$Description, 1, 70)
    # ego_MF@result$Description <- substring(ego_MF@result$Description, 1, 70)
    # p_CC <- barplot(ego_BP, showCategory = 10, title = "Cellular Component")
    # p_BP <- barplot(ego_CC, showCategory = 10, title = "Biological Process")
    # p_MF <- barplot(ego_MF, showCategory = 10, title = "Molecular Function")
    # pc <- p_BP | p_CC | p_MF
    # ggsave("GO.pdf",plot = pc, device = "pdf",width = 12, height = 10)
    return(ego_BP)
    # ego_results <- list(ego_ALL = ego_ALL, ego_BP = ego_BP, ego_CC = ego_CC, ego_MF = ego_MF)
    return(ego_results)
  } else if(enrich.method == "KEGG") {
    genelist <- bitr(DEgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    genelist <- pull(genelist, ENTREZID)
    ekegg <- enrichKEGG(gene = genelist, organism = "hsa")
    p <- dotplot(ekegg, showCategory = 20)
    ggsave("KEGG.pdf",plot = p, device = "pdf", width = 12, height = 10)
    return(ekegg)
  }
}

# granulocytes
load('wclc/rds/markers/markers.gra.rda')
up <- rownames(markers.gra[markers.gra$avg_logFC > 0,])
dn <- rownames(markers.gra[markers.gra$avg_logFC < 0,])
gra.up <- EnrichmentAnalysis(DEgenes = up, enrich.method = 'GO')
openxlsx::write.xlsx(gra.up$ego_BP@result, file = 'wclc/gra.rcc.up.xlsx')
gra.dn <- EnrichmentAnalysis(DEgenes = dn, enrich.method = 'GO')
edo <- enrichGO(gene = dn, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
p <- cnetplot(edo)
ggsave(filename = "wclc/gra/lcc.network.pdf",width = 6,height = 6)
openxlsx::write.xlsx(gra.dn$ego_BP@result, file = 'wclc/gra.rcc.dn.xlsx')
# monocytes
load('wclc/markers.mono.rda')
up <- rownames(markers.mono[markers.mono$avg_logFC > 0,])
mono.up <- EnrichmentAnalysis(DEgenes = up, enrich.method = 'GO')
openxlsx::write.xlsx(mono.up$ego_BP@result, file = 'wclc/mono.rcc.up.xlsx')
dn <- rownames(markers.mono[markers.mono$avg_logFC < 0,])
mono.dn <- EnrichmentAnalysis(DEgenes = dn, enrich.method = 'GO')
openxlsx::write.xlsx(mono.dn$ego_BP@result, file = 'wclc/mono.rcc.dn.xlsx')
# marcophages
load('wclc/markers.macro.rda')
up <- rownames(markers.macro[markers.macro$avg_logFC > 0,])
macro.up <- EnrichmentAnalysis(DEgenes = up, enrich.method = 'GO')
openxlsx::write.xlsx(macro.up$ego_BP@result, file = 'wclc/macro.rcc.up.xlsx')
dn <- rownames(markers.macro[markers.macro$avg_logFC < 0,])
macro.dn <- EnrichmentAnalysis(DEgenes = dn, enrich.method = 'GO')
openxlsx::write.xlsx(macro.dn$ego_BP@result, file = 'wclc/macro.rcc.dn.xlsx')
# DC
load('wclc/markers.DC.rda')
up <- rownames(markers.DC[markers.DC$avg_logFC > 0,])
DC.up <- EnrichmentAnalysis(DEgenes = up, enrich.method = 'GO')
openxlsx::write.xlsx(DC.up$ego_BP@result, file = 'wclc/DC.rcc.up.xlsx')
dn <- rownames(markers.DC[markers.DC$avg_logFC < 0,])
DC.dn <- EnrichmentAnalysis(DEgenes = dn, enrich.method = 'GO')
openxlsx::write.xlsx(DC.dn$ego_BP@result, file = 'wclc/DC.rcc.dn.xlsx')


# ---------------- 4. Extra dataset ----------------



# ---------------- 5. LRAnalysis ----------------

library(Seurat)
mye <- readRDS('wclc/rds/myeloid.tumor.sub.rds')
VlnPlot(mye, features = 'DYSF',group.by = "clMidwayPr")

LCC.mye <- subset(mye, cells = rownames(subset(mye@meta.data,TissueSiteSimple=="left")))
RCC.mye <- subset(mye, cells = rownames(subset(mye@meta.data,TissueSiteSimple=="right")))
saveRDS(LCC.mye,file = 'wclc/rds/LCC.mye.rds')
saveRDS(RCC.mye,file = 'wclc/rds/RCC.mye.rds')

rm(list = ls())

cellCommunication <- function(seurat.obj, group.by, out.dir, sampleName){
  cellchat <- createCellChat(object = seurat.obj, group.by = group.by)
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = T)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  df.net <- subsetCommunication(cellchat, slot.name = 'netP')
  openxlsx::write.xlsx(df.net, file = paste0(out.dir, sampleName, '.lr.net.xlsx'))
  cellchat <- aggregateNet(cellchat)
  df.net.2 <- subsetCommunication(cellchat)
  openxlsx::write.xlsx(df.net.2, file = paste0(out.dir, sampleName, '.lr.xlsx'))
  groupsize <- as.numeric(table(cellchat@idents))
  pdf(file = paste0(out.dir, sampleName, '.net.pdf'))
  netVisual_circle(cellchat@net$count, vertex.weight = groupsize, weight.scale = T,
                   label.edge = F, title.name = 'Number of interactions')
  netVisual_circle(cellchat@net$weight, vertex.weight = groupsize,weight.scale = T,
                   label.edge = F, title.name = 'Interaction weights/strength')
  dev.off()
  save(cellchat, file = paste0(out.dir, sampleName, '.cellchat.rda'))
  
  return(cellchat)
}

# right-sided crc
RCC.mye <- readRDS('wclc/RCC.mye.rds')
rcc.cellchat <- cellCommunication(seurat.obj = RCC.mye, group.by = 'clMidwayPr', out.dir = 'wclc/',sampleName = 'rcc')

mat <- rcc.cellchat@net$count 
par(mfrow=c(3,3),xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = as.numeric(table(rcc.cellchat@idents)),weight.scale = T,arrow.width = 0.2,arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

mat <- rcc.cellchat@net$weight
par(mfrow=c(3,3),xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = as.numeric(table(rcc.cellchat@idents)),weight.scale = T,arrow.width = 0.2,arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
rcc.cellchat@netP$pathways
pathways.show <- c("CCL")
levels(rcc.cellchat@idents)
pdf(file = 'wclc/cellchat/rcc.vis.pdf',width = 6,height = 6)
par(mfrow = c(1,1))
netVisual_aggregate(rcc.cellchat,signaling = pathways.show,layout = "circle")
netVisual_aggregate(rcc.cellchat,signaling = pathways.show,layout = "chord")
netVisual_heatmap(rcc.cellchat,signaling = pathways.show, color.heatmap = "Reds")
dev.off()

pairLR.CCL <- extractEnrichedLR(rcc.cellchat,signaling = pathways.show,geneLR.return = FALSE)
#vertex.receiver = c(0,7,12,16)
pdf(file = 'wclc/cellchat/rcc.CCL23.pdf',width = 5,height = 5)
netVisual_individual(rcc.cellchat,signaling = pathways.show,pairLR.use = pairLR.CCL[2,],layout = "circle")
dev.off()

rm(list = ls());gc()

# left-sided crc
LCC.mye <- readRDS('wclc/LCC.mye.rds')
lcc.cellchat <- cellCommunication(seurat.obj = LCC.mye, group.by = 'clMidwayPr', out.dir = 'wclc/',sampleName = 'lcc')
mat <- lcc.cellchat@net$count 
#mat <- cellchat2@net$weight #如果想展示互作的强度/概率图
par(mfrow=c(3,3),xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = as.numeric(table(lcc.cellchat@idents)),weight.scale = T,arrow.width = 0.2,arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

mat <- lcc.cellchat@net$weight
par(mfrow=c(3,3),xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = as.numeric(table(lcc.cellchat@idents)),weight.scale = T,arrow.width = 0.2,arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
lcc.cellchat@netP$pathways
pathways.show <- c("CCL")
levels(lcc.cellchat@idents)
pdf(file = 'wclc/cellchat/lcc.vis.pdf',width = 6,height = 6)
par(mfrow = c(1,1))
netVisual_aggregate(lcc.cellchat,signaling = pathways.show,layout = "circle")
netVisual_aggregate(lcc.cellchat,signaling = pathways.show,layout = "chord")
netVisual_heatmap(lcc.cellchat,signaling = pathways.show, color.heatmap = "Reds")
dev.off()

pairLR.CCL <- extractEnrichedLR(lcc.cellchat,signaling = pathways.show,geneLR.return = FALSE)
#vertex.receiver = c(0,7,12,16)
pdf(file = 'wclc/cellchat/lcc.CCL23.pdf',width = 5,height = 5)
netVisual_individual(lcc.cellchat,signaling = pathways.show,pairLR.use = pairLR.CCL[1,],layout = "circle")
dev.off()

rm(list = ls());gc()

load("wclc/cellchat/rcc.cellchat.rda")
rcc.cellchat <- cellchat
p <- netVisual_bubble(rcc.cellchat,signaling = "CCL")
ggsave(filename = "wclc/cellchat/rcc.CCL.bubble.pdf",p,width = 8,height = 3)

load("wclc/cellchat/lcc.cellchat.rda")
lcc.cellchat <- cellchat
p <- netVisual_bubble(lcc.cellchat,signaling = "CCL")
ggsave(filename = "wclc/cellchat/lcc.CCL.bubble.pdf",p,width = 8,height = 3)

pairs <- extractEnrichedLR(rcc.cellchat, signaling = "CCL", geneLR.return = F)
p <- netVisual_individual(rcc.cellchat, signaling = "CCL",pairLR.use = pairs[2,], layout = 'circle')



# ---------------- 6. ReallocatedALLCelltypes ----------------

mye <- readRDS('wclc/rds/myeloid.tumor.sub.rds')
gra <- readRDS('wclc/rds/gra.rds')
DC <- readRDS('wclc/rds/DC.rds')
macro <- readRDS('wclc/rds/macro.rds')
mono <- readRDS('wclc/rds/mono.rds')

ct.gra <- gra@meta.data[,c("orig.ident","celltype")]
ct.DC <- DC@meta.data[,c("orig.ident","celltype")]
ct.macro <- macro@meta.data[,c("orig.ident","celltype")]
ct.mono <- mono@meta.data[,c("orig.ident","celltype")]

celltype <- rbind(ct.gra,ct.DC,ct.macro,ct.mono)
rm(gra,DC,macro,mono,ct.gra,ct.DC,ct.macro,ct.mono);gc()

mye$celltype <- plyr::mapvalues(colnames(mye),from = rownames(celltype),to = celltype$celltype)
saveRDS(mye,file = 'wclc/rds/myeloid.tumor.sub.rds')

mypal <- c('#0C7D47','#6D9BB2','#E5498D','#E24132','#13AFCE',
           '#39568A','#E9967F','#8091B4','#86BCB0','#87684B',
           '#EA6C24','#04899E','#7F3B83')
p <- DimPlot(mye,cols = mypal,reduction = 'tsne',group.by = 'celltype',pt.size = 0.8) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank())
ggsave("wclc/myeloid/celltype.pdf",p,width = 7,height = 4.8)


library(Seurat)
crc <- readRDS('wclc/rds/crc.dr.rds')

mye <- readRDS('wclc/rds/myeloid.tumor.sub.rds')
mye.meta <- mye@meta.data[,c('clTopLevel','clMidwayPr','cl295v11SubFull','celltype')]
rm(mye)
B <- readRDS('wclc/rds/B.tumor.sub.rds')
old.ident <- names(table(as.character(B$cl295v11SubFull)))
new.ident <- c('IGD+IgM+B',"GC-like B","CD40+ GC-like B")
B$celltype <- plyr::mapvalues(B$cl295v11SubFull, from = old.ident, to = new.ident)
B.meta <- B@meta.data[,c('clTopLevel','clMidwayPr','cl295v11SubFull','celltype')]

mast <- readRDS('wclc/rds/Mast.tumor.sub.rds')
old.ident <- names(table(as.character(mast$cl295v11SubFull)))
new.ident <- c('Mast')
mast$celltype <- plyr::mapvalues(mast$cl295v11SubFull, from = old.ident, to = new.ident)
mast.meta <- mast@meta.data[,c('clTopLevel','clMidwayPr','cl295v11SubFull','celltype')]

plasma <- readRDS('wclc/rds/plasma.tumor.sub.rds')
old.ident <- names(table(as.character(plasma$cl295v11SubFull)))
new.ident <- c('IgA+Plasma',"IgG+Plasma","IgG+ prolif Plasma")
plasma$celltype <- plyr::mapvalues(plasma$cl295v11SubFull, from = old.ident, to = new.ident)
plasma.meta <- plasma@meta.data[,c('clTopLevel','clMidwayPr','cl295v11SubFull','celltype')]

strom <- readRDS('wclc/rds/stromal.tumor.sub.rds')
old.ident <- names(table(as.character(strom$cl295v11SubFull)))
new.ident <- c('Arterial Endothelial',"Capillary Endothelial","Capillary Endothelial",
               "Endothelial","Venous Endothelial", "Lymphatic Endothelial",
               "Capillary-like Endothelial", 'Arterial-like Endothelial', 'Endothelial',
               'Endothelial tip cell',"prolif Endothelial","Endothelial",
               "Venous-like Endothelial","Endothelial","Pericyte",
               "Pericyte", "Pericyte", "Pericyte",
               "Pericyte", "prolif Pericyte", "stem-cell-niche Fibroblast",
               'stem-cell-niche Fibroblast',"BMP-producing Fibroblast","BMP-producing Fibroblast",
               "CCL8+ Fibroblast", "Myofibroblast","CXCL14+CAF",
               "GREM1+CAF", "MMP3+CAF", "CCL8+CAF",
               "stem-cell-niche Fibroblast", "Smooth Muscle", "Schwann")
strom$celltype <- plyr::mapvalues(strom$cl295v11SubFull, from = old.ident, to = new.ident)
strom.meta <- strom@meta.data[,c('clTopLevel','clMidwayPr','cl295v11SubFull','celltype')]

ILC <- readRDS('wclc/rds/TNKILC.tumor.sub.rds')
old.ident <- names(table(as.character(ILC$cl295v11SubFull)))
new.ident <- c('CD4+IL7R+T',"CD4+IL7R+SELL+T","CD4+IL7R+HSP+T",
               "CD4+IL7R+CCL5+T", "CD4+IL17+T", "CD4+Tfh",
               "CD4+CXCL13+T", "CD4+Treg", "prolif CD4+Treg",
               "CD8+IL7R+T", "CD8+GZMK+T", "CD8+IL7R+T",
               "CD8+IL17+T", "CD8+CXCL13+T", "CD8+CXCL13+HSP+T",
               "prolif CD8+CXCL13+T", "gd-like T", "gd-like PDCD1+T",
               "prolif gd-like T", "PLZF+T", "prolif PLZF+T",
               "unkT", "CD16A+NK", "GZMK+NK",
               "XCL1+NK", "ILC3")
ILC$celltype <- plyr::mapvalues(ILC$cl295v11SubFull, from = old.ident, to = new.ident)
ILC.meta <- ILC@meta.data[,c('clTopLevel','clMidwayPr','cl295v11SubFull','celltype')]

Epi <- readRDS('wclc/rds/Epi.tumor.sub.rds')
old.ident <- names(table(as.character(Epi$cl295v11SubFull)))
new.ident <- c('TA-like Stem',"TA-like Immature Goblet Stem",
               "prolif TA-like Stem", "Enterocyte 1",
               "Enterocyte 2", "Immature Goblet",
               "Goblet Enterocyte", "Goblet",
               'Best4', "Tuft", "Enteroendocrine")
Epi$celltype <- plyr::mapvalues(Epi$cl295v11SubFull, from = old.ident, to = new.ident)
Epi.meta <- Epi@meta.data[,c('clTopLevel','clMidwayPr','cl295v11SubFull','celltype')]

crc.tumor.meta <- rbind(B.meta,Epi.meta,ILC.meta,mast.meta,mye.meta,plasma.meta,strom.meta)
crc.tumor.meta$celltype <- as.character(crc.tumor.meta$celltype)
crc.tumor.meta$celltype <- gsub("unkT","unkNK",crc.tumor.meta$celltype)
crc.tumor.meta$clMidwayPr <- as.character(crc.tumor.meta$clMidwayPr)
old.ident <- names(table(crc.tumor.meta$clMidwayPr))
new.ident <- c('B','DC','Endo','Epi','Fibro','Granulo','TNKILC','Macro','Mast','Mono','TNKILC',
               'Peri','Plasma','Schwann','SmoothMuscle','TNKILC','TNKILC','TNKILC','TNKILC')
crc.tumor.meta$celltype2 <- plyr::mapvalues(crc.tumor.meta$clMidwayPr,from = old.ident,to = new.ident)
saveRDS(crc.tumor.meta,file = 'wclc/rds/crc.tumor.meta.rds')

## modify metadata
crc.tumor.meta <- readRDS('wclc/rds/crc.tumor.meta.rds')
crc.tumor.meta.full <- readRDS('wclc/rds/crc.tumor.meta.full.rds')
crc.tumor.meta.full$TissueSiteSimple <- as.character(crc.tumor.meta.full$TissueSiteSimple)
crc.tumor.meta$TissueSiteSimple <- plyr::mapvalues(rownames(crc.tumor.meta),
                                                   from = rownames(crc.tumor.meta.full),
                                                   to = crc.tumor.meta.full$TissueSiteSimple)
crc.tumor.meta$Tissue <- 'Tumor'
saveRDS(crc.tumor.meta, file = 'wclc/rds/crc.tumor.meta.modified.rds')

crc.normal.meta.full <- readRDS('wclc/rds/crc.normal.meta.full.rds')
crc.normal.meta <- crc.normal.meta.full[,c('clTopLevel','clMidwayPr','cl295v11SubFull','TissueSiteSimple')]
crc.normal.meta$celltype <- plyr::mapvalues(crc.normal.meta$cl295v11SubFull, 
                                            from = crc.tumor.meta$cl295v11SubFull,
                                            to = crc.tumor.meta$celltype)
crc.normal.meta$celltype <- as.character(crc.normal.meta$celltype)
crc.normal.meta$celltype2 <- plyr::mapvalues(crc.normal.meta$clMidwayPr,
                                             from = crc.tumor.meta$clMidwayPr,
                                             to = crc.tumor.meta$celltype2)
crc.normal.meta$Tissue <- 'Normal'
saveRDS(crc.normal.meta, file = 'wclc/rds/crc.normal.meta.modified.rds')

# crc.meta <- rbind(crc.tumor.meta, crc.normal.meta)
# crc.meta$celltype <- as.character(crc.meta$celltype)
# saveRDS(crc.meta, file = 'wclc/rds/crc.meta.modified.rds')

library(ggpubr)
# process normal data
normal.tbl <- table(crc.normal.meta$celltype, crc.normal.meta$TissueSiteSimple)
write.table(normal.tbl,file = 'tmp.txt',col.names = T,row.names = T,sep = '\t',quote = F)
normal.tbl <- read.table('tmp.txt',header = T,row.names = 1,sep = '\t')
normal.tbl$left <- normal.tbl$left/nrow(crc.normal.meta)*100
normal.tbl$right <- normal.tbl$right/nrow(crc.normal.meta)*100
normal.tbl$celltype <- rownames(normal.tbl)
norm.dat <- reshape2::melt(normal.tbl)
# process tumor data
tumor.tbl <- table(crc.tumor.meta$celltype, crc.tumor.meta$TissueSiteSimple)
write.table(tumor.tbl,file = 'tmp.txt',col.names = T,row.names = T,sep = '\t',quote = F)
tumor.tbl <- read.table('tmp.txt',header = T,row.names = 1,sep = '\t')
tumor.tbl$left <- tumor.tbl$left/nrow(crc.tumor.meta)*100
tumor.tbl$right <- tumor.tbl$right/nrow(crc.tumor.meta)*100
tumor.tbl$celltype <- rownames(tumor.tbl)
tumor.dat <- reshape2::melt(tumor.tbl)

norm.dat$celltype <- factor(
  norm.dat$celltype,
  levels = c(
    # B
    'CD40+ GC-like B','GC-like B','IGD+IgM+B',
    # Epi
    'Best4','Enterocyte 1','Enterocyte 2','Enteroendocrine','Goblet','Goblet Enterocyte',
    'Immature Goblet','prolif TA-like Stem','TA-like Immature Goblet Stem','TA-like Stem','Tuft',
    # Mast
    'Mast',
    # DC
    'AS-DC','C1Q+DC2','DC1','DC2','IL22RA2+DC','mregDC','pDC',
    # Granulo
    'N1 TANs',
    # Macro
    'C1QC+TAM',
    # Mono
    'CD14+ monocytes',
    # Plasma
    'IgA+Plasma','IgG+ prolif Plasma','IgG+Plasma',
    # Endo
    'Arterial-like Endothelial','Arterial Endothelial','Capillary-like Endothelial',
    'Capillary Endothelial','Endothelial','Endothelial tip cell','Lymphatic Endothelial',
    'prolif Endothelial','Venous-like Endothelial','Venous Endothelial',
    # Fibro
    'BMP-producing Fibroblast','CCL8+ Fibroblast','CCL8+CAF','CXCL14+CAF','Myofibroblast',
    'stem-cell-niche Fibroblast',
    # Peri
    'Pericyte','prolif Pericyte',
    # Schwann
    'Schwann',
    # SmoothMuscle
    'Smooth Muscle',
    # ILC
    'ILC3',
    # NK
    'CD16A+NK','GZMK+NK','unkNK','XCL1+NK',
    # TCD4
    'CD4+CXCL13+T','CD4+IL17+T','CD4+IL7R+CCL5+T','CD4+IL7R+HSP+T','CD4+IL7R+SELL+T',
    'CD4+IL7R+T','CD4+Tfh','CD4+Treg','prolif CD4+Treg',
    # TCD8
    'CD8+CXCL13+HSP+T','CD8+CXCL13+T','CD8+GZMK+T','CD8+IL17+T','CD8+IL7R+T','prolif CD8+CXCL13+T',
    # Tgd
    'gd-like PDCD1+T','gd-like T','prolif gd-like T',
    # TZBTB16
    'PLZF+T','prolif PLZF+T'
  )
)
p1 <- ggbarplot(norm.dat,x = 'celltype',y = 'value',color = 'black',merge = T,fill = 'variable') +
  ylab('Cell Proportion (%)') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank()) +
  stat_compare_means(ref.group = '.all.',label = '..p.signif..',hide.ns = T)
ggsave('wclc/total/norm.cell.prop.pdf',p1,width = 18,height = 5)

tumor.dat$celltype <- factor(
  tumor.dat$celltype,
  levels = c(
    # B
    'CD40+ GC-like B','GC-like B','IGD+IgM+B',
    # Epi
    'Best4','Enterocyte 1','Enterocyte 2','Enteroendocrine','Goblet','Goblet Enterocyte',
    'Immature Goblet','prolif TA-like Stem','TA-like Immature Goblet Stem','TA-like Stem','Tuft',
    # Mast
    'Mast',
    # DC
    'AS-DC','C1Q+DC2','DC1','DC2','IL22RA2+DC','mregDC','pDC',
    # Granulo
    'N1 TANs','N2 TANs',
    # Macro
    'C1QC+TAM',
    # Mono
    'CD14+ monocytes','CD14+CD16+ monocytes','CD16+ monocytes',
    # Plasma
    'IgA+Plasma','IgG+ prolif Plasma','IgG+Plasma',
    # Endo
    'Arterial-like Endothelial','Arterial Endothelial','Capillary-like Endothelial',
    'Capillary Endothelial','Endothelial','Endothelial tip cell','Lymphatic Endothelial',
    'prolif Endothelial','Venous-like Endothelial','Venous Endothelial',
    # Fibro
    'BMP-producing Fibroblast','CCL8+ Fibroblast','CCL8+CAF','CXCL14+CAF','GREM1+CAF',
    'MMP3+CAF','Myofibroblast','stem-cell-niche Fibroblast',
    # Peri
    'Pericyte','prolif Pericyte',
    # Schwann
    'Schwann',
    # SmoothMuscle
    'Smooth Muscle',
    # ILC
    'ILC3',
    # NK
    'CD16A+NK','GZMK+NK','unkNK','XCL1+NK',
    # TCD4
    'CD4+CXCL13+T','CD4+IL17+T','CD4+IL7R+CCL5+T','CD4+IL7R+HSP+T','CD4+IL7R+SELL+T',
    'CD4+IL7R+T','CD4+Tfh','CD4+Treg','prolif CD4+Treg',
    # TCD8
    'CD8+CXCL13+HSP+T','CD8+CXCL13+T','CD8+GZMK+T','CD8+IL17+T','CD8+IL7R+T','prolif CD8+CXCL13+T',
    # Tgd
    'gd-like PDCD1+T','gd-like T','prolif gd-like T',
    # TZBTB16
    'PLZF+T','prolif PLZF+T'
  )
)
p2 <- ggbarplot(tumor.dat,x = 'celltype',y = 'value',color = 'black',merge = T,fill = 'variable') +
  ylab('Cell Proportion (%)') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank()) +
  stat_compare_means(ref.group = '.all.',label = '..p.signif..',hide.ns = T,method = 't.test')
ggsave('wclc/total/tumor.cell.prop.pdf',p2,width = 18,height = 5)

save(norm.dat, tumor.dat, file = 'wclc/rds/total.cell.prop.rda')

# mycols <- c(
#   "#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6",
#   "#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762",
#   "#004D43","#8FB0FF","#997D87","#5A0007","#809693","#6A3A4C",
#   "#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#61615A",
#   "#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA",
#   "#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018",
#   "#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED",
#   "#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C",
#   "#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1",
#   "#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459",
#   "#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA",
#   "#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329",
#   "#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1",
#   "#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C",
#   "#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625",
#   "#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534",
#   "#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72"
# )
# p3 <- ggplot(crc.tumor.meta, aes(x = TissueSiteSimple, fill = celltype2)) +
#   geom_bar(position = "fill") + theme_bw() + 
#   scale_fill_manual(values = mycols[1:length(table(crc.tumor.meta$celltype2))]) +
#   mytheme +
#   guides(fill = guide_legend(ncol = 3)) +
#   labs(y = "Cell Fraction")
# crc.tumor <- readRDS('wclc/rds/crc.tumor.rds')
# a <- table(crc.tumor$orig.ident, crc.tumor$TissueSiteSimple)
# b <- table(crc.tumor$orig.ident, crc.tumor$celltype2)
# c <- cbind(a,b)
# c <- read.table("wclc/cell.type.count.txt",sep = '\t',header = T,row.names = 1)
# c$LMR <- c$TNKILC/c$Mono
# c$NLR <- c$Granulo/c$TNKILC
# c$Tissue <- ifelse(c$left > 0, "LCC","RCC")
# c$Patient <- rownames(c)

### Cell communication between myeloid subsets
library(Seurat)
library(CellChat)
mye <- readRDS('wclc/rds/myeloid.tumor.sub.rds')

cellCommunication <- function(seurat.obj, group.by, out.dir, sampleName){
  cellchat <- createCellChat(object = seurat.obj, group.by = group.by)
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = T)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  df.net <- subsetCommunication(cellchat, slot.name = 'netP')
  openxlsx::write.xlsx(df.net, file = paste0(out.dir, sampleName, '.lr.net.xlsx'))
  cellchat <- aggregateNet(cellchat)
  df.net.2 <- subsetCommunication(cellchat)
  openxlsx::write.xlsx(df.net.2, file = paste0(out.dir, sampleName, '.lr.xlsx'))
  groupsize <- as.numeric(table(cellchat@idents))
  pdf(file = paste0(out.dir, sampleName, '.net.pdf'))
  netVisual_circle(cellchat@net$count, vertex.weight = groupsize, weight.scale = T,
                   label.edge = F, title.name = 'Number of interactions')
  netVisual_circle(cellchat@net$weight, vertex.weight = groupsize,weight.scale = T,
                   label.edge = F, title.name = 'Interaction weights/strength')
  dev.off()
  save(cellchat, file = paste0(out.dir, sampleName, '.cellchat.rda'))
  return(cellchat)
}

LCC.mye <- subset(mye, cells = rownames(subset(mye@meta.data,TissueSiteSimple=="left")))
RCC.mye <- subset(mye, cells = rownames(subset(mye@meta.data,TissueSiteSimple=="right")))
rcc.cellchat <- cellCommunication(
  seurat.obj = RCC.mye, group.by = 'celltype', 
  out.dir = 'wclc/myeloid/cellcommu/rcc/',sampleName = 'rcc'
)
lcc.cellchat <- cellCommunication(
  seurat.obj = LCC.mye, group.by = 'celltype', 
  out.dir = 'wclc/myeloid/cellcommu/lcc/',sampleName = 'lcc'
)

## **************** Append  ****************

library(CellChat)
load('wclc/myeloid/cellcommu/lcc/lcc.cellchat.rda')
groupsize <- as.numeric(table(cellchat@idents))
pdf(file = paste0('wclc/myeloid/cellcommu/lcc/lcc.weighted.net.pdf'))
netVisual_circle(cellchat@net$weight, vertex.weight = groupsize,weight.scale = T,
                 label.edge = F, title.name = 'Interaction weights/strength')
dev.off()

load('wclc/myeloid/cellcommu/rcc/rcc.cellchat.rda')
groupsize <- as.numeric(table(cellchat@idents))
pdf(file = paste0('wclc/myeloid/cellcommu/rcc/rcc.weighted.net.pdf'))
netVisual_circle(cellchat@net$weight, vertex.weight = groupsize,weight.scale = T,
                 label.edge = F, title.name = 'Interaction weights/strength')
dev.off()

## ******************************************


lcc.path <- xlsx::read.xlsx("wclc/myeloid/cellcommu/lcc/lcc.lr.net.xlsx",sheetIndex = 1)
lcc.path$overall <- paste(lcc.path$source,lcc.path$target,lcc.path$pathway_name,sep = '-')
lcc.path$pathway_name <- as.character(lcc.path$pathway_name)
rcc.path <- xlsx::read.xlsx("wclc/myeloid/cellcommu/rcc/rcc.lr.net.xlsx",sheetIndex = 1)
rcc.path$overall <- paste(rcc.path$source,rcc.path$target,rcc.path$pathway_name,sep = '-')
rcc.path$pathway_name <- as.character(rcc.path$pathway_name)

lcc.path.uniq <- lcc.path[lcc.path$pathway_name==setdiff(lcc.path$pathway_name,rcc.path$pathway_name),]
lcc.path.uniq$`log10(prob)` <- log10(lcc.path.uniq$prob)

uniq.path <- setdiff(rcc.path$pathway_name,lcc.path$pathway_name)
rcc.path.uniq <- rcc.path[rcc.path$pathway_name==setdiff(rcc.path$pathway_name,lcc.path$pathway_name),]
rcc.path.uniq <- c()
for(i in 1:length(uniq.path)){
  path.uniq.i <- rcc.path[rcc.path$pathway_name==uniq.path[i],]
  rcc.path.uniq <- rbind(rcc.path.uniq,path.uniq.i)
}

rcc.path.uniq$pairs <- paste(rcc.path.uniq$source,rcc.path.uniq$target,sep = '->')
rcc.path.uniq$`log10(prob)` <- log10(rcc.path.uniq$prob)

mycol <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(6))
p <- ggplot(lcc.path.uniq, aes(x=source,y=target,size=`log10(prob)`,color=`log10(prob)`)) +
  geom_point() + scale_color_gradientn(colors = mycol) +
  theme_minimal() + ggtitle('LCC specific pathway - BTLA') +
  theme(axis.text = element_text(size = 14,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5,size = 16,face = 'bold'))
ggsave('wclc/myeloid/cellcommu/lcc.uniq.BTLA.path.pdf',p,width = 8,height = 6)
# EGF
rcc.EGF <- rcc.path.uniq[rcc.path.uniq$pathway_name=='EGF',]
p <- ggplot(rcc.EGF, aes(x=source,y=target,size=`log10(prob)`,color=`log10(prob)`)) +
  geom_point() + scale_color_gradientn(colors = mycol) +
  theme_minimal() + ggtitle('RCC specific pathway - EGF') +
  theme(axis.text = element_text(size = 14,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5,size = 16,face = 'bold'))
ggsave('wclc/myeloid/cellcommu/rcc.uniq.EGF.path.pdf',p,width = 5,height = 5.6)
# EPO
rcc.EPO <- rcc.path.uniq[rcc.path.uniq$pathway_name=='EPO',]
p <- ggplot(rcc.EPO, aes(x=source,y=target,size=`log10(prob)`,color=`log10(prob)`)) +
  geom_point() + scale_color_gradientn(colors = mycol) +
  theme_minimal() + ggtitle('RCC specific pathway - EPO') +
  theme(axis.text = element_text(size = 14,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5,size = 16,face = 'bold'))
ggsave('wclc/myeloid/cellcommu/rcc.uniq.EPO.path.pdf',p,width = 6.5,height = 5.8)
# IGF
rcc.IGF <- rcc.path.uniq[rcc.path.uniq$pathway_name=='IGF',]
p <- ggplot(rcc.IGF, aes(x=source,y=target,size=`log10(prob)`,color=`log10(prob)`)) +
  geom_point() + scale_color_gradientn(colors = mycol) +
  theme_minimal() + ggtitle('RCC specific pathway - IGF') +
  theme(axis.text = element_text(size = 14,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5,size = 16,face = 'bold'))
ggsave('wclc/myeloid/cellcommu/rcc.uniq.IGF.path.pdf',p,width = 7.2,height = 4)

# IL4
rcc.IL4 <- rcc.path.uniq[rcc.path.uniq$pathway_name=='IL4',]
p <- ggplot(rcc.IL4, aes(x=source,y=target,size=`log10(prob)`,color=`log10(prob)`)) +
  geom_point() + scale_color_gradientn(colors = mycol) +
  theme_minimal() + ggtitle('RCC specific pathway - IL4') +
  theme(axis.text = element_text(size = 14,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5,size = 16,face = 'bold'))
ggsave('wclc/myeloid/cellcommu/rcc.uniq.IL4.path.pdf',p,width = 7.4,height = 6.8)

# NT
rcc.NT <- rcc.path.uniq[rcc.path.uniq$pathway_name=='NT',]
p <- ggplot(rcc.NT, aes(x=source,y=target,size=`log10(prob)`,color=`log10(prob)`)) +
  geom_point() + scale_color_gradientn(colors = mycol) +
  theme_minimal() + ggtitle('RCC specific pathway - NT') +
  theme(axis.text = element_text(size = 14,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5,size = 16,face = 'bold'))
ggsave('wclc/myeloid/cellcommu/rcc.uniq.NT.path.pdf',p,width = 6,height = 4)

# OSM
rcc.OSM <- rcc.path.uniq[rcc.path.uniq$pathway_name=='OSM',]
p <- ggplot(rcc.OSM, aes(x=source,y=target,size=`log10(prob)`,color=`log10(prob)`)) +
  geom_point() + scale_color_gradientn(colors = mycol) +
  theme_minimal() + ggtitle('RCC specific pathway - OSM') +
  theme(axis.text = element_text(size = 14,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5,size = 16,face = 'bold'))
ggsave('wclc/myeloid/cellcommu/rcc.uniq.OSM.path.pdf',p,width = 7.4,height = 6.4)

# SEMATOSTATIN
rcc.SEMATOSTATIN <- rcc.path.uniq[rcc.path.uniq$pathway_name=='SEMATOSTATIN',]
p <- ggplot(rcc.SEMATOSTATIN, aes(x=source,y=target,size=`log10(prob)`,color=`log10(prob)`)) +
  geom_point() + scale_color_gradientn(colors = mycol) +
  theme_minimal() + ggtitle('RCC specific pathway - SEMATOSTATIN') +
  theme(axis.text = element_text(size = 14,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5,size = 16,face = 'bold'))
ggsave('wclc/myeloid/cellcommu/rcc.uniq.SEMATOSTATIN.path.pdf',p,width = 4,height = 6.4)

# TNF
rcc.TNF <- rcc.path.uniq[rcc.path.uniq$pathway_name=='TNF',]
p <- ggplot(rcc.TNF, aes(x=source,y=target,size=`log10(prob)`,color=`log10(prob)`)) +
  geom_point() + scale_color_gradientn(colors = mycol) +
  theme_minimal() + ggtitle('RCC specific pathway - TNF') +
  theme(axis.text = element_text(size = 14,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5,size = 16,face = 'bold'))
ggsave('wclc/myeloid/cellcommu/rcc.uniq.TNF.path.pdf',p,width = 6,height = 6.4)


# p <- ggplot(rcc.path.uniq, aes(x=pathway_name,y=pairs,size=`log10(prob)`,color=`log10(prob)`)) +
#   geom_point() + scale_color_gradientn(colors = mycol) +
#   theme_minimal() + ggtitle('RCC specific pathway') +
#   coord_flip() +
#   theme(axis.text = element_text(size = 14,colour = 'black'),
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
#         plot.title = element_text(hjust = .5,size = 16,face = 'bold'))
# ggsave('wclc/myeloid/cellcommu/rcc.uniq.path.pdf',p,width = 16,height = 6)

save(lcc.path.uniq,rcc.path.uniq,file = 'wclc/rds/mye.lcc_rcc.uniq.pathway.rda')


lcc.inte <- xlsx::read.xlsx('wclc/myeloid/cellcommu/lcc/lcc.lr.xlsx',sheetIndex = 1)
rcc.inte <- xlsx::read.xlsx('wclc/myeloid/cellcommu/rcc/rcc.lr.xlsx',sheetIndex = 1)

# lcc.uniq.lr.pair
lcc.inte$pathway_name <- as.character(lcc.inte$pathway_name)
lcc.uniq.lr <- lcc.inte[lcc.inte$pathway_name=='BTLA',]
rcc.inte$pathway_name <- as.character(rcc.inte$pathway_name)
rcc.uniq.path <- setdiff(rcc.inte$pathway_name,lcc.inte$pathway_name)
rcc.uniq.lr <- NULL
for(i in 1:length(rcc.uniq.path)){
  rcc.uniq.i <- rcc.inte[rcc.inte$pathway_name==rcc.uniq.path[i],]
  rcc.uniq.lr <- rbind(rcc.uniq.lr,rcc.uniq.i)
}

lcc.uniq.lr$interaction_name_3 <- paste(lcc.uniq.lr$source,lcc.uniq.lr$target,sep = ' | ')
rcc.uniq.lr$interaction_name_3 <- paste(rcc.uniq.lr$source,rcc.uniq.lr$target,sep = ' | ')
lcc.uniq.lr$tissue <- "LCC"
rcc.uniq.lr$tissue <- "RCC"
save(lcc.uniq.lr,rcc.uniq.lr,file = 'wclc/rds/mye.lcc_rcc.uniq.lr.rda')

lcc.dat <- lcc.uniq.lr[,c('tissue','interaction_name_2','interaction_name_3','prob')]
lcc.dat$logProb <- log10(lcc.dat$prob)
lcc.dat$prob <- NULL
rcc.dat <- rcc.uniq.lr[,c('tissue','interaction_name_2','interaction_name_3','prob')]
rcc.dat$logProb <- log10(rcc.dat$prob)
rcc.dat$prob <- NULL

dat <- rbind(lcc.dat, rcc.dat)
mycol <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(6))
p <- ggplot(dat,aes(x=interaction_name_3,y=interaction_name_2,fill=logProb)) +
  geom_raster() + scale_fill_gradientn(colours = mycol) +
  theme(panel.grid = element_line(size = 0.1,colour = 'grey40'),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.1,colour = 'grey40',fill = NA),
        axis.title = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 15,colour = 'black'))
ggsave('wclc/myeloid/cellcommu/uniq.lr.pdf',p,width = 48,height = 10)


load('wclc/rds/mye.lcc_rcc.uniq.lr.rda')
lcc.dat <- lcc.uniq.lr[,c('tissue','interaction_name_2','source','target','prob')]
lcc.dat$logProb <- log10(lcc.dat$prob)
lcc.dat$prob <- NULL
rcc.dat <- rcc.uniq.lr[,c('tissue','interaction_name_2','source','target','prob','pathway_name')]
rcc.dat$logProb <- log10(rcc.dat$prob)
rcc.dat$prob <- NULL

mycol <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(6))
p1 <- ggplot(lcc.dat,aes(x=source,y=target,fill=logProb,size=logProb,color=logProb)) +
  geom_point() + scale_color_gradientn(colours = mycol) +
  ggtitle('BTLA signaling (BTLA - TNFRSF14)') +
  theme(panel.grid = element_line(size = 0.1,colour = 'grey40'),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.1,colour = 'grey40',fill = NA),
        axis.title = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = .5,size = 16,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 15,colour = 'black'))
ggsave('wclc/myeloid/cellcommu/lr/lcc.uniq.lr.pdf',p1,width = 6,height = 5.5)

rcc.EGF <- rcc.dat[rcc.dat$pathway_name == 'EGF',]
rcc.EPO <- rcc.dat[rcc.dat$pathway_name == 'EPO',]
rcc.IGF <- rcc.dat[rcc.dat$pathway_name == 'IGF',]
rcc.IL4 <- rcc.dat[rcc.dat$pathway_name == 'IL4',]
rcc.NT <- rcc.dat[rcc.dat$pathway_name == 'NT',]
rcc.OSM <- rcc.dat[rcc.dat$pathway_name == 'OSM',]
rcc.SEMATOSTATIN <- rcc.dat[rcc.dat$pathway_name == 'SEMATOSTATIN',]
rcc.TNF <- rcc.dat[rcc.dat$pathway_name == 'TNF',]

p2.1 <- ggplot(rcc.EGF,aes(x=source,y=target,fill=logProb,size=logProb,color=logProb)) +
  geom_point() + scale_color_gradientn(colours = mycol) +
  ggtitle('EGF signaling (BTC - EGFR)') +
  theme(panel.grid = element_line(size = 0.1,colour = 'grey40'),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.1,colour = 'grey40',fill = NA),
        axis.title = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = .5,size = 16,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 15,colour = 'black'))
ggsave('wclc/myeloid/cellcommu/lr/rcc.uniq.BTC-EGFR.pdf',p2.1,width = 3.5,height = 5.5)
p2.2 <- ggplot(rcc.EPO,aes(x=source,y=target,fill=logProb,size=logProb,color=logProb)) +
  geom_point() + scale_color_gradientn(colours = mycol) +
  ggtitle('EPO signaling (EPO - EPOR)') +
  theme(panel.grid = element_line(size = 0.1,colour = 'grey40'),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.1,colour = 'grey40',fill = NA),
        axis.title = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = .5,size = 16,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 15,colour = 'black'))
ggsave('wclc/myeloid/cellcommu/lr/rcc.uniq.EPO-EPOR.pdf',p2.2,width = 5.5,height = 5.5)
p2.3 <- ggplot(rcc.IGF,aes(x=source,y=target,fill=logProb,size=logProb,color=logProb)) +
  geom_point() + scale_color_gradientn(colours = mycol) +
  ggtitle('IGF signaling (IGF1 - (ITGA6+ITGB4))') +
  theme(panel.grid = element_line(size = 0.1,colour = 'grey40'),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.1,colour = 'grey40',fill = NA),
        axis.title = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = .5,size = 16,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 15,colour = 'black'))
ggsave('wclc/myeloid/cellcommu/lr/rcc.uniq.IGF1-ITGA6_ITGB4.pdf',p2.3,width = 6,height = 4)
p2.4 <- ggplot(rcc.IL4,aes(x=source,y=target,fill=logProb,size=logProb,color=logProb)) +
  geom_point() + scale_color_gradientn(colours = mycol) +
  ggtitle('IL4 signaling (IL4 - IL4R)') +
  theme(panel.grid = element_line(size = 0.1,colour = 'grey40'),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.1,colour = 'grey40',fill = NA),
        axis.title = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = .5,size = 16,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 15,colour = 'black'))
ggsave('wclc/myeloid/cellcommu/lr/rcc.uniq.IL4-IL4R.pdf',p2.4,width = 6,height = 5.5)
p2.5 <- ggplot(rcc.NT,aes(x=source,y=target,fill=logProb,size=logProb,color=logProb)) +
  geom_point() + scale_color_gradientn(colours = mycol) +
  ggtitle('NT signaling (BDNF - NTRK2)') +
  theme(panel.grid = element_line(size = 0.1,colour = 'grey40'),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.1,colour = 'grey40',fill = NA),
        axis.title = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = .5,size = 16,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 15,colour = 'black'))
ggsave('wclc/myeloid/cellcommu/lr/rcc.uniq.BDNF-NTRK2.pdf',p2.5,width = 5,height = 5)
p2.6 <- ggplot(rcc.OSM,aes(x=source,y=target,fill=logProb,size=logProb,color=logProb)) +
  geom_point() + scale_color_gradientn(colours = mycol) +
  ggtitle('OSM signaling (OSM - (LIFR+IL6ST))') +
  theme(panel.grid = element_line(size = 0.1,colour = 'grey40'),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.1,colour = 'grey40',fill = NA),
        axis.title = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = .5,size = 16,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 15,colour = 'black'))
ggsave('wclc/myeloid/cellcommu/lr/rcc.uniq.OSM-LIFR_IL6ST.pdf',p2.6,width = 6.5,height = 6)
p2.7 <- ggplot(rcc.SEMATOSTATIN,aes(x=source,y=target,fill=logProb,size=logProb,color=logProb)) +
  geom_point() + scale_color_gradientn(colours = mycol) +
  ggtitle('SEMATOSTATIN signaling (CORT - SSTR3)') +
  theme(panel.grid = element_line(size = 0.1,colour = 'grey40'),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.1,colour = 'grey40',fill = NA),
        axis.title = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = .5,size = 16,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 15,colour = 'black'))
ggsave('wclc/myeloid/cellcommu/lr/rcc.uniq.CORT-SSTR3.pdf',p2.7,width = 3,height = 5.5)
p2.8 <- ggplot(rcc.TNF,aes(x=source,y=target,fill=logProb,size=logProb,color=logProb)) +
  geom_point() + scale_color_gradientn(colours = mycol) +
  ggtitle('TNF signaling (TNF - TNFRSF1A)') +
  theme(panel.grid = element_line(size = 0.1,colour = 'grey40'),
        panel.background = element_blank(),
        panel.border = element_rect(size = 0.1,colour = 'grey40',fill = NA),
        axis.title = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = .5,size = 16,colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 15,colour = 'black'))
ggsave('wclc/myeloid/cellcommu/lr/rcc.uniq.TNF-TNFRSF1A.pdf',p2.8,width = 5.5,height = 5.5)


