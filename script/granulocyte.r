# ---------------- granulocyte ----------------

mye <- readRDS('wclc/rds/myeloid.tumor.sub.rds')
gra <- subset(mye, cells = rownames(subset(mye@meta.data, clMidwayPr == 'Granulo')))

# N1 TAN: TNF,CCL3,ICAM1,EPAS1(HIF2A)
# N2 Neutrophils: CCL2, CCL3, CCL4, CCL8, CCL17, CXCL1, CXCL2, CXCL8(IL8), CXCL16
VlnPlot(gra, features = 'TNF', group.by = 'TissueSiteSimple')

features <- c("TNF","CCL3","ICAM1","EPAS1",
              "CCL2","CCL4","CCL8","CCL17","CXCL1",'CXCL2','CXCL8','CXCL16')

p.features <- DotPlot(gra, features = features,group.by = "TissueSiteSimple") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 14,face = "bold"),
        axis.text.y = element_text(size = 14,face = "bold"),
        panel.border = element_rect(size = 1.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdYlBu"))) +
  guides(size = guide_legend(title = "pct.exp"), color = guide_colorbar(title = "ave.exp"))

ggsave(filename = 'wclc/gra.lcc.rcc.pdf',p.features,width = 8,height = 3)

saveRDS(gra, file = 'wclc/rds/gra.rds')

gra <- NormalizeData(object = gra, normalization.method = "LogNormalize", scale.factor = 10000)
gra <- FindVariableFeatures(gra,selection.method="vst",nfeatures=2000)
gra <- ScaleData(gra, features = VariableFeatures(gra))
gra <- RunPCA(object = gra, features = VariableFeatures(gra))
# Find neighbors and clusters with harmony batch correction
gra <- FindNeighbors(gra, dims = 1:30, reduction = "pca")
gra <- FindClusters(gra, resolution = 0.2)
gra <- RunUMAP(gra, dims = 1:30, reduction = "pca")
gra <- RunTSNE(gra,dims = 1:30)

saveRDS(gra,file = 'wclc/rds/gra.rds')

mypal <- c('#1F78B4','#E31A1C','#FF7F00','#FB9A99')
p <- DimPlot(gra,cols = mypal,label = T,label.size = 5) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1,colour = 'black'))
ggsave(filename = 'wclc/gra/seurat.cluster.pdf',p,width = 5,height = 5)

p <- ggplot(gra@meta.data, aes(x = TissueSiteSimple, fill = seurat_clusters)) +
  geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values = mypal) +
  mytheme +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(filename = 'wclc/gra/cell.prop.pdf',p,width = 5,height = 8)


pathway <- xlsx::read.xlsx("wclc/xlsx/pathway.tmp.xlsx",sheetIndex = 3)
p <- ggplot(pathway, aes(x=Tissue,y=Description,color=Genes,size=-log10(pvalue))) +
  geom_point() + theme_bw() + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 14,colour = "black"),
        legend.position = "left",
        axis.title = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20),limits=c(1,4))
ggsave(filename = "wclc/gra/pathway.lcc.rcc.dotplot.pdf",p,width = 8,height = 6)


features <- c("TNF","CCL3")
p.features <- DotPlot(gra, features = features,group.by = "TissueSiteSimple") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 14,face = "bold"),
        axis.text.y = element_text(size = 14,face = "bold"),
        panel.border = element_rect(size = 1.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdYlBu"))) +
  guides(size = guide_legend(title = "pct.exp"), color = guide_colorbar(title = "ave.exp"))

ggsave(filename = "wclc/gra/TNF.CCL3.pdf",p.features,width = 5,height = 4)


load('wclc/rds/markers/markers.gra.rda')

# load expression data
coad.expr <- read.table('data/TCGA/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt', sep = "\t", header = T)
coad.expr <- coad.expr[!duplicated(coad.expr$Hugo_Symbol),]
rownames(coad.expr) <- coad.expr$Hugo_Symbol
coad.expr <- coad.expr[, -c(1,2)]
colnames(coad.expr) <- gsub("\\.", "-", colnames(coad.expr))
coad.expr <- data.frame(t(coad.expr))
gra.up <- markers.gra[markers.gra$avg_logFC > 0,]
top20 <- head(rownames(gra.up),20)
up.comm <- intersect(top20, colnames(coad.expr))

gra.dn <- markers.gra[markers.gra$avg_logFC < 0,]
tail20 <- tail(rownames(gra.dn),20)
dn.comm <- intersect(tail20, colnames(coad.expr))

features <- c(up.comm, dn.comm)

Idents(gra) <- "TissueSiteSimple"
exp.mat <- AverageExpression(gra,features = features)$RNA
tmp <- log2(exp.mat)
tmp.2 <- t(tmp)
pheatmap::pheatmap(tmp.2,cluster_rows = F,
                   color = colorRampPalette(RColorBrewer::brewer.pal(8,"Reds"))(8),
                   border=NA, filename = "wclc/gra/LCC.RCC.features.pdf",
                   width = 8,height = 3)

old.ident <- c("left","right")
new.ident <- c("N1 TANs", "N2 TANs")
gra$celltype <- ifelse(gra$TissueSiteSimple == "left","N1 TANs","N2 TANs")
saveRDS(gra, file = "wclc/rds/gra.rds")

DimPlot(gra,group.by = "celltype",split.by = "TissueSiteSimple",cols = c("#4E9BD2","#DB4740"))

lcc.gra <- subset(gra,cells = rownames(subset(gra@meta.data, TissueSiteSimple == "left")))
rcc.gra <- subset(gra,cells = rownames(subset(gra@meta.data, TissueSiteSimple == "right")))

p <- DimPlot(lcc.gra,group.by = "celltype",cols = c("#4E9BD2","#DB4740"))
l1 <- cowplot::get_legend(p)
p1 <- DimPlot(lcc.gra,group.by = "celltype",cols = c("#4E9BD2","#DB4740")) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1,colour = "black"))
ggsave(filename = "wclc/gra/lcc.N1.TAN.pdf",p1,width = 4,height = 4)

p <- DimPlot(rcc.gra,group.by = "celltype",cols = c("#DB4740","#4E9BD2"))
l2 <- cowplot::get_legend(p)
p2 <- DimPlot(rcc.gra,group.by = "celltype",cols = c("#DB4740","#4E9BD2")) +
  ggtitle('') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(size = 1,colour = "black"))
ggsave(filename = "wclc/gra/lcc.N1.TAN.pdf",p1,width = 4,height = 4)

pc <- cowplot::plot_grid(plotlist = list(p1,l1,p2,l2),rel_widths = 2, rel_heights = 4)
ggsave(filename = "wclc/gra/LCC.RCC.comb.pdf",pc,width = 5,height = 5)


### survival

mytheme.3 <- theme(legend.title = element_text(size = 14), 
                   legend.text = element_text(size = 14),
                   axis.text.x = element_text(size = 14),
                   axis.text.y = element_text(size = 14),
                   axis.line.x = element_line(size = 1),
                   axis.line.y = element_line(size = 1),
                   axis.title.x = element_text(size = 16),
                   axis.title.y = element_text(size = 16))

# clin
clin <- read.table('data/TCGA/data_clinical_patient.txt',skip = 4, header = T,sep = '\t')
small_clin <- data.frame(clin$PATIENT_ID, clin$OS_MONTHS, clin$OS_STATUS, clin$PRIMARY_SITE_PATIENT)
small_clin <- small_clin[small_clin$clin.OS_STATUS != "[Not Available]",]
small_clin <- small_clin[small_clin$clin.PRIMARY_SITE_PATIENT != "[Not Available]",]
small_clin <- small_clin[small_clin$clin.PRIMARY_SITE_PATIENT != "[Discrepancy]",]
small_clin$clin.PRIMARY_SITE_PATIENT <- as.character(small_clin$clin.PRIMARY_SITE_PATIENT)

colon <- data.frame(TissueSite = names(table(small_clin$clin.PRIMARY_SITE_PATIENT)),
                    TissueSiteSimple = c('right','right','left','right','left','left',
                                         'left','left','left'),stringsAsFactors = F)
small_clin$TissueSiteSimple <- plyr::mapvalues(small_clin$clin.PRIMARY_SITE_PATIENT, 
                                               from = colon$TissueSite, to = colon$TissueSiteSimple)

# create binary variable if dead or alive
Dead <- rep(NA, dim(small_clin)[1])
for(i in 1:length(small_clin$clin.OS_STATUS)){
  if(stringr::str_split(small_clin$clin.OS_STATUS,":",simplify = T)[,2][i] == "DECEASED"){
    Dead[i] <- 1
  }else if(stringr::str_split(small_clin$clin.OS_STATUS,":",simplify = T)[,2][i] == "LIVING"){
    Dead[i] <- 0
  }
}

small_clin$Dead <- Dead
small_clin$clin.OS_MONTHS <- as.integer(small_clin$clin.OS_MONTHS)
save(small_clin, file = 'wclc/rds/small.clin.rda')

rcc.clin <- small_clin[small_clin$TissueSiteSimple == 'right',]
lcc.clin <- small_clin[small_clin$TissueSiteSimple == 'left',]

save(lcc.clin, rcc.clin, small_clin, file = 'wclc/rds/small.clin.rda')
rm(small_clin,clin);gc()

# load expression data
coad.expr <- read.table('data/TCGA/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt', sep = "\t", header = T)
coad.expr <- coad.expr[!duplicated(coad.expr$Hugo_Symbol),]
rownames(coad.expr) <- coad.expr$Hugo_Symbol
coad.expr <- coad.expr[, -c(1,2)]
colnames(coad.expr) <- gsub("\\.", "-", colnames(coad.expr))
coad.expr <- data.frame(t(coad.expr))

# just grab genes in signatures (RCC)
load('wclc/rds/markers.gra.rda')
gra.up <- markers.gra[markers.gra$avg_logFC > 0,]
top20 <- head(rownames(gra.up),20)
up.comm <- intersect(top20, colnames(coad.expr))
coad.expr.2 <- data.frame(coad.expr[,up.comm])
rownames(coad.expr.2) <- rownames(coad.expr)
colnames(coad.expr.2) <- gsub("expr.","",colnames(coad.expr.2))
coad.expr.3 <- data.frame(scale(coad.expr.2))
coad.expr.3$PATIENT_ID <- gsub("-01","",rownames(coad.expr.3))
# average over those with multiple measurements
ave.exp <- aggregate(.~PATIENT_ID, FUN = mean, data = coad.expr.3)
# RCC
merged.data <- merge(rcc.clin, ave.exp, by.x = "clin.PATIENT_ID", by.y = "PATIENT_ID")
# survival analysis
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS, Dead) ~., data = merged.data[,c(2,6:ncol(merged.data))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix((merged.data[,7:ncol(merged.data)] - colMeans(merged.data[,7:ncol(merged.data)]))) %*% coefs
merged.data$risk_score <- risk.score
group <- rep(NA, dim(merged.data)[1])
for(i in 1:dim(merged.data)[1]){
  if(risk.score[i]>0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data$group <- group
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)
p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 9, y = 0.15, size = 5) +
  annotate("text", label = paste0("RCC(n = 145)"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "RCC.Granulo\nsignatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave(filename = 'RCC.Granulo.surv.pdf',p.surv,width = 5.5,height = 4)

# LCC
gra.dn <- markers.gra[markers.gra$avg_logFC < 0,]
tail20 <- tail(rownames(gra.dn),20)
dn.comm <- intersect(tail20, colnames(coad.expr))
coad.expr.2 <- data.frame(coad.expr[,dn.comm])
rownames(coad.expr.2) <- rownames(coad.expr)
colnames(coad.expr.2) <- gsub("expr.","",colnames(coad.expr.2))
coad.expr.3 <- data.frame(scale(coad.expr.2))
coad.expr.3$ELANE <- NULL
coad.expr.3$PATIENT_ID <- gsub("-01","",rownames(coad.expr.3))
# average over those with multiple measurements
ave.exp <- aggregate(.~PATIENT_ID, FUN = mean, data = coad.expr.3)
# merge data
merged.data <- merge(lcc.clin, ave.exp, by.x = "clin.PATIENT_ID", by.y = "PATIENT_ID")
# survival analysis
coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS, Dead) ~., data = merged.data[,c(2,6:ncol(merged.data))])
coefs <- as.matrix(coxfit$coefficients)
risk.score <- as.matrix((merged.data[,7:ncol(merged.data)] - colMeans(merged.data[,7:ncol(merged.data)]))) %*% coefs
merged.data$risk_score <- risk.score
group <- rep(NA, dim(merged.data)[1])
for(i in 1:dim(merged.data)[1]){
  if(risk.score[i]>0){
    group[i] <- "High"
  }else{
    group[i] <- "Low"
  }
}
merged.data$group <- group
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)
p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 9, y = 0.15, size = 5) +
  annotate("text", label = paste0("LCC(n = 213)"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "LCC.Granulo\nsignatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
ggsave(filename = 'LCC.Granulo.surv.pdf',p.surv,width = 5.5,height = 4)



# ---------------------------------------------------------------------------------------------

library(Seurat)
gra <- readRDS("wclc/rds/gra.rds")
p <- DimPlot(gra,group.by = 'celltype',pt.size = 0.8,cols = c("#4E9BD2","#DB4740")) + 
  theme(legend.position = 'top',
        axis.line = element_blank(),
        plot.title = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black',fill = NA))
ggsave('wclc/gra/celltype.pdf',p,width = 4,height = 4.2)

## feature
p1 <- FeaturePlot(gra,features = 'TNF') +
  theme(axis.line = element_blank(),
        plot.title = element_text(hjust = .5),
        panel.border = element_rect(size = 1,colour = 'black',fill = NA))
ggsave('wclc/gra/features/TNF.pdf',p1,width = 4.2,height = 4)
p2 <- FeaturePlot(gra,features = 'ICAM1') +
  theme(axis.line = element_blank(),
        plot.title = element_text(hjust = .5),
        panel.border = element_rect(size = 1,colour = 'black',fill = NA))
ggsave('wclc/gra/features/ICAM1.pdf',p2,width = 4.2,height = 4)

genes <- c('CCL2','CCL3','CCL4','CCL8','CCL17','CXCL1','CXCL2','CXCL8','CXCL16')
for(i in 1:length(genes)){
  p <- FeaturePlot(gra,features = genes[i]) +
    theme(axis.line = element_blank(),
          plot.title = element_text(hjust = .5),
          panel.border = element_rect(size = 1,colour = 'black',fill = NA))
  ggsave(paste0('wclc/gra/features/',genes[i],'.pdf'),p,width = 4.2,height = 4)
}

# fea1 <- c('TNF','ICAM1','IFNB1','CCL3','TNFSF10','IL21')
# fea2 <- c('IL10','TGFB1','TGFB2','TGFB3','IL22','MIF','IL6')
# fea2 <- c('CCL2','CCL3','CCL4','CCL8','CCL17','CXCL1','CXCL2','CXCL8','CXCL16')

gra <- AddModuleScore(gra,features = list(c('TNF','ICAM1')),name = 'Antitumor')
p <- VlnPlot(gra,features = 'Antitumor1',group.by = 'TissueSiteSimple',
        y.max = 2, pt.size = 0, cols = c("#4E9BD2","#DB4740")) + 
  NoLegend() + ylab('Antitumor Score') +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black',fill = NA),
        axis.text.x = element_text(angle = 0,hjust = .5)) +
  ggpubr::stat_compare_means(comparisons = list(c('left','right')),label = 'p.signif')
ggsave('wclc/gra/AntitumorScore.pdf',p,width = 4,height = 4.2)



library(garnett)
library(org.Hs.eg.db)
gra <- readRDS("wclc/rds/gra.rds")
dat <- as(as.matrix(gra@assays$RNA@counts),"sparseMatrix")
pd <- new('AnnotatedDataFrame', data = gra@meta.data)
fData <- data.frame(gene_short_name = rownames(dat), row.names = rownames(dat))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(dat, phenoData = pd, featureData = fd)
mycds <- estimateSizeFactors(mycds)

marker_file_path <- "data/neutrophils.markers.txt"
marker_check <- check_markers(
  mycds,marker_file_path,db = org.Hs.eg.db,
  cds_gene_id_type = "SYMBOL",
  marker_file_gene_id_type = "SYMBOL"
)
marker_check <- marker_check[-12,]
plot_markers(marker_check)

mycds_classifier <- train_cell_classifier(
  cds = mycds,
  marker_file = marker_file_path,
  db = org.Hs.eg.db,
  cds_gene_id_type = "SYMBOL",
  num_unknown = 10,
  marker_file_gene_id_type = "SYMBOL"
)

mycds <- classify_cells(
  mycds,mycds_classifier,
  db = org.Hs.eg.db,
  cluster_extend = T,
  cds_gene_id_type = "SYMBOL"
)

gra$CellType <- pData(mycds)$cell_type
gra$CellType <- gsub("Unknown","Granulocytes",gra$CellType)
saveRDS(gra, file = "wclc/rds/gra.rds")

p <- DimPlot(gra,group.by = 'CellType',pt.size = 0.8,cols = c('grey',"#4E9BD2","#DB4740")) + 
  theme(legend.position = 'top',
        axis.line = element_blank(),
        plot.title = element_blank(),
        panel.border = element_rect(size = 1,colour = 'black',fill = NA))
ggsave('wclc/gra/celltype2.pdf',p,width = 5,height = 5.2)


gra <- readRDS('wclc/rds/gra.rds')

# CD66b+, CD11b+, CD170+, PDL1+
gene.1 <- 'CD274'
# CD66b+, CD11b+, CD101+, CD177+, CD170-, CD4+, HLA-DR+, CD86+, CD15+
gene.2 <- c('CD101','CD177','CD4','HLA-DRA','HLA-DRB5','HLA-DRB1','CD86','FUT4')

gene <- c(gene.1,gene.2)

gra <- AddModuleScore(gra,features = list(gene.1),name = 'Pro.tumor')
gra <- AddModuleScore(gra,features = list(gene.2),name = 'Anti.tumor')

p1 <- VlnPlot(gra,features = 'Pro.tumor1',group.by = 'TissueSiteSimple',y.max = 3.5) +
  ylab('Pro-tumor Score') + NoLegend() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 20,colour = 'black'),
        axis.text = element_text(size = 18,colour = 'black'),
        axis.text.x = element_text(angle = 0,hjust = .5,vjust = .5)) +
  ggpubr::stat_compare_means(comparisons = list(c('left','right')))
p2 <- VlnPlot(gra,features = 'Anti.tumor1',group.by = 'TissueSiteSimple',y.max = 0.75) +
  ylab('Anti-tumor Score') + NoLegend() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 20,colour = 'black'),
        axis.text = element_text(size = 18,colour = 'black'),
        axis.text.x = element_text(angle = 0,hjust = .5,vjust = .5)) +
  ggpubr::stat_compare_means(comparisons = list(c('left','right')))

cowplot::plot_grid(plotlist = list(p1,p2),ncol = 1)

p1 <-FeaturePlot(gra,features = 'CD274',split.by = 'TissueSiteSimple',pt.size = 1)
p2 <- FeaturePlot(gra,features = 'HLA-DRA',split.by = 'TissueSiteSimple',pt.size = 1)
p3 <- FeaturePlot(gra,features = 'HLA-DRB1',split.by = 'TissueSiteSimple',pt.size = 1)

cowplot::plot_grid(plotlist = list(p1,p2,p3),nrow = 3)



library(Seurat)
gra <- readRDS('wclc/rds/gra.rds')
Phagocytosis <- list(c('MRC1','CD163','MERTK','C1QB'))
Angiogenesis <- list(c('CCND2','CCNE1','CD44','CXCR4','E2F3','EDN1','EZH2',
                       'FGF18','FGFR1','FYN','HEY1','ITGAV','JAG1','JAG2',
                       'MMP9','NOTCH1','PDGFA','PTK2','SPP1','STC1','TNFAIP6',
                       'TYMP','VAV2','VCAN','VEGFA'))
Proliferation <- list(c('AURKA','BUB1','CCNB1','CCND1','CCNE1','DEK','E2F1','FEN1','FOXM1',
                        'H2AFZ','HMGB2','MCM2','MCM3','MCM4','MCM5','MCM6','MKI67','MYBL2',
                        'PCNA','PLK1','TOP2A','TYMS','ZWINT'))
gra <- AddModuleScore(gra,features = Phagocytosis, name = 'Phagocytosis', assay = 'RNA')
gra <- AddModuleScore(gra,features = Angiogenesis, name = 'Angiogenesis', assay = 'RNA')
gra <- AddModuleScore(gra,features = Proliferation, name = 'Proliferation', assay = 'RNA')

meta <- gra@meta.data[,c('orig.ident','TissueSiteSimple','Phagocytosis1','Angiogenesis1','Proliferation1')]
meta$orig.ident <- as.character(meta$orig.ident)
meta$TissueSiteSimple <- as.character(meta$TissueSiteSimple)
PhagocytosisDa <- unname(tapply(meta$Phagocytosis1, meta$orig.ident, mean))
AngiogenesisDa <- unname(tapply(meta$Angiogenesis1, meta$orig.ident, mean))
ProliferationDa <- unname(tapply(meta$Proliferation1, meta$orig.ident, mean))

df <- cbind(PhagocytosisDa,AngiogenesisDa,ProliferationDa)
df <- as.data.frame(df)
df$sample <- names(table(meta$orig.ident))
df$group <- plyr::mapvalues(df$sample,from = meta$orig.ident,to = meta$TissueSiteSimple)


