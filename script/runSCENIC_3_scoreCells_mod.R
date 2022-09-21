#' OpenDev
#' 
openDev <- function(fileName, devType, ...)
{
  if(devType == "pdf")
    pdf(paste0(fileName, ".pdf"), ...)
  
  if(devType == "png")
    png(paste0(fileName, ".png", type = "cairo"),...)
  
  if(devType == "cairo_pdf")
    grDevices::cairo_pdf(paste0(fileName, ".pdf"),...)
}

openDevHeatmap <- function(fileName, devType)
{
  if(devType != "pdf")
  {
    if(devType == "png") openDev(fileName = fileName, devType = devType, width = 1200, height = 1200)
    if(devType == "png") openDev(fileName = fileName, devType = devType)
    fileName <- NA
  }else{
    fileName <- paste0(fileName, ".pdf")
  }
  return(fileName)
}

closeDevHeatmap <- function(devType){
  if(devType != "pdf") dev.off()
}


#' @param scenicOptions - Fields used: Intermediate file
#' @param exp.mat.norm - Log transformed expression matrix
#' @param skipBinaryThresholds - 
#' @param skipHeatmap
#' @param skipTsne
#' @param server - boolean, Determine if run on server
#' @param dist.method - the distance measure to be used. This must be one of "`euclidean`", "`maximum`",
#' "`manhattan`", "`canberra`", "`binary`" or "`minkowski`". Any unambiguous substring can be given.
#' @param hclust.method - the agglomeration method to be used. This should be (an unambiguous abbreviation of)
#' one of "`ward.D`", "`ward.D2`", "`single`", "`complete`", "`average`"(=UPGMA), "`mcquitty`"(=WPGMA),
#' "`median`"(=WPGMC) or "`centroid`"(=UPGMC).

runSCENIC_3_scoreCells_mod <- function(scenicOptions, exp.mat.norm,
                                       skipBinaryThresholds = FALSE, 
                                       skipHeatmap = FALSE, 
                                       skipTsne = FALSE,
                                       server = TRUE,
                                       dist.method = "spear",
                                       hclust.method = "ward.D2") {
  
  pacman::p_load(AUCell,Seurat,SCENIC)

  message(paste0(format(Sys.time(), "%H:%M"), "\tstep1. regulons filtering..."))
  
  nCores <- getSettings(scenicOptions, "nCores")
  regulons <- loadInt(scenicOptions, "regulons")
  regulons <- regulons[order(lengths(regulons), decreasing = TRUE)]
  regulons <- regulons[lengths(regulons) >= 10]
  
  if(length(regulons) < 2)
    stop("Not enough regulons with at least 10 genes!!!")

  message(paste0(format(Sys.time(), "%H:%M"), "\tstep2. Renaming regulons..."))
  
  regulons <- setNames(lapply(names(regulons), function(tf){
    sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))
  }), names(regulons))
  
  names(regulons) <- paste0(names(regulons), " (", lengths(regulons), "g)")
  saveRDS(regulons, file = getIntName(scenicOptions, "aucell_regulons"))

  message(paste0(format(Sys.time(), "%H:%M"), "\tStep3. Analyzing the network activity in each individual cell"))
  
  msg <- paste0("\tNumber of regulons to evaluate on cells: ", length(regulons),
                "\nBiggest (non-extended) regulons: \n",
                paste("\t", grep("_extended",names(regulons),invert = T,value = T)[1:min(length(regulons),10)],collapse = "\n"))
  message(msg)
  
  set.seed(getSettings(scenicOptions, "seed"))
  
  if(server) {
    exp.mat.log <- readRDS(exp.mat.norm)

    openDev(fileName = getIntName(scenicOptions, "aucell_genesStatsPlot"),
            devType = getSettings(scenicOptions, "devType"))
    AUCellRankings <- AUCell_buildRankings(
      exprMat = exp.mat.log, nCores = getSettings(scenicOptions, "nCores"),
      plotStats = TRUE, verbose = getSettings(scenicOptions,"verbose")
    )
    abline(v = AUCellRankings@nGenesDetected["1%"], col = "skyblue3", lwd = 5, lty = 3)
    dev.off()
    saveRDS(AUCellRankings, file = getIntName(scenicOptions, "aucell_rankings"))
  }
  
  regulonAUC <- AUCell_calcAUC(regulons, AUCellRankings, aucMaxRank = AUCellRankings@nGenesDetected["1%"], nCores = nCores)
  variableRegulons <- names(which(apply(getAUC(regulonAUC),1,sd) > 0))
  reguDist <- as.dist(1 - cor(t(getAUC(regulonAUC)[variableRegulons,]), method = dist.method))
  reguClust <- hclust(reguDist, method = hclust.method)
  regulonClusters <- setNames(dynamicTreeCut::cutreeDynamic(reguClust,distM = as.matrix(reguDist),verbose = F),reguClust$labels)
  
  regulonOrder <- reguClust$labels[reguClust$order]
  regulonOrder <- regulonOrder[order(regulonClusters[regulonOrder], decreasing = TRUE)]
  regulonAUC <- regulonAUC[regulonOrder]
  saveRDS(regulonAUC, file = getIntName(scenicOptions, "aucell_regulonAUC"))
  
  cells_AUCellThresholds <- NULL
  
  if(!skipBinaryThresholds) {
    cells_AUCellThresholds <- AUCell_exploreThresholds(regulonAUC, smallestPopPercent = getSettings(scenicOptions,"aucell/smallestPopPercent"),
                                                       assignCells = TRUE, plotHist = FALSE, verbose = FALSE, nCores = nCores)
    saveRDS(cells_AUCellThresholds, file = getIntName(scenicOptions, "aucell_thresholds"))
    
    regulonsCell <- getAssignments(cells_AUCellThresholds)
    trhAssignment <- signif(getThresholdSelected(cells_AUCellThresholds),3)
    commentThresholds <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$comment))
    table2edit <- cbind(regulon = names(cells_AUCellThresholds), threshold = trhAssignment[names(cells_AUCellThresholds)],
                        nCellsAssigned = lengths(regulonsCell)[names(cells_AUCellThresholds)],
                        AUCellComment = commentThresholds[names(cells_AUCellThresholds)],
                        nGenes = gsub("[\\(g\\)]", "", regmatches(names(cells_AUCellThresholds),gregexpr("\\(.*?\\)",names(cells_AUCellThresholds)))),
                        clusteringOrder = seq_along(cells_AUCellThresholds),
                        clusterGroup = regulonClusters[names(cells_AUCellThresholds)],
                        onlyNonDuplicatedExtended = (names(cells_AUCellThresholds) %in% onlyNonDuplicatedExtended(names(cells_AUCellThresholds))),
                        personalNotes = "")
    write.table(table2edit, file = getIntName(scenicOptions, "aucell_thresholdsTxt"), row.names = F, quote = F)
    rm(trhAssignment,regulonsCell,table2edit,cells_AUCellThresholds,commentThresholds);gc()
  }
  
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tFinished running AUCell.")
  message(msg)
  
  if(!skipHeatmap) {
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tPlotting heatmap...")
    message(msg)
    
    nCellsHeatmap <- min(500, ncol(regulonAUC))
    cells2plot <- sample(colnames(regulonAUC), nCellsHeatmap)
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists = "null")
    if(!is.null(cellInfo)) {
      cellInfo <- data.frame(cellInfo)[cells2plot, , drop = F]
    }
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists = "null")
    fileName <- getOutName(scenicOptions, "s3_AUCheatmap")
    fileName <- openDevHeatmap(fileName = fileName, devType = getSettings(scenicOptions, "devType"))
    NMF::aheatmap(getAUC(regulonAUC)[,cells2plot], annCol = cellInfo, annColors = colVars, main = "AUC",
                  sub = paste("Subset of", nCellsHeatmap, " random cells"), filename = fileName)
    closeDevHeatmap(devType = getSettings(scenicOptions, "devType"))
  }
  
  if(!skipTsne) {
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tPlotting t-SNE...")
    message(msg)
    
    tSNE_fileName <- tsneAUC(scenicOptions, aucType = "AUC", onlyHighConf = FALSE)
    tSNE <- readRDS(tSNE_fileName)
    fileName <- getOutName(scenicOptions, "s3_AUCtSNE_colAct")
    plotTsne_AUCellHtml(scenicOptions = scenicOptions, exprMat = exp.mat.log, fileName = fileName, tSNE = tSNE)
    
    sub = ""
    if("type" %in% names(tSNE)) sub <- paste0("t-SNE on", tSNE$type)
    
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions,"cellInfo"),ifNotExists="null")
    colVars <- loadFile(scenicOptions,getDatasetInfo(scenicOptions,"colVars"), ifNotExists="null")
    pdf(paste0(getOutName(scenicOptions, "s3_AUCtSNE_colProps"), ".pdf"))
    plotTsne_cellProps(tSNE$Y, cellInfo = cellInfo, colVars = colVars, cex = 1, sub = sub)
    dev.off()
  }
  
  if (F) {
    message(paste0(format(Sys.time(), "%H:%M"), "\tRegulon binary conversion and visualization..."))
    
    # 查看调整阈值的代码(Optional)
    # Using shiny to interactively optimize thresholds
    aucellApp <- plotTsne_AUCellApp(scenicOptions = scenicOptions, exprMat = exp.mat.log)
    savedSelections <- shiny::runApp(aucellApp)
    # save optimized thresholds
    newThresholds <- savedSelections$thresholds
    scenicOptions@fileNames$int['aucell_thresholds',1] <- "SCENIC/int/newThresholds.rds"
    saveRDS(newThresholds, file = getIntName(scenicOptions, "aucell_thresholds"))
    saveRDS(scenicOptions, file = "SCENIC/Optimized.scenicOptions.rds")
  }
}