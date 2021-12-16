#' A dimUMAPPlot Function
#'
#' This function allows you to perform UMAP visualization of gene of interest list.
#' @param ann Annotation table. Table must consist column Sample (Participant 
#' sample name), PTID (Participant), Time (longitudinal time points), group,
#' name of the group, group_donor (combined string using group:Sample)
#' @param rnaObj The seurat scRNA object in case of single cell RNA data (optional).
#' @param countMat Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column),
#' in case count matrix for expression data (optional).
#' @param gene_oi Genes of interest to explore, required
#' @param nPC Number of PCAs to be used for UMAP, Default is 30
#' @param groupName User-defined group name column from annotation table or seurat
#' annotation column, required
#' @param plotname User-defined output file name, required
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords dimUMAPPlot
#' @export
#' @examples
#' ##Count/genescore matrix data
#' #dimUMAPPlot(ann=annotation, countMat=countData, nPC=15, gene_oi=var_gene, 
#' #groupName="celltype", plotname="variable", filePATH=filePATH, fileName="ATAC")
#' 
#' ##Single cell RNA data
#' #dimUMAPPlot(rnaObj=SeuratObj, nPC=15, gene_oi=var_gene, groupName="celltype",
#' #plotname="variable", filePATH=filePATH, fileName="scRNA")

dimUMAPPlot <- function(ann, rnaObj=NULL, countMat=NULL, nPC=30, gene_oi=NULL, groupName=NULL, plotname=NULL, filePATH=NULL, fileName=NULL) {
    
  cat(date(),": Visualizing UMAP\n")
  if(!is.null(rnaObj)) {
      #Remove scaled data if any
      rnaObj <- DietSeurat(rnaObj, counts = TRUE, data = TRUE,
                     scale.data = FALSE, features = NULL, assays = NULL, dimreducs = NULL, graphs = NULL)
  } else if(!is.null(countMat)) {
      #create seurat object
      rnaObj <- CreateSeuratObject(counts=countMat, project = "proj", assay = "RNA",
                                 min.cells = 0, min.features = 0, names.field = 1,
                                 names.delim = "_", meta.data = NULL)
      metaData <- rnaObj@meta.data
      rnaObj@meta.data <- data.frame(metaData, ann, stringsAsFactors = F)
      metaData <- rnaObj@meta.data
      rnaObj <- FindVariableFeatures(rnaObj, selection.method = "vst", nfeatures = 3000)
      save(rnaObj, file=paste(filePATH,"/",fileName,"-genematatrix-seuratobj.Rda", sep=""))
  } else {
      cat(now(),": Please provide Seurat object as input for scRNA or archR object for scATAC\n")
      stop()
  }

  if(is.null(nPC)) {
    nPC=30 #Number of PCA to be used, default 30
  }

  if(is.null(groupName)) {
    cat(now(),": Please provide avgGroup\n")
  }
    
  #Get the variable features
  seurat_vargenes <- VariableFeatures(rnaObj)
  top_seurat <- seurat_vargenes[1:length(gene_oi)]
  gene_oi <- intersect(gene_oi, row.names(rnaObj))
    
  #Plot
  rnaObj <- ScaleData(rnaObj, features = gene_oi)
  rnaObj <- suppressMessages(RunPCA(rnaObj, features=gene_oi, npcs = nPC, approx=TRUE), classes = "message")
  rnaObj <- suppressMessages(RunUMAP(rnaObj, reduction = "pca", dims = 1:nPC), classes = "message")
  #To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation
  p1 <- DimPlot(object = rnaObj, reduction = "pca", group.by = groupName, label = T)
  p2 <- DimPlot(object = rnaObj, reduction = 'umap', group.by = groupName, label = T)
    
  png(paste(filePATH,"/",fileName,"-UMAP-",plotname,"-Genes.png", sep=""), width=16, height=5, res=200, units = "in")
  print(plot_grid(p1, p2, align="hv"))
  dev.off()
    
  return(rnaObj)
    
}
