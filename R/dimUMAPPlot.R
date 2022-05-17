#' dimUMAPPlot Function
#'
#' This function allows to perform UMAP visualization of gene of interest list.
#' @param data_object Input \emph{PALMO} S4 object. Contains annotation table and
#' single cell data stored as Seurat scRNA object.
#' @param gene_oi Genes of interest to explore, required
#' @param nPC Number of PCAs to be used for UMAP, Default is 30
#' @param group_column User-defined group name column from annotation table or
#' Seurat annotation column. Example, group_column='celltype' (required)
#' @param plotname User-defined output file name (required)
#' @param repel UMAP plot with labels repel=TRUE. Default FALSE
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory \emph{PATH} Default, current
#' directory
#' @return UMAP plot
#' @keywords dimUMAPPlot
#' @export
#' @examples
#' \dontrun{
#' dimUMAPPlot(data_object=pamo_obj, nPC=15, gene_oi=stable_gene,
#' group_column='celltype', plotname='stable')
#' }

dimUMAPPlot <- function(data_object, nPC = 30, gene_oi, group_column, plotname = NULL, repel = FALSE,
    filePATH = NULL, fileName = NULL) {

    message(date(), ": Visualizing UMAP")
    ## If filename or filepath null
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    if (is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }

    ## Number of PCs
    if (is.null(nPC)) {
        nPC = 30  #Number of PCA to be used, default 30
    }

    ## group_column is required
    if (is.null(group_column)) {
        stop(date(), ": Please provide avgGroup like 'group', 'celltype'.")
    }
    ## Plotname is required
    if (is.null(plotname)) {
        plotname <- "dimred_output"
    }

    ## Get the variable features
    if (!is.null(data_object@curated$SeuratObj)) {
        rnaObj <- data_object@curated$SeuratObj
        seurat_vargenes <- VariableFeatures(rnaObj)
        top_seurat <- seurat_vargenes[1:length(gene_oi)]
        gene_oi <- intersect(gene_oi, row.names(rnaObj))

        ## Plot
        rnaObj <- ScaleData(rnaObj, features = gene_oi)
        rnaObj <- suppressMessages(RunPCA(rnaObj, features = gene_oi,
                                          npcs = nPC, approx = TRUE),
                                   classes = "message")
        rnaObj <- suppressMessages(RunUMAP(rnaObj, reduction = "pca",
                                           dims = 1:nPC),
                                   classes = "message")

        ## UMAP plot
        p1 <- DimPlot(object = rnaObj, reduction = "pca",
                      group.by = group_column, label = TRUE, repel = repel)
        p2 <- DimPlot(object = rnaObj, reduction = "umap",
                      group.by = group_column, label = TRUE, repel = repel)

        png(paste(filePATH, "/", fileName, "-UMAP-", plotname, "-Genes.png",
                  sep = ""), width = 16, height = 5, res = 200, units = "in")
        print(plot_grid(p1, p2, align = "hv"))
        dev.off()

        message(date(), ": Please check output directory for PC, UMAP plot")
        ## plot result
        print(plot_grid(p1, p2, align = "hv"))
    } else {
        stop(date(), ": PALMO object do not contain Seurat object. Please
        create PALMO object.")
    }
}
