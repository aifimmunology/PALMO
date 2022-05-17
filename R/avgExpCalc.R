#' avgExpCalc Function
#'
#' This function allows you to calculate average gene expression on
#' log-normalized data by group defined by user. This function uses Seurat
#' function AverageExpression
#' (https://satijalab.org/seurat/reference/averageexpression)
#' @param data_object Input \emph{PALMO} S4 object. Contains annotation table
#' and expression matrix or data frame. Rows represent gene/proteins column
#' represents participant samples (same as annotation table Sample column)
#' @param assay Single cell data Assay type ('RNA', 'SCT'). Default 'RNA'
#' @param group_column Calculate average expression by given group
#' like 'celltype' or 'cluster'
#' @return PALMO object with avg expression
#' @keywords avgExpCalc
#' @export
#' @examples
#' \dontrun{
#' palmo_obj=avgExpCalc(data_object=palmo_obj, assay='RNA',
#' group_column='celltype')
#' }

avgExpCalc <- function(data_object, assay = "RNA", group_column) {

    message(date(), ": Calculating scRNA Average expression")

    ## Get the data
    if (!is.null(data_object@curated$SeuratObj)) {
        anndata <- data_object@curated$anndata
        dataObj <- data_object@curated$SeuratObj

        ## Add sample group to metadata Define sample group and Calculate
        ## average expression
        dataObj@meta.data$group <- dataObj@meta.data[, group_column]
        dataObj@meta.data$Sample_group <- paste(dataObj@meta.data$Sample,
                                                dataObj@meta.data$group,
                                                sep = ":")
        dataObj@meta.data$Sample_group <- gsub(" ", "_", dataObj@meta.data$Sample_group)
        metaData <- dataObj@meta.data

        ## Check assay
        DefaultAssay(dataObj) <- assay
        ## Average expression on log-scaled data
        dataObj@assays[[assay]]@counts <- dataObj@assays[[assay]]@data
        scrna_avgmat <- AverageExpression(object = dataObj, assays = assay,
                                          slot = "counts",
                                          group.by = "Sample_group",
                                          verbose = TRUE)
        cn <- data.frame(colnames(scrna_avgmat[[assay]]))
        cn <- gsub("data\\[, 1]", "", as.character(row.names(cn)))
        scrna_avgmat <- data.frame(scrna_avgmat[[assay]], check.names = FALSE,
                                   stringsAsFactors = FALSE)
        colnames(scrna_avgmat) <- cn
        message(date(), ": scRNA Average expression finished")

        ## Keep genes with avgExpression > zero
        rowDF <- rowSums(scrna_avgmat)
        rowDF <- rowDF[rowDF > 0]
        mat <- scrna_avgmat[names(rowDF), ]
        message(date(), ": Keeping genes with avg expression >0")

        ## Create annotation
        cn <- data.frame(Sample_group = colnames(mat))
        temp <- data.frame(do.call(rbind, strsplit(cn$Sample_group,
                                                   split = ":")),
                           stringsAsFactors = FALSE)
        cn <- data.frame(cn, Sample = temp$X1, group = temp$X2,
                         stringsAsFactors = FALSE)
        row.names(cn) <- cn$Sample_group
        cn <- merge(cn, anndata, by = "Sample", all = TRUE)
        cn <- cn[!is.na(cn$Sample_group), ]
        row.names(cn) <- cn$Sample_group
        ann <- cn
        ann[, ncol(ann) + 1] <- ann$group
        colnames(ann) <- c(colnames(ann)[1:ncol(ann) - 1], group_column)
        ann$Sample_group_i <- paste(ann$group, ann$PTID, sep = ":")
        rm(cn)

        ## Add CV result
        data_object@curated$anndata <- ann
        data_object@curated$data <- mat
        data_object@rownames <- row.names(mat)
        data_object@colnames <- colnames(mat)
        return(data_object)

    } else {
        stop(date(), ": Seurat object is absent")
    }


}
