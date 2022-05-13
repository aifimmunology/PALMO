#' mergePALMOdata Function
#'
#' This function allows to merge expression data from bulk or single cell data
#' with annotation file provided by user. It will remove the annotations with
#' missing information from Sample name, donor/participant and time variable.
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information and expression data from Bulk or single cell data.
#' @param datatype Input data type 'bulk' or 'singlecell'
#' @return PALMO object with merged annotation and data matrix
#' @keywords mergePALMOdata
#' @export
#' @examples
#' \dontrun{
#' palmo_obj <- mergePALMOdata(data_object=palmo_obj)
#' }

mergePALMOdata <- function(data_object, datatype) {

    ## Get the annotation data and data matrix
    anndata <- data_object@curated$ann
    data <- data_object@raw$data

    if (datatype == "bulk") {
        data_object@datatype <- "bulk"
        colnames(data) <- gsub("-", "_", colnames(data))
        overlap <- intersect(row.names(anndata), colnames(data))
        if (length(overlap) > 2) {
            anndata <- anndata[overlap, ]
            data <- data.frame(data, check.names=FALSE, stringsAsFactors=FALSE)
            data <- data[, overlap]
            data_object@curated$anndata <- anndata
            data_object@curated$data <- data
            return(data_object)
        } else {
            stop(date(), ": Please provide appropriate sample_column field")
        }

    } else if (datatype == "singlecell") {
        data_object@datatype <- "singlecell"
        # Sample overlap for data.frame type input data
        class_data <- class(data)
        if (class(data)[1] %in% c("data.frame", "matrix")) {
            # for data.frame or Matrix format single cell/omics data
            cn <- data.frame(Sample_group = colnames(data))
            temp <- data.frame(do.call(rbind, strsplit(cn$Sample_group,
                                split = ":")), stringsAsFactors = FALSE)
            cn <- data.frame(cn, Sample = temp$X1, group = temp$X2,
                             stringsAsFactors = FALSE)
            cn$Sample <- gsub("-", "_", cn$Sample)
            cn <- merge(anndata, cn, by = "Sample")
            cn <- cn[!is.na(cn$Sample_group), ]
            row.names(cn) <- cn$Sample_group
            ann <- cn
            rm(cn)

            data <- data.frame(data, check.names=FALSE, stringsAsFactors=FALSE)
            data <- data[, ann$Sample_group]

            ## Keep genes with avgExpression > zero
            rowDF <- rowSums(data)
            rowDF <- rowDF[rowDF > 0]
            data <- data[names(rowDF), ]

            ## Assign to data object
            data_object@raw$data <- NULL  #reduce object size
            data_object@curated$anndata <- ann
            data_object@curated$data <- data
            return(data_object)

        } else if (class(data)[1] %in% c("Seurat")) {
            # For single cell RNA Seurat object
            checkSample <- intersect(colnames(data@meta.data), "Sample")
            if (length(checkSample) == 1) {
                data@meta.data$Sample <- gsub("-", "_", data@meta.data$Sample)
                overlap <- intersect(anndata$Sample, data@meta.data$Sample)
                anndata <- anndata[overlap, ]
                # in-case subset of samples only
                data <- subset(x = data, subset = Sample %in% overlap)

                # Aggregate data at celltypes (psuedo-bulk) For single cell
                # data merge annotation and single cell metadata
                metaData <- data@meta.data
                anndata_temp <- anndata[metaData$Sample, ]
                metaData <- cbind(metaData, anndata_temp)
                colnames(metaData) <- make.names(colnames(metaData),unique=TRUE)
                data@meta.data <- metaData

                # Final object
                data_object@raw$data <- NULL  #reduce object size
                data_object@curated$anndata <- anndata
                data_object@curated$SeuratObj <- data
                return(data_object)

            } else {
                stop(date(), ": Please provide appropriate sample_column field")
            }
        }
    } else {
        stop(date(), ": Please provide datatype 'bulk' or 'singlecell'")
    }
}
