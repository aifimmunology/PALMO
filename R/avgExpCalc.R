#' A avgExpCalc Function
#'
#' This function allows you to calculate average gene expression on
#' long-normalized data by group defined by user
#' @param dataObj scRNA object with log-normalized data
#' @param assay Single cell data Assay type ("RNA", "SCT"). Default "RNA"
#' @param group.by Calculate average expression by given group
#' @keywords avgExpCalc
#' @export
#' @examples
#' ##Input Expression data
#' #avgExpCalc(dataObj, group.by)

avgExpCalc <- function(dataObj, assay="RNA", group.by) {

  cat(date(),": Calculating scRNA Average expression\n")
  DefaultAssay(dataObj) <- assay
  dataObj@assays[[assay]]@counts <- dataObj@assays[[assay]]@data #average expression on log-scaled data
  scrna_avgmat <- AverageExpression(object=dataObj, assays=assay, slot="counts", group.by="Sample_group", verbose=TRUE)
  cn <- data.frame(colnames(scrna_avgmat[[assay]]))
  cn <- gsub("data\\[, 1]", "", as.character(row.names(cn)))
  scrna_avgmat <- data.frame(scrna_avgmat[[assay]], check.names = F, stringsAsFactors = F)
  colnames(scrna_avgmat) <- cn
  cat(date(),": scRNA Average expression finished\n")

  return(scrna_avgmat)
}
