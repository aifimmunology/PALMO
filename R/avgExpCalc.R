#' A avgExpCalc Function
#'
#' This function allows you to calculate average gene expression on 
#' long-normalized data by group defined by user
#' @param dataObj scRNA object with log-normalized data
#' @param group.by Calculate average expression by given group
#' @keywords avgExpCalc
#' @export
#' @examples
#' ##Input Expression data
#' #avgExpCalc(dataObj, group.by)

avgExpCalc <- function(dataObj, group.by) {
  
  cat(date(),": Calculating scRNA Average expression\n")
  DefaultAssay(dataObj) <- "RNA"
  dataObj@assays$RNA@counts <- dataObj@assays$RNA@data #average expression on log-scaled data
  scrna_avgmat <- AverageExpression(object=dataObj, assays="RNA", slot="counts", group.by="Sample_group", verbose=TRUE)
  cn <- data.frame(colnames(scrna_avgmat$RNA))
  cn <- gsub("data\\[, 1]", "", as.character(row.names(cn)))
  scrna_avgmat <- data.frame(scrna_avgmat$RNA, check.names = F, stringsAsFactors = F)
  colnames(scrna_avgmat) <- cn
  cat(date(),": scRNA Average expression finished\n")
  
  return(scrna_avgmat)
}
