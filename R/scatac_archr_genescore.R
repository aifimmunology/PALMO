#' A scatac_archr_genescore Function
#'
#' This function allows you to calculate genescore matrix from scATAC archR object.
#' This function requires archR package installed and scATAC object created.
#' @param ArchRProj archR scATAC object for input single cell ATAC longitudinal data
#' @param groupBy Group label to be used to calculate average gene expression by
#' group label, Eg. "celltype"
#' @keywords scatac_archr_genescore
#' @export
#' @examples
#' ##Input scATAC data
#' #genescore <- scatac_archr_genescore(ArchRProj=proj, groupBy="celltype")

scatac_archr_genescore <- function(ArchRProj, groupBy) {
  
  #if(!require("ArchR")){ install.packages("ArchR"); library("ArchR"); }
  #Define sample group and Calculate Genescore matrix
  #Average expression: Get the aggregate data for celltype from scATAC
  cat(date(), ": Creating Group-based Genescore Matrix\n")
  res <- getGroupSE(ArchRProj = ArchRProj,
                    useMatrix = "GeneScoreMatrix", groupBy = groupBy,
                    divideN = TRUE, scaleTo = NULL,
                    threads = getArchRThreads(),
                    verbose = FALSE,
                    logFile = createLogFile("getGroupSE")
  )
  scatac_avgmat <- res@assays@data@listData[["GeneScoreMatrix"]]
  cat(date(), ": Group-based Genescore Matrix finished\n")
  #row.names(scatac_avgmat)[1:5]
  ele_metdata <- data.frame(res@elementMetadata@listData, row.names(scatac_avgmat))
  row.names(scatac_avgmat) <- ele_metdata$name
  scatac_avgmat <- data.frame(scatac_avgmat, check.names = F, stringsAsFactors = F)
  #colSums(scatac_avgmat)
  return(scatac_avgmat)
  cat(date(), ": Group-based Genescore Matrix done\n")
}
