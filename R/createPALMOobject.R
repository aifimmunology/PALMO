#' createPALMOobject Function
#'
#' This function allows to create PALMO object using Annotation dataframe and
#' Data dataframe. The Data can be bulk data or single cell data. The bulk input
#' data should consists of rows as genes/proteins/features and column as Sample
#' name (same as user-defined Samples in Annotation dataframe). The single cell
#' data should be Seurat object (please check
#' https://search.r-project.org/CRAN/refmans/SeuratObject/html/CreateSeuratObject.html).
#' In case Seurat object not available then user can use function
#' createPALMOfromsinglecellmatrix to create PALMO object. The Seurat
#' object/metadata should have Sample column curesponding to Annotation
#' dataframe.
#' @param anndata Annotation dataframe. It consist of information such as
#' \emph{Sample} (sample name), \emph{PTID} (donor/participant), \emph{Time}
#' (longitudinal timepoints)
#' @param data Data can be bulk data or single cell data
#' @return PALMO S4 object
#' @keywords createPALMOobject
#' @export
#' @examples
#' \dontrun{
#' palmo_obj=createPALMOobject(anndata, data)
#' }

createPALMOobject <- function(anndata, data) {

    ## Create PALMO object
    obj <- palmo(nDonors = 0)

    ## an object from the class PALMO
    obj@raw$ann <- anndata
    obj@raw$data <- data
    obj@rownames <- row.names(data)
    obj@colnames <- colnames(data)
    message(date(), ": The PALMO object is created\n")

    ## Creating output directory
    mainDir <- getwd()
    subDir <- "output"
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    obj@filePATH <- file.path(mainDir, subDir)
    message(date(), ": For outout files, the output directory is created\n")

    return(obj)
}
