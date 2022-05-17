#' createPALMOfromsinglecellmatrix Function
#'
#' This function allows to create Seurat object from counts and metadata as
#' mentioned in
#' https://search.r-project.org/CRAN/refmans/SeuratObject/html/CreateSeuratObject.html.
#' The seurat object then stored in a newly created PALMO object.
#' @param metadata Metadata associated with singe cell information. For example
#' rownames are unique cell_barcode and columns are information on each
#' cell_barcode like Sample (source of cell_barcode)
#' @param data Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample
#' column)
#' @param anndata Annotation dataframe. It consist of information such as
#' \emph{Sample} (sample name), \emph{PTID} (donor/participant), \emph{Time}
#' (longitudinal timepoints)
#' @return PALMO object with scRNA
#' @keywords createPALMOfromsinglecellmatrix
#' @export
#' @examples
#' \dontrun{
#' palmo_obj=createPALMOfromsinglecellmatrix(counts, metadata, annotation)
#' }

createPALMOfromsinglecellmatrix <- function(data, metadata, anndata = NULL) {

    ## Get the overlap between single cell annotations and expression dataframe
    ov <- intersect(row.names(metadata), colnames(data))
    metadata <- metadata[ov, ]
    data <- data[, ov]
    if (length(ov) > 1) {
        message(date(), "Number of cells overlapped=", length(ov))
        # Create seurat object
        scrna_obj <- CreateSeuratObject(counts = data, meta.data = metadata)
        # Create PALMO object
        palmoobj <- createPALMOobject(anndata = anndata, data = scrna_obj)
        return(palmoobj)
    } else {
        warning(date(), "Sample of matrix (column names) and annotation (row
                names) are not matching. Please check input data.")
    }

}
