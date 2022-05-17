#' naFilter Function
#'
#' This function allows users to filter genes/features having expression NAs
#' with user-defined threshold.
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information and expression data from Bulk or single cell data.
#' @param na_cutoff Threshold to remove NAs. Range is 0-1 (Default 0.4 means
#' >40% NAs are not allowed).
#' @return PALMO object with filtered NAs
#' @keywords naFilter
#' @export
#' @examples
#' \dontrun{
#' palmo_obj <- naFilter(data_object=palmo_obj, na_cutoff=0.4)
#' }

naFilter <- function(data_object, na_cutoff = 0.4) {

    data <- data_object@curated$data
    ## Remove genes with >40%NAs (optional)
    row.na <- apply(data, 1, function(x) {
        sum(is.na(x))
    })
    ## Select features with NAs <na_cutoff
    row.non.na <- row.na[row.na <= (na_cutoff * ncol(data))]
    data <- data[names(row.non.na), ]
    data_object@curated$data <- data
    return(data_object)
}
