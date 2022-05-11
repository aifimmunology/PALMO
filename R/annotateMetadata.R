#' annotateMetadata Function
#'
#' This function allows to add user-defined sample, participant, and time column
#' to a PALMO object in standard format.
#' @param data_object Input \emph{PALMO} S4 object. Contains annotation table
#' and expression matrix or data frame. Rows represent gene/proteins column
#' represents participant samples (same as annotation table Sample column)
#' @param sample_column Name of Sample column in user input annotation data
#' frame. Default 'Sample'
#' @param donor_column Name of Donor/participant column in user input annotation
#' data frame. Default 'PTID'
#' @param time_column Name of Time column in user input annotation data
#' frame. Default 'Time'
#' @param group_column Optional. Calculate average expression by given group
#' like 'celltype' or 'cluster'
#' @return PALMO object
#' @keywords annotateMetadata
#' @export
#' @examples
#' \dontrun{
#' annotateMetadata(data_object=palmo_obj, sample_column='Sample',
#' donor_column='PTID', time_column='Time')
#' }

annotateMetadata <- function(data_object, sample_column = "Sample",
                             donor_column = "PTID",
                             time_column = "Time",
                             group_column=NULL) {

    ## Get annotation data
    anndata <- data_object@raw$ann
    ## Add Sample, PTID and Time parameters
    anndata$Sample <- gsub("-", "_", anndata[, sample_column])
    anndata$PTID <- anndata[, donor_column]
    anndata$Time <- anndata[, time_column]

    ## Is NA
    if (sum(is.na(anndata$Sample)) > 0 | sum(is.na(anndata$PTID)) > 0 |
        sum(is.na(anndata$Time)) > 0) {
        warning(date(), ": Caution -> The Sample_column, Donor_column or
        Time_column contains missing value or NA. Missing data is removed.\n")
        anndata <- anndata[!is.na(anndata$Sample), ]
        anndata <- anndata[!is.na(anndata$PTID), ]
        anndata <- anndata[!is.na(anndata$Time), ]
    }

    ## Get the number of timepoints for each donor
    donor_tp <- data.frame(table(anndata$PTID))
    donor_tp <- donor_tp[donor_tp$Freq > 1, ]
    if (nrow(donor_tp) == 0) {
        stop(date(), ": Error -> Number of timepoints are less than 2.\n")
    }
    data_object@nDonors <- nrow(donor_tp)

    ## assign rownames with sample name
    row.names(anndata) <- make.names(anndata$Sample, unique = TRUE)

    ## Assign sample_group
    if(!is.null(group_column)) {
      anndata$group <- gsub("-", "_", anndata[, group_column])
      anndata$Sample_group <- paste(anndata$Sample, anndata$group, sep=":")
    }

    ## Assign annotation data to object
    data_object@curated$ann <- anndata

    return(data_object)
}
