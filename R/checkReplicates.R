#' checkReplicates Function
#'
#' This function allows you to check for any replicates in data. If present then
#' merge expression of samples by median provided mergeReplicates=TRUE
#' @param data_object Input \emph{PALMO} S4 object. Contains annotation table
#' and expression matrix or data frame. Rows represent gene/proteins column
#' represents participant samples (same as annotation table Sample column)
#' @param mergeReplicates Merge replicates expression data by Median.
#' Default FALSE
#' @return PALMO object with merged replicates
#' @keywords checkReplicates
#' @export
#' @examples
#' \dontrun{
#' palmo_obj=checkReplicates(data_object=palmo_obj, mergeReplicates=TRUE)
#' }

checkReplicates <- function(data_object, mergeReplicates = FALSE) {

    ## Check whether replicates are TRUE
    anndata <- data_object@curated$anndata
    mat <- data_object@curated$data

    ## Check Replicates
    ann <- anndata
    ann$PTID_Time <- paste(ann$PTID, ann$Time, sep = "_")
    ann_temp <- ann[!duplicated(ann$Sample), ]
    rep <- data.frame(table(ann_temp$PTID, ann_temp$Time))
    rep <- rep[rep$Freq > 1, ]
    if (nrow(rep) > 0) {
        message(date(), ": Replicates for data found")
        print(rep)
        rep_PTID_Time <- unique(paste(rep$Var1, rep$Var2, sep = "_"))

        ## Is Group information available
        group_info <- intersect(colnames(ann), "group")
        if (length(group_info) == 0) {
            ann$group <- ann$Time
            ann$Sample_group <- ann$PTID_Time
        }

        if (mergeReplicates == TRUE) {
            message(date(), ": Merging replicates by Median")
            ann_rep <- ann[ann$PTID_Time %in% rep_PTID_Time, ]
            PTID_Time <- unique(ann_rep$PTID_Time)
            ann$Sample_old <- ann$Sample
            ann$Sample_group_old <- ann$Sample_group
            ## Create unique annotations
            for (i in 1:length(PTID_Time)) {
                pt <- PTID_Time[i]
                ann_rep_sub <- ann_rep[ann_rep$PTID_Time %in% pt, ]
                ann[ann$PTID_Time %in% pt, ]$Sample <- pt
                if (length(group_info) == 0) {
                  ann[ann$PTID_Time %in% pt, ]$Sample_group <- ann[ann$PTID_Time %in% pt, ]$Sample
                } else {
                  ann[ann$PTID_Time %in% pt, ]$Sample_group <- paste(pt, ann[ann$PTID_Time %in%
                    pt, ]$group, sep = ":")
                }
            }

            ## Check for columns
            check <- all.equal(row.names(ann), colnames(mat))
            if (check == FALSE) {
                mat <- mat[, row.names(ann)]
            }

            ## Aggregate (replicates) expression using median
            uni_sample_group <- unique(ann$Sample_group)
            med_agg_res <- lapply(uni_sample_group, function(sg) {
                ann_sub <- ann[ann$Sample_group %in% sg, ]
                mat_sub <- mat[, row.names(ann_sub)]
                if (nrow(ann_sub) > 1) {
                  mat_agg <- apply(mat_sub, 1, function(x) {
                    median(x, na.rm = TRUE)
                  })
                } else {
                  mat_agg <- mat_sub
                }
                mat_agg <- data.frame(mat_agg, check.names = FALSE,
                                      stringsAsFactors = FALSE)
                colnames(mat_agg) <- sg
                return(mat_agg)
            })
            med_agg_res <- do.call(cbind, med_agg_res)
            ## Aggregated result
            med_agg_res <- data.frame(med_agg_res, check.names = FALSE,
                                      stringsAsFactors = FALSE)

            ## Check annotation data
            ann <- ann[!duplicated(ann$Sample_group), ]
            row.names(ann) <- ann$Sample_group
            check <- all.equal(row.names(ann), colnames(med_agg_res))
            if (check == FALSE) {
                med_agg_res <- med_agg_res[, row.names(ann)]
            }

            ## Add data to PALMO object
            data_object@curated$ann_old <- anndata
            data_object@curated$data_old <- mat
            data_object@curated$anndata <- ann
            data_object@curated$data <- med_agg_res
            message(date(), ": Check PALMO object (anndata, data). To ignore
                    replicates use mrgereplicates=FALSE.")
        } else {
            warning(date(), ": Merging Replicates ignored.")
        }
    } else {
        message(date(), ": No Replicates found")
    }

    return(data_object)
}
