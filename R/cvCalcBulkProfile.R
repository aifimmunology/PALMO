#' cvCalcBulkProfile Function
#'
#' This function allows to calculate Intra-donor variations in bulk data over
#' longitudinal timepoints and visualize in a CV vs Mean plot. Plots stored in
#' output directory.
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information and expression data from Bulk or single cell data.
#' @param fileName User-defined filename, Default outputFile
#' @param filePATH User-defined output directory \emph{PATH} Default, current directory
#' @param cl Number of clusters. Use nCores-1 to run parallel. Default 2
#' @return PALMO object with CV profile cv_all
#' @keywords cvCalcBulkProfile
#' @export
#' @examples
#' \dontrun{
#' cvCalcBulkProfile(data_object=palmo_obj)
#' }

cvCalcBulkProfile <- function(data_object, cl = 2, fileName = NULL, filePATH = NULL) {

    message(date(), ": Performing Coefficient of variance analysis\n")
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    if (is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }

    ## get the data
    ann <- data_object@curated$anndata
    mat <- data_object@curated$data
    check_data <- all.equal(row.names(ann), colnames(mat))
    if (check_data == FALSE) {
        stop(date(), ": Annotation of samples (rows) and datamatrix columns do not match")
    }

    # CV vs Mean
    unigene <- row.names(mat)
    uniSample <- sort(unique(ann$PTID))

    message(date(), ": Performing CV calculations\n")
    op <- pboptions(type = "timer")  # default
    res <- pblapply(uniSample, cl = cl, function(uS) {
        # print(uS)
        meta_df <- ann[ann$PTID %in% uS, ]
        if (nrow(meta_df) > 1) {
            df <- mat[unigene, meta_df$Sample]
            df <- data.frame(df, NAs = apply(df, 1, function(x) {
                sum(is.na(x))
            }), zeros = apply(df, 1, function(x) {
                sum(x == 0)
            }), mean = rowMeans(df, na.rm = TRUE), sd = apply(df, 1, sd, na.rm = TRUE), var = apply(df,
                1, var, na.rm = TRUE), stringsAsFactors = FALSE)
            df$CV <- 100 * df$sd/df$mean
            df <- df[, c("mean", "sd", "var", "CV", "NAs", "zeros")]
            df$feature <- row.names(df)
            df$group <- uS
            return(df)
        }
    })
    pboptions(op)
    cv_all <- do.call(rbind, res)
    cv_all <- data.frame(cv_all, check.names = FALSE, stringsAsFactors = FALSE)

    ## Add CV result
    data_object@result$cv_all <- cv_all

    ## histogram of CV
    plot1 <- ggplot(cv_all, aes(x = mean, y = CV)) + geom_point(size = 0.5, color = "grey") +
        scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") + facet_wrap(~group) +
        theme_classic()
    print(plot1)

    message(date(), ": Done\n")
    return(data_object)
}
