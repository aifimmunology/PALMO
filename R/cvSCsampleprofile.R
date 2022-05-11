#' cvSCsampleprofile Function
#'
#' This function allows to calculate Intra-donor variations in single cell data
#' at sample level over longitudinal timepoints and visualize in a CV vs Mean
#' plot. Plots stored in output directory.
#' @param data_object Input \emph{PALMO} S4 object. Contains annotation table and
#' expression matrix or data frame. Rows represent gene/proteins column
#' represents participant samples (same as annotation table Sample column)
#' @param meanThreshold Average expression threshold to filter lowly expressed
#' genes Default is 0.1 (log2 scale)
#' @param cvThreshold Coefficient of variation threshold to select variable and
#' stable genes Default is 10 for single cell RNA (100*SD/mean)
#' @param cl Number of clusters. Use nCores-1 to run parallel. Default 2
#' @param plot_log10 Optional, Plot CV vs Mean on log10 scale. Default FALSE
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory \emph{PATH} Default, current
#' directory
#' @return PALMO object with CV list
#' @keywords cvSCsampleprofile
#' @examples
#' \dontrun{
#' palmo_obj <- cvSCsampleprofile(data_object=palmo_obj,
#' housekeeping_genes=c('GAPDH', 'ACTB'), fileName='scrna')
#' }

cvSCsampleprofile <- function(data_object, meanThreshold = NULL, cvThreshold = NULL, cl = 2,
    plot_log10 = FALSE, fileName = NULL, filePATH = NULL) {

    message(date(), ": Performing Sample-wise Coefficient of variance analysis\n")

    ## If filename or filepath null
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    if (is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }

    ## meanThrehold
    if (is.null(meanThreshold)) {
        meanThreshold <- 0
        message(date(), ": Using mean threshold >= 0\n")
    }
    data_object@meanThreshold <- meanThreshold
    ## cvThrehold
    if (is.null(cvThreshold)) {
        cvThreshold <- 10
        message(date(), ": Using CV threshold 10\n")
    }
    data_object@cvThreshold <- cvThreshold

    ## Get the data
    ann <- data_object@curated$anndata
    mat <- data_object@curated$data
    check_data <- all.equal(row.names(ann), colnames(mat))
    if (check_data == FALSE) {
        stop(date(), ": Annotation of samples (rows) and datamatrix columns do not match")
    }

    ## Calculate CV vs Mean for all genes per celltype
    unigene <- row.names(mat)
    uniSample <- sort(unique(ann$PTID))
    ann$group_donor <- paste(ann$group, ann$PTID, sep = ":")
    uniSamplegroup <- as.character(unique(ann$group_donor))

    ## CV-mean plot
    message(date(), ": Plotting Sample wise CV analysis\n")
    pdf(paste(filePATH, "/", fileName, "-CV-SampleGroup-Plot.pdf", sep = ""), width = 5, height = 5)
    op <- pboptions(type = "timer")  # default
    res1 <- pblapply(uniSamplegroup, function(uS) {
        # print(uS)
        ann_df <- ann[ann$group_donor %in% uS, ]
        if (nrow(ann_df) > 1) {
            df <- mat[unigene, ann_df$Sample_group]
            df <- data.frame(df, nonZero = apply(df, 1, function(x) {
                sum(x != 0)
            }), mean = rowMeans(df, na.rm = TRUE), sd = apply(df, 1, sd, na.rm = TRUE), var = apply(df,
                1, var, na.rm = TRUE), stringsAsFactors = FALSE)
            df$CV <- 100 * df$sd/df$mean
            # the CV becomes very high for data with 0
            df <- df[df$mean >= meanThreshold, ]  #minimum expression >2^0.1=1
            dp2a <- df[df$mean >= meanThreshold & df$CV > cvThreshold, c("mean", "sd", "var",
                "CV")]
            dp2a <- dp2a[order(dp2a$CV, dp2a$mean, decreasing = TRUE), ]
            dp2b <- df[df$mean >= meanThreshold & df$CV < cvThreshold, ]
            dp2b <- dp2b[order(-dp2b$mean, dp2b$CV, decreasing = FALSE), ]
            plot1 <- ggplot(df, aes(x = mean, y = CV)) + geom_point(size = 0.5, color = "grey") +
                labs(title = paste(uS, " timepoints=", nrow(ann_df), sep = "")) + geom_text_repel(data = dp2a[1:10,
                ], aes(x = mean, y = CV, label = row.names(dp2a[1:10, ])), col = "red", size = 2,
                max.overlaps = 20) + geom_text_repel(data = dp2b[1:10, ], aes(x = mean, y = CV,
                label = row.names(dp2b[1:10, ])), col = "blue", size = 2, max.overlaps = 20) +
                theme_classic()
            if (plot_log10 == TRUE) {
                plot1 <- plot1 + scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10")
            }
            print(plot1)
        }
        return(NULL)
    })
    pboptions(op)
    dev.off()

    message(date(), ": Done. Please check output directory for results.\n")

}
