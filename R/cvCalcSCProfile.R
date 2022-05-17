#' cvCalcSCProfile Function
#'
#' This function allows to calculate Intra-donor variations in single cell data
#' over longitudinal timepoints and visualize in a CV vs Mean plot. Plots stored
#' in output directory.
#' @param data_object Input \emph{PALMO} S4 object. Contains annotation table and
#' expression matrix or data frame. Rows represent gene/proteins column
#' represents participant samples (same as annotation table Sample column)
#' @param meanThreshold Average expression threshold to filter lowly expressed
#' genes Default is 0.1 (log2 scale)
#' @param housekeeping_genes Optional, vector of housekeeping genes. Default is
#' c('ACTB', 'GAPDH')
#' @param cl Number of clusters. Use nCores-1 to run parallel. Default 2
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory \emph{PATH} Default, current
#' directory
#' @return PALMO object with CV profile list
#' @keywords cvCalcSCProfile
#' @examples
#' \dontrun{
#' palmo_obj <- cvCalcSCProfile(data_object=palmo_obj,
#' housekeeping_genes=c('GAPDH', 'ACTB'), fileName='scrna')
#' }

cvCalcSCProfile <- function(data_object, meanThreshold = NULL,
                            housekeeping_genes = NULL, cl = 2,
                            fileName = NULL, filePATH = NULL) {

    message(date(), ": Performing Coefficient of variance analysis")

    ## If filename or filepath null
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    if (is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }
    ## Assign housekeeping_genes
    if (is.null(housekeeping_genes)) {
        housekeeping_genes <- c("ACTB", "GAPDH")
        data_object@housekeeping_genes <- housekeeping_genes
    }

    ## meanThrehold
    if (is.null(meanThreshold)) {
        meanThreshold <- 0
        message(date(), ": Using mean threshold >= 0")
    }
    data_object@meanThreshold <- meanThreshold

    ## Get the data
    ann <- data_object@curated$anndata
    mat <- data_object@curated$data
    check_data <- all.equal(row.names(ann), colnames(mat))
    if (check_data == FALSE) {
        stop(date(), ": Annotation of samples (rows) and datamatrix columns
             do not match")
    }

    ## Calculate CV vs Mean for all genes per celltype
    unigene <- row.names(mat)
    uniSample <- sort(unique(ann$PTID))
    ann$group_donor <- paste(ann$group, ann$PTID, sep = ":")
    uniSamplegroup <- as.character(unique(ann$group_donor))

    ## All genes CV calculations
    message(date(), ": Performing CV calculations")
    op <- pboptions(type = "timer")  # default
    res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
        # print(uS)
        ann_df <- ann[ann$group_donor %in% uS, ]
        if (nrow(ann_df) > 1) {
            df <- mat[unigene, ann_df$Sample_group]
            df <- data.frame(df,
                    zeros = apply(df, 1, function(x) { sum(x != 0) }),
                    mean = rowMeans(df, na.rm = TRUE),
                    sd = apply(df, 1, sd, na.rm = TRUE),
                    var = apply(df, 1, var, na.rm = TRUE),
                    stringsAsFactors = FALSE)
            df$cv <- 100 * df$sd/df$mean
            df$cv <- ifelse(df$mean >= meanThreshold, df$cv, NA)
            df <- df[, c("mean", "sd", "var", "cv")]
            df$gene <- row.names(df)
            df$group <- uS
            return(df)
        }
    })
    pboptions(op)
    cv_all <- do.call(rbind, res)
    cv_all <- data.frame(cv_all, check.names = FALSE, stringsAsFactors = FALSE)
    cv_all$select <- ifelse(cv_all$mean >= meanThreshold, "Y", "N")

    ## Plot results
    df <- cv_all[cv_all$mean >= meanThreshold, ]
    p1 <- ggplot(cv_all, aes(x = mean)) +
        geom_histogram(aes(color = select), fill = "white", binwidth = 0.1) +
        labs(title = "Mean expression (log10)") +
        scale_x_continuous(trans = "log10")
    p2 <- ggplot(df, aes(x = cv)) +
        labs(title = "CV (mean/SD %)") +
        geom_histogram(binwidth = 1, color = "black", fill = "white")
    ## Housekeeping genes data
    df <- df[df$gene %in% housekeeping_genes, ]
    if (nrow(df) > 0) {
        p3 <- ggplot(df, aes(x = mean, y = cv)) +
            geom_point() + labs(title = "Housekeeping genes") +
            facet_wrap(~gene)

        pdf(paste(filePATH, "/", fileName, "-CVPlot.pdf", sep = ""),
            width = 12, height = 5)
        print(plot_grid(p1, p2, p3, ncol = 3))
        dev.off()
        print(plot_grid(p1, p2, p3, ncol = 3))
    } else {
        pdf(paste(filePATH, "/", fileName, "-CVPlot.pdf", sep = ""),
            width = 8, height = 5)
        print(plot_grid(p1, p2, ncol = 2))
        dev.off()
        print(plot_grid(p1, p2, ncol = 2))
    }

    ## Add CV result
    data_object@result$cv_all <- cv_all

    message(date(), ": Done. Please check output directory for Plots/results.")
    return(data_object)
}
