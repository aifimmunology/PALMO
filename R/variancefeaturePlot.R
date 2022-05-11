#' variancefeaturePlot Function
#'
#' This function allows user to plot variance observed in the data by provided
#' featureList
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information
#' and expression data from Bulk or single cell data.
#' @param vardata Variance result obtained from lmeVariance function
#' @param featureSet Column of interest to focus on, Default is 'PTID'
#' @param Residual Add residual in plot, Default FALSE
#' @param top_n Number of top features to show. Default is 10.
#' @param cols The colors associated with features. Default is NULL.
#' @return variance plot list
#' @keywords variancefeaturePlot
#' @export
#' @examples
#' \dontrun{
#' variancefeaturePlot(data_object=palmo_obj, top_n=15)
#' }

variancefeaturePlot <- function(data_object = NULL, vardata = NULL, featureSet = "PTID", Residual = FALSE,
    top_n = 15, cols = NULL) {

    ## Check for input data
    if (!is.null(data_object)) {
        data <- data_object@result$variance_decomposition
    } else if (!is.null(vardata)) {
        data <- vardata
    } else {
        stop(date(), ": Please enter decomposition data frame or PALMO object.\n")
    }

    # If column of interest do not match with parameters
    check <- intersect(colnames(data), featureSet)
    if (length(check) != length(featureSet)) {
        stop(date(), ": Input featureSet do not match with data provided.\n")
    }

    if (Residual == TRUE) {
        featureList <- unique(c(featureSet, "Residual"))
    } else {
        featureList <- featureSet
    }

    # Result dataframe
    data <- data[, featureList]
    splots <- list()
    for (i in 1:length(featureList)) {
        column_oi <- featureList[i]
        data_sub <- data[order(data[, column_oi], decreasing = TRUE), ]

        # Top variables
        if (top_n > nrow(data_sub)) {
            data_sub <- melt(data.matrix(data_sub))
        } else {
            data_sub <- melt(data.matrix(data_sub[1:top_n, ]))
        }

        # Assign column names
        colnames(data_sub) <- c("feature", "variable", "value")

        # Orderby
        data_sub$variable <- factor(data_sub$variable, levels = rev(unique(c(column_oi, as.character(data_sub$variable)))))
        data_sub$feature <- factor(data_sub$feature, levels = rev(unique(as.character(data_sub$feature))))

        # color
        if (!is.null(cols)) {
            cols_list <- cols
        } else {
            cols_list <- rainbow(length(featureSet))
        }
        if (Residual == TRUE) {
            cols_list <- c(cols_list, "grey")  #Add grey color for Residual
        }
        names(cols_list) <- featureList

        # Column of interest specific plot
        plot <- ggplot(data_sub, aes(x = feature, y = value, fill = variable)) + geom_bar(stat = "identity",
            position = "stack") + labs(x = "Features", y = "% Variance explained", title = column_oi) +
            scale_fill_manual(values = cols_list) + theme_bw() + theme(axis.text.x = element_text(angle = 90,
            hjust = 0.5, vjust = 1), legend.position = "right") + coord_flip()
        print(plot)
        splots[[i]] <- plot
    }
    names(splots) <- featureSet

    # Return plots
    return(splots)
}
