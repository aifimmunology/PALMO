#' gene_featureplot Function
#'
#' This function allows to output the user-defined input features expression
#' in graphical format. Users can select x-axis as donor/participant
#' (x_group_by='PTID') and expression on y-axis organized by variable time
#' (var_oi='Time'). Add group facet feature like facet_by='celltype'.
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information and expression data from Bulk or single cell data.
#' @param anndata Optional, Annotation dataframe consist of information such as
#' Sample (sample name), PTID (donor/participant), Time (longitudinal timepoints)
#' @param data Optional, Data can be bulk data or single cell data
#' @param featureList User-defined feature or genelist as a vector
#' @param x_group_by x-axis grouping variable like 'PTID'
#' @param var_oi x-axis subgrouping variable like 'Time'
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ncol Number of columns in the plot grid
#' @param facet_by A set of variables or expressions
#' @param compare_means Add mean comparison p-value in a plot (for more
#' information refer http://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html)
#' @param x_text_angle xaxis text angle on ggplot
#' @param text_font font size on ggplot
#' @return gene plot
#' @keywords gene_featureplot
#' @export
#' @examples
#' \dontrun{
#' plots <- gene_featureplot(data_object=palmo_obj,
#' featureList=c('LILRA4', 'CLEC9A'))
#' }

gene_featureplot <- function(data_object = NULL, data = NULL, anndata = NULL, featureList, x_group_by = "PTID",
    var_oi = "Time", xlab = "group_by", ylab = "Value/Expression", ncol = 2, facet_by = NULL,
    compare_means = FALSE, x_text_angle = NULL, text_font = NULL) {

    ## Check for input data
    if (!is.null(data_object)) {
        data <- data_object@curated$data
        anndata <- data_object@curated$anndata
    }

    if (is.null(data) & is.null(data_object)) {
        stop(date(), ": Please enter decomposition data frame or PALMO object.\n")
    }

    ## get the data
    ann <- anndata
    mat <- data
    check_data <- all.equal(row.names(ann), colnames(mat))
    if (check_data == FALSE) {
        stop(date(), ": Annotation of samples (rows) and datamatrix columns do not match")
    }

    ## Check genes
    featureList <- intersect(featureList, row.names(mat))
    uniSample <- as.character(unique(ann[, x_group_by]))
    splots <- list()
    for (i in 1:length(featureList)) {
        featureName <- featureList[i]
        df <- data.frame(exp = as.numeric(mat[featureName, ]), ann, stringsAsFactors = FALSE)
        df$group_oi <- df[, x_group_by]
        df$group_oi <- factor(df$group_oi, levels = uniSample)
        df$var_oi <- factor(df[, var_oi], levels = unique(df[, var_oi]))

        p1 <- ggplot(df, aes(x = group_oi, y = exp)) + geom_boxplot(outlier.shape = NA) + labs(x = xlab,
            y = ylab, title = featureName) + geom_jitter(width = 0.2, aes(shape = var_oi)) +
            scale_shape_manual(values = 1:nlevels(df$var_oi)) + theme_classic()  #+ theme(legend.title = element_blank())

        ## number of group
        nGroup <- length(unique(df$group_oi))
        if (nGroup == 2 & compare_means == TRUE) {
            p1 <- p1 + ggpubr::stat_compare_means(method = "wilcox.test")
        } else if (nGroup > 2 & compare_means == TRUE) {
            p1 <- p1 + ggpubr::stat_compare_means(method = "anova")
        }

        if (!is.null(facet_by)) {
            p1 <- p1 + facet_wrap(~df[, facet_by])
        }
        if (!is.null(x_text_angle)) {
            p1 <- p1 + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1),
                legend.position = "right")
        }
        if (!is.null(text_font)) {
            p1 <- p1 + theme(text = element_text(size = text_font))
        }
        splots[[i]] <- p1
        print(p1)
    }
    names(splots) <- featureList

    return(splots)
}
