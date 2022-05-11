#' sample_correlation Function
#'
#' This function allows to perform sample correlation (by group like
#' celltype, or by donor).
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation information
#' and expression data from Bulk or single cell data.
#' @param donor_sep Sample and group(celltype) separator like (donor_sep=':')
#' @param cor_method (Optional) Correlation method 'pearson' or 'spearman'. Default
#' is 'spearman'
#' @param group_position (Optional) Position of group 1st or 2nd like (group_position=2, in
#' Donor:group)
#' @param group_by Cluster correlation heatplot by 'donor' or 'group'
#' @param col_break Value between 0-1.
#' @param col_max Maximum color limit (Default 1)
#' @param column_names_fontsize Font size of the column names, Default 4
#' @param row_names_fontsize Font size of the row names, Default 4
#' @param column_title_fontsize Font size of the column title, Default 6
#' @param row_title_fontsize Font size of the row title, Default 6
#' @param plotHeight Height of the plot in inch, Default 20 in
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory \emph{PATH} Default, current
#' directory
#' @return PALMO object with correlation cor_res dataframe
#' @keywords sample_correlation
#' @export
#' @examples
#' \dontrun{
#' palmo_obj <- sample_correlation(data_object=palmo_obj, donor_sep=':')
#' }

sample_correlation <- function(data_object, donor_sep = ":", cor_method = "spearman", group_position = 2,
    group_by = "donor", col_break = NULL, col_max = 1, column_names_fontsize = 4, row_names_fontsize = 4,
    row_title_fontsize = 6, column_title_fontsize = 6, plotHeight = 20, fileName = NULL, filePATH = NULL) {

    message(date(), ": Performing Sample correlation anlaysis\n")
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    if (is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }
    data_object@donor_sep <- donor_sep
    data_object@cor_method <- cor_method

    # Create annotation
    data <- data_object@curated$data
    data_annotation <- data.frame(id = colnames(data))
    temp <- data.frame(do.call(rbind, strsplit(data_annotation$id, split = donor_sep)), stringsAsFactors = FALSE)

    if (ncol(temp) > 1) {
        data_annotation <- data.frame(data_annotation, Sample = temp$X1, group = temp$X2, stringsAsFactors = FALSE)
        row.names(data_annotation) <- data_annotation$id
        if (group_position == 1) {
            colnames(data_annotation) <- c("sample", "group", "donor")
        } else {
            colnames(data_annotation) <- c("sample", "donor", "group")
        }

        # Sample Correlation map
        cor_mat <- rcorr(as.matrix(data), type = cor_method)
        res <- cor_mat$r

        # color break
        if (is.null(col_break)) {
            col_bar <- colorRamp2(c(-col_max, 0, col_max), c("blue", "white", "red"))
        } else {
            col_break <- abs(col_break)
            col_bar <- colorRamp2(c(-col_max, -col_break, 0, col_break, col_max), c("blue", "skyblue",
                "white", "pink", "red"))
        }
        # Plot heatmap
        if (group_by == "donor") {
            ha_col <- HeatmapAnnotation(df = data.frame(donor = data_annotation$donor))
            ht1 <- Heatmap(data.matrix(res), cluster_rows = FALSE, cluster_columns = FALSE, row_split = as.character(data_annotation$donor),
                column_split = as.character(data_annotation$donor), na_col = "grey", col = col_bar,
                row_names_max_width = unit(10, "cm"), column_names_gp = gpar(fontsize = column_names_fontsize),
                row_names_gp = gpar(fontsize = row_names_fontsize), top_annotation = ha_col,
                row_title_gp = gpar(fontsize = row_title_fontsize), column_title_gp = gpar(fontsize = column_title_fontsize),
                heatmap_legend_param = list(title = cor_method, heatmap_legend_side = "right"))
        } else {
            ha_col <- HeatmapAnnotation(df = data.frame(Group = data_annotation$group))
            ht1 <- Heatmap(data.matrix(res), cluster_rows = FALSE, cluster_columns = FALSE, row_split = as.character(data_annotation$group),
                column_split = as.character(data_annotation$group), na_col = "grey", col = col_bar,
                row_names_max_width = unit(10, "cm"), column_names_gp = gpar(fontsize = column_names_fontsize),
                row_names_gp = gpar(fontsize = row_names_fontsize), top_annotation = ha_col,
                row_title_gp = gpar(fontsize = row_title_fontsize), column_title_gp = gpar(fontsize = column_title_fontsize),
                heatmap_legend_param = list(title = cor_method, heatmap_legend_side = "right"))
        }
    } else {

        row.names(data_annotation) <- data_annotation$id
        # Sample Correlation map
        cor_mat <- rcorr(as.matrix(data), type = cor_method)
        res <- cor_mat$r
        # Plot heatmap
        ht1 <- Heatmap(data.matrix(res), cluster_rows = FALSE, cluster_columns = FALSE, na_col = "grey",
            col = colorRamp2(c(-1, -0.9, 0, 0.9, 1), c("blue", "skyblue", "white", "pink", "red")),
            row_names_max_width = unit(10, "cm"), column_names_gp = gpar(fontsize = column_names_fontsize),
            row_names_gp = gpar(fontsize = row_names_fontsize), heatmap_legend_param = list(title = cor_method,
                heatmap_legend_side = "right"))

    }

    pdf(paste(filePATH, "/", fileName, "-SampleCorrelationplot.pdf", sep = ""), width = plotHeight,
        height = plotHeight)
    draw(ht1)
    dev.off()
    print(draw(ht1))

    data_object@result$sample_cor <- res

    message(date(), ": Please check output directory for results\n")
    return(data_object)
}
