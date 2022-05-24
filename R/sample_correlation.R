#' sample_correlation Function
#'
#' This function allows to perform sample correlation (by group like
#' celltype, or by donor).
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation information
#' and expression data from Bulk or single cell data.
#' @param cor_method (Optional) Correlation method 'pearson' or 'spearman'. Default
#' is 'spearman'
#' @param group_by Cluster correlation heat plot by 'donor' or 'group'
#' @param col_break Value between 0 and 1
#' @param col_max Maximum color limit (Default 1)
#' @param cluster_rows \emph{ComplexHeatmap} cluster rows, Default FALSE
#' @param cluster_columns \emph{ComplexHeatmap} cluster columns, Default FALSE
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
#' palmo_obj <- sample_correlation(data_object=palmo_obj, group_by="Time")
#' }

sample_correlation <- function(data_object,
                               cor_method = "spearman",
                               group_by = NULL,
                               col_break = NULL,
                               col_max = 1,
                               cluster_rows = FALSE,
                               cluster_columns = FALSE,
                               column_names_fontsize = 4,
                               row_names_fontsize = 4,
                               row_title_fontsize = 6,
                               column_title_fontsize = 6,
                               plotHeight = 20,
                               fileName = NULL, filePATH = NULL) {

    message(date(), ": Performing Sample correlation anlaysis")
    if(is.null(fileName)) {
        fileName <- "outputFile"
    }
    if(is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }
    data_object@cor_method <- cor_method

    #Create annotation
    ann <- data_object@curated$anndata
    data <- data_object@curated$data

    #Sample Correlation map
    cor_mat <- rcorr(as.matrix(data), type = cor_method)
    cor_res <- cor_mat$r

    # color break
    if(is.null(col_break)) {
      col_bar <- colorRamp2(c(-col_max, 0, col_max),
                            c("blue", "white", "red"))
    } else {
      col_break <- abs(col_break)
      col_bar <- colorRamp2(c(-col_max, -col_break, 0,
                              col_break, col_max),
                            c("blue", "skyblue", "white", "pink", "red"))
    }

    if(!is.null(group_by)) {

      # Plot heatmap
      ha_col <- HeatmapAnnotation(df=data.frame(groupBy=ann[,group_by]),
                                  show_legend = FALSE)
      ht1 <- Heatmap(data.matrix(cor_res),
                     cluster_rows = cluster_rows,
                     cluster_columns = cluster_columns,
                     row_split = as.character(ann[,group_by]),
                     column_split = as.character(ann[,group_by]),
                     na_col = "grey", col = col_bar,
                     row_names_max_width = unit(10, "cm"),
                     column_names_gp=gpar(fontsize=column_names_fontsize),
                     row_names_gp = gpar(fontsize = row_names_fontsize),
                     top_annotation = ha_col,
                     row_title_gp = gpar(fontsize = row_title_fontsize),
                     column_title_gp=gpar(fontsize=column_title_fontsize),
                     heatmap_legend_param = list(title = cor_method,
                            heatmap_legend_side = "right"))
      } else {
        # Plot heatmap
        ha_col <- HeatmapAnnotation(df=data.frame(Donor=ann$PTID,
                                                  Time=ann$Time),
                                    show_legend = FALSE)
        ht1 <- Heatmap(data.matrix(cor_res),
                       cluster_rows = cluster_rows,
                       cluster_columns = cluster_columns,
                       row_split = as.character(ann$PTID),
                       column_split = as.character(ann$PTID),
                       na_col = "grey", col = col_bar,
                       row_names_max_width = unit(10, "cm"),
                       column_names_gp=gpar(fontsize=column_names_fontsize),
                       row_names_gp = gpar(fontsize = row_names_fontsize),
                       top_annotation = ha_col,
                       row_title_gp = gpar(fontsize = row_title_fontsize),
                       column_title_gp=gpar(fontsize=column_title_fontsize),
                       heatmap_legend_param = list(title = cor_method,
                                  heatmap_legend_side = "right"))
    }

    pdf(paste(filePATH, "/", fileName, "-SampleCorrelationplot.pdf", sep = ""),
        width = plotHeight, height = plotHeight)
    draw(ht1)
    dev.off()
    print(draw(ht1))

    data_object@result$sample_cor <- cor_mat

    message(date(), ": Please check output directory for results")
    return(data_object)
}
