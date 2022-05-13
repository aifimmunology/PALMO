#' genecircosPlot Function
#'
#' This function allows to Circos Plot for gene list of interest by group
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information and expression data from Bulk or single cell data. Rows
#' represents gene/proteins column represents group:donor (or donor:group)
#' @param data Expression matrix or data frame. Rows represents gene/proteins
#' column represents group:donor (or donor:group)
#' @param geneList User-defined genes of interest
#' @param group_position Default 1, use 2 when columns are donor:group format
#' @param group_oi Optional, User-defined groups to consider and order plot
#' @param titleName Title of the plot
#' @param colorThreshold User-defined color threshold (same as cvThreshold,
#' like 5)
#' @param colorMax Maximum CV value in heatplot ("max", numeric or NULL)
#' @param colorscale Show color scale, TRUE or FALSE (default).
#' @return Circos plots and dataframe
#' @keywords genecircosPlot
#' @export
#' @examples
#' \dontrun{
#' genecircosPlot(data_object=palmo_obj, geneList=c('IL32','CCL5','TCF7'))
#' }

genecircosPlot <- function(data = NULL, data_object = NULL, geneList,
                           group_position = 1, group_oi = NULL,
                           titleName = "", colorThreshold = 10,
                           colorMax=NULL, colorscale=FALSE) {

    ## If PALMO object not available
    if (!is.null(data_object)) {
        data <- data_object@result$cv_all
    }

    ## Split the group and donor
    data_annotation <- data.frame(do.call(rbind, strsplit(colnames(data),
                                split = ":")), stringsAsFactors = FALSE)
    data_annotation$id <- colnames(data)
    if (group_position == 1) {
        colnames(data_annotation) <- c("group", "donor", "sample")
    } else {
        colnames(data_annotation) <- c("donor", "group", "sample")
    }

    ## Group of interest
    if (is.null(group_oi)) {
        group_oi <- sort(unique(data_annotation$group))
        Sample_group <- unique(data_annotation$sample)
    } else {
        group_oi <- group_oi
        Sample_group <- data_annotation$sample
    }

    ## Filter group_oi
    data_annotation <- data_annotation[data_annotation$group %in% group_oi, ]
    data_annotation <- data_annotation[order(data_annotation$sample), ]

    ## Max color scale (CV)
    if(is.null(colorMax)) {
      colorMax <- 50
    } else if(colorMax == "max") {
      colorMax <- max(data, na.rm=TRUE)
    } else {
      colorMax <- as.numeric(colorMax)
    }

    ## Show heatplot
    data_mat <- data[geneList, data_annotation$sample]
    ht1 <- Heatmap(data.matrix(data_mat),
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   column_split=factor(data_annotation$group, levels=group_oi),
                   na_col = "grey",
                   col = colorRamp2(c(0, colorThreshold, (colorThreshold + 0.1),
                            colorMax), c("blue", "white", "pink", "red")),
                   row_names_max_width = unit(10, "cm"),
                   column_names_gp = gpar(fontsize = 5),
                   row_names_gp = gpar(fontsize = 6),
                   heatmap_legend_param = list(title = "CV",
                                heatmap_legend_side = "right")
                   )
    print(ht1)

    ## Define order
    column_ht <- unlist(column_order(ht1))
    data_mat <- data_mat[, column_ht]
    data_annotation <- data_annotation[column_ht, ]
    output_data_mat <- data_mat
    # For visualization change NAs with 0
    data_mat[is.na(data_mat)] <- 0

    data_mat <- t(data_mat)
    col_fun1 = colorRamp2(c(0, 0.001, colorThreshold, (colorThreshold + 0.1),
                    colorMax), c("grey", "blue", "white", "pink", "red"))
    lastSector <- length(unique(data_annotation$group))

    ## Show circos plot
    circos.clear()
    circos.par(start.degree = 90, gap.after = c(rep(2, lastSector - 1), 40))
    # circos.par(start.degree = 90, gap.degree = 5)
    circos.heatmap(data_mat,
                   col = col_fun1,
                   split = data_annotation$group,
                   rownames.side = "outside",
                   na.col = "grey",
                   bg.border = "black")  #, track.height = 0.4)
    circos.track(track.index = get.current.track.index(),
                 panel.fun = function(x, y) {
        if (CELL_META$sector.numeric.index == lastSector) {
            # the last sector
            cn = rev(colnames(data_mat))
            n = length(cn)
            circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                        1:n - 0.5, cn, cex = 0.5, adj = c(0, 0.5),
                        facing = "inside")
        }
    }, bg.border = NA)
    circos.clear()
    text(0, 1, titleName)

    if(colorscale == TRUE) {
      lgd1 = Legend(title = "Exp", col_fun = col_fun1)
      grid.draw(lgd1)
    }

    return(output_data_mat)

}
