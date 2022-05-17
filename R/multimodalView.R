#' multimodalView Function
#'
#' This function allows to visualize the multimodal view genes of interest
#' by celltypes/ groups defined by use
#' @param modality1 Variation or Expression matrix/data frame. Rows represents
#' gene/proteins column represents group:donor (group and donor separated by :)
#' @param modality2 Variation or Expression matrix/data frame. Rows represents
#' gene/proteins column represents group:donor (group and donor separated by :)
#' @param geneList Genes of interest to explore
#' @param group_position Default 1, use 2 when columns are donor:group format
#' @param group_oi Optional, User-defined groups to consider and order
#' @param colorThreshold User-defined color threshold in color space
#' @param colorMax Maximum CV value in heatplot ("max", numeric or NULL)
#' @param colorscale Show color scale, TRUE or FALSE (default).
#' @param plotHeight User-defined Plot size (in)
#' @param titleName Title of the plot
#' @param fileName User defined filename
#' @param filePATH User-defined output directory \emph{PATH} to save result
#' @return Multimodal plot and data list
#' @keywords multimodalView
#' @export
#' @examples
#' \dontrun{
#' multimodalView(modality1=scrna_cv_res, modality2=scatac_cv_res, geneList)
#' }

multimodalView <- function(modality1, modality2, group_oi = NULL,
                           geneList, colorThreshold = 10,
                           group_position = NULL, plotHeight = 10,
                           titleName = "",
                           colorMax=NULL, colorscale=FALSE,
                           fileName = NULL, filePATH = NULL) {

    ## Define filename
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    ## Define output path
    if (is.null(filePATH)) {
        filePATH <- paste(getwd(), "/output", sep = "")
        dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
    }
    ## Define color threshold
    if (is.null(colorThreshold)) {
        colorThreshold <- 10
    }
    ## Define group column
    if (is.null(group_position)) {
        group_position <- 1
    }

    ## Overlap
    ov <- intersect(colnames(modality1), colnames(modality2))
    if(length(ov) == 0) {
      stop(date(), "Error: Column names do not match")
    }

    ## Intersect column and add which were missing in modality2
    cn <- colnames(modality1)
    not_in_cn <- setdiff(cn, colnames(modality2))
    mat <- data.frame(matrix(NA, ncol = length(not_in_cn),
                             nrow = length(row.names(modality2))))
    colnames(mat) <- not_in_cn
    modality2 <- data.frame(modality2, mat, check.names = FALSE,
                            stringsAsFactors = FALSE)
    modality2 <- modality2[, cn]

    ## Split column names by celltypes (:)
    data_annotation <- data.frame(do.call(rbind, strsplit(cn, split = ":")),
                                  stringsAsFactors = FALSE)
    data_annotation$sample <- cn
    if (group_position == 1) {
        colnames(data_annotation) <- c("group", "donor", "sample")
    } else {
        colnames(data_annotation) <- c("donor", "group", "sample")
    }

    ## IF celltype of interest are not given
    if (is.null(group_oi)) {
        group_oi <- unique(data_annotation$group)
    }

    ## Select group
    data_annotation <- data_annotation[data_annotation$group %in% group_oi, ]
    Samplecelltype <- data_annotation$sample

    ## Select genes
    geneList1 <- intersect(geneList, row.names(modality1))
    geneList2 <- intersect(geneList, row.names(modality2))
    geneList <- unique(geneList1, geneList2)

    ## Get the data
    mat1 <- t(modality1[geneList, Samplecelltype])
    mat2 <- t(modality2[geneList, Samplecelltype])
    colnames(mat1) <- geneList
    colnames(mat2) <- geneList
    output_mat1 <- mat1
    output_mat2 <- mat2
    ## For visualization change NAs with 0
    mat1[is.na(mat1)] <- 0
    mat2[is.na(mat2)] <- 0

    ## Max color scale (CV)
    if(is.null(colorMax)) {
      colorMax <- 50
    } else if(colorMax == "max") {
      colorMax <- max(c(max(mat1, na.rm=TRUE),
                        max(mat2, na.rm=TRUE)),
                      na.rm=TRUE)
    } else {
      colorMax <- as.numeric(colorMax)
    }

    ## Define color function
    col_fun1 = colorRamp2(c(0, 1e-04, colorThreshold, (colorThreshold + 0.1),
              colorMax), c("grey", "blue", "white", "pink", "red"))
    lastSector <- length(unique(data_annotation$group))

    ## Print Plot
    circos.clear()
    pdf(file = paste(filePATH, "/", fileName, "-Multimodal-circularplot.pdf",
                     sep = ""), width = plotHeight, height = plotHeight)
    circos.par(start.degree = 90, gap.after = c(rep(2, lastSector - 1), 40))
    # circos.par(start.degree = 90, gap.degree = 5)
    circos.heatmap(mat1, col = col_fun1,
                   split = factor(data_annotation$group, levels = group_oi),
                   rownames.side = "outside",
                   na.col = "grey", bg.border = "black")
    circos.heatmap(mat2, col = col_fun1,
                   split = factor(data_annotation$group, levels = group_oi),
                   na.col = "grey", bg.border = "black")
    circos.track(track.index = get.current.track.index(),
                 panel.fun = function(x, y) {
        if (CELL_META$sector.numeric.index == lastSector) {
            # the last sector
            cn = rev(colnames(mat1))
            n = length(cn)
            circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                        1:n - 0.5, cn, cex = 0.5, adj = c(0, 0.5),
                        facing = "inside")
        }
    }, bg.border = NA)
    circos.clear()
    text(0, 1, titleName)
    lgd1 = Legend(title = "Exp", col_fun = col_fun1)
    grid.draw(lgd1)
    dev.off()
    message(date(), ": Check output directory for results")

    ## Vizualize
    circos.clear()
    circos.par(start.degree = 90, gap.after = c(rep(2, lastSector - 1), 40))
    # circos.par(start.degree = 90, gap.degree = 5)
    circos.heatmap(mat1, col = col_fun1,
                   split = factor(data_annotation$group, levels = group_oi),
                   rownames.side = "outside", na.col = "grey",
                   bg.border = "black")
    circos.heatmap(mat2, col = col_fun1,
                   split = factor(data_annotation$group, levels = group_oi),
                   na.col = "grey", bg.border = "black")
    circos.track(track.index = get.current.track.index(),
                 panel.fun = function(x, y) {
        if (CELL_META$sector.numeric.index == lastSector) {
            # the last sector
            cn = rev(colnames(mat1))
            n = length(cn)
            circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                        1:n - 0.5, cn, cex = 0.5, adj = c(0, 0.5),
                        facing = "inside")
        }
    }, bg.border = NA)
    circos.clear()
    text(0, 0, "")

    if(colorscale == TRUE) {
      lgd1 = Legend(title = "Exp", col_fun = col_fun1)
      grid.draw(lgd1)
    }

    return(list(mod1 = output_mat1, mod2 = output_mat2))
}
