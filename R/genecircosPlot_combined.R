#' A genecircosPlot_combined Function
#'
#' This function allows you to Circos Plot for gene list of interest by group
#' @param data Expression matrix or data frame. Rows represents gene/proteins
#' column represents group:donor (group and donor separated by :)
#' @param geneList Genes of interest to explore
#' @param groupColumn Default 1, use 2 when columns are donor:group format
#' @param group_oi Optional, User-defined groups to consider
#' @param colorThreshold User-defined color threshold in colorspace
#' @keywords genecircosPlot
#' @export
#' @examples
#' ##Circos Plot for genes expression in a group
#' #geneList <- c("IL32","CCL5","TCF7","IL7R","LEF1")
#' #res <- genecircosPlot(data=cv_res, geneList=geneList)

genecircosPlot_combined <- function(data, geneList, groupColumn=NULL, group_oi=NULL, colorThreshold=NULL) {
    
    if(is.null(colorThreshold)) {
      colorThreshold <- 10
    }
    if(is.null(groupColumn)) {
      groupColumn <- 1
    }
    
    #Split the group and donor
    data_annotation <- data.frame(do.call(rbind, strsplit(colnames(data), split = ":")), stringsAsFactors = F)
    data_annotation$id <- colnames(data)
    if(groupColumn==1) {
      colnames(data_annotation) <- c("group", "donor", "sample")
    } else {
      colnames(data_annotation) <- c("donor","group", "sample")
    }
    data_annotation <- data_annotation[data_annotation$group %in% group_oi,]
    data_annotation <- data_annotation[order(data_annotation$sample),]
    
    if(is.null(group_oi)) {
      group_oi <- sort(unique(data_annotation$group))
      Sample_group <- unique(data_annotation$sample)
    } else {
      Sample_group <- data_annotation$sample
    }
    
    data_mat <- data[geneList,data_annotation$sample]
    ht1 <- Heatmap(data.matrix(data_mat), cluster_rows =F,  cluster_columns = F, 
                   column_split = factor(data_annotation$group, levels = group_oi),
                   na_col = "grey", col = colorRamp2(c(0,colorThreshold,(colorThreshold+0.1),(colorThreshold+50)), c("blue","white","pink","red")),
                   row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                   column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 6),
                   column_title = "Celltype",
                   heatmap_legend_param = list(title = "CV",heatmap_legend_side = "right") )
    #ht1
    
    #Define order
    column_ht <- unlist(column_order(ht1))
    data_mat <- data_mat[,column_ht]
    data_annotation <- data_annotation[column_ht,]
    
    data_mat <- t(data_mat)
    col_fun1 = colorRamp2(c(0,colorThreshold,(colorThreshold+0.1),(colorThreshold+50)), c("blue","white","pink","red"))
    lastSector <- length(unique(data_annotation$group))
    
    circos.clear()
    circos.par(start.degree = 90, gap.after = c(rep(2, lastSector-1), 40))
    #circos.par(start.degree = 90, gap.degree = 5)
    circos.heatmap(data_mat, col = col_fun1, split=data_annotation$group, rownames.side = "outside", na.col = "grey", bg.border="black")
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
      if(CELL_META$sector.numeric.index == lastSector) { # the last sector
        cn = rev(colnames(data_mat))
        n = length(cn)
        circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                    1:n - 0.5, cn, 
                    cex = 0.5, adj = c(0, 0.5), facing = "inside")
      }
    }, bg.border = NA)
    circos.clear()
    
    return(data_mat)
    
}
