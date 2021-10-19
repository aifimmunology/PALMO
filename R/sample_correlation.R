#' A sample_correlation Function
#'
#' This function allows to perform sample correlation (by group like
#' celltype, ot by donor).
#' @param data Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (if celltype with in donor then sample:
#' celltype, separated by :)
#' @param column_sep Sample and celltype seperator like (:)
#' @param method Correlation method "pearson" or "spearman"
#' @param groupColumn Data column names consists group (Donor-group) at 2nd place
#' or 1st place(like PTIDxGroupX, 2 or GroupXPTIDx, 1)
#' @param clusterBy Cluster correlation result by "donor" or "group". Default
#' donor
#' @param max Maximum color limit (Default, 0.9 correlation)
#' @param column_names_fontsize Font size of the column names, default 4
#' @param row_names_fontsize Font size of the row names, default 4
#' @param column_title_fontsize Font size of the column title, default 6
#' @param row_title_fontsize Font size of the row title, default 6
#' @param plotHeight Height of the plot (in), deafult 20in
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords sample_correlation
#' @export
#' @examples
#' #res <- sample_correlation(data=datamatrix, column_sep=":", method="spearman")

sample_correlation <- function(data, column_sep=":",
                       method="spearman", groupColumn=2,
                       clusterBy="donor", max=0.9,
                       column_names_fontsize=4, row_names_fontsize=4,
                       row_title_fontsize=6, column_title_fontsize=6,
                       plotHeight=20, fileName=NULL, filePATH=NULL) {

  cat(date(),": Performing Sample correlation anlaysis\n")
  if(is.null(fileName)) {
    fileName <- "outputFile"
  }
  if(is.null(filePATH)) {
    filePATH <- paste(getwd(), "/output", sep="")
    dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  }

  #Create annotation
  data_annotation <- data.frame(id=colnames(data))
  temp <- data.frame(do.call(rbind, strsplit(data_annotation$id,
             split = column_sep)), stringsAsFactors = F)

  if(ncol(temp)>1) {
    data_annotation <- data.frame(data_annotation, Sample=temp$X1, group=temp$X2, stringsAsFactors = F)
    row.names(data_annotation) <- data_annotation$id
    if(groupColumn==1) {
      colnames(data_annotation) <- c("sample","group", "donor")
    } else {
      colnames(data_annotation) <- c("sample","donor","group")
    }

    #Sample Correlation map
    cor_mat <- rcorr(as.matrix(data), type=method)
    res <- cor_mat$r

    #Plot heatmap
    if(clusterBy == "donor") {
      ha_col <- HeatmapAnnotation(df=data.frame(donor=data_annotation$donor))
      ht1 <- Heatmap(data.matrix(res),
                     cluster_rows =F, cluster_columns = F,
                     row_split = as.character(data_annotation$donor),
                     column_split = as.character(data_annotation$donor),
                     na_col = "grey",
                     col = colorRamp2(c(-1,-max,0,max,1),
                                      c("blue","skyblue","white","pink","red")),
                     row_names_max_width=unit(10, "cm"),
                     column_names_gp = gpar(fontsize = column_names_fontsize), row_names_gp = gpar(fontsize =row_names_fontsize),
                     top_annotation = ha_col,
                     row_title_gp = gpar(fontsize = row_title_fontsize),
                     column_title_gp = gpar(fontsize = column_title_fontsize),
                     heatmap_legend_param = list(title = method,
                                                 heatmap_legend_side = "right")
      )
    } else {
      ha_col <- HeatmapAnnotation(df=data.frame(Group=data_annotation$group))
      ht1 <- Heatmap(data.matrix(res),
                   cluster_rows =F, cluster_columns = F,
                   row_split = as.character(data_annotation$group),
                   column_split = as.character(data_annotation$group),
                   na_col = "grey",
                   col = colorRamp2(c(-1,-max,0,max,1),
                          c("blue","skyblue","white","pink","red")),
                   row_names_max_width=unit(10, "cm"),
                   column_names_gp = gpar(fontsize = column_names_fontsize), row_names_gp = gpar(fontsize =row_names_fontsize),
                   top_annotation = ha_col,
                   row_title_gp = gpar(fontsize = row_title_fontsize),
                   column_title_gp = gpar(fontsize = column_title_fontsize),
                   heatmap_legend_param = list(title = method,
                                        heatmap_legend_side = "right")
                   )
    }
  } else {

    row.names(data_annotation) <- data_annotation$id
    #Sample Correlation map
    cor_mat <- rcorr(as.matrix(data), type=method)
    res <- cor_mat$r
    #Plot heatmap
    ht1 <- Heatmap(data.matrix(res),
                   cluster_rows =F, cluster_columns = F,
                   na_col = "grey", col = colorRamp2(c(-1,-0.9,0,0.9,1),
                                                     c("blue","skyblue","white","pink","red")),
                   row_names_max_width=unit(10, "cm"),
                   column_names_gp = gpar(fontsize = column_names_fontsize),
                   row_names_gp = gpar(fontsize =row_names_fontsize),
                   heatmap_legend_param = list(title = method,
                         heatmap_legend_side = "right")
                   )

  }

  pdf(paste(filePATH,"/",fileName,"-SampleCorrelationplot.pdf", sep=""),
      width=plotHeight, height=plotHeight)
  draw(ht1)
  dev.off()
  print(draw(ht1))

  cat(date(),": Please check output directory for results\n")
  return(res)
}
