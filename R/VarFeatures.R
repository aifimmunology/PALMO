#' VarFeatures Function
#'
#' This function allows user to identify variable genes in participants across
#' longitudinal timepoints using single cell expression data. The coefficient of
#' variation (CV) calculated using \code{cvCalcSC} function. Users can identify
#' cvThreshold in different datasets using housekeeping genes CV distribution.
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information and expression data from Bulk or single cell data.
#' @param group_oi Group of interest to focus on. Example among celltypes focus
#' on selected ones. Default is NULL.
#' @param cvThreshold Coefficient of variation threshold to select variable and
#' stable genes Default is 10 for single cell RNA \code{(100*SD/mean)}
#' @param donorThreshold Donor threshold number to be used, Default is number of
#' participants
#' @param groupThreshold Group label threshold number to be used, Default is
#' \code{(number of participants x group labels)/2}
#' @param topFeatures Number of features to be selected from each group, Default
#' is 25
#' @param fileName User defined filename
#' @param filePATH User-defined output directory \emph{PATH} to load \emph{CV}
#' result obtained from \code{cvCalcSC} function
#' @return PALMO object with variable (var_genes) features
#' @keywords VarFeatures
#' @export
#' @examples
#' \dontrun{
#' palmo_obj <- VarFeatures(data_object=palmo_obj, cvThreshold=10)
#' }

VarFeatures <- function(data_object, group_oi = NULL, cvThreshold = NULL,
                        donorThreshold = NULL, groupThreshold = NULL,
                        topFeatures = 25, filePATH = NULL, fileName = NULL) {

    message(date(), ": Identifying Variable features")

    ## If filename or filepath null
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    if (is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }

    ## Load CV result
    if (!is.null(data_object@result$cv_meanthreshold)) {
        cv_res <- data_object@result$cv_meanthreshold
    } else {
        stop(date(), ": Please run cvCalcSC function before VarFeatures
             function")
    }

    variable_gene <- data_object@result$variable_gene
    non_variable_gene <- data_object@result$non_variable_gene

    ## get the annotation data
    data_object@curated$anndata$Sample_group_i <- paste(data_object@curated$anndata$group,
                                data_object@curated$anndata$PTID, sep = ":")
    ann <- data_object@curated$anndata

    ## group list
    if (is.null(group_oi)) {
        group_oi <- unique(as.character(ann$group))
    }

    ## Select group of interest
    ann_sub <- ann[ann$group %in% group_oi, ]
    if (nrow(ann_sub) < 1) {
        stop(date(), ": Group of interest features do not match with annotation
             group eg.", unique(ann$group)[1:3])
    }
    Sample_group <- unique(ann_sub$Sample_group_i)
    Sample_group <- intersect(Sample_group, colnames(cv_res))

    if (is.null(donorThreshold)) {
        donorThreshold <- length(unique(ann_sub$PTID))
        message(date(), ": Donor threshold defined ", donorThreshold)
    } else if (donorThreshold > length(unique(ann$PTID))) {
        donorThreshold <- length(unique(ann_sub$PTID))
        message(date(), ": Donors were larger than unique donors. Donor
        threshold defined ", donorThreshold)
    }

    ## Calculate group threshold limit as 90% of num of donor * num of groups
    gThr <- round(length(unique(ann_sub$PTID)) *
                      length(unique(ann_sub$group)) * 0.9)
    if (is.null(groupThreshold)) {
        groupThreshold <- round(length(unique(ann_sub$PTID)) *
                                    length(unique(ann_sub$group)) * 0.5)
        message(date(), ": Groupwise threshold defined ", groupThreshold)
    } else if (groupThreshold > gThr) {
        groupThreshold <- round(length(unique(ann_sub$PTID)) *
                                    length(unique(ann_sub$group)) * 0.9)
        message(date(), ": Number of groups were larger than unique
        donors x groups. Groupwise threshold defined ", groupThreshold)
    }

    ## Summary of variable genes
    temp <- data.frame(do.call(rbind,strsplit(as.character(variable_gene$donor),
                            split = ":")), stringsAsFactors=FALSE)
    variable_gene$PTID <- temp$X2
    variable_gene$group <- temp$X1
    var_genelist <- variable_gene[variable_gene$group %in% group_oi, ]
    var_genelist <- data.frame(table(var_genelist$gene), stringsAsFactors=FALSE)
    var_genelist <- var_genelist[order(var_genelist$Freq, decreasing=TRUE), ]
    var_geneName <- as.character(var_genelist$Var1)
    # create Variable matrix
    var_geneName <- intersect(var_geneName, row.names(cv_res))
    var_mat <- cv_res[var_geneName, ]
    var_mat[var_mat <= cvThreshold] <- NA
    save(var_mat, file = paste(filePATH, "/", fileName, "-variableMatrix.Rda",
                               sep = ""))

    ## Define the super-variable genes
    plot1 <- ggplot(var_genelist, aes(x = Freq)) +
        geom_histogram(binwidth = 1) +
        labs(title = "Variable genes occurance from each sample")
    ## Atleast in all donor x group
    super_variable1 <- var_genelist[var_genelist$Freq >= groupThreshold, ]
    super_variable2 <- var_genelist[var_genelist$Freq >= donorThreshold &
                                var_genelist$Freq < groupThreshold, ]
    ## Top 25
    gn <- as.character(super_variable1$Var1[1:25])
    data_mat <- var_mat[gn, Sample_group]
    rn <- data.frame(do.call(rbind, strsplit(colnames(data_mat), split = ":")),
                     stringsAsFactors = FALSE)
    # set.seed(2020)
    ha_col <- HeatmapAnnotation(df = data.frame(group = rn$X1),
                                annotation_name_gp = gpar(fontsize = 6),
                                simple_anno_size = unit(0.3, "cm"))
    ht1 <- Heatmap(data.matrix(data_mat),
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   column_split = factor(rn$X1, levels = group_oi),
                   na_col = "grey",
                   col = colorRamp2(c(0, cvThreshold, (cvThreshold + 0.001),
                    (cvThreshold + 50)), c("blue", "white", "pink", "red")),
                   row_names_max_width = unit(10, "cm"),
                   column_names_gp = gpar(fontsize = 5),
                   row_names_gp = gpar(fontsize = 6),
                   row_title = "Super variable 25",
                   column_title_gp = gpar(fontsize = 4),
                   top_annotation = ha_col,
                   heatmap_legend_param = list(title = "CV",
                                    heatmap_legend_side = "right"))
    pdf(paste(filePATH, "/", fileName, "-Super-variable-Top25.pdf", sep = ""),
        width = 10, height = 3.5)
    print(plot1)
    print(ht1)
    dev.off()
    print(ht1)
    write.csv(data_mat, paste(filePATH, "/", fileName,
                "-Super-variable-Top25.csv", sep = ""))

    ## Define the variable genes
    var_list <- as.character(super_variable2$Var1)
    dfx <- c()
    for (i in 1:length(group_oi)) {
        groupName <- group_oi[i]
        df <- variable_gene[variable_gene$gene %in% var_list &
                                variable_gene$group %in% groupName,]
        df1 <- data.frame(table(df$gene, df$group))
        df1 <- df1[df1$Freq >= donorThreshold &
                       order(df1$Freq, decreasing = TRUE), ]
        df <- df[df$gene %in% df1$Var1, ]
        dfx <- rbind(dfx, df[1:topFeatures, ])
    }
    dfx <- dfx[!is.na(dfx$mean), ]
    var_gene <- unique(dfx$gene)
    write.csv(dfx, file = paste(filePATH, "/", fileName,
                                "-variable-genelist.csv", sep = ""))
    ## Plot heatmap
    data_mat <- var_mat[var_gene, Sample_group]
    ht2 <- Heatmap(data.matrix(data_mat),
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   column_split = factor(rn$X1, levels = group_oi),
                   na_col = "grey",
                   col = colorRamp2(c(0, cvThreshold, (cvThreshold + 0.001),
                        (cvThreshold + 50)), c("white", "blue", "pink", "red")),
                   row_names_max_width = unit(10, "cm"),
                   column_names_gp = gpar(fontsize = 5),
                   row_names_gp = gpar(fontsize = 5),
                   row_title = paste("Variable genes:", length(var_gene)),
                   column_title_gp = gpar(fontsize = 4),
                   top_annotation = ha_col,
                   heatmap_legend_param = list(title = "CV",
                                heatmap_legend_side = "right"))

    pdf(paste(filePATH, "/", fileName, "-variable-Features.pdf", sep = ""),
        width = 10, height = 10)
    print(ht2)
    dev.off()
    print(ht2)
    write.csv(data_mat, paste(filePATH, "/", fileName, "-variable-Features-",
                              topFeatures, ".csv", sep = ""))

    ## Add the result
    data_object@result$var_genes <- dfx
    message(date(), ": Check output directory for results")

    return(data_object)
}
