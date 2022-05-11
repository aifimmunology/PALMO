#' outlierDetect Function
#'
#' This function allows users to perform outlier analysis on bulk data by
#' calculating z-score. Outlier genes defined as mean/SD = |Z| > z_cutoff.
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information and expression data from Bulk or single cell data.
#' @param z_cutoff |Z| cutoff threshold to find potential outliers (Eg.
#' \emph{z_cutoff}=2, equals to \code{Mean/SD} 2)
#' @param plotWidth User-defined plot width, Default 10 in
#' @param plotHeight User-defined plot height, Default 5 in
#' @param group_column Include group by outlier analysis (celltype, cluster)
#' @param cl Number of clusters. Use nCores-1 to run parallel. Default 2
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory \emph{PATH} Default, current
#' directory
#' @return PALMO object with outlier_res dataframe
#' @keywords outlierDetect
#' @export
#' @examples
#' \dontrun{
#' palmo_obj <- outlierDetect(data_object=palmo_obj, z_cutoff=2)
#' }

outlierDetect <- function(data_object, z_cutoff = NULL, plotWidth = 10,
                          plotHeight = 5, group_column = NULL, cl = 2,
                          fileName = NULL, filePATH = NULL) {

    message(date(), ": Performing Outlier anlaysis\n")
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    if (is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }
    ## Z-cutoff
    if (is.null(z_cutoff)) {
        z_cutoff <- 2
    }
    data_object@z_cutoff <- z_cutoff

    ## get the data
    ann <- data_object@curated$anndata
    mat <- data_object@curated$data
    check_data <- all.equal(row.names(ann), colnames(mat))
    if (check_data == FALSE) {
        stop(date(), ": Annotation of samples (rows) and datamatrix columns do not match")
    }

    ## Input
    rowN <- row.names(mat)
    uniTime <- as.character(unique(ann$Time))
    uniSample <- sort(unique(ann$PTID))

    if (is.null(group_column)) {

        ## Calculate Z-score
        op <- pboptions(type = "timer")  # default
        outlier_res <- pblapply(rowN, cl = cl, function(geneName) {
            df <- data.frame(exp = as.numeric(mat[geneName, ]), ann,
                             stringsAsFactors = FALSE)
            dfx <- lapply(uniSample, function(y) {
                temp <- df[df$PTID %in% y, ]
                if (nrow(temp) > 0) {
                  temp$gene <- geneName
                  temp$meanDev <- temp$exp - mean(temp$exp, na.rm = TRUE)
                  temp$z <- (temp$exp - mean(temp$exp, na.rm = TRUE))/(sd(temp$exp, na.rm = TRUE))
                  temp$outlier <- ifelse(abs(temp$z) >= z_cutoff, temp$z, 0)
                  temp <- temp[temp$outlier != 0, ]
                  return(data.frame(temp))
                }
            })
            dfx <- do.call(rbind, dfx)
            return(dfx)
        })
        pboptions(op)
    } else {

        ann$bygroup <- ann[, group_column]
        ann$Sample_bygroup <- paste(ann$Sample, ann$bygroup, sep = ":")
        ## Input group
        uniGroup <- sort(unique(ann$bygroup))

        ann <- ann[order(ann$PTID, ann$bygroup), ]
        uniSample_group <- unique(ann$Sample_bygroup)
        mat <- mat[, row.names(ann)]
        # all.equal(row.names(ann), colnames(mat))

        ## Calculate Z-score (outlier analysis)
        op <- pboptions(type = "timer")  # default
        outlier_res <- pblapply(rowN, cl = cl, function(geneName) {
            df <- data.frame(exp = as.numeric(mat[geneName, ]), ann, stringsAsFactors = FALSE)
            dfx <- lapply(uniSample, function(uS) {
                temp <- df[df$PTID %in% uS, ]
                res <- lapply(uniGroup, function(uG) {
                  temp <- temp[temp$bygroup %in% uG, ]
                  if (nrow(temp) > 0) {
                    temp$gene <- geneName
                    temp$meanDev <- temp$exp - mean(temp$exp, na.rm = TRUE)
                    temp$z <- (temp$exp - mean(temp$exp, na.rm = TRUE))/(sd(temp$exp, na.rm = TRUE))
                    temp$outlier <- ifelse(abs(temp$z) >= z_cutoff, temp$z, 0)
                    temp <- temp[temp$outlier != 0, ]
                    return(temp)
                  }
                })
                res <- do.call(rbind, res)
                return(res)
            })
            dfx <- do.call(rbind, dfx)
            return(dfx)
        })
        pboptions(op)
    }

    ## Combine data
    outlier_res <- do.call(rbind, outlier_res)
    outlier_res <- outlier_res[!is.na(outlier_res$Sample), ]
    outlier_res <- outlier_res[, !colnames(outlier_res) %in% "outlier"]
    write.csv(outlier_res, file = paste(filePATH, "/", fileName, "-Outlier-result.csv", sep = ""),
        row.names = FALSE)

    ## Plot
    df <- outlier_res
    if (nrow(df) > 1) {
        df$Sample <- factor(df$Sample, levels = unique(ann$Sample))
        df$direction <- ifelse(df$z > 0, "> Z", "< -Z")
        df$direction <- factor(df$direction, levels = c("> Z", "< -Z"))
        # Z-plot
        plot1 <- ggplot(df, aes(x = Sample, y = z, fill = direction)) +
            geom_violin(scale = "width") +
            geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white",
                         position = position_dodge(preserve = "single")) +
            ggforce::geom_sina(size = 0.1) +
            labs(x = "", y = "Z-score", title = paste("Outlier events |Z| >", z_cutoff), fill = "") +
            facet_wrap(~direction, scales = "free_y", ncol = 1) +
            scale_fill_manual(values = c(`> Z` = "red", `< -Z` = "blue")) +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6),
                  axis.text.y = element_text(size = 6),
                  legend.position = "right")

        ## Mean-deviation plot
        plot2 <- ggplot(df, aes(x = Sample, y = meanDev, fill = direction)) +
            geom_violin(scale = "width") +
            geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", position = position_dodge(preserve = "single")) +
            ggforce::geom_sina(size = 0.1) +
            labs(x = "", y = "Mean deviation",
                 title = paste("Outlier events |Z| >", z_cutoff), fill = "") +
            facet_wrap(~direction, scales = "free_y", ncol = 1) +
            scale_fill_manual(values = c(`> Z` = "red", `< -Z` = "blue")) +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6),
                  axis.text.y = element_text(size = 6),
                  legend.position = "right")

        ## Count plot
        df_up <- df[df$z > 0, ]
        df_up <- data.frame(table(df_up$PTID, df_up$Time))
        df_down <- df[df$z < 0, ]
        df_down <- data.frame(table(df_down$PTID, df_down$Time))
        df1 <- rbind(data.frame(df_up, direction = "> Z"), data.frame(df_down, direction = "< -Z"))
        df1$id <- paste(df1$Var1, df1$Var2, sep = "")
        df1$direction <- factor(df1$direction, levels = c("> Z", "< -Z"))
        df1$label <- ifelse(abs(df1$Freq) > 0, df1$Freq, NA)
        plot3 <- ggplot(df1, aes(x = id, y = Freq, fill = direction)) +
            geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
            labs(x = "", y = "# Features", title = paste("Outlier events |Z| >", z_cutoff), fill = "") +
            scale_fill_manual(values = c(`> Z` = "red", `< -Z` = "blue")) +
            geom_text(aes(label = label), position = position_dodge(width = 0.9), size = 2, vjust = 0) +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6),
                  axis.text.y = element_text(size = 6),
                  legend.position = "right")

        message(date(), ": Results are in ...Outlier-result.csv sheet. Now generating plots.\n")
        if (is.null(group_column)) {
            ## Plot
            pdf(paste(filePATH, "/", fileName, "-Outlier-Boxplot.pdf", sep = ""), width = plotWidth,
                height = plotHeight)
            print(plot1)
            print(plot2)
            print(plot3)
            dev.off()
        } else {
            plot1b <- ggplot(df, aes(x = Sample, y = z, color = bygroup)) +
                geom_violin(scale = "width") +
                # geom_boxplot(width=0.1, fill='white') +
                labs(x = "", y = "Z-score", title = "Mean Deviation") +
                ggforce::geom_sina(size = 0.5) +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6),
                      axis.text.y = element_text(size = 6),
                      legend.position = "right")

            plot2b <- ggplot(df, aes(x = Sample, y = meanDev, color = bygroup)) +
                geom_violin(scale = "width") +
                # geom_boxplot(width=0.1, fill='white') +
                labs(x = "", y = "Mean Deviation", title = "Mean Deviation") +
                ggforce::geom_sina(size = 0.5) +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6),
                      axis.text.y = element_text(size = 6),
                      legend.position = "right")

            ## Plot
            pdf(paste(filePATH, "/", fileName, "-Outlier-Boxplot.pdf", sep = ""), width = plotWidth,
                height = plotHeight)
            print(plot1)
            print(plot2)
            print(plot1b)
            print(plot2b)
            dev.off()
        }

        ## Visualize
        print(plot2)
        print(plot1)
        print(plot3)

        ## Print result
        df <- outlier_res[order(outlier_res$z, decreasing = TRUE), ]
        data_object@result$outlier_res <- df

        message(date(), ": Please check output directory for results\n")
    } else {
        message(date(), ": Did not see events with given z cutoff\n")
    }
    return(data_object)

}
