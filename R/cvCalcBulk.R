#' cvCalcBulk Function
#'
#' This function allows to calculate Intra-donor variations in bulk data over
#' longitudinal timepoints. The coefficient of variation (CV=SD/mean) is
#' calculated in Bulk data in same donor/participant across timepoints.
#' @param data_object Input \emph{PALMO} S4 object. Contains annotation table
#' and expression matrix or data frame. Rows represent gene/proteins column
#' represents participant samples (same as annotation table Sample column)
#' @param meanThreshold Average expression threshold to filter lowly expressed
#' genes Default is 1 (log2 scale)
#' @param cvThreshold Coefficient of variation threshold to select variable and
#' stable genes. Default is 5 for bulk data. Users can use 10-20 for single cell
#' average expression data.
#' @param median_cvThreshold Optional, median of CVs from each donor/participant
#' calculated. Threshold used to differentiate variable and stable features
#' across donors/participants. Default, same as \emph{cvThreshold}.
#' @param donorThreshold Donor threshold number to be used, Default is number of
#' participants
#' @param naThreshold Optional, For a give feature % of donors/participants
#' showing non-NA CVs (NAs appear due to expression ~0 or absent). Default is 1
#' means all donors/participants to consider. 0.5 means from 4 donors atleast 2
#' donors should have non-NA CVs for a given feature.
#' @param housekeeping_genes Optional, vector of housekeeping genes. Default is
#' c("ACTB", "GAPDH")
#' @param plot_log10 Optional, Plot CV vs Mean on log10 scale. Default FALSE
#' @param selectedFeatures Optional, focus on selected genes/features.
#' @param median_cv_max Optional, Remove features with greater than median CV
#' Default is NULL
#' @param plotWidth Optional, heat plot width 5 in
#' @param plotHeight Optional, heat plot height 8 in
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @return PALMO object with CV list
#' @keywords cvCalcBulk
#' @export
#' @examples
#' \dontrun{
#' palmo_obj=cvCalcBulk(data_object=palmo_obj, meanThreshold=0.1, cvThreshold=5)
#' }

cvCalcBulk <- function(data_object, meanThreshold=1,
                       cvThreshold=5,
                       median_cvThreshold=NULL,
                       donorThreshold=NULL,
                       housekeeping_genes=NULL,
                       naThreshold=1,
                       plot_log10=FALSE,
                       selectedFeatures=NULL,
                       median_cv_max=NULL,
                       plotWidth=5, plotHeight=8,
                       fileName=NULL, filePATH=NULL) {

    message(date(),": Performing Coefficient of variance analysis")
    if(is.null(fileName)) {
        fileName <- "outputFile"
    }
    if(is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }

    ## meanThrehold
    if(is.null(meanThreshold)) {
        meanThreshold <- 0
        message(date(),": Using mean threshold >= 0")
    }
    data_object@meanThreshold <- meanThreshold
    data_object@cvThreshold <- cvThreshold

    ## Cumulative cvThrehold (across participants)
    if(is.null(median_cvThreshold)) {
        median_cvThreshold <- cvThreshold
        message(date(),": Using median CV threshold (across donors) same as CV
          threshold at single donor")
        data_object@median_cvThreshold <- median_cvThreshold
    }

    ## Assign housekeeping_genes
    if(!is.null(housekeeping_genes)) {
        data_object@housekeeping_genes <- housekeeping_genes
    }

    ## get the data
    ann <- data_object@curated$anndata
    mat <- data_object@curated$data
    check_data <- all.equal(row.names(ann), colnames(mat))
    if(check_data == FALSE) {
        stop(date(),": Annotation of samples (rows) and datamatrix columns do
             not match")
    }

    ## If features selected before hand
    if(!is.null(selectedFeatures)) {
        message(date(),": Filtering for selected features")
        mat <- mat[selectedFeatures,]
    }

    ## CV vs Mean
    unigene <- row.names(mat)
    uniSample <- sort(unique(ann$PTID))

    variable_gene <- NA
    stable_gene <- NA
    pdf(paste(filePATH,"/",fileName,"-CV-Sample-Plot.pdf", sep=""),
        width=5, height=5)
    for(i in 1:length(uniSample)) {
        uS <- uniSample[i]
        #print(uS)
        meta_df <- ann[ann$PTID %in% uS,]
        if(nrow(meta_df)>1) {
            df <- mat[unigene, meta_df$Sample]
            df <- data.frame(df, NAs=apply(df,1,function(x){sum(is.na(x))}),
                             Zeros=apply(df,1,function(x){sum(x==0)}),
                             mean=rowMeans(df, na.rm=TRUE),
                             sd=apply(df,1,sd, na.rm=TRUE),
                             var=apply(df,1,var, na.rm=TRUE),
                             stringsAsFactors = FALSE)
            df$CV <- 100*df$sd/df$mean
            dp1 <- df
            dp1 <- dp1[!is.na(dp1$mean),]
            dp1 <- dp1[abs(dp1$mean) >= meanThreshold,]

            ## Plot
            plot1 <- ggplot(dp1, aes(x=mean, y=CV)) +
                geom_point(size=0.5, color="grey") +
                labs(title=paste(uniSample[i], " (g=", nrow(dp1),"/",
                                 nrow(df),")", sep="")) +
                theme_classic()
            if(plot_log10 ==TRUE) {
                plot1 <- plot1 + scale_x_continuous(trans='log10') +
                    scale_y_continuous(trans='log10')
            }

            ## Variable genes
            dp2a <- dp1[abs(dp1$CV)> cvThreshold,]
            dp2a <- dp2a[!is.na(dp2a$CV),]
            #dp2a <- dp1 %>% filter(abs(CV)> cvThreshold)
            if(nrow(dp2a)>0) {
                dp2a <- dp2a[order(abs(dp2a$CV), abs(dp2a$mean),
                                   decreasing = TRUE),]
                if(nrow(dp2a)>10) {
                    dp2a_label <- dp2a[1:10,]
                } else {
                    dp2a_label <- dp2a
                }
                plot1 <- plot1 +geom_text_repel(data=dp2a_label,
                            aes(x=mean, y=CV, label=row.names(dp2a_label)),
                            col="red", size=2, max.overlaps=20)

                ## Add variable genes
                variable_gene <- rbind(variable_gene,
                                       data.frame(donor=uS,
                                        feature=row.names(dp2a),
                                        dp2a[,c("mean","sd","var", "CV", "NAs",
                                                "Zeros")],
                                        stringsAsFactors = FALSE))
            }

            ## Stable features
            dp2b <- dp1[abs(dp1$CV) <= cvThreshold,]
            #dp2b <- dp1 %>% filter(abs(CV) <= cvThreshold)
            dp2b <- dp2b[!is.na(dp2b$CV),]
            if(nrow(dp2b)>0) {
                dp2b <- dp2b[order(-abs(dp2b$mean), abs(dp2b$CV),
                                   decreasing = FALSE),]
                if(nrow(dp2b)>10) {
                    dp2b_label <- dp2b[1:10,]
                } else {
                    dp2b_label <- dp2b
                }
                plot1 <- plot1 + geom_text_repel(data=dp2b_label,
                                        aes(x=mean, y=CV,
                                            label=row.names(dp2b_label)),
                                        col="blue", size=2, max.overlaps=20)

                #Add stable genes
                stable_gene <- rbind(stable_gene,
                                     data.frame(donor=uS,
                                                feature=row.names(dp2b),
                                                dp2b[,c("mean","sd","var","CV",
                                                        "NAs", "Zeros")],
                                                stringsAsFactors = FALSE))

            }
            print(plot1)

            ## Heatmap for variance
            #df <- df[abs(df$mean) >= meanThreshold,]
            temp <- data.frame(feature=row.names(df), var=df$CV,
                               stringsAsFactors = FALSE)
            colnames(temp)[-1] <- paste(uS, "_", colnames(temp)[-1], sep="")
            if(i==1) {
                res <- temp
            } else {
                res <- merge(res, temp, by="feature", all=TRUE)
            }
        }

    }
    dev.off()

    ## Remove blank
    if(!is.null(nrow(variable_gene))) {
        variable_gene <- variable_gene[!is.na(variable_gene$donor),]
    }
    if(!is.null(nrow(stable_gene))) {
        stable_gene <- stable_gene[!is.na(stable_gene$donor),]
    }

    ## Decompose result into mean and CV (variance)
    if(ncol(res) == 2) {
        res_var <- data.frame(res[,-1], stringsAsFactors = FALSE)
        colnames(res_var) <- colnames(res)[2]
        row.names(res_var) <- res$feature
        colnames(res_var) <- gsub("_var", "", colnames(res_var))
        var_mat <- NA
        stable_mat <- NA

    } else if(ncol(res) > 2){
        res_var <- res[,-1]
        row.names(res_var) <- res$feature
        colnames(res_var) <- gsub("_var", "", colnames(res_var))

        #May be some samples are bad
        uniSample <- intersect(uniSample, colnames(res_var))

        #If donor cutoff NULL then use all donors
        if(is.null(donorThreshold)) {
          #Calculate feature Mean of CV by all donors
          res_var$Mean <- apply(res_var[,uniSample],1,function(x){
            mean(abs(x),na.rm=TRUE)
          })
          #Calculate feature Median of CV by all donors
          res_var$Median <- apply(res_var[,uniSample],1,function(x){
            median(abs(x),na.rm=TRUE)
          })
        } else {

          if(donorThreshold > length(uniSample)) {
            message(date(),": donorThreshold greater than number of donors")
          }
          #Calculate feature Mean of CV by donor threshold
          res_var$Mean <- apply(res_var[,uniSample],1,function(x){
            x <- sort(x)[1:donorThreshold]
            mean(abs(x),na.rm=TRUE)
          })
          #Calculate feature Median of CV by donor threshold
          res_var$Median <- apply(res_var[,uniSample],1,function(x){
            x <- sort(x)[1:donorThreshold]
            median(abs(x),na.rm=TRUE)
          })
        }

        #Calculate number of NAs in all donors
        res_var$NAs <- apply(res_var[,uniSample],1,function(x){
          sum(is.na(x),na.rm=TRUE)
        })

        ## Remove genes with CV cutoff (high noise)
        if(!is.null(median_cv_max)) {
            res_var <- res_var[res_var$Median < median_cv_max,]
        }

        ## Calculate CV range and check whether greater than median
        max_cv <- max(abs(res_var$Median), na.rm=TRUE)
        if(median_cvThreshold> max_cv) {
            median_cvThreshold <- max_cv
            message(date(),": Input median_cv_max higher than maximum CV
                    range.")
        }

        #Remove NAs
        unisample_naThreshold <- ceiling(length(uniSample)*naThreshold)
        res_var <- res_var[res_var$NAs <= unisample_naThreshold,]

        #histogram of CV
        suppressMessages(dx <- melt(res_var[,uniSample]), classes = "message")
        plot1 <- ggplot(dx, aes(x=value)) +
            geom_histogram(aes(y=..density..), colour="black", fill="skyblue",
                           bins=50)+
            labs(x="CV") +
            theme_classic()
        pdf(paste(filePATH,"/",fileName,"-CV-distribution.pdf", sep=""),
            width=5, height=5)
        print(plot1)
        dev.off()
        print(plot1)

        #color list
        if(min(res_var$Median,na.rm=TRUE)<0) {
            col_list <- colorRamp2(c(-max_cv, -cvThreshold-0.01, -cvThreshold,
                                     -cvThreshold/2, 0, cvThreshold/2,
                                     cvThreshold, cvThreshold+0.01, max_cv),
                                   c("red", "pink", "blue", "skyblue", "black",
                                     "skyblue", "blue", "pink", "red"))
        } else {
            col_list <- colorRamp2(c(0, cvThreshold/2, cvThreshold,
                                     cvThreshold+0.01, max_cv),
                                   c("black", "skyblue", "blue", "pink", "red"))
        }

        ## Plot variable genes CV
        mat <- res_var[order(-(res_var$NAs), abs(res_var$Median),
                             decreasing = TRUE),]
        mat <- mat[abs(mat$Median) > median_cvThreshold,]
        if(nrow(mat)>1) {
            var_mat <- mat
            mat <- mat[,uniSample]
            if(nrow(mat)>50) {
                mat <- mat[1:50,]
                write.csv(mat, file=paste(filePATH,"/",fileName,
                                          "-CV-Variable-Matrix.csv", sep=""))
            }
            #Avoid boxplot error
            column_mat <- mat[,uniSample]
            column_mat[is.na(column_mat)] <- 0
            check_0 <- unique(as.numeric(unlist(apply(column_mat,2,
                                                      function(x) {
                                                          unique(x,na.rm=TRUE)
                                                          }) )))
            if(length(check_0) != 1) {
                column_ha <- HeatmapAnnotation(CV = anno_boxplot(column_mat))
                ht1 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                        cluster_columns = FALSE,
                        na_col = "grey", col = col_list,
                        row_names_max_width=unit(10, "cm"),
                        column_names_gp = gpar(fontsize = 5),
                        row_names_gp = gpar(fontsize = 6),
                        top_annotation = column_ha,
                        column_title="Variable features",
                        heatmap_legend_param = list(title = "median CV(%)",
                                            heatmap_legend_side = "right") )
            } else {
                ht1 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                        cluster_columns = FALSE,
                        na_col = "grey", col = col_list,
                        row_names_max_width=unit(10, "cm"),
                        column_names_gp = gpar(fontsize = 5),
                        row_names_gp = gpar(fontsize = 6),
                        column_title="Stable features",
                        heatmap_legend_param = list(title = "median CV(%)",
                                            heatmap_legend_side = "right") )
            }
            print(ht1)
            pdf(paste(filePATH,"/",fileName,"-CV-Variable-Heatplot.pdf",
                      sep=""), width=plotWidth, height=plotHeight)
            draw(ht1)
            dev.off()
        } else {
            var_mat <- NA
            message(date(),": Variable features do not found. Check data
            for missing values or change CV parameters")
        }

        ## Plot stable genes CV
        mat <- res_var[order(res_var$NAs, abs(res_var$Median),
                             decreasing = FALSE),]
        mat <- mat[abs(mat$Median) <= median_cvThreshold,]
        if(nrow(mat)>1) {
            stable_mat <- mat
            mat <- mat[,uniSample]
            if(nrow(mat)>50) {
                mat <- mat[1:50,]
                write.csv(mat, file=paste(filePATH,"/",fileName,
                                          "-CV-Stable-Matrix.csv", sep=""))
            }
            #Avoid boxplot error
            column_mat <- mat[,uniSample]
            column_mat[is.na(column_mat)] <- 0
            check_0 <- unique(as.numeric(unlist(apply(column_mat,2,
                                                      function(x) {
                                                          unique(x,na.rm=TRUE)
                                                          }) )))
            if(length(check_0) != 1) {
                column_ha <- HeatmapAnnotation(CV = anno_boxplot(column_mat))
                ht2 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                        cluster_columns = FALSE,
                        na_col = "grey", col = col_list,
                        row_names_max_width=unit(10, "cm"),
                        column_names_gp = gpar(fontsize = 5),
                        row_names_gp = gpar(fontsize = 6),
                        top_annotation = column_ha,
                        column_title="Stable features",
                        heatmap_legend_param = list(title = "median CV(%)",
                                            heatmap_legend_side = "right") )
            } else {
                ht2 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                        cluster_columns = FALSE,
                        na_col = "grey", col = col_list,
                        row_names_max_width=unit(10, "cm"),
                        column_names_gp = gpar(fontsize = 5),
                        row_names_gp = gpar(fontsize = 6),
                        column_title="Stable features",
                        heatmap_legend_param = list(title = "median CV(%)",
                                            heatmap_legend_side = "right") )
            }

            print(ht2)
            pdf(paste(filePATH,"/",fileName,"-CV-Stable-Heatplot.pdf", sep=""),
                width=plotWidth, height=plotHeight)
            draw(ht2)
            dev.off()
        } else {
            stable_mat <- NA
            message(date(),": Stable features do not found. Check data for
            missing values or change CV parameters. Summary of CV range
                    in given data:")
            message(summary(res_var$Median))
        }

        ## Housekeeping genes
        if(!is.null(housekeeping_genes)) {
            housekeeping_genes <- intersect(housekeeping_genes,
                                            row.names(res_var))
            mat <- res_var[housekeeping_genes,uniSample]
            #Avoid boxplot error
            column_mat <- mat[,uniSample]
            column_mat[is.na(column_mat)] <- 0
            check_0 <- unique(as.numeric(unlist(apply(column_mat,2,
                                 function(x) {unique(x,na.rm=TRUE) }) )))
            if(length(check_0) != 1) {
                column_ha <- HeatmapAnnotation(CV = anno_boxplot(column_mat))
                ht3 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                        cluster_columns = FALSE,
                        na_col = "grey", col = col_list,
                        row_names_max_width=unit(10, "cm"),
                        column_names_gp = gpar(fontsize = 5),
                        row_names_gp = gpar(fontsize = 6),
                        top_annotation = column_ha,
                        column_title="Housekepping features",
                        heatmap_legend_param = list(title = "median CV(%)",
                                heatmap_legend_side = "right"))
            } else {
                ht3 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                        cluster_columns = FALSE,
                        na_col = "grey", col = col_list,
                        row_names_max_width=unit(10, "cm"),
                        column_names_gp = gpar(fontsize = 5),
                        row_names_gp = gpar(fontsize = 6),
                        column_title="Housekepping features",
                        heatmap_legend_param = list(title = "median CV(%)",
                                            heatmap_legend_side = "right"))
            }

            print(ht3)
            pdf(paste(filePATH,"/",fileName,
                      "-CV-housekeeping_genes-Heatplot.pdf", sep=""),
                width=plotWidth, height=plotHeight)
            draw(ht3)
            dev.off()
        }

    }

    #All CV
    write.csv(res_var, file=paste(filePATH,"/",fileName,"-CV-result.csv",
                                  sep=""))
    rm(res, temp)
    #Variable genes
    write.csv(variable_gene, file=paste(filePATH,"/",
                    fileName,"-CV-VariableFeatures-result.csv", sep=""))
    #Stable genes
    write.csv(stable_gene, file=paste(filePATH,"/",
                    fileName,"-CV-StableFeatures-result.csv", sep=""))

    ## Add CV result
    data_object@result$cv_all <- res_var
    data_object@result$variable_gene <- variable_gene
    data_object@result$non_variable_gene <- stable_gene
    data_object@result$var_mat <- var_mat
    data_object@result$stable_mat <- stable_mat

    message(date(),": Done. Please check output directory for Plots/results.")
    return(data_object)
}
