#' A cvCalcBulk Function
#'
#' This function allows to calculate Intra-donor variations over longitudinal
#' timepoints. The coefficient of variation (CV) is calculated in Bulk
#' data without group information. CV calculated across samples. It requires
#' longitudinal data matrix/data frame and annotation file.
#' @param ann Annotation table. Table must consist column Sample (Participant
#' sample name), PTID (Participant), Time (longitudinal time points)
#' @param mat Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column)
#' @param meanThreshold Average expression threshold to filter lowly expressed
#' genes Default is 0.1 (log2 scale)
#' @param cvThreshold Coefficient of variation threshold to select variable and
#' stable genes Default is 5 for bulk olink data
#' @param housekeeping_genes Optional list of housekeeping genes to focus on.
#' Default is ACTB, GAPDH
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords cvCalcBulk
#' @export

cvCalcBulk <- function(ann, mat, meanThreshold=NULL, cvThreshold,
                       housekeeping_genes=NULL, fileName=NULL, filePATH=NULL) {

    cat(date(),": Performing Coefficient of variance analysis\n")
    if(is.null(fileName)) {
      fileName <- "outputFile"
    }
    if(is.null(filePATH)) {
      filePATH <- paste(getwd(), "/output", sep="")
      dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
    }
    if(is.null(meanThreshold)) {
      cat(date(),": Mean threshold not defined\n")
      meanThreshold <- 0 #absolute mean value
    }

    #CV vs Mean
    unigene <- row.names(mat)
    uniSample <- sort(unique(ann$PTID))

    variable_gene <- NA
    stable_gene <- NA
    pdf(paste(filePATH,"/",fileName,"-CV-Sample-Plot.pdf", sep=""), width=5, height=5)
    for(i in 1:length(uniSample)) {
      uS <- uniSample[i]
      #print(uS)
      meta_df <- ann[ann$PTID %in% uS,]
      if(nrow(meta_df)>1) {
        df <- mat[unigene, meta_df$Sample]
        df <- data.frame(df, NAs=apply(df,1,function(x){sum(is.na(x))}),
                         mean=rowMeans(df, na.rm=T), sd=apply(df,1,sd, na.rm=T),
                         var=apply(df,1,var, na.rm=T), stringsAsFactors = F)
        df$CV <- 100*df$sd/df$mean
        #summary(df$mean)
        #summary(df$var)
        dp1 <- df
        dp1 <- dp1[dp1$NAs <= (0.5*nrow(meta_df)),] #missing value threshold 50%NAs
        dp1 <- dp1[abs(dp1$mean) >= meanThreshold,]
        #Plot
        plot1 <- ggplot(dp1, aes(x=mean, y=CV)) +
          geom_point(size=0.5, color="grey") +
          labs(title=paste(uniSample[i], " (g=", nrow(dp1),"/",nrow(df),")", sep="")) +
          theme_classic()

        #Variable genes
        dp2a <- dp1[abs(dp1$CV)> cvThreshold,]
        if(nrow(dp2a)>0) {
          dp2a <- dp2a[order(abs(dp2a$CV), abs(dp2a$mean), decreasing = T),]
          if(nrow(dp2a)>10) {
            dp2a_label <- dp2a[1:10,]
          } else {
            dp2a_label <- dp2a
          }
          plot1 <- plot1 + geom_text_repel(data=dp2a_label,
                           aes(x=mean, y=CV, label=row.names(dp2a_label)),
                           col="red", size=2, max.overlaps=20)

          #Add variable genes
          variable_gene <- rbind(variable_gene,
                                 data.frame(donor=uS, gene=row.names(dp2a),
                                            dp2a[,c("mean","sd","var","CV")],
                                            stringsAsFactors = F))
        }

        dp2b <- dp1[abs(dp1$CV) <= cvThreshold,]
        if(nrow(dp2b)>0) {
          dp2b <- dp2b[order(-abs(dp2b$mean), abs(dp2b$CV), decreasing = F),]
          if(nrow(dp2b)>10) {
            dp2b_label <- dp2b[1:10,]
          } else {
            dp2b_label <- dp2b
          }
          plot1 <- plot1 + geom_text_repel(data=dp2b_label,
                                           aes(x=mean, y=CV,
                                               label=row.names(dp2b_label)),
                                           col="blue", size=2,
                                           max.overlaps=20)

          #Add stable genes
          stable_gene <- rbind(stable_gene,
                               data.frame(donor=uS, gene=row.names(dp2b),
                                          dp2b[,c("mean","sd","var","CV")],
                                          stringsAsFactors = F))

        }
        print(plot1)

        #Heatmap for variance
        temp <- data.frame(gene=row.names(df), var=df$CV, stringsAsFactors = F)
        colnames(temp)[-1] <- paste(uS, "_", colnames(temp)[-1], sep="")
        if(i==1) {
          res <- temp
        } else {
          res <- merge(res, temp, by="gene", all=T)
        }
      }

    }
    dev.off()

    #Remove blank
    variable_gene <- variable_gene[!is.na(variable_gene$donor),]
    stable_gene <- stable_gene[!is.na(stable_gene$donor),]

    #decompose result into mean and CV (variance)
    res_var <- res[,-1]
    row.names(res_var) <- res$gene
    colnames(res_var) <- gsub("_var", "", colnames(res_var))

    #May be some samples are bad
    uniSample <- intersect(uniSample, colnames(res_var))

    res_var$Mean <- apply(res_var[,uniSample],1,function(x){
      mean(abs(x),na.rm=T)
    })
    res_var$Median <- apply(res_var[,uniSample],1,function(x){
      median(abs(x),na.rm=T)
    })
    write.csv(res_var, file=paste(filePATH,"/",fileName,"-CV-result.csv", sep=""))
    rm(res, temp)

    #histogram of CV
    suppressMessages(dx <- melt(res_var[,uniSample]), classes = "message")
    plot1 <- ggplot(dx, aes(x=value)) +
      #geom_histogram(color="black", fill="skyblue", bins = 50) +
      geom_histogram(aes(y=..density..), colour="black", fill="skyblue", bins=50)+
      labs(x="CV") +
      #geom_density(alpha=.2, fill="#FF6666") +
      #geom_vline(xintercept = 5, color="red", linetype="dashed") +
      theme_classic()
    pdf(paste(filePATH,"/",fileName,"-CV-distribution.pdf", sep=""), width=5, height=5)
    print(plot1)
    dev.off()
    print(plot1)

    #Variable genes
    write.csv(variable_gene, file=paste(filePATH,"/",fileName,"-CV-VariableFeatures-result.csv", sep=""))
    #Stable genes
    write.csv(stable_gene, file=paste(filePATH,"/",fileName,"-CV-StableFeatures-result.csv", sep=""))

    #color list
    if(length(res_var$Median<0)>1) {
      col_list <- colorRamp2(c(-(cvThreshold*10),-(cvThreshold*4),-(cvThreshold*2),cvThreshold,0,cvThreshold,(cvThreshold*2),(cvThreshold*4),(cvThreshold*10)),
                             c("brown","red","yellow","blue","white","blue","yellow","red","brown"))
    } else {
      col_list <- colorRamp2(c(0,cvThreshold,(cvThreshold*2),(cvThreshold*4),(cvThreshold*10)),
                             c("white","blue","yellow","red","brown"))
    }

    #Plot variable genes CV
    mat <- res_var[order(abs(res_var$Median), decreasing = T),]
    mat <- mat[abs(mat$Median) > cvThreshold,]
    if(nrow(mat)>50) {
      mat <- mat[1:50,]
      write.csv(mat, file=paste(filePATH,"/",fileName,"-CV-Variable-Matrix.csv", sep=""))
    }
    mat <- mat[,uniSample]
    column_ha <- HeatmapAnnotation(CV = anno_boxplot(mat[,uniSample]))
    ht1 <- Heatmap(data.matrix(mat), cluster_rows =F,  cluster_columns = F,
                   na_col = "grey", col = col_list,
                   row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                   column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 6),
                   top_annotation = column_ha,
                   heatmap_legend_param = list(title = "CV(%)",heatmap_legend_side = "right") )

    #Plot stable genes CV
    mat <- res_var[order(abs(res_var$Median), decreasing = F),]
    mat <- mat[abs(mat$Median) <= cvThreshold,]
    if(nrow(mat)>50) {
      mat <- mat[1:50,]
      write.csv(mat, file=paste(filePATH,"/",fileName,"-CV-Stable-Matrix.csv", sep=""))
    }
    mat <- mat[,uniSample]
    column_ha <- HeatmapAnnotation(CV = anno_boxplot(mat[,uniSample]))
    ht2 <- Heatmap(data.matrix(mat), cluster_rows =F,  cluster_columns = F,
                   na_col = "grey", col = col_list,
                   row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                   column_names_gp = gpar(fontsize = 5),
                   row_names_gp = gpar(fontsize = 6),
                   top_annotation = column_ha,
                   heatmap_legend_param = list(title = "CV(%)",
                                               heatmap_legend_side = "right") )

    cat(date(),": Saving CV plots in output directory\n")
    pdf(paste(filePATH,"/",fileName,"-CV-Stable-Variable-Heatplot.pdf", sep=""),
        width=5, height=8)
    draw(ht1)
    draw(ht2)
    dev.off()
    draw(ht1)
    draw(ht2)

    #Housekeeping genes
    if(!is.null(housekeeping_genes)) {
      housekeeping_genes <- intersect(housekeeping_genes, row.names(res_var))
      mat <- res_var[housekeeping_genes,uniSample]
      column_ha <- HeatmapAnnotation(CV = anno_boxplot(mat[,uniSample]))
      ht3 <- Heatmap(data.matrix(mat), cluster_rows =F,  cluster_columns = F,
                     na_col = "grey", col = col_list,
                     row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                     column_names_gp = gpar(fontsize = 5),
                     row_names_gp = gpar(fontsize = 6),
                     top_annotation = column_ha,
                     heatmap_legend_param = list(title = "CV(%)",
                                                 heatmap_legend_side = "right")
                     )
      pdf(paste(filePATH,"/",fileName,"-CV-housekeeping_genes-Heatplot.pdf", sep=""),
          width=5, height=5)
      draw(ht3)
      dev.off()
    }

    cat(date(),": Done. Please check output directory for results.\n")
    res <- list(CV=res_var, variable_genes=variable_gene,
                stable_genes=stable_gene)
    return(res)
}
