#' A cvCalcBulk Function
#'
#' This function allows to calculate Intra-donor variations over longitudinal
#' timepoints. The coefficient of variation (CV) is calculated in Bulk
#' data without group information. CV calculated across samples. It requires 
#' longitudinal data matrix/data frame and annotation file.
#' @param ann Annotation table. the must consist column Sample, 
#' Participant sample name; PTID, participant; Time, longitudinal time frame
#' @param mat Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column)
#' @param meanThreshold Average expression threshold to filter lowly expressed 
#' genes Default is 0.1 (log2 scale)
#' @param cvThreshold Coefficient of variation threshold to select variable and 
#' stable genes Default is 5 for bulk olink data
#' @param housekeeping_genes Optional list of housekeeping genes to focus on. 
#' Default is NULL
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords cvCalcBulk
#' @export

cvCalcBulk <- function(mat=NULL, ann=NULL, meanThreshold=NULL, cvThreshold=NULL, housekeeping_genes=NULL, fileName=NULL, filePATH=NULL) {
    
    cat(date(),": Performing Coefficient of variance analysis\n")
    if(is.null(fileName)) {
      fileName <- "outputFile"
    }
    if(is.null(filePATH)) {
      filePATH <- paste(getwd(), "/output", sep="")
      dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
    }
    
    #CV vs Mean
    unigene <- row.names(mat)
    uniSample <- sort(unique(ann$PTID))
    pdf(paste(filePATH,"/",fileName,"-CV-Sample-Plot.pdf", sep=""), width=5, height=5)
    for(i in 1:length(uniSample)) {
      uS <- uniSample[i]
      meta_df <- ann[ann$PTID %in% uS,]
      df <- mat[unigene, meta_df$Sample]
      df <- data.frame(df, NAs=apply(df,1,function(x){sum(is.na(x))}), mean=rowMeans(df, na.rm=T), sd=apply(df,1,sd, na.rm=T), var=apply(df,1,var, na.rm=T), stringsAsFactors = F)
      df$CV <- 100*df$sd/df$mean
      summary(df$mean)
      summary(df$var)
      dp1 <- df
      dp1 <- dp1[dp1$NAs <= (0.5*nrow(meta_df)),] #missing value threshold 50%NAs
      dp2a <- dp1[dp1$mean >= meanThreshold & dp1$CV> cvThreshold,]
      dp2b <- dp1[dp1$mean >= meanThreshold & dp1$CV< cvThreshold,]
      dp2a <- dp2a[order(dp2a$CV, dp2a$mean, decreasing = T),]
      dp2b <- dp2b[order(-dp2b$mean, dp2b$CV, decreasing = F),]
      plot1 <- ggplot(dp1, aes(x=mean, y=CV)) + geom_point(size=0.5, color="grey") +
        #scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
        labs(title=uniSample[i]) +
        geom_text_repel(data=dp2a[1:10,], aes(x=mean, y=CV, label=row.names(dp2a[1:10,])), col="red", size=2, max.overlaps=20) +
        geom_text_repel(data=dp2b[1:10,], aes(x=mean, y=CV, label=row.names(dp2b[1:10,])), col="blue", size=2, max.overlaps=20) +
        theme_classic()
      print(plot1)
      
      #Find stable and variable genes
      if(i==1) {
        variable_gene <- data.frame(donor=uS, gene=row.names(dp2a), dp2a[,c("mean","sd","var","CV")], stringsAsFactors = F)
        stable_gene <- data.frame(donor=uS, gene=row.names(dp2b), dp2b[,c("mean","sd","var","CV")], stringsAsFactors = F)
      } else {
        variable_gene <- rbind(variable_gene, data.frame(donor=uS, gene=row.names(dp2a), dp2a[,c("mean","sd","var","CV")], stringsAsFactors = F))
        stable_gene <- rbind(stable_gene, data.frame(donor=uS, gene=row.names(dp2b), dp2b[,c("mean","sd","var","CV")], stringsAsFactors = F))
      }
      
      #Heatmap for variance
      temp <- data.frame(gene=row.names(df), var=df$CV, stringsAsFactors = F)
      colnames(temp)[-1] <- paste(uS, "_", colnames(temp)[-1], sep="")
      if(i==1) {
        res <- temp
      } else {
        res <- merge(res, temp, by="gene", all=T)
      }
    }
    dev.off()
    
    #decompose result into mean and CV (variance)
    res_var <- res[,-1]
    row.names(res_var) <- res$gene
    colnames(res_var) <- gsub("_var", "", colnames(res_var))
    res_var$Mean <- apply(res_var[,uniSample],1,function(x){mean(x,na.rm=T)})
    res_var$Median <- apply(res_var[,uniSample],1,function(x){median(x,na.rm=T)})
    write.csv(res_var, file=paste(filePATH,"/",fileName,"-CV-result.csv", sep=""))
    rm(res, temp)
    
    #histogram of CV
    suppressMessages(dx <- melt(res_var[,uniSample]), classes = "message")
    plot1 <- ggplot(dx, aes(x=value)) +
      #geom_histogram(color="black", fill="skyblue", bins = 50) +
      geom_histogram(aes(y=..density..), colour="black", fill="skyblue", bins = 50)+
      #geom_density(alpha=.2, fill="#FF6666") +
      geom_vline(xintercept = 5, color="red", linetype="dashed") +
      theme_classic()
    pdf(paste(filePATH,"/",fileName,"-CV-distribution.pdf", sep=""), width=5, height=5)
    print(plot1)
    dev.off()
    
    #Variable genes
    write.csv(variable_gene, file=paste(filePATH,"/",fileName,"-CV-VariableFeatures-result.csv", sep=""))
    #Stable genes
    write.csv(stable_gene, file=paste(filePATH,"/",fileName,"-CV-StableFeatures-result.csv", sep=""))
    
    #Plot variable genes CV
    mat <- res_var[order(res_var$Median, decreasing = T),]
    mat <- mat[mat$Median>cvThreshold,]
    if(nrow(mat)>50) {
      mat <- mat[1:50,]
      write.csv(mat, file=paste(filePATH,"/",fileName,"-CV-Variable-Matrix.csv", sep=""))
    }
    mat <- mat[,uniSample]
    column_ha <- HeatmapAnnotation(CV = anno_boxplot(mat[,uniSample]))
    ht1 <- Heatmap(data.matrix(mat), cluster_rows =F,  cluster_columns = F,
                   na_col = "grey", col = colorRamp2(c(0,5,10,20,50), c("white","blue","yellow","red","brown")),
                   row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                   column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 6),
                   top_annotation = column_ha,
                   heatmap_legend_param = list(title = "CV(%)",heatmap_legend_side = "right") )
    
    #Plot stable genes CV
    mat <- res_var[order(res_var$Median, decreasing = F),]
    mat <- mat[mat$Median<cvThreshold,]
    if(nrow(mat)>50) {
      mat <- mat[1:50,]
      write.csv(mat, file=paste(filePATH,"/",fileName,"-CV-Stable-Matrix.csv", sep=""))
    }
    mat <- mat[,uniSample]
    column_ha <- HeatmapAnnotation(CV = anno_boxplot(mat[,uniSample]))
    ht2 <- Heatmap(data.matrix(mat), cluster_rows =F,  cluster_columns = F,
                   na_col = "grey", col = colorRamp2(c(0,5,10,20,50), c("white","blue","yellow","red","brown")),
                   row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "white"),
                   column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 6),
                   top_annotation = column_ha,
                   heatmap_legend_param = list(title = "CV(%)",heatmap_legend_side = "right") )
    pdf(paste(filePATH,"/",fileName,"-CV-Stable-Variable-Heatplot.pdf", sep=""), width=5, height=8)
    draw(ht1)
    draw(ht2)
    dev.off()

    res <- list(CV=res_var, variable_genes=variable_gene, stable_genes=stable_gene)
    return(res)
}