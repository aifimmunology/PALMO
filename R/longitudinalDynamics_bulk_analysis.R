#' A longitudinalDynamics_bulk_analysis Function
#'
#' This function allows you to perform analysis of longitudinal data from bulk
#' dataset. It requires longitudinal bulk data matrix/data frame at sample level
#' and annotation file.
#' @param metadata Annotation table must consist column Sample, 
#' Participant sample name; PTID, participant; Time, longitudinal time frame
#' @param datamatrix Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column)
#' @param datatype Data input can be bulk or singlecell
#' @param omics User defined name like RNA, ATAC, Proteomics, FLOW
#' @param fileName User defined filename
#' @param features Variance analysis carried out for the features provided such 
#' as c("PTID", "Time", "Sex")
#' @param meanThreshold Average expression threshold to filter lowly expressed 
#' genes Default is 0.1 (log2 scale)
#' @param cvThreshold Coefficient of variation threshold to select variable and 
#' stable genes Default is 10 for single cell RNA (100*SD/mean)
#' @param method Type of variations to be analyzed, c("intra-donor", "inter-donor",
#' "iqr"), Default is all
#' @param filePATH User-defined output directory path
#' @keywords longitudinalDynamics_bulk_analysis
#' @export
#' @examples
#' ##Bulk data
#' #longitudinalDynamics(metadata=ann, datamatrix=data, datatype="bulk", omics="Protein",
#' #fileName="olink", features=c("PTID", "Time", "Sex"), meanThreshold=1, cvThreshold=5,
#' #outputDirectory="output")

longitudinalDynamics_bulk_analysis <- function(metadata=NULL, datamatrix=NULL, datatype=NULL, omics=NULL,
    fileName=NULL, filePATH=NULL, features=NULL, method=NULL, meanThreshold=NULL, cvThreshold=NULL) {
    
    #------------------------------------------------
    if(is.null(metadata)) {
        cat(date(),": Please provide the annotation dataframe")
        stop()
    }
    if(is.null(datamatrix)) {
        cat(date(),": Please provide the expression datamatrix/dataframe")
        stop()
    }
    #------------------------------------------------
        
    #------------------------------------------------
    #PCA plot
    row.has.na <- apply(datamatrix,1,function(x){any(is.na(x))})
    datamatrix_nonNA <- datamatrix[!row.has.na,]
    if(nrow(datamatrix_nonNA)>5) {
        res.pca <- suppressMessages(prcomp(t(datamatrix_nonNA),  center=T, scale = TRUE), classes = "message")
        titleName <- paste(omics," PCA: sample=",ncol(datamatrix_nonNA), " gene=",nrow(datamatrix_nonNA), sep="")
        plot1 <- fviz_pca_ind(res.pca, col.ind = metadata$PTID, geom.ind =c("point", "text"),  labelsize = 3, addEllipses=FALSE, ellipse.level=0.95, title=titleName) + theme(axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black"), legend.position = "none") + theme_classic()
        pdf(paste(filePATH,"/",fileName,"-PCA-Plot.pdf", sep=""), width=6, height=6)
        print(plot1)
        dev.off()

        #TSNE plot
        set.seed(1)
        er1 <- try(Rtsne(t(datamatrix_nonNA), dims = 2, verbose=FALSE, max_iter = 500), silent=TRUE)
        if(inherits(er1, "try-error")) {
            tsne_model <- Rtsne(t(datamatrix_nonNA), dims = 2, perplexity=1, verbose=FALSE, max_iter = 500)
        } else {
            tsne_model <- Rtsne(t(datamatrix_nonNA), dims = 2, verbose=FALSE, max_iter = 500)
        }
        tsne_df <- data.frame(sample=colnames(datamatrix_nonNA), tSNE1=tsne_model$Y[,1], tSNE2=tsne_model$Y[,2], metadata, stringsAsFactors = F) ## getting the two dimension matrix
        plot2 <- ggplot(tsne_df, aes(x=tSNE1, y=tSNE2, color=factor(PTID))) +
            geom_point(size=2) +
            ggtitle("t-SNE") +
            geom_text_repel(data=tsne_df, aes(x=tSNE1, y=tSNE2), label=tsne_df$sample, size=2, segment.size = 0.1, segment.alpha=0.9, max.overlaps = 20, color="grey") +
            theme_classic() + theme(legend.position = "bottom")
        pdf(paste(filePATH,"/",fileName,"-UMAP-Plot.pdf", sep=""), width=6, height=6)
        print(plot2)
        dev.off()

        #Sample variability (Correlation)
        cor_mat <- rcorr(as.matrix(datamatrix_nonNA), type="pearson")
        res <- cor_mat$r
        #Plot heatmap
        ha_col <- HeatmapAnnotation(df=data.frame(PTID=metadata$PTID))
        ht1 <- Heatmap(data.matrix(res), cluster_rows =F,  cluster_columns = F,
                       row_split = as.character(metadata$PTID), column_split = as.character(metadata$PTID),
                       na_col = "grey", col = colorRamp2(c(-1,0,0.9,1), c("black","white","pink","red")),
                       row_names_max_width=unit(10, "cm"), #rect_gp = gpar(col = "black"),
                       column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize =6),
                       top_annotation = ha_col,
                       heatmap_legend_param = list(title = "Pearson",heatmap_legend_side = "right") )
        pdf(paste(filePATH,"/",fileName,"-SampleCorrelation-Heatplot.pdf", sep=""), width=10, height=10)
        draw(ht1)
        dev.off()
    } else {
        cat("All features have missing values.\n")
    }
    #------------------------------------------------

    #------------------------------------------------
    #Remove genes with >40%NAs
    row.na <- apply(datamatrix,1,function(x) { sum(is.na(x)) })
    if(length(row.na[row.na>0])>0) {
        row.non.na <- row.na[row.na < (0.4*ncol(datamatrix))] #select features with NAs <40%
        if(length(row.non.na)>1) {
            datamatrix <- datamatrix[names(row.non.na),]
        } else {
            stop("More than 40% missing values in data. Try to impute data and run again.\n")
        }
    }
    rowN <- data.frame(row.names(datamatrix), stringsAsFactors = F)
    #------------------------------------------------
    
    #------------------------------------------------
    #Variance decomposition
    if("inter-donor" %in% method) {
        lmem_res <- lmeVariance(ann=metadata, mat=datamatrix, features=features, meanThreshold=meanThreshold, fileName=fileName, filePATH=filePATH)
    } else {
        lmem_res <- NULL
    }
    #------------------------------------------------

    #------------------------------------------------
    #CV vs Mean
    if("intra-donor" %in% method) {
        cv_res <- cvCalcBulk(mat=datamatrix, ann=metadata, meanThreshold=meanThreshold, cvThreshold=cvThreshold, fileName=fileName, filePATH=filePATH)
        CV <- cv_res$CV
        variable_genes <- cv_res$variable_genes
        stable_genes <- cv_res$stable_genes
    } else {
        CV <- NULL
        variable_genes <- NULL
        stable_genes <- NULL
    }
    #------------------------------------------------

    #------------------------------------------------
    #Calculate IQR (outlier analysis)
    if("iqr" %in% method) {
        IQR_res <- iqrBulk(ann=metadata, mat=datamatrix, fileName=fileName, filePATH=filePATH)
    } else {
        IQR_res <- NULL
    }
    #----------------------------------
    
    #------------------------------------------------
    #return result object
    long_res <- list(
        CV=CV,
        variable_genes=variable_genes,
        stable_genes=stable_genes,
        varianceDecomposition=lmem_res,
        iqr=IQR_res
    )
    return(long_res)
    #------------------------------------------------
}
