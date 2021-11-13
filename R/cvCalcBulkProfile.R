#' A cvCalcBulkProfile Function
#'
#' This function allows to calculate Intra-donor variations over longitudinal
#' timepoints. The coefficient of variation (CV) is calculated in Bulk
#' data without group information. CV calculated across samples. It requires 
#' longitudinal data matrix/data frame and annotation file.
#' @param ann Annotation table. Table must consist column Sample (Participant 
#' sample name), PTID (Participant), Time (longitudinal time points)
#' @param mat Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column)
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords cvCalcBulkProfile
#' @export

cvCalcBulkProfile <- function(ann, mat, fileName=NULL, filePATH=NULL) {
    
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
    
    cat(date(),": Performing CV calculations\n")
    op <- pboptions(type = "timer") # default
    res <- pblapply(uniSample,function(uS) {
      #print(uS)
      meta_df <- ann[ann$PTID %in% uS,]
      if(nrow(meta_df)>1) {
        df <- mat[unigene, meta_df$Sample]
        df <- data.frame(df, NAs=apply(df,1,function(x){sum(is.na(x))}), mean=rowMeans(df, na.rm=T), sd=apply(df,1,sd, na.rm=T), var=apply(df,1,var, na.rm=T), stringsAsFactors = F)
        df$CV <- 100*df$sd/df$mean
        df <- df[,c("mean","sd","var","CV")]
        df$gene <- row.names(df)
        df$group <- uS
        return(df)
        
        # #histogram of CV
        # plot1 <- ggplot(df, aes(x=mean, y=CV)) + geom_point(size=0.5, color="grey") +
        #   #scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
        #   labs(title=paste(uS, " timepoints=",nrow(meta_df), sep="")) +
        #   theme_classic()
        # print(plot1)
      }
    })
    pboptions(op)
    cv_all <- do.call(rbind, res)
    cv_all <- data.frame(cv_all, check.names=F, stringsAsFactors = F)
    return(cv_all)
}