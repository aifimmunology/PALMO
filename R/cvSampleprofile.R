#' A cvSampleprofile Function
#'
#' This function allows to calculate Intra-donor variations over longitudinal
#' timepoints. The coefficient of variation is calculated in single cell
#' data. It requires longitudinal data matrix/data frame and annotation file.
#' @param ann Annotation table. Table must consist column Sample (Participant 
#' sample name), PTID (Participant), Time (longitudinal time points)
#' @param mat Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column)
#' @param meanThreshold Average expression threshold to filter lowly expressed 
#' genes Default is 0.1 (log2 scale)
#' @param cvThreshold Coefficient of variation threshold to select variable and 
#' stable genes Default is 10 for single cell RNA (100*SD/mean)
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords cvSampleprofile
#' @export

cvSampleprofile <- function(mat, ann, meanThreshold=NULL, cvThreshold=NULL, fileName=NULL, filePATH=NULL) {
  
  cat(date(),": Performing Sample-wise Coefficient of variance analysis\n")
  
  #If filename or filepath null
  if(is.null(fileName)) {
    fileName <- "outputFile"
  }
  if(is.null(filePATH)) {
    filePATH <- paste(getwd(), "/output", sep="")
    dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  }
  
  #meanThrehold and cvThreshold
  if(is.null(meanThreshold)) {
    meanThreshold <- 0.1
    cat(date(),": Using mean threshold 0.1\n")
  }
  if(is.null(cvThreshold)) {
    cvThreshold <- 10
    cat(date(),": Using cv threshold 10\n")
  }
  
  #Calculate CV vs Mean for all genes per celltype
  unigene <- row.names(mat)
  uniSample <- sort(unique(ann$PTID))
  ann$group_donor <- paste(ann$group, ann$PTID, sep=":") 
  uniSamplegroup <- as.character(unique(ann$group_donor))
  
  #Variable genes
  cat(date(),": Plotting Sample wise CV analysis\n")
  pdf(paste(filePATH,"/",fileName,"-CV-SampleGroup-Plot.pdf", sep=""), width=5, height=5)
  op <- pboptions(type = "timer") # default
  res1 <- pblapply(uniSamplegroup,function(uS) {
    #print(uS)
    ann_df <- ann[ann$group_donor %in% uS,]
    if(nrow(ann_df) > 1 ) {
      df <- mat[unigene, ann_df$Sample_group]
      df <- data.frame(df, nonZero=apply(df,1,function(x){sum(x!=0)}), mean=rowMeans(df, na.rm=T), sd=apply(df,1,sd, na.rm=T), var=apply(df,1,var, na.rm=T), stringsAsFactors = F)
      df$CV <- 100*df$sd/df$mean
      #the CV becomes very high for data with 0
      df <- df[df$mean >= meanThreshold,] #minimum expression >2^0.1=1
      dp2a <- df[df$mean >= meanThreshold & df$CV > cvThreshold, c("mean", "sd", "var", "CV")]
      dp2a <- dp2a[order(dp2a$CV, dp2a$mean, decreasing = T),]
      dp2b <- df[df$mean >= meanThreshold & df$CV< cvThreshold,]
      dp2b <- dp2b[order(-dp2b$mean, dp2b$CV, decreasing = F),]
      plot1 <- ggplot(df, aes(x=mean, y=CV)) + geom_point(size=0.5, color="grey") +
        #scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
        labs(title=paste(uS, " timepoints=",nrow(ann_df), sep="")) +
        geom_text_repel(data=dp2a[1:10,], aes(x=mean, y=CV, label=row.names(dp2a[1:10,])), col="red", size=2, max.overlaps=20) +
        geom_text_repel(data=dp2b[1:10,], aes(x=mean, y=CV, label=row.names(dp2b[1:10,])), col="blue", size=2, max.overlaps=20) +
        theme_classic()
      print(plot1)
    }
      return(NULL)
  })
  pboptions(op)
  dev.off()
  
  cat(date(),": Done. Please check output directory for results.\n")
  
}
