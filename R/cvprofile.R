#' A cvprofile Function
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
#' @param housekeeping_genes Optional list of housekeeping genes to focus on.
#' Default is ACTB, GAPDH
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords cvprofile
#' @export

cvprofile <- function(mat, ann, meanThreshold=NULL, housekeeping_genes=NULL, fileName=NULL, filePATH=NULL) {

  cat(date(),": Performing Coefficient of variance analysis\n")

  #If filename or filepath null
  if(is.null(fileName)) {
    fileName <- "outputFile"
  }
  if(is.null(filePATH)) {
    filePATH <- paste(getwd(), "/output", sep="")
    dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  }
  if(is.null(housekeeping_genes)) {
    housekeeping_genes <- c("ACTB","GAPDH")
  }


  #meanThrehold and cvThreshold
  if(is.null(meanThreshold)) {
    meanThreshold <- 0
    cat(date(),": Using mean threshold >= 0\n")
  }

  #Calculate CV vs Mean for all genes per celltype
  unigene <- row.names(mat)
  uniSample <- sort(unique(ann$PTID))
  ann$group_donor <- paste(ann$group, ann$PTID, sep=":")
  uniSamplegroup <- as.character(unique(ann$group_donor))

  #All genes CV calculations
  cat(date(),": Performing CV calculations\n")
  op <- pboptions(type = "timer") # default
  res <- pblapply(uniSamplegroup,function(uS) {
    #print(uS)
    ann_df <- ann[ann$group_donor %in% uS,]
    if(nrow(ann_df) > 1 ) {
      df <- mat[unigene, ann_df$Sample_group]
      df <- data.frame(df, zeros=apply(df,1,function(x){sum(x!=0)}), mean=rowMeans(df, na.rm=T), sd=apply(df,1,sd, na.rm=T), var=apply(df,1,var, na.rm=T), stringsAsFactors = F)
      df$cv <- 100*df$sd/df$mean
      df$cv <- ifelse(df$mean >= meanThreshold, df$cv, NA)
      df <- df[,c("mean","sd","var","cv")]
      df$gene <- row.names(df)
      df$group <- uS
      return(df)
    }
  })
  pboptions(op)
  cv_all <- do.call(rbind, res)
  cv_all <- data.frame(cv_all, check.names=F, stringsAsFactors = F)
  cv_all$select <- ifelse(cv_all$mean >= meanThreshold, "Y", "N")

  df <- cv_all[cv_all$mean >= meanThreshold, ]
  p1 <- ggplot(cv_all, aes(x=mean)) +
      geom_histogram(aes(color=select), fill="white", binwidth=0.1) +
      labs(title="Mean expression (log10)") + scale_x_continuous(trans='log10')
  p2 <- ggplot(df, aes(x=cv)) +
      labs(title="CV (mean/SD %)") + geom_histogram(binwidth=1, color="black", fill="white")
  # Housekeeping genes data
  df <- df[df$gene %in% housekeeping_genes,]
  if(nrow(df)>0) {
    p3 <- ggplot(df, aes(x=mean, y=cv)) + geom_point() +
      labs(title="Housekeeping genes") + facet_wrap(~gene)

    pdf(paste(filePATH,"/",fileName,"-CVPlot.pdf", sep=""), width=12, height=5)
    print(plot_grid(p1, p2, p3, ncol=3))
    dev.off()
    print(plot_grid(p1, p2, p3, ncol=3))
  } else {
    pdf(paste(filePATH,"/",fileName,"-CVPlot.pdf", sep=""), width=8, height=5)
    print(plot_grid(p1, p2,ncol=2))
    dev.off()
    print(plot_grid(p1, p2, ncol=2))
  }

  cat(date(),": Done. Please check output directory for results.\n")
  return(cv_all)
}
