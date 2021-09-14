#' A iqrSC Function
#'
#' This function allows you to perform outlier analysis on single cell based
#' average expression matrix by defining IQR. Outlier genes defined as exp >
#' mean +/- 2SD.
#' @param ann Annotation table. the must consist column Sample, 
#' Participant sample name; PTID, participant; Time, longitudinal time frame
#' @param mat Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column)
#' @param groupName Group label used to calculate average gene expression
#' and check for outliers
#' @param SD_threshold Standard devaition limit to find outlliers (Eg. SD_threshold=
#' 2, equals to Mean+/-2SD)
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords iqrSC
#' @export
#' @examples
#' #filePATH <- getwd()
#' #IQR_res <- iqrSC(ann=metadata, mat=datamatrix, fileName="olink", filePATH=filePATH)

iqrSC <- function(ann, mat, groupName=NULL, SD_threshold=NULL, fileName=NULL, filePATH=NULL) {

    cat(date(),": Performing Outlier anlaysis\n")
    if(is.null(fileName)) {
      fileName <- "outputFile"
    }
    if(is.null(filePATH)) {
      filePATH <- getwd()
    }
    if(!is.null(groupName)) {
      ann$group <- ann[,groupName]
      ann$Sample_group <- paste(ann$Sample, ann$grop, sep=":")
    }
    if(!is.null(SD_threshold)) {
      SD_threshold <- 2 #+/- 2SD threshold
    }
  
    #Check overlap
    overlap <- intersect(row.names(ann), colnames(mat))
    ann <- ann[overlap,]
    mat <- mat[,overlap]
    
    #Calculate IQR (outlier analysis)
    rowN <- data.frame(row.names(mat), stringsAsFactors = F)
    uniTime <- as.character(unique(ann$Time))
    uniSample <- sort(unique(ann$PTID))
    uniGroup <- sort(unique(ann$group))
    
    ann <- ann[order(ann$PTID, ann$group),]
    uniSample_group <- unique(ann$Sample_group)
    mat <- mat[,row.names(ann)]
    #all.equal(row.names(ann), colnames(mat))
    
    op <- pboptions(type = "timer") # default
    IQR_res <- pbapply(rowN,1,function(geneName) {
      df <- data.frame(exp=as.numeric(mat[geneName,]), ann, stringsAsFactors = F)
      #df$PTID <- factor(df$PTID, levels = uniSample)
      #plot1 <- ggpubr::ggline(df, x = "PTID", y = "exp", add.params = list(shape="Time"), add = c("mean_se", "jitter"), ylab = "NPX", xlab = "Donor", title = geneName, legend = "right", outlier.shape = NA) + scale_shape_manual(values = 0:10)
      #plot1 <- ggpubr::ggline(df, x = "Time", y = "exp", color="PTID", add.params = list(shape="PTID"), add = c("mean_se", "jitter"), ylab = "NPX", xlab = "Donor", title = geneName, legend = "right", outlier.shape = NA) + scale_shape_manual(values = 0:10)
      #plot1
      
      dfx <- lapply(uniSample, function(y1) {
        temp <- df[df$PTID %in% y1,]
        res <- lapply(uniGroup, function(y2) {
          temp <- temp[temp$group %in% y2,]
          temp$z <- temp$exp - mean(temp$exp, na.rm=T)
          upper <- mean(temp$exp, na.rm=T) + (SD_threshold*sd(temp$exp, na.rm=T))
          lower <- mean(temp$exp, na.rm=T) - (SD_threshold*sd(temp$exp, na.rm=T))
          temp$outlier <- ifelse(temp$exp >= upper | temp$exp <= lower, temp$z, 0)
          return(temp)
        })
        res <- do.call(rbind, res)
        return(res)
      })
      dfx <- do.call(rbind, dfx)
      return(dfx)
      pb$tick()
      Sys.sleep(1 / totolIter)
    })
    pboptions(op)
    names(IQR_res) <- as.character(rowN[,1])
    
    res <- apply(rowN,1,function(geneName) {
      temp <- IQR_res[[geneName]]
      temp <- temp[uniSample_group,]
      return(temp$outlier)
    })
    res <- data.frame(t(res), stringsAsFactors = F)
    row.names(res) <- as.character(rowN[,1])
    colnames(res) <- uniSample_group
    
    #Plot
    df <- melt(data.matrix(res))
    df <- df[df$value != 0,]
    df <- df[!is.na(df$Var1),]
    cn <- data.frame(do.call(rbind, strsplit(as.character(df$Var2), split = ":")), stringsAsFactors = F)
    df$group <- as.character(cn$X2)
    df$Sample <- as.character(cn$X1)
    plot1 <- ggplot(df, aes(x=group, y=value, color=Sample)) +
      geom_violin(scale="width") +
      #geom_boxplot(width=0.1, fill="white") +
      labs(x="", y="Z-score (>2SD)") +
      ggforce::geom_sina(size=0.5) +
      #facet_wrap(~Sample, ncol = 1, scales = "free") +
      theme_classic() + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1, size=6), axis.text.y = element_text(size=6), legend.position = "right")
    
    img_width <- length(unique(df$group))+length(unique(df$Sample))
    pdf(paste(filePATH,"/",fileName,"-IQR-Boxplot.pdf", sep=""), width=img_width, height=5)
    print(plot1)
    dev.off()

    return(list(IQR_res=IQR_res, res=res))
}
