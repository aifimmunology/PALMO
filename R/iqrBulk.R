#' A iqrBulk Function
#'
#' This function allows you to perform outlier analysis using defining IQR in 
#' bulk dataset. Outlier genes defined as exp > mean +/- 2SD.
#' @param ann Annotation table. the must consist column Sample, 
#' Participant sample name; PTID, participant; Time, longitudinal time frame
#' @param mat Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column)
#' @param SD_threshold Standard devaition limit to find outlliers (Eg. SD_threshold=
#' 2, equals to Mean+/-2SD)
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords iqrBulk
#' @export
#' @examples
#' #filePATH <- getwd()
#' #IQR_res <- iqrBulk(ann=metadata, mat=datamatrix, fileName="olink", filePATH=filePATH)

iqrBulk <- function(ann, mat, SD_threshold=NULL, fileName=NULL, filePATH=NULL) {

    cat(date(),": Performing Outlier anlaysis\n")
    if(is.null(fileName)) {
      fileName <- "outputFile"
    }
    if(is.null(filePATH)) {
      filePATH <- paste(getwd(), "/output", sep="")
      dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
    }
    if(is.null(SD_threshold)) {
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
    op <- pboptions(type = "timer") # default
    IQR_res <- pbapply(rowN,1,function(geneName) {
      df <- data.frame(exp=as.numeric(mat[geneName,]), ann, stringsAsFactors = F)
      #df$PTID <- factor(df$PTID, levels = uniSample)
      #plot1 <- ggpubr::ggline(df, x = "PTID", y = "exp", add.params = list(shape="Time"), add = c("mean_se", "jitter"), ylab = "NPX", xlab = "Donor", title = geneName, legend = "right", outlier.shape = NA) + scale_shape_manual(values = 0:10)
      #plot1 <- ggpubr::ggline(df, x = "Time", y = "exp", color="PTID", add.params = list(shape="PTID"), add = c("mean_se", "jitter"), ylab = "NPX", xlab = "Donor", title = geneName, legend = "right", outlier.shape = NA) + scale_shape_manual(values = 0:10)
      #plot1
      
      dfx <- lapply(uniSample, function(y) {
        temp <- df[df$PTID %in% y,]
        temp$z <- temp$exp - mean(temp$exp, na.rm=T)
        upper <- mean(temp$exp, na.rm=T) + (SD_threshold*sd(temp$exp, na.rm=T))
        lower <- mean(temp$exp, na.rm=T) - (SD_threshold*sd(temp$exp, na.rm=T))
        temp$outlier <- ifelse(temp$exp >= upper | temp$exp <= lower, temp$z, 0)
        return(data.frame(temp))
      })
      dfx <- do.call(rbind, dfx)
      return(dfx$outlier)
    })
    pboptions(op)
    IQR_res <- data.frame(t(IQR_res), stringsAsFactors = F)
    row.names(IQR_res) <- rowN[,1]
    colnames(IQR_res) <- colnames(mat)
    
    #Plot
    df <- melt(data.matrix(IQR_res))
    df <- df[df$value != 0,]
    df <- df[!is.na(df$Var1),]
    #df$Var2 <- factor(df$Var2, levels = sort(unique(ann$Sample)))
    plot1 <- ggplot(df, aes(x=Var2, y=value)) +
      geom_violin(scale="width") +
      #geom_boxplot(width=0.1, fill="white") +
      labs(x="", y="Z-score (>2SD)") +
      ggforce::geom_sina(size=0.5) +
      theme_classic() + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1, size=6), axis.text.y = element_text(size=6), legend.position = "right")
    pdf(paste(filePATH,"/",fileName,"-IQR-Boxplot.pdf", sep=""), width=10, height=5)
    print(plot1)
    dev.off()
    print(plot1)
    df <- df[order(df$value, decreasing = T),]
    colnames(df) <- c("Gene", "sample","Z")
    write.csv(df, file=paste(filePATH,"/",fileName,"-IQR-result.csv", sep=""), row.names = F)
    
    # genelist <- as.character(unique(df$Gene)[1:12])
    # splots <- list()
    # for(i in 1:length(genelist)) {
    #   geneName <- genelist[i]
    #   df <- data.frame(exp=as.numeric(mat[geneName,]), ann, stringsAsFactors = F)
    #   df$PTID <- factor(df$PTID, levels = uniSample)
    #   plot1 <- ggline(df, x = "PTID", y = "exp", add.params = list(shape="Time"), add = c("mean_se", "jitter"), ylab = "NPX", xlab = "Donor", title = geneName, legend = "right", outlier.shape = NA) +
    #     scale_shape_manual(values = 0:10) + theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
    #   splots[[i]] <- plot1
    # }
    # pdf(paste(filePATH,"/",fileName,"-IQR-Outliers.pdf", sep=""), width=15, height=15)
    # print(plot_grid(plotlist=splots, ncol= 3, align="hv"))
    # dev.off()

    return(IQR_res)
}