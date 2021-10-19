#' A outlierDetect Function
#'
#' This function allows you to perform outlier analysis on bulk data by
#' calculating z-score. Outlier genes defined as
#' (exp-avgExp/SD) > mean + 2SD or (exp-avgExp/SD) < mean - 2SD.
#' @param ann Annotation table. Table must consist column Sample (Participant
#' sample name), PTID (Participant), Time (longitudinal time points)
#' @param mat Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column)
#' @param SD_threshold Standard deviation limit to find outlliers (Eg. SD_threshold=
#' 2, equals to Mean+/-2SD)
#' @param plotWidth User-defined plot width, Default 10 in
#' @param plotHeight User-defined plot height, Default 5 in
#' @param groupBy Include groupwise outlier analysis (TRUE or FALSE)
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords outlierDetect
#' @export
#' @examples
#' #filePATH <- getwd()
#' #outlier_res <- outlierDetect(ann=metadata, mat=datamatrix)

outlierDetect <- function(ann, mat, SD_threshold=2,
                          plotWidth=10, plotHeight=5,
                          groupBy=FALSE,
                          fileName=NULL, filePATH=NULL) {

    cat(date(),": Performing Outlier anlaysis\n")
    if(is.null(fileName)) {
      fileName <- "outputFile"
    }
    if(is.null(filePATH)) {
      filePATH <- paste(getwd(), "/output", sep="")
      dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
    }
    

    #Check overlap
    overlap <- intersect(row.names(ann), colnames(mat))
    ann <- ann[overlap,]
    mat <- mat[,overlap]

    #Input
    rowN <- data.frame(row.names(mat), stringsAsFactors = F)
    uniTime <- as.character(unique(ann$Time))
    uniSample <- sort(unique(ann$PTID))

    if(groupBy == FALSE) {

      #Calculate Z-score
      op <- pboptions(type = "timer") # default
      outlier_res <- pbapply(rowN,1,function(geneName) {
        df <- data.frame(exp=as.numeric(mat[geneName,]), ann, stringsAsFactors = F)
        dfx <- lapply(uniSample, function(y) {
          temp <- df[df$PTID %in% y,]
          if(nrow(temp)>0) {
            temp$gene <- geneName
            temp$meanDev <- temp$exp - mean(temp$exp, na.rm=T)
            temp$z <- (temp$exp - mean(temp$exp, na.rm=T))/(sd(temp$exp, na.rm=T))
            temp$outlier <- ifelse(abs(temp$z) >= SD_threshold, temp$z, 0)
            temp <- temp[temp$outlier != 0,]
            # #temp$z <- temp$exp - mean(temp$exp, na.rm=T)
            # temp$z <- (temp$exp - mean(temp$exp, na.rm=T))/(sd(temp$exp, na.rm=T))
            # upper <- mean(temp$exp, na.rm=T) + (SD_threshold*sd(temp$exp, na.rm=T))
            # lower <- mean(temp$exp, na.rm=T) - (SD_threshold*sd(temp$exp, na.rm=T))
            # temp$outlier <- ifelse(temp$exp >= upper | temp$exp <= lower, temp$z, 0)
            return(data.frame(temp))
          }
        })
        dfx <- do.call(rbind, dfx)
        return(dfx)
      })
      pboptions(op)
    } else {

      #Input group
      uniGroup <- sort(unique(ann$group))

      ann <- ann[order(ann$PTID, ann$group),]
      uniSample_group <- unique(ann$Sample_group)
      mat <- mat[,row.names(ann)]
      #all.equal(row.names(ann), colnames(mat))

      #Calculate Z-score (outlier analysis)
      op <- pboptions(type = "timer") # default
      outlier_res <- pbapply(rowN,1,function(geneName) {
        df <- data.frame(exp=as.numeric(mat[geneName,]), ann, stringsAsFactors = F)
        dfx <- lapply(uniSample, function(uS) {
          temp <- df[df$PTID %in% uS,]
          res <- lapply(uniGroup, function(uG) {
            temp <- temp[temp$group %in% uG,]
            if(nrow(temp)>0) {
              temp$gene <- geneName
              temp$meanDev <- temp$exp - mean(temp$exp, na.rm=T)
              temp$z <- (temp$exp - mean(temp$exp, na.rm=T))/(sd(temp$exp, na.rm=T))
              temp$outlier <- ifelse(abs(temp$z) >= SD_threshold, temp$z, 0)
              temp <- temp[temp$outlier != 0,]
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

    #Combine data
    outlier_res <- do.call(rbind, outlier_res)
    outlier_res <- outlier_res[!is.na(outlier_res$Sample),]
    write.csv(outlier_res, file=paste(filePATH,"/",fileName,"-Outlier-result.csv", sep=""), row.names = F)

    #Plot
    df <- outlier_res
    df$Sample <- factor(df$Sample, levels = unique(ann$Sample))
    #Z-plot
    plot1 <- ggplot(df, aes(x=Sample, y=outlier)) +
      geom_violin(scale="width") +
      #geom_boxplot(width=0.1, fill="white") +
      labs(x="", y=paste("Z-score (>",SD_threshold,"SD)"), title="Z (>SD)") +
      ggforce::geom_sina(size=0.5) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1, size=6),
            axis.text.y = element_text(size=6), legend.position = "right")

    #Mean-deviation plot
    plot2 <- ggplot(df, aes(x=Sample, y=meanDev)) +
      geom_violin(scale="width") +
      #geom_boxplot(width=0.1, fill="white") +
      labs(x="", y="Mean Deviation", title="Mean Deviation") +
      ggforce::geom_sina(size=0.5) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1, size=6),
            axis.text.y = element_text(size=6), legend.position = "right")

    if(groupBy == FALSE) {
      #Plot
      pdf(paste(filePATH,"/",fileName,"-Outlier-Boxplot.pdf", sep=""),
        width=plotWidth, height=plotHeight)
      print(plot1)
      print(plot2)
      dev.off()
    } else {
      plot1b <- ggplot(df, aes(x=Sample, y=outlier, color=group)) +
        geom_violin(scale="width") +
        #geom_boxplot(width=0.1, fill="white") +
        labs(x="", y=paste("Z-score (>",SD_threshold,"SD)"), title="Z (>SD)") +
        ggforce::geom_sina(size=0.5) +
        theme_classic() +
        theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1, size=6),
              axis.text.y = element_text(size=6), legend.position = "right")

      plot2b <- ggplot(df, aes(x=Sample, y=meanDev, color=group)) +
        geom_violin(scale="width") +
        #geom_boxplot(width=0.1, fill="white") +
        labs(x="", y="Mean Deviation", title="Mean Deviation") +
        ggforce::geom_sina(size=0.5) +
        theme_classic() +
        theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1, size=6),
              axis.text.y = element_text(size=6), legend.position = "right")

      #Plot
      pdf(paste(filePATH,"/",fileName,"-Outlier-Boxplot.pdf", sep=""),
          width=plotWidth, height=plotHeight)
      print(plot1)
      print(plot2)
      print(plot1b)
      print(plot2b)
      dev.off()
    }

    #Visualize
    print(plot2)
    print(plot1)

    # #Calculate MAD
    # op <- pboptions(type = "timer") # default
    # outlier_res <- pbapply(rowN,1,function(geneName) {
    #   df <- data.frame(exp=as.numeric(mat[geneName,]), ann, stringsAsFactors = F)
    #   dfx <- lapply(uniSample, function(y) {
    #     temp <- df[df$PTID %in% y,]
    #     temp$gene <- geneName
    #     MAD <- median(abs(temp$exp - mean(temp$exp, na.rm=T)))
    #     thresold <- SD_threshold*MAD
    #     temp$outlier <- ifelse(abs(temp$exp) >= SD_threshold, temp$exp, 0)
    #     temp <- temp[temp$outlier != 0,]
    #     # #temp$z <- temp$exp - mean(temp$exp, na.rm=T)
    #     # temp$z <- (temp$exp - mean(temp$exp, na.rm=T))/(sd(temp$exp, na.rm=T))
    #     # upper <- mean(temp$exp, na.rm=T) + (SD_threshold*sd(temp$exp, na.rm=T))
    #     # lower <- mean(temp$exp, na.rm=T) - (SD_threshold*sd(temp$exp, na.rm=T))
    #     # temp$outlier <- ifelse(temp$exp >= upper | temp$exp <= lower, temp$z, 0)
    #     return(data.frame(temp))
    #   })
    #   dfx <- do.call(rbind, dfx)
    #   return(dfx)
    # })
    # pboptions(op)
    # #Combine data
    # outlier_res <- do.call(rbind, outlier_res)

    df <- outlier_res[order(outlier_res$outlier, decreasing = T),]
    cat(date(),": Please check output directory for results\n")
    return(df)
}
