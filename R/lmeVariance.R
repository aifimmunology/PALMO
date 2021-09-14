#' A lmeVariance Function
#'
#' This function allows you to calculate inter-donor variation between participants
#' over longitudinal time points. It uses linear mixed model to calculate variance
#' contribution from each given feature list
#' @param ann Annotation table. the must consist column Sample, 
#' Participant sample name; PTID, participant; Time, longitudinal time frame
#' @param mat Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column)
#' @param features Variance analysis carried out for the features provided such 
#' as c("PTID", "Time", "Sex")
#' @param meanThreshold Average expression threshold to filter lowly expressed 
#' genes/features Default is 0
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory PATH Default, current directory
#' @keywords lmeVariance
#' @export
#' @examples
#' ##Input Expression data
#' #filePATH <- getwd()
#' #lmem_res <- lmeVariance(ann=metadata, mat=datamatrix,
#' #features=c("PTID", "Time", "Sex"),
#' #meanThreshold=0.1, fileName="RNA", filePATH=filePATH)

lmeVariance <- function(ann, mat, features, meanThreshold=NULL, fileName=NULL, filePATH=NULL) {
    
    cat(date(),": Performing variance decomposition\n")
    if(is.null(fileName)) {
      fileName <- "outputFile"
    }
    if(is.null(filePATH)) {
      filePATH <- paste(getwd(), "/output", sep="")
      dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
    }
    if(is.null(meanThreshold)) {
      meanThreshold <- 0
    }
    
    #Define featureList
    featureList <- c(features, "Residual")
    rowN <- data.frame(row.names(mat))
    op <- pboptions(type = "timer") # default
    suppressMessages(lmem_res <- pbapply(rowN,1,function(geneName) {
      df <- data.frame(exp=as.numeric(mat[geneName,]), ann, stringsAsFactors = F)
      
      #formula
      form <- paste(paste("(1|", features,")", sep=""), collapse=" + ")
      form <- as.formula(paste("exp ~ ", form, sep=""))
      #linear mixed effect model
      lmem <- lmer(formula=form, data = df)
      lmem_re <- as.data.frame(VarCorr(lmem))
      row.names(lmem_re) <- lmem_re$grp
      lmem_re <- lmem_re[featureList,]
      #lmem_re$vcov/sum(lmem_re$vcov)
      fix_effect <- fixef(lmem) #get fixed effect
      lmem_re$CV <- lmem_re$sdcor/fix_effect ##Calculate CV
      return(c(mean(df$exp, na.rm=T), median(df$exp, na.rm=T), sd(df$exp, na.rm=T), max(df$exp, na.rm=T),
               (lmem_re$vcov)/sum(lmem_re$vcov)) )
    }) , classes = "message")
    pboptions(op)
    lmem_res <- data.frame(t(lmem_res), stringsAsFactors = F)
    colnames(lmem_res) <- c("mean", "median", "sd", "max", featureList)
    lmem_res <- data.frame(Gene=rowN[,1], lmem_res, stringsAsFactors = F)
    row.names(lmem_res) <- lmem_res$Gene
    lmem_res <- lmem_res[order(lmem_res[,featureList[1]], decreasing = T),]
    write.csv(lmem_res, file=paste(filePATH,"/",fileName,"-Variance.csv", sep=""))
    
    #Overall gene mean result
    plot1 <- ggplot(lmem_res, aes(x=mean)) + geom_histogram() + labs(title="Mean Expression")
    plot2 <- ggplot(lmem_res, aes(x=log10(mean+1))) + geom_histogram() + labs(title="Mean Expression (log10)")
    pdf(paste(filePATH,"/",fileName,"-Meanplot.pdf", sep=""), width=7, height=3.5)
    suppressMessages(print(plot_grid(plot1, plot2, align="hv")), classes = "message")
    dev.off()
    
    #create data matrix
    res <- lmem_res[lmem_res$max > meanThreshold,featureList]
    df <- melt(data.matrix(res))
    df$value <- df$value * 100
    df$feature <- paste(df$Var1, df$Var2, sep="_")
    
    sigFeatures <- c()
    for(i in 1:length(features)) {
      dfx <- res[order(res[,features[i]], decreasing = T),]
      sigFeatures <- c(sigFeatures, paste(row.names(dfx[1:5,]),features[i],sep="_"))
    }
    df1 <- df[df$feature %in% sigFeatures,]
    plot1 <- ggplot(df, aes(x=Var2, y=value, fill=Var2)) +
      geom_violin(scale="width") +
      geom_boxplot(width=0.1, fill="white") +
      #ggpubr::stat_compare_means(label.y = 140) + # Add global p-value
      geom_text_repel(data=df1, aes(x=Var2, y=value), label=df1$Var1, size=2, segment.size = 0.1, segment.alpha=0.9, max.overlaps = 20, color="black") +
      labs(x="FeatureList", y="Variance Explained (%)") + theme_classic()
    #Plot
    pdf(paste(filePATH,"/",fileName,"-VarianceExplained-Boxplot.pdf", sep=""), width=5, height=5)
    print(plot1)
    dev.off()
    print(plot1)
    
    #Plot variable genes
    genelist <- unique(as.character(df1$Var1))
    uniSample <- sort(unique(ann$PTID))
    pdf(paste(filePATH,"/",fileName,"-VarianceExplained-Geneplot.pdf", sep=""), width=5, height=5)
    for(i in 1:length(genelist)) {
      geneName <- genelist[i]
      df <- data.frame(exp=as.numeric(mat[geneName,]), ann, stringsAsFactors = F)
      df$PTID <- factor(df$PTID, levels = uniSample)
      plot1 <- ggline(df, x = "PTID", y = "exp", add.params = list(shape="Time"), add = c("mean_se", "jitter", "boxplot"), ylab = "Expression", xlab = "Donor", title = geneName, legend = "right", outlier.shape = NA) +
        scale_shape_manual(values = 0:10) +
        theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1, color="black"), axis.text.y = element_text(color="black"))
      #ggplot(df, aes(x=PTID, y=exp)) + geom_boxplot() + geom_jitter(aes(shape=Time), position=position_jitter(0.2)) + labs(y="Expression") + theme_classic()
      print(plot1)
    }
    dev.off()
    
    #return data frame
    return(lmem_res)
}
