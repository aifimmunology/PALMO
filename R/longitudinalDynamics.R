#' A longitudinalDynamics Function
#'
#' This function allows you to perform analysis of longitudinal dataset.
#' It requires longitudinal data matrix/data frame and annotation file.
#' @param metadata Annotation table. Table must consist column Sample (Participant
#' sample name), PTID (Participant), Time (longitudinal time points)
#' @param data Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column).
#' For single cell, Single cell RNA Seurat object, if datatype is single cell RNA
#' and Single cell ATAC genescore matrix or data frame
#' @param datatype Data input can be bulk or singlecell
#' @param omics User defined name like RNA, ATAC, Proteomics, FLOW
#' @param fileName User defined filename
#' @param featureSet Variance analysis carried out on the featureSet provided such
#' as c("PTID", "Time", "Sex")
#' @param meanThreshold Average expression threshold to filter lowly expressed
#' genes Default is 0.1 (log2 scale)
#' @param cvThreshold Coefficient of variation threshold to select variable and
#' stable genes Default is 10 for single cell RNA (100*SD/mean)
#' @param NA_threshold Number of NAs in data (numeric value or NULL). Default,
#' 40% * number of columns.
#' @param column_sep Separator of "PTID" and "Time" in "Sample" column of
#' Annotation table like column_sep="W" for PTID1W1, column_sep=":" for
#' PTID1W1:Tcell
#' @param clusterBy for sample correlation cluster columns by ("donor", "group")
#' @param coding_genes Selecting protein coding/user-defined gene list only
#' @param avgGroup Group label to be used to calculate average gene expression by
#' group label
#' @param housekeeping_genes Optional list of housekeeping genes to focus on
#' Default is NULL
#' @param group_oi Group of interest to focus on, Default is NULL
#' @param nPC Number of PCAs to be used for UMAP, Default is 15
#' @param donorThreshold Donor threshold number to be used, Default is number of
#' participants
#' @param groupThreshold Group label threshold number to be used, Default is
#' (number of participants x group labels)/2
#' @param topFeatures Number of features to be selected from each group, Default
#' is 25
#' @param method Sample correlation analysis ("pearson", "spearman"). Default is
#' "spearman"
#' @param doOutlier To perform outlier analysis (TRUE or FALSE). Default FALSE
#' @param z_cutoff |Z| cutoff threshold to find potential outliers (Eg. z_cutoff=
#' 2, equals to Mean/SD 2)
#' @param outputDirectory User-defined output directory Default, output
#' @keywords longitudinalDynamics
#' @export

longitudinalDynamics <- function(metadata=NULL, data=NULL, datatype=NULL, omics=NULL,
    featureSet=NULL, meanThreshold=1, cvThreshold=5,
    NA_threshold=0.4, column_sep=NULL,
    coding_genes=NULL, avgGroup=NULL,
    housekeeping_genes=c("ACTB", "GAPDH"), group_oi=NULL, nPC=15,
    donorThreshold=NULL, groupThreshold=NULL, topFeatures=25,
    method="spearman", clusterBy="donor", z_cutoff=2,
    doOutlier=FALSE, fileName=NULL, outputDirectory=NULL) {

    #------------------------------------------------
    #Create output directory
    mainDir <- getwd()
    if(is.null(outputDirectory)) {
        outputDirectory <- "output"
    }
    subDir <- outputDirectory
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    filePATH <- file.path(mainDir, subDir)
    cat(date(),": The output directory is created.\n")

    #Filename
    if(is.null(fileName)) {
      fileName <- "outputFile"
    }
    #------------------------------------------------

    #------------------------------------------------
    if(is.null(metadata)) {
        cat(date(),": Please provide the annotation dataframe")
        stop()
    } else {
        metadata <- data.frame(metadata, check.names = F, stringsAsFactors = F)
    }

    #check Metadata
    check_features <- intersect(c("Sample", "PTID", "Time"), colnames(metadata))
    if(length(check_features) != 3) {
        stop("Please include columns: Sample, PTID, Time in your metdata")
    }
    #------------------------------------------------

    #------------------------------------------------
    #Check enough samples are available
    row.names(metadata) <- metadata$Sample
    uniSample <- as.character(sort(unique(metadata$PTID)))
    nPTID <- unique(metadata$PTID)
    nSample <- table(metadata$PTID)

    #Is dataset is single cell or Bulk data
    if(datatype == "bulk") {
        cat(date(),": Bulk data Analysis started\n")
        datamatrix <- data
        overlap <- intersect(metadata$Sample, colnames(datamatrix))
        if(length(overlap)==0) {
          cat(date(),": No overlap between the data columnn and annotation samples")
          stop()
        } else {
          cat(date(),": Number of samples matched:",length(overlap),"\n")
        }

        metadata <- metadata[overlap,]
        datamatrix <- data.frame(datamatrix, check.names = F, stringsAsFactors = F)
        datamatrix <- datamatrix[,overlap]

        #Remove genes with >40%NAs (optional)
        if(!is.null(NA_threshold)) {
          row.na <- apply(datamatrix,1,function(x) { sum(is.na(x)) })
          row.non.na <- row.na[row.na < (NA_threshold*ncol(datamatrix))] #select features with NAs say <40%
          datamatrix <- datamatrix[names(row.non.na),]
        }
        rowN <- data.frame(row.names(datamatrix), stringsAsFactors = F)

        #Features contributing towards donor variations
        lmem_res <- lmeVariance(ann=metadata, mat=datamatrix,
                                featureSet=featureSet,
                                meanThreshold=meanThreshold,
                                fileName = fileName,
                                filePATH = filePATH)


        #Intra-donor variations over time
        #CV vs Mean
        cv_res <- cvCalcBulk(mat=datamatrix,
                             ann=metadata,
                             meanThreshold=meanThreshold,
                             cvThreshold=cvThreshold,
                             fileName=fileName,
                             filePATH=filePATH)

        #Detect outliers (if any)
        if(doOutlier == TRUE) {
          #Sample variability (Correlation)
          cor_res <- sample_correlation(data=datamatrix,
                                      column_sep = column_sep,
                                      method="spearman",
                                      clusterBy=clusterBy,
                                      fileName=fileName,
                                      filePATH=filePATH)

          outlier_res <- outlierDetect(ann=metadata,
                                     mat=datamatrix,
                                     z_cutoff=z_cutoff,
                                     fileName=fileName,
                                     filePATH=filePATH)

          #Final result
          final_res <- list(varDecomposition=lmem_res,
                            intraDonorCV=cv_res,
                            sample_cor_res=cor_res,
                            outlier_res=outlier_res)
        } else {
          #Final result
          final_res <- list(varDecomposition=lmem_res,
                            intraDonorCV=cv_res)
        }

    } else if(datatype == "singlecell" & omics == "rna") {
        cat(date(),": scRNA Analysis started\n")
        dataObj <- data
        # Sample overlap
        overlap <- intersect(metadata$Sample, dataObj@meta.data$Sample)
        cat(date(),": Number of samples matched:",length(overlap),"\n")

        metadata <- metadata[overlap,]
        #in-case subset of samples only
        dataObj <- subset(x = dataObj, subset = Sample %in% overlap)

        #Aggregate data at celltypes (psuedo-bulk)
        metaData <- dataObj@meta.data
        metadata1 <- metadata[metaData$Sample,]
        metaData <- cbind(metaData, metadata1)
        dataObj@meta.data <- metaData

        if(is.null(avgGroup)) {
          stop("avgGroup parameter need to be provided. For example Sex or Celltype.")
        } else if(length(intersect(avgGroup, colnames(metaData)))==0) {
          stop("avgGroup parameter need to be provided. For example Sex or Celltype.")
        }

        #Define sample group and Calculate average expression
        dataObj@meta.data$Sample_group <- paste(dataObj@meta.data$Sample,
                                                dataObj@meta.data[,avgGroup],
                                                sep=":")
        dataObj@meta.data$Sample_group <- gsub(" ", "_", dataObj@meta.data$Sample_group)
        metaData <- dataObj@meta.data

        #Calculate average expression across group/celltype
        Idents(dataObj) <- "Sample_group"
        #table(Idents(dataObj))
        scrna_avgmat <- avgExpCalc(dataObj=dataObj, group.by="Sample_group")

        #Keep genes with avgExpression > zero
        rowDF <- rowSums(scrna_avgmat)
        rowDF <- rowDF[rowDF > 0]
        mat <- scrna_avgmat[names(rowDF),]

        #Create data annotation
        cn <- data.frame(Sample_group=colnames(mat))
        temp <- data.frame(do.call(rbind, strsplit(cn$Sample_group, split = ":")), stringsAsFactors = F)
        cn <- data.frame(cn, Sample=temp$X1, group=temp$X2, stringsAsFactors = F)
        row.names(cn) <- cn$Sample_group
        cn <- merge(cn, metadata, by="Sample", all=TRUE)
        cn <- cn[!is.na(cn$Sample_group),]
        row.names(cn) <- cn$Sample_group
        ann <- cn
        ann$Sample_group_i <- paste(ann$group, ann$PTID, sep=":")
        rm(cn)

        #Final matrix
        Overlap <- intersect(colnames(mat), row.names(ann))
        ann <- ann[Overlap,]
        mat <- mat[,Overlap]

        #Check the mean expression and CV cross groups (celltypes)
        cv_profile <- cvprofile(mat=mat,
                                ann=ann,
                                meanThreshold=meanThreshold,
                                fileName=fileName,
                                filePATH=filePATH)

        #Sample Celltype Mean-CV plot (check output folder)
        cvSampleprofile(mat=mat, ann=ann,
                        meanThreshold=meanThreshold,
                        cvThreshold = cvThreshold,
                        fileName=fileName,
                        filePATH=filePATH)

        #Variance decomposition
        lmem_res <- lmeVariance(ann=ann, mat=mat,
                                featureSet=c(featureSet,"group"),
                                meanThreshold=meanThreshold,
                                fileName=fileName,
                                filePATH=filePATH)

        #Calculate CV
        cv_res <- cvCalcSC(mat=mat,
                           ann=ann,
                           meanThreshold=meanThreshold,
                           cvThreshold = cvThreshold,
                           housekeeping_genes=housekeeping_genes,
                           fileName=fileName,
                           filePATH=filePATH)

        #Find stable and variable features in longitudinal data
        stable_gene <- StableFeatures(ann=ann,
                                group_oi=group_oi,
                                meanThreshold=meanThreshold,
                                cvThreshold = cvThreshold,
                                donorThreshold=donorThreshold,
                                groupThreshold=groupThreshold,
                                topFeatures=topFeatures,
                                fileName=fileName,
                                filePATH=filePATH)

        var_gene <- VarFeatures(ann=ann,
                                group_oi=group_oi,
                                meanThreshold=meanThreshold,
                                cvThreshold = cvThreshold,
                                donorThreshold=donorThreshold,
                                groupThreshold=groupThreshold,
                                topFeatures=topFeatures,
                                fileName=fileName,
                                filePATH=filePATH)

        #UMAP Plot
        #Top variable and stable features used for UMAP
        dimUMAPPlot(rnaObj=dataObj,
                    nPC=nPC,
                    gene_oi=unique(stable_gene$gene),
                    groupName=avgGroup,
                    plotname="stable",
                    fileName=fileName, filePATH=filePATH)
        dimUMAPPlot(rnaObj=dataObj,
                    nPC=nPC,
                    gene_oi=unique(var_gene$gene),
                    groupName=avgGroup,
                    plotname="variable",
                    fileName=fileName, filePATH=filePATH)

        #Sample correlation and outlier plot
        if(doOutlier == TRUE) {
          cor_res <- sample_correlation(data=mat,
                                        column_sep = column_sep,
                                        method=method,
                                        clusterBy="group",
                                        fileName=fileName,
                                        filePATH=filePATH)
          outlier_res <- outlierDetect(ann=ann, mat=mat,
                                       z_cutoff=z_cutoff,
                                       groupBy=TRUE,
                                       fileName=fileName,
                                       filePATH=filePATH)

          #Final result
          final_res <- list(varDecomposition=lmem_res,
                            intraDonorCV=cv_res,
                            stable_gene=stable_gene,
                            var_gene=var_gene,
                            sample_cor_res=cor_res,
                            outlier_res=outlier_res)
        } else {
          #Final result
          final_res <- list(varDecomposition=lmem_res,
                            intraDonorCV=cv_res,
                            stable_gene=stable_gene,
                            var_gene=var_gene)
        }

    } else if(datatype == "singlecell" & omics == "atac") {
        cat(date(),": scATAC Analysis started\n")
        atacObj <- data

        #Get the annotation
        temp <- data.frame(do.call(rbind, strsplit(colnames(atacObj),split = ":")), stringsAsFactors = F)
        cn <- data.frame(id=colnames(atacObj), Sample=temp$X1, group=temp$X2, check.names=F, stringsAsFactors = F)
        cn <- cn[cn$Sample %in% metadata$Sample,]

        #Sample overlap
        overlap <- as.character(unique(cn$Sample))
        cat(date(),": Number of samples matched:",length(overlap),"\n")
        atac_overlap <- cn$id

        #Sample overlap
        metadata <- metadata[overlap,]
        atacObj <- atacObj[,atac_overlap]

        #Keep genes with avgExpression > zero
        rowDF <- rowSums(atacObj)
        rowDF <- rowDF[rowDF > 0]
        mat <- atacObj[names(rowDF),]

        #Create annotation
        cn <- data.frame(Sample_group=colnames(mat))
        temp <- data.frame(do.call(rbind, strsplit(cn$Sample_group,
                           split = ":")), stringsAsFactors = F)
        cn <- data.frame(cn, Sample=temp$X1, group=temp$X2,
                         stringsAsFactors = F)
        row.names(cn) <- cn$Sample_group
        cn <- merge(cn, metadata, by="Sample", all=TRUE)
        cn <- cn[!is.na(cn$Sample_group),]
        row.names(cn) <- cn$Sample_group
        ann <- cn
        ann$Sample_group_i <- paste(ann$group, ann$PTID, sep=":")
        rm(cn)

        #Final matrix
        Overlap <- intersect(colnames(mat), row.names(ann))
        ann <- ann[Overlap,]
        mat <- mat[,Overlap]

        #CV profile
        #Check the mean expression and CV cross groups (celltypes)
        cv_profile <- cvprofile(mat=mat,
                                ann=ann,
                                meanThreshold=meanThreshold,
                                fileName=fileName,
                                filePATH=filePATH)

        #Sample Celltype Mean-CV plot (check output folder)
        cvSampleprofile(mat=mat, ann=ann,
                        meanThreshold=meanThreshold,
                        cvThreshold = cvThreshold,
                        fileName=fileName,
                        filePATH=filePATH)

        #Variance decomposition
        lmem_res <- lmeVariance(ann=ann, mat=mat,
                                featureSet=c(featureSet,"group"),
                                meanThreshold=meanThreshold,
                                fileName=fileName,
                                filePATH=filePATH)

        #Calculate CV
        cv_res <- cvCalcSC(mat=mat,
                           ann=ann,
                           meanThreshold=meanThreshold,
                           cvThreshold = cvThreshold,
                           housekeeping_genes=housekeeping_genes,
                           fileName=fileName,
                           filePATH=filePATH)

        #Find stable and variable features in longitudinal data
        stable_gene <- StableFeatures(ann=ann,
                                      group_oi=group_oi,
                                      meanThreshold=meanThreshold,
                                      cvThreshold = cvThreshold,
                                      donorThreshold=donorThreshold,
                                      groupThreshold=groupThreshold,
                                      topFeatures=topFeatures,
                                      fileName=fileName,
                                      filePATH=filePATH)

        var_gene <- VarFeatures(ann=ann,
                                group_oi=group_oi,
                                meanThreshold=meanThreshold,
                                cvThreshold = cvThreshold,
                                donorThreshold=donorThreshold,
                                groupThreshold=groupThreshold,
                                topFeatures=topFeatures,
                                fileName=fileName,
                                filePATH=filePATH)

        #Sample correlation and outlier plot
        if(doOutlier == TRUE) {
          cor_res <- sample_correlation(data=mat,
                                      column_sep = column_sep,
                                      method="spearman",
                                      clusterBy="group",
                                      fileName=fileName,
                                      filePATH=filePATH)
          outlier_res <- outlierDetect(ann=ann, mat=mat,
                                     z_cutoff=z_cutoff,
                                     groupBy=TRUE,
                                     fileName=fileName,
                                     filePATH=filePATH)

          #Final result
          final_res <- list(varDecomposition=lmem_res,
                            intraDonorCV=cv_res,
                            stable_gene=stable_gene,
                            var_gene=var_gene,
                            sample_cor_res=cor_res,
                            outlier_res=outlier_res)
        } else {
          #Final result
          final_res <- list(varDecomposition=lmem_res,
                            intraDonorCV=cv_res,
                            stable_gene=stable_gene,
                            var_gene=var_gene)
        }



    } else {
        stop("ERROR: Please provide datatype and omics name")
    }
    #------------------------------------------------

    #------------------------------------------------
    cat(date(),": Done\n")
    return(final_res)
    #------------------------------------------------

}
