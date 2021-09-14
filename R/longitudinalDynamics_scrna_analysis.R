#' A longitudinalDynamics_scrna_analysis Function
#'
#' This function allows you to perform analysis of longitudinal dataset from
#' single cell RNA.
#' @param metadata Annotation table must consist column Sample, 
#' Participant sample name; PTID, participant; Time, longitudinal time frame
#' @param dataObj Single cell Expression data. The scRNA data should be seurat
#' object. To get more information visit https://cran.r-project.org/web/packages/SeuratObject/SeuratObject.pdf
#' @param datatype Default singlecell
#' @param omics Default scRNA
#' @param fileName User defined filename
#' @param features Variance analysis carried out for the features provided such 
#' as c("PTID", "Time", "Sex")
#' @param meanThreshold Average expression threshold to filter lowly expressed 
#' genes Default is 0.1 (log2 scale)
#' @param cvThreshold Coefficient of variation threshold to select variable and 
#' stable genes Default is 10 for single cell RNA (100*SD/mean)
#' @param coding_genes Selecting protein coding/user-defined gene list only  
#' @param avgGroup Group label to be used to calculate average gene expression by
#' group label
#' @param housekeeping_genes Optional list of housekeeping genes to focus on
#' Default is NULL
#' @param group_oi Group of interest to focus on, Default is NULL
#' @param nPC Number of PCAs to be used for UMAP, Default is 30
#' @param donorThreshold Donor threshold number to be used, Default is number of
#' participants
#' @param groupThreshold Group label threshold number to be used, Default is 
#' (number of participants x group labels)/2
#' @param topFeatures Number of features to be selected from each group, Default
#' is 25
#' @param method Type of variations to be analyzed, c("intra-donor", "inter-donor",
#' "iqr"), Default is all
#' @param filePATH User-defined output directory path
#' @keywords longitudinalDynamics_scrna_analysis
#' @export
#' @examples
#' ##Single cell RNA data
#' #load("data/pbmc1-subset.Rda")
#' #load("data/data_Annotation.Rda")
#' #metadata=ann
#' #long_res <- longitudinalDynamics(metadata=ann, scrnaObj=pbmc,
#' #datatype="singlecell", omics="rna",
#' #fileName="scrna", features=c("PTID", "Time"),
#' #meanThreshold=0.1, cvThreshold=10,
#' #coding_genes=FALSE, avgGroup="celltype",
#' #housekeeping_genes=housekeeping_genes,
#' #nPC=15,
#' #donorThreshold=4, groupThreshold=40, topFeatures=25,
#' #method=c("intra-donor","inter-donor","iqr"),
#' #outputDirectory="output")

longitudinalDynamics_scrna_analysis <- function(metadata=NULL, dataObj=NULL, datatype="singlecell", omics="scRNA",
                            fileName=NULL, filePATH=NULL, features=NULL, meanThreshold=0.1, cvThreshold=10,
                            coding_genes=NULL, avgGroup=NULL,
                            housekeeping_genes=NULL, group_oi=NULL, nPC=NULL,
                            method=NULL, 
                            donorThreshold=NULL, groupThreshold=NULL, topFeatures=25) {
    
    #------------------------------------------------
    # #Label celltypes if not available
    # if(!is.null(celllabel)) {
    #     if(celllabel == TRUE) {
    #         cat(date(),": Applying cell labels\n")
    #         #Label celltypes
    #         dataObj <- labelingCelltypes(srnaObj=dataObj)
    #         dataObj@meta.data$celltype <- dataObj@meta.data$predicted.celltype.l2
    #         dataObj@meta.data$celltype_score <- dataObj@meta.data$predicted.celltype.l2.score
    #     }
    # }
    #------------------------------------------------

    #------------------------------------------------
    #For single cell data merge annotation and single cell metadata
    metaData <- dataObj@meta.data
    row.names(metadata) <- metadata$Sample
    metadata1 <- metadata[metaData$Sample,]
    metaData <- cbind(metaData, metadata1)
    dataObj@meta.data <- metaData
    #------------------------------------------------

    #------------------------------------------------
    cat(date(),": Defining sample group [Sample -",avgGroup,"]\n")
    #Define sample group and Calculate average expression
    if(is.null(avgGroup)) {
        stop("avgGroup parameter need to be provided. For example Sex or Celltype.")
    } else if(length(intersect(avgGroup, colnames(metaData)))==0) {
        stop("avgGroup parameter need to be provided. For example Sex or Celltype.")
    }
    dataObj@meta.data$Sample_group <- paste(dataObj@meta.data$Sample, dataObj@meta.data[,avgGroup], sep=":")
    dataObj@meta.data$Sample_group <- gsub(" ", "_", dataObj@meta.data$Sample_group)
    metaData <- dataObj@meta.data
    #Save the metadata dataframe
    save(metaData, file=paste(filePATH,"/",datatype,"-",fileName,"-metaData.Rda", sep=""))
    
    # #cat(unique(dataObj@meta.data$Sample_group))
    # scrna_avgmat <- averageExpression(rnaObject=dataObj, type="expression")
    # #Save the result dataframe
    # save(scrna_avgmat, file=paste("output/",datatype,"-",fileName,"-averageExpression.Rda", sep=""))
    # ov1 <- intersect(colnames(scrna_avgmat1), colnames(scrna_avgmat))
    # ge1 <- intersect(row.names(scrna_avgmat1), row.names(scrna_avgmat))
    # scrna_avgmat1 <- scrna_avgmat1[ge1,ov1]
    # scrna_avgmat <- scrna_avgmat[ge1,ov1]
    # plot(as.numeric(scrna_avgmat1[,10]), as.numeric(scrna_avgmat[,10]))

    cat(date(),": Performing scRNA average Expression\n")
    #Calculate average expression
    Idents(dataObj) <- "Sample_group"
    #table(Idents(dataObj))
    dataObj@assays$RNA@counts <- dataObj@assays$RNA@data #average expression on log-scaled data
    scrna_avgmat <- AverageExpression(object=dataObj, assays="RNA", slot="counts", group.by="Sample_group", verbose=TRUE)
    cn <- data.frame(colnames(scrna_avgmat$RNA))
    cn <- gsub("data\\[, 1]", "", as.character(row.names(cn)))
    scrna_avgmat <- data.frame(scrna_avgmat$RNA, check.names = F, stringsAsFactors = F)
    colnames(scrna_avgmat) <- cn
    #Save the result dataframe
    save(scrna_avgmat, file=paste(filePATH,"/",datatype,"-",fileName,"-averageExpression.Rda", sep=""))
    cat(date(),": scRNA Average expression finished\n")
    
    # #find zeroes
    # findnonZero <- apply(scrna_avgmat,1,function(x){ sum(x==0)})
    # findnonZero <- findnonZero[findnonZero < (ncol(scrna_avgmat)-6)]
    # scrna_avgmat <- scrna_avgmat[names(findnonZero),]
    #------------------------------------------------

    #------------------------------------------------
    # Select coding-genes only
    if(!is.null(coding_genes)) {
        goi <- intersect(row.names(scrna_avgmat), coding_genes)
        cat(date(),": Selecting protein coding genes only [",nrow(scrna_avgmat),"->", length(goi),"]\n")
        if(length(goi)>0) {
            scrna_avgmat <- scrna_avgmat[goi,]
        } else {
            cat(date(),": No overlap with coding genes provided\n")
            stop()
        }
    }
    #------------------------------------------------

    #------------------------------------------------
    #Keep genes with avgExpression > zero
    rowDF <- rowSums(scrna_avgmat)
    rowDF <- rowDF[rowDF > 0]
    mat <- scrna_avgmat[names(rowDF),]
    #Create annotation
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

    Overlap <- intersect(colnames(mat), row.names(ann))
    ann <- ann[Overlap,]
    mat <- mat[,Overlap]
    write.table(sort(unique(ann$group)), file=paste(filePATH,"/",fileName,"-group.txt", sep=""), row.names = F, col.names=F, quote=F)
    #------------------------------------------------

    #------------------------------------------------
    if(is.null(housekeeping_genes)) {
        housekeeping_genes <- c("GAPDH", "ACTB", "C1ORF43", "CHMP2A", "EMC7", "GPI", "PSMB2", "PSMB4", "RAB7A", "REEP5", "SNRPD3", "VCP", "VPS29")
    }
    
    #CV vs Mean
    if("intra-donor" %in% method) {
        cv_res <- cvCalcSC(mat=mat, ann=ann, meanThreshold=meanThreshold, cvThreshold=cvThreshold, housekeeping_genes=housekeeping_genes, filePATH=filePATH, fileName=fileName)
        cat(date(),": CV analysis done. Results in output directory.\n")

        #Find stable and variable features in longitudinal data
        var_gene <- VarFeatures(ann=ann, group_oi=group_oi, meanThreshold=meanThreshold, cvThreshold=cvThreshold, donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=topFeatures, filePATH=filePATH, fileName=fileName)
        stable_gene <- StableFeatures(ann=ann, group_oi=group_oi, meanThreshold=meanThreshold, cvThreshold=cvThreshold, donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=topFeatures, filePATH=filePATH, fileName=fileName)
        cat(date(),": Variable and stable features are identified. Now making UMAP.\n")
        
        #UMAP Plot
        dimUMAPPlot(rnaObj=dataObj, nPC=nPC, gene_oi=var_gene, groupName=avgGroup, plotname="variable", filePATH=filePATH, fileName=fileName)
        dimUMAPPlot(rnaObj=dataObj, nPC=nPC, gene_oi=stable_gene, groupName=avgGroup, plotname="stable", filePATH=filePATH, fileName=fileName)
        dimUMAPPlot(rnaObj=dataObj, nPC=nPC, gene_oi=unique(c(stable_gene,var_gene)), groupName=avgGroup, plotname="stable-variable", filePATH=filePATH, fileName=fileName)
        
        CV <- cv_res
        variable_gene <- var_gene
        stable_gene <- stable_gene
        cat(date(),": UMAP done.\n") 
    } else {
        CV <- NULL
        variable_gene <- NULL
        stable_gene <- NULL
    }
    #------------------------------------------------

    #------------------------------------------------
    #Variance decomposition
    if("inter-donor" %in% method) {
        lmem_res <- lmeVariance(ann=ann, mat=mat, features=c(features,"group"), meanThreshold=meanThreshold, fileName=fileName, filePATH=filePATH)
        cat(date(),": Variance decomposition analysis done. Results are in output directory.\n")
    } else {
        lmem_res <- NULL
    }
    #------------------------------------------------
    
    #------------------------------------------------
    #Calculate IQR (Outlier analysis)
    if("iqr" %in% method) {
        IQR <- iqrSC(ann=ann, mat=mat, fileName=fileName, filePATH=filePATH)
    } else {
        IQR <- NULL
    }
    #------------------------------------------------

    #------------------------------------------------
    #return result object
    long_res <- list(
        ann=ann,
        scrnaMetadata=metaData,
        scrna_avgdataframe=scrna_avgmat,
        CV=CV,
        variable_gene=variable_gene,
        stable_gene=stable_gene,
        varianceDecomposition=lmem_res,
        iqr=IQR$IQR_res,
        outlier=IQR$res
    )
    return(long_res)
    #------------------------------------------------

}
