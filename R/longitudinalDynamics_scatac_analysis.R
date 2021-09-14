#' A longitudinalDynamics_scatac_analysis Function
#'
#' This function allows you to perform analysis of longitudinal dataset from
#' single cell atac based genescore matrix.
#' @param metadata Annotation table must consist column Sample, 
#' Participant sample name; PTID, participant; Time, longitudinal time frame
#' @param genescore Genescore Expression matrix or data frame derived from scATAC.
#' Rows represents genes and column represents participant samples (same as 
#' annotation table Sample column)
#' @param datatype Data input can be bulk or singlecell
#' @param omics User defined name like RNA, ATAC, Proteomics, FLOW
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
#' @keywords longitudinalDynamics_scatac_analysis
#' @export
#' @examples
#' ##Single cell ATAC data
#' #longitudinalDynamics(metadata=ann, atacObj=atacObj, 
#' #datatype="singlecell", omics="atac",fileName="scatac",
#' #features=c("PTID", "Time", "group"),
#' #meanThreshold=0.1, cvThreshold=10,
#' #housekeeping_genes=NULL, group_oi=NULL,
#' #nPC=30, donorThreshold=4, groupThreshold=40, topFeatures=25,
#' #method=c("intra-donor","inter-donor","iqr"),
#' #outputDirectory="output")

longitudinalDynamics_scatac_analysis <- function(metadata=NULL, genescore=NULL, datatype=NULL, omics=NULL,
                            fileName=NULL, filePATH=NULL, features=NULL, meanThreshold=0.1, cvThreshold=10,
                            coding_genes=NULL, avgGroup=NULL,
                            housekeeping_genes=NULL, group_oi=NULL, nPC=NULL,
                            method=NULL, 
                            donorThreshold=NULL, groupThreshold=NULL, topFeatures=25) {

    #------------------------------------------------
    cat(date(), ": Running scATAC longitudinal analysis\n")
    scatac_avgmat <- data.frame(genescore, check.names=F, stringsAsFactors = F)
    #Keep genes with avgExpression > zero
    rowDF <- rowSums(scatac_avgmat)
    rowDF <- rowDF[rowDF > 0]
    scatac_avgmat <- scatac_avgmat[names(rowDF),]
    #---------------------------------------
    
    #------------------------------------------------
    # Select coding-genes only
    if(!is.null(coding_genes)) {
        if(coding_genes != FALSE) {
            cat(date(),": Selecting protein coding genes only\n")
            goi <- intersect(row.names(scatac_avgmat), coding_genes)
            scatac_avgmat <- scatac_avgmat[goi,]
        }
    }
    #------------------------------------------------
        
    #------------------------------------------------
    #Variance decomposition
    mat <- scatac_avgmat
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
        
        # #UMAP Plot
        # dimUMAPPlot(ann=ann, countMat=mat, nPC=nPC, gene_oi=var_gene, groupName=avgGroup, plotname="variable", filePATH=filePATH, fileName=fileName)
        # dimUMAPPlot(ann=ann, countMat=mat, nPC=nPC, gene_oi=stable_gene, groupName=avgGroup, plotname="stable", filePATH=filePATH, fileName=fileName)
        # dimUMAPPlot(ann=ann, countMat=mat, nPC=nPC, gene_oi=unique(stable_gene,var_gene), groupName=avgGroup, plotname="stable-variable", filePATH=filePATH, fileName=fileName)
        # cat(date(),": UMAP done.\n")

        CV <- cv_res
        variable_gene <- var_gene
        stable_gene <- stable_gene
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
        IQR_res <- iqrSC(mat=mat, ann=ann, filePATH=filePATH, fileName=fileName)
    } else {
        IQR_res <- NULL
    }
    #------------------------------------------------
    
    #------------------------------------------------
    #return result object
    long_res <- list(
        ann=ann,
        scatacMetadata=metaData,
        scatac_avgdataframe=scatac_avgmat,
        CV=CV,
        variable_gene=variable_gene,
        stable_gene=stable_gene,
        varianceDecomposition=lmem_res,
        iqr=IQR_res
    )
    return(long_res)
    #------------------------------------------------
}
