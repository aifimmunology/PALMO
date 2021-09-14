#' A longitudinalDynamics Function
#'
#' This function allows you to perform analysis of longitudinal dataset.
#' It requires longitudinal data matrix/data frame and annotation file.
#' @param metadata Annotation table must consist column Sample, 
#' Participant sample name; PTID, participant; Time, longitudinal time frame
#' @param datamatrix Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples (same as annotation table Sample column)
#' @param scrnaObj Single cell RNA Seurat object, if datatype is single cell RNA
#' @param atacObj Single cell ATAC genematrix or data frame, if datatype is 
#' single cell ATAC
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
#' @param outputDirectory User-defined output directory Default, output
#' @keywords longitudinalDynamics
#' @export
#' @examples
#' ##Bulk data
#' #load("data/data_Annotation.Rda")
#' #load("data/Olink_NPX_log2_Protein.Rda")
#' #longitudinalDynamics(metadata=ann, datamatrix=data, datatype="bulk", omics="Protein",
#' #fileName="olink", features=c("PTID", "Time", "Sex"), meanThreshold=1, cvThreshold=5,
#' #outputDirectory="output")
#' 
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

longitudinalDynamics <- function(metadata=NULL, datamatrix=NULL,
    scrnaObj=NULL, atacObj=NULL, datatype=NULL, omics=NULL,
    fileName, features=NULL, meanThreshold=1, cvThreshold=5,
    coding_genes=NULL, avgGroup=NULL,
    housekeeping_genes=NULL, group_oi=NULL, nPC=30,
    donorThreshold=NULL, groupThreshold=NULL, topFeatures=25,
    method=NULL,
    outputDirectory=NULL) {
    
    #------------------------------------------------
    #load packages
    #loadLibrary()
    #------------------------------------------------

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
    #setwd(file.path(mainDir, subDir))
    #------------------------------------------------
    
    #------------------------------------------------
    if(is.null(metadata)) {
        cat(date(),": Please provide the annotation dataframe")
        stop()
    } else {
        metadata <- data.frame(metadata, check.names = F, stringsAsFactors = F)
    }
    #------------------------------------------------
    
    #------------------------------------------------
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
    nID <- unique(metadata$PTID)
    nSample <- table(metadata$PTID)
    
    #Check single cell or bulk data
    if(datatype =="singlecell" & !is.null(scrnaObj)) {
        overlap <- intersect(metadata$Sample, scrnaObj@meta.data$Sample)
    } else if(datatype =="singlecell" & !is.null(atacObj)) {
        temp <- data.frame(do.call(rbind, strsplit(colnames(atacObj), split = ":")), stringsAsFactors = F)
        cn <- data.frame(id=colnames(atacObj), Sample=temp$X1, group=temp$X2, check.names=F, stringsAsFactors = F)
        cn <- cn[cn$Sample %in% metadata$Sample,]
        overlap <- as.character(unique(cn$Sample))
        atac_overlap <- cn$id
    } else if(datatype =="bulk" & !is.null(datamatrix)) {
        overlap <- intersect(colnames(datamatrix), metadata$Sample)
    } else {
        cat(date(),": Please provide dataframe for bulk data and scRNA Seurat Object for scrna
         with Sample column matched with metadta")
        stop()
    }
    
    if(length(overlap)==0) {
        cat(date(),": No overlap between the columnnames and annotation samples")
        stop()
    } else if( (length(nID)<2) & (sum(as.numeric(nSample<2))>0) ) {
        cat(date(),": Please include multiple donors (min=2) and samples (min=2)")
        stop()
    } else {
        cat(date(),": Number of samples matched:",length(overlap),"\n")
    }
    
    if(datatype =="singlecell" & !is.null(scrnaObj)) {
        metadata <- metadata[overlap,]
        scrnaObj <- subset(x = scrnaObj, subset = Sample %in% overlap)
    } else if(datatype =="singlecell" & !is.null(atacObj)) {
        metadata <- metadata[overlap,]
        atacObj <- data.frame(atacObj, check.names = F, stringsAsFactors = F)
        atacObj <- atacObj[,atac_overlap]
    } else if(datatype =="bulk" & !is.null(datamatrix)) {
        metadata <- metadata[overlap,]
        datamatrix <- data.frame(datamatrix, check.names = F, stringsAsFactors = F)
        datamatrix <- datamatrix[,overlap]
    } else {
        cat(date(),": Please provide the expression datamatrix/dataframe")
        stop()
    }
    #------------------------------------------------

    #------------------------------------------------
    #Analysis to be performed
    if(is.null(method)) {
        method <- c("inter-donor", "intra-donor", "iqr")
    } else if(!method %in% c("inter-donor", "intra-donor", "iqr")) {
        method <- c("inter-donor", "intra-donor", "iqr")
    }
    #------------------------------------------------

    #------------------------------------------------
    #Is dataset is single cell or Bulk data
    if(datatype == "bulk") {
        cat(date(),": Bulk data Analysis started\n")
        long_res <- longitudinalDynamics_bulk_analysis(metadata=metadata, datamatrix=datamatrix, datatype=datatype, omics=omics, 
            fileName=fileName, filePATH=filePATH, features=features, meanThreshold=meanThreshold, cvThreshold=cvThreshold, method=method) 
    } else if(datatype == "singlecell") {
        if(omics == "rna") {
        cat(date(),": scRNA Analysis started\n")
        long_res <- longitudinalDynamics_scrna_analysis(metadata=metadata, dataObj=scrnaObj, datatype=datatype, omics=omics,
            fileName=fileName, filePATH=filePATH, features=features, meanThreshold=meanThreshold, cvThreshold=cvThreshold,
            coding_genes=coding_genes, avgGroup=avgGroup,
            housekeeping_genes=housekeeping_genes, group_oi=group_oi, nPC=nPC,
            donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=topFeatures, method=method) 
        } else if(omics == "atac") {
        cat(date(),": scATAC Analysis started\n")
        long_res <- longitudinalDynamics_scatac_analysis(metadata=metadata, genescore=atacObj, datatype=datatype, omics=omics,
            fileName=fileName, filePATH=filePATH, features=features, meanThreshold=meanThreshold, cvThreshold=cvThreshold,
            coding_genes=coding_genes, avgGroup=avgGroup,
            housekeeping_genes=housekeeping_genes, group_oi=group_oi, nPC=nPC,
            donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=topFeatures, method=method) 
        } else {
        cat(date(),"ERROR: Please provide datatype and omics name\n")
        stop()
        }
    } else {
        stop("ERROR: Please provide datatype and omics name") 
    }
    #------------------------------------------------
    
    #------------------------------------------------
    cat(date(),": Done\n")
    return(long_res)
    #------------------------------------------------
    
}
