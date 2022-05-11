#' palmo class
#'
#' This function creates \emph{PALMO} class object. All the raw data and results from
#' PALMO are stored in this object.
#' @field raw list, contains user entered annotation and expression dataframe or
#' object
#' @field curated list, contains curated input data
#' @field result list, output from \emph{PALMO} stored in result list
#' @field nDonors numeric, number of donors in the input data
#' @field rownames character, row names of the expression data
#' @field colnames character, column names of the expression data
#' @field housekeeping_genes character, user-defined housekeeping genes listed
#' @field datatype character, datatype used like bulk or singlecell
#' @field omics character, omics such as RNA, scRNA, scATAC
#' @field featureSet character, parameters used for variance analysis
#' @field meanThreshold numeric, Average expression threshold
#' @field cvThreshold numeric, CV threshold
#' @field median_cvThreshold numeric, median of CV threshold (from inter-donor)
#' @field groupName character, group defined by user like celltype, cluster
#' @field group_oi character, selected types from user-defined group list
#' @field donorThreshold numeric, minimum donors to explore
#' @field groupThreshold numeric, minimum group types to explore
#' @field topFeatures numeric, number of top features to retrieve
#' @field donor_sep character, donor and group separator such as ':'
#' @field cor_method character, correelation method 'pearson', 'spearman'
#' @field clusterBy character, cluster by a group (celltype or cluster)
#' @field z_cutoff numeric, z-cutoff value for outlier analysis
#' @field filePATH character, PATH of outout directory
#' @exportClass palmo
#' @return PALMO S4 class

palmo <- setClass(Class = "palmo", slots = c(raw = "list", curated = "list", result = "list", 
    nDonors = "numeric", rownames = "character", colnames = "character", housekeeping_genes = "character", 
    datatype = "character", omics = "character", featureSet = "character", meanThreshold = "numeric", 
    cvThreshold = "numeric", median_cvThreshold = "numeric", groupName = "character", group_oi = "character", 
    donorThreshold = "numeric", groupThreshold = "numeric", topFeatures = "numeric", donor_sep = "character", 
    cor_method = "character", clusterBy = "character", z_cutoff = "numeric", filePATH = "character"))
