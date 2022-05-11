#' StableFeatures Function
#'
#' This function allows user to identify stable genes in participants across
#' longitudinal timepoints using single cell expression data. The coefficient of
#' variation (CV) calculated using \code{cvCalcSC} function. Users can identify
#' cvThreshold in different datasets using housekeeping genes CV distribution.
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information and expression data from Bulk or single cell data.
#' @param group_oi Group of interest to focus on. Example among celltypes focus
#' on selected ones. Default is NULL.
#' @param cvThreshold Coefficient of variation threshold to select variable and
#' stable genes Default is 10 for single cell RNA \code{(100*SD/mean)}
#' @param donorThreshold Donor threshold number to be used, Default is number of
#' participants
#' @param groupThreshold Group label threshold number to be used, Default is
#' \code{(number of participants x group labels)/2}
#' @param housekeeping_genes Optional list of housekeeping genes to focus on.
#' Default is ACTB, GAPDH
#' @param topFeatures Number of features to be selected from each group, Default
#' is 25
#' @param fileName User defined filename
#' @param filePATH User-defined output directory path to load the CV result obtained
#' from cvCalcSC function
#' @return PALMO object with stable (stable_genes) features
#' @keywords StableFeatures
#' @export
#' @examples
#' \dontrun{
#' palmo_obj <- StableFeatures(data_object=palmo_obj, cvThreshold=10)
#' }

StableFeatures <- function(data_object, group_oi=NULL,
                           cvThreshold=NULL, donorThreshold=NULL,
                           housekeeping_genes=NULL,
                           groupThreshold=NULL, topFeatures=25,
                           filePATH=NULL, fileName=NULL) {

    message(date(),": Identifying Stable features\n")

    ## If filename or filepath null
    if(is.null(fileName)) {
        fileName <- "outputFile"
    }
    if(is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }

    #Load data
    if(!is.null(data_object@result$cv_meanthreshold)) {
        cv_res <- data_object@result$cv_meanthreshold
    } else {
        stop(date(),": Please run cvCalcSC function before StableFeatures function\n")
    }
    variable_gene <- data_object@result$variable_gene
    non_variable_gene <- data_object@result$non_variable_gene

    ## get the annotation data
    data_object@curated$anndata$Sample_group_i <- paste(data_object@curated$anndata$group, data_object@curated$anndata$PTID, sep=":")
    ann <- data_object@curated$anndata

    #group list
    if(is.null(group_oi)) {
        group_oi <- unique(as.character(ann$group))
    }

    #Select group of interest
    ann_sub <- ann[ann$group %in% group_oi,]
    if(nrow(ann_sub)<1) {
        stop(date(), ": Group of interest features do not match with annotation group eg.",unique(ann$group)[1:3],"\n")
    }
    Sample_group <- unique(ann_sub$Sample_group_i)
    Sample_group <- intersect(Sample_group, colnames(cv_res))

    if(is.null(donorThreshold)) {
        donorThreshold <- length(unique(ann_sub$PTID))
        message(date(),": Donor threshold defined ",donorThreshold,"\n")
    } else if(donorThreshold > length(unique(ann$PTID))) {
        donorThreshold <- length(unique(ann_sub$PTID))
        message(date(),": Donors were larger than unique donors. Donor threshold defined ",donorThreshold,"\n")
    }

    gThr <- round(length(unique(ann_sub$PTID)) * length(unique(ann_sub$group)) * 0.9)
    if(is.null(groupThreshold)) {
        groupThreshold <- round(length(unique(ann_sub$PTID)) * length(unique(ann_sub$group)) * 0.5)
        message(date(),": Groupwise threshold defined ",groupThreshold,"\n")
    } else if(groupThreshold > gThr) {
        groupThreshold <- round(length(unique(ann_sub$PTID)) * length(unique(ann_sub$group)) * 0.9)
        message(date(),": Number of groups were larger than unique groupsxsample. Groupwise threshold defined ",groupThreshold,"\n")
    }

    #Summary of non-variable genes
    temp <- data.frame(do.call(rbind, strsplit(as.character(non_variable_gene$donor), split = ":")), stringsAsFactors = FALSE)
    non_variable_gene$PTID <- temp$X2
    non_variable_gene$group <- temp$X1
    stable_genelist <- non_variable_gene[non_variable_gene$group %in% group_oi, ]
    stable_genelist <- data.frame(table(stable_genelist$gene))
    stable_genelist <- stable_genelist[order(stable_genelist$Freq, decreasing = TRUE),]
    stable_gene <- as.character(stable_genelist$Var1)
    #create Stable matrix
    stable_gene <- intersect(stable_gene, row.names(cv_res))
    stable_mat <- cv_res[stable_gene,]
    stable_mat[stable_mat > cvThreshold] <- NA
    save(stable_mat, file=paste(filePATH,"/",fileName,"-stableMatrix.Rda", sep=""))

    #Define the super-stable genes
    plot1 <- ggplot(stable_genelist, aes(x=Freq)) + geom_histogram(binwidth=1) +
        labs(title="Stable genes occurance from each sample")
    super_stable1 <- stable_genelist[stable_genelist$Freq >= groupThreshold,] #atleast in donor x group x4
    super_stable2 <- stable_genelist[stable_genelist$Freq >= donorThreshold & stable_genelist$Freq < groupThreshold,]
    #plot heatmap (super-stable)
    gn <- as.character(super_stable1$Var1[1:25])
    data_mat <- stable_mat[gn,Sample_group]
    rn <- data.frame(do.call(rbind, strsplit(colnames(data_mat), split = ":")), stringsAsFactors = FALSE)
    #set.seed(2020)
    ha_col <- HeatmapAnnotation(df=data.frame(group=rn$X1),
                                annotation_name_gp = gpar(fontsize = 6),
                                simple_anno_size = unit(0.3, "cm"))
    ht1 <- Heatmap(data.matrix(data_mat), cluster_rows =FALSE,
                   cluster_columns = FALSE,
                   column_split = factor(rn$X1, levels = group_oi),
                   na_col = "grey",
                   col = colorRamp2(c(0,cvThreshold,(cvThreshold+0.001),(cvThreshold+10)), c("blue","white","pink","red")),
                   row_names_max_width=unit(10, "cm"),
                   column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 6),
                   row_title = "Super stable 25",
                   column_title_gp = gpar(fontsize = 4),
                   top_annotation = ha_col,
                   heatmap_legend_param = list(title = "CV",heatmap_legend_side = "right") )
    pdf(paste(filePATH,"/",fileName,"-Super-stable-Top25.pdf", sep=""), width=10, height=3.5)
    print(plot1)
    print(ht1)
    dev.off()
    print(ht1)
    write.csv(data_mat, paste(filePATH,"/",fileName,"-Super-stable-Top25.csv", sep=""))

    #Housekeeping genes
    if(is.null(housekeeping_genes)) { housekeeping_genes <- c("ACTB", "GAPDH") }
    housekeeping_genes <- intersect(housekeeping_genes, row.names(cv_res))
    if(length(housekeeping_genes)>0) {
        data_mat_hg <- cv_res[housekeeping_genes,Sample_group]
        rn_hg <- data.frame(do.call(rbind, strsplit(colnames(data_mat_hg), split = ":")), stringsAsFactors = FALSE)
        #set.seed(2020)
        ha_col <- HeatmapAnnotation(df=data.frame(group=rn_hg$X1),
                                    annotation_name_gp = gpar(fontsize = 6),
                                    simple_anno_size = unit(0.3, "cm"))
        ht1_hg <- Heatmap(data.matrix(data_mat_hg), cluster_rows =FALSE,
                          cluster_columns = FALSE,
                          column_split = factor(rn_hg$X1, levels = group_oi),
                          na_col = "grey", col = colorRamp2(c(0,cvThreshold,(cvThreshold+0.001),(cvThreshold+10)), c("blue","white","pink","red")),
                          row_names_max_width=unit(10, "cm"),
                          column_names_gp = gpar(fontsize = 5),
                          row_names_gp = gpar(fontsize = 6),
                          row_title = "Stable (housekeeping) genes",
                          column_title_gp = gpar(fontsize = 4),
                          top_annotation = ha_col,
                          heatmap_legend_param = list(title = "CV",heatmap_legend_side = "right") )
        pdf(paste(filePATH,"/",fileName,"-housekeepingGenes.pdf", sep=""), width=10, height=3.5)
        print(ht1_hg)
        dev.off()
    }

    #Define the stable genes
    stable_list <- as.character(super_stable2$Var1)
    dfx <- c()
    for(i in 1:length(group_oi)) {
        groupName <- group_oi[i]
        df <- non_variable_gene[non_variable_gene$gene %in% stable_list & non_variable_gene$group %in% groupName,]
        df1 <- data.frame(table(df$gene, df$group))
        df1 <- df1[df1$Freq >= donorThreshold & order(df1$Freq, decreasing = TRUE),]
        df <- df[df$gene %in% df1$Var1,]
        dfx <- rbind(dfx, df[1:topFeatures,])
    }
    dfx <- dfx[!is.na(dfx$mean),]
    stable_gene <- unique(dfx$gene)
    write.csv(dfx, file=paste(filePATH,"/",fileName,"-stable-genelist.csv", sep=""))
    #plot heatmap
    data_mat <- stable_mat[stable_gene,Sample_group]
    ht2 <- Heatmap(data.matrix(data_mat), cluster_rows =FALSE,
                   cluster_columns = FALSE,
                   column_split = factor(rn$X1, levels = group_oi),
                   na_col = "grey",
                   col = colorRamp2(c(0,cvThreshold,(cvThreshold+0.001),(cvThreshold+10)), c("white","blue","pink","red")),
                   row_names_max_width=unit(10, "cm"),
                   column_names_gp = gpar(fontsize = 5),
                   row_names_gp = gpar(fontsize = 5),
                   row_title = paste("Stable genes:", length(stable_gene)),
                   column_title_gp = gpar(fontsize = 4),
                   top_annotation = ha_col,
                   heatmap_legend_param = list(title = "CV",heatmap_legend_side = "right") )

    pdf(paste(filePATH,"/",fileName,"-stable-Features.pdf", sep=""), width=10, height=10)
    print(ht2)
    dev.off()
    print(ht2)
    write.csv(data_mat, paste(filePATH,"/",fileName,"-stable-Features-",topFeatures,".csv", sep=""))

    ## Add the result
    data_object@result$stable_genes <- dfx
    message(date(),": Check output directory for results\n")

    return(data_object)

}
