#' sclongitudinalDEG Function
#'
#' This function allows ser to calculate differential expressed genes in the
#' direction of given time points (if timepoints>3 otherwise DEGs between two
#' timepoints). A hurdle model was fit to each participant independently in
#' order to identify participant-specific longitudinal transcriptomic changes.
#' Genes that were expressed in at least 10% of cells per participant were
#' considered for this analysis. The models were fit on the input normalized
#' data, modeling the timepoints as a continuous variable within each cell type
#' and adjusting for the batch only if any timepoints from the same participant
#' were run across multiple batches.
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information and expression data from Bulk or single cell data.
#' @param scassay Single cell assay from scRNA seurat object (Default "RNA")
#' @param group_column Column of interest such as "celltype" to analyze DEGs in
#' participant over time
#' @param group_oi Features of interest such as  specific celltypes
#' c("CD4_Naive", "CD4_TEM")
#' @param mincellsexpressed Average expression threshold to filter lowly
#' expressed genes/features Default is 0.1
#' @param removelnc Remove lincRNAs, mitochondrial and ribosomal genes from
#' analysis incldes (^RP|^MT-|^LINC|orf) (TRUE/FALSE). Default is TRUE
#' @param adjfac Factors to be adjusted for such as batch, sex
#' @param baseline Donors (PTID) to be considered as baseline. Deafult NULL
#' @param CDR_column (Optional) cellular detection rate column name
#' @param addCDR (Optional) Add CDR while performing differential analysis.
#' Default is FALSE
#' @param plotWidth User-defined plot width, Default 10 in
#' @param plotHeight User-defined plot height, Default 10 in
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory \emph{PATH} Default, current
#' directory
#' @keywords sclongitudinalDEG
#' @import Seurat, MAST, tidyverse, pheatmap
#' @export
#' @examples
#' \dontrun{
#' palmo_obj <- sclongitudinalDEG(ann=metadata, dataObj=pbmc, scassay="RNA",
#' group_column="celltype")
#' }

sclongitudinalDEG <- function(data_object, scassay="RNA",
                              group_column, group_oi=NULL,
                              mincellsexpressed = 0.1, removelnc = "TRUE",
                              adjfac = NULL,
                              baseline = NULL, addCDR=FALSE, CDR_column=NULL,
                              plotWidth=10, plotHeight=10,
                              fileName=NULL, filePATH=NULL) {

    ## Filename
    if(is.null(fileName)) {
        fileName <- "outputFile"
    }
    ## Filepath
    if(is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }

    ## Define donor column
    donorIDcol <- "PTID"
    ## loading single cell data
    scobj <- data_object@curated$SeuratObj
    ann <- data_object@curated$anndata

    ## Add CDR value (cellular detection rate)
    if(addCDR == TRUE) {
        scobj@meta.data$CDR <- scale(scobj@meta.data[,CDR_column])
    }

    ## Any baseline
    if(!is.null(baseline)) {
        ann[ann$PTID %in% baseline,]$PTID <- "Baseline"
        baseline <- TRUE
    } else {
        baseline <- FALSE
    }

    ## fetch single cell counts
    sccounts <- scobj[[scassay]]@data

    ## fetch single cell meta data
    scmetadata <- scobj@meta.data %>%
        rownames_to_column(var = "barcodes")
    filtercolumn <- c("Sample", setdiff(colnames(scmetadata), colnames(ann)))
    if(length(filtercolumn)>1) {
        ## Merge clinical metadata with scobject metadata
        scmetadata <- merge(scmetadata[,filtercolumn], ann, by = "Sample", all.x = TRUE) %>%
            column_to_rownames(var = "barcodes")
        scmetadata <- scmetadata[colnames(sccounts), , drop = FALSE]
        #table(rownames(scmetadata) == colnames(sccounts))
    }

    ## fetch PTIDs and group/celltypes
    ptids <- sort(as.character(unique(scmetadata[[donorIDcol]])))
    ptids <- ptids[!is.na(ptids)]
    group_types <- sort(as.character(unique(scmetadata[[group_column]])))
    group_types <- group_types[!is.na(group_types)]

    ##
    if(!is.null(group_oi)) {
        group_types <- intersect(group_oi, group_types)
    }

    ## initiating loop for per donor analysis
    lmLS <- lapply(ptids, function(ptid) {
      message(date(),": Fitting a ZLM model for donorID: ", ptid)
      ## Consider baseline
      if(baseline==TRUE) {
        message(date(),": Performing baseline analysis")
        scmetadata_sub <- subset(scmetadata, scmetadata[[donorIDcol]] == ptid |
                              scmetadata[[donorIDcol]] == "Baseline")
        scmetadata_sub[[donorIDcol]] <- ifelse(scmetadata_sub[[donorIDcol]] %in% "Baseline", ptid,
                              scmetadata_sub[[donorIDcol]])
      } else {
        # select single cell barcodes that belong to the ptid
        scmetadata_sub <- subset(scmetadata, scmetadata[[donorIDcol]] == ptid)
      }

      sccounts_sub <- sccounts[, rownames(scmetadata_sub), drop = FALSE]
      # initiating loop for per cell type within a donor analysis
      lapply(group_types, function(ct) {
        message(paste("Fitting a ZLM model for donorID: ", ptid,
                ", within group: ", ct, sep = ""))

        # subset counts on selected cell type barcodes
        selcells <- subset(scmetadata_sub, scmetadata_sub[[group_column]] == ct)
        selcells <- rownames(selcells)
        sccounts_sub_ct <- sccounts_sub[, selcells, drop = FALSE]

        # filter on genes expressed in atleast 10% of cells
        min_expr = mincellsexpressed
        selgenes1 <- data.frame(num_cells_expressed = rowSums(as.matrix(sccounts_sub_ct) > min_expr)) %>%
                    rownames_to_column(var = "Gene") %>%
                    filter(num_cells_expressed >= min_expr*length(selcells))

        # remove long non-coding rnas (conditional)
        if(removelnc == "TRUE") {
            rmgenes <- selgenes1$Gene[grep("^RP|^MT-|^LINC|orf", selgenes1$Gene)]
            selgenes2 <- selgenes1 %>% filter(!Gene %in% rmgenes) %>% .$Gene
            sccounts_sub_ct <- sccounts_sub_ct[selgenes2, ]
        } else {
            sccounts_sub_ct <- sccounts_sub_ct[selgenes1$Gene, ]
        }

        # make fdata for single cell assay (SCA) object for input to MAST
        fdat <- data.frame(rownames(x = sccounts_sub_ct))
        colnames(x = fdat)[1] <- "primerid"
        rownames(x = fdat) <- fdat[, 1]

        # changing the column name of the adjusting factor (if provided) to "adjfactor"
        scmetadata_sub <- scmetadata_sub[colnames(sccounts_sub_ct), , drop = FALSE]

        # define heatmap color palatte
        colorPalette <- c("purple", "black", "yellow")
        colorPalette <- colorRampPalette(colors = colorPalette)(100)

        # fetch timepoints of PTID
        tps1 <- unique(scmetadata_sub$Time)

        # if number of timepoints >= 3, treat as continuous
        suppressMessages(
        if(length(tps1) >= 3) {
          #message(paste("modeling time continuous(t=",length(tps1),
          #"): ",ptid,"-", ct, sep = ""))

          # make sca object with time as numeric continuous
          cdat <- scmetadata_sub
          cdat$Time <- as.numeric(gsub("[[:alpha:]]", "",cdat$Time))

          # make SCA object
          sca <- FromMatrix(exprsArray = as.matrix(sccounts_sub_ct),
                              cData = cdat, fData = fdat)

          # Fit hurdle model
          # add adjusting factor and CDR if provided
          if(is.null(adjfac) & addCDR == FALSE) {
            form <- as.formula(paste("~Time", sep=""))
          } else if(is.null(adjfac) & addCDR == TRUE) {
            form <- as.formula(paste("~Time + CDR", sep=""))
          } else if(!is.null(adjfac) & addCDR == FALSE) {
            form <- paste(adjfac, collapse=" + ")
            form <- as.formula(paste("~Time + ", form, sep=""))
          } else if(!is.null(adjfac) & addCDR == TRUE) {
            form <- paste(adjfac, collapse=" + ")
            form <- as.formula(paste("~Time + CDR + ", form, sep=""))
          }

          zlmCond <- zlm(formula = form, sca = sca, method = "bayesglm",
                         ebayes = TRUE, parallel = TRUE)

          # perform Likelihood Ratio Test (LRT)
          lrt = "Time"
          summaryGroup <- MAST::summary(object = zlmCond, doLRT = lrt)

          # get summary table
          summaryDt <- summaryGroup$datatable

          # extract p-values and coefficients
          pval <- summaryDt[summaryDt$component == "H" & summaryDt$contrast == lrt,
                            c("primerid", "contrast", "Pr(>Chisq)")]
          coefDF <- summaryDt[summaryDt$component == "logFC" & summaryDt$contrast == lrt,
                            c("primerid", "coef")]
          lmdf <- merge(pval, coefDF, by = "primerid")
          colnames(lmdf)[3] <- "nomP"
          lmdf$adjP <- p.adjust(lmdf$nomP, method = "BH")
          dir = ifelse(lmdf$coef > 0, "increasing over time",
                       "decreasing over time")
          lmdf <- data.frame(lmdf, donorID = ptid, group_type = ct,
                                       dir = dir, stringsAsFactors = FALSE)
          lmdf <- lmdf[order(lmdf$coef, decreasing = TRUE),]

          outFile1 <- paste(filePATH,"/",fileName,"_",ptid,"_", ct,
                                      "_DEGs.txt", sep = "")
          write.table(lmdf, file = outFile1, sep = "\t", quote = FALSE,
                                row.names = FALSE)

          # take top 50 genes by direction
          topG <- lmdf %>% group_by(dir) %>%
            top_n(min(50, length(primerid)), abs(coef)) %>%
            as.data.frame() %>% .$primerid

          # plot heatmap of significant genes: taking mean gene expression per timepoint
          sigDF <- as.data.frame(sccounts_sub_ct[topG, ])
          mCounts <- sigDF %>% rownames_to_column(var = "Gene") %>%
                        gather(barcodes, value, -Gene) %>%
                        mutate(Sample = scmetadata_sub$Sample[match(barcodes,
                                        table = rownames(scmetadata_sub))]) %>%
                        group_by(Sample, Gene) %>%
                        dplyr::summarize(m = mean(value)) %>%
                        as.data.frame() %>%
                        spread(Sample, m) %>%
                        column_to_rownames(var = "Gene")

          # scale matrix and define breaks
          mCounts2 <- t(scale(t(mCounts))) %>% as.data.frame()
          limit <- range(mCounts2)
          limit <- c(ceiling(limit[1]), floor(limit[2]))
          limit <- min(abs(limit))
          if(limit == 0) { limit = 1 }

          # column annotation of heatmap
          colannot1 <- data.frame(Sample = colnames(mCounts2)) %>%
                        mutate(donorID = ptid,
                               group_type = ct,
                               Time = scmetadata_sub$Time[match(Sample,
                                             table = scmetadata_sub$Sample)],
                               Time = gsub("[[:alpha:]]", "",Time),
                               #Time = gsub("\\D", "", Time),
                               Time = as.numeric(Time)) %>%
                        column_to_rownames(var = "Sample")
                    colannot1 <- colannot1[colnames(mCounts2), , drop = FALSE]

          # order cols by Time
          orderCol1 <- colannot1 %>% rownames_to_column() %>%
                        arrange(Time) %>% .$rowname

          # order rows by gene directions
          orderRow <- lmdf %>%
                        filter(primerid %in% topG) %>%
                        as.data.frame() %>%
                        select(primerid, coef, dir) %>%
                        arrange(coef)
          orderRow1 <- orderRow$primerid

          # gaps row
          gaprow <- as.numeric(table(orderRow$dir)[1])

          # annotation colors
          ann_colors <- list(Time = c("aliceblue", "yellow", "orange", "red"))

          # plot
          outFile2 <- gsub(".txt", "_avgExpHeatmap.pdf", outFile1)
          pdf(file = outFile2, width = plotWidth, height = plotHeight)
          pheatmap(mat = mCounts2[orderRow1, orderCol1],
                             breaks = c(min(mCounts2),
                                        seq(from = -1*limit,
                                            to = limit,
                                            length.out = 99),
                                        max(mCounts2)),
                             color = colorPalette,
                             cellwidth = 10,
                             cellheight = 3,
                             cluster_cols = FALSE,
                             cluster_rows = FALSE,
                             show_rownames = TRUE,
                             show_colnames = TRUE,
                             gaps_row = gaprow,
                             annotation = colannot1,
                             annotation_colors = ann_colors,
                             fontsize_row = 3,
                             fontsize_col = 3,
                             fontsize_number = 20,
                             border_color = NA)
          dev.off()

          ## Plotting single cells heatmap
          ## downsampling to a random 50 or min cells per sample:
          ## for plotting purpose only
          #set.seed(seed = 123)
          colannot2 <- cdat %>% rownames_to_column(var = "barcodes") %>%
              select(barcodes, Sample, Time, !!donorIDcol, !!group_column) %>%
              group_by(Sample) %>%
              sample_n(min(50, length(barcodes))) %>%
              as.data.frame() %>%
              select(-Sample) %>%
              column_to_rownames(var = "barcodes")

          sccounts_sub_ct_v2 <- sccounts_sub_ct[topG, rownames(colannot2)]
          # order col by Time
          orderCol2 <- colannot2 %>% rownames_to_column() %>%
            arrange(Time) %>% .$rowname

          # scale matrix and define breaks
          sccounts_sub_ct_v2 <- sccounts_sub_ct_v2[apply(sccounts_sub_ct_v2, 1, function(x)  length(unique(x)) != 1), ]
          sccounts_sub_ct_v2 <- t(scale(t(as.matrix(sccounts_sub_ct_v2)))) %>% as.data.frame()

          # set limit for color scale on heatmap
          limit <- range(sccounts_sub_ct_v2)
          limit <- c(ceiling(limit[1]), floor(limit[2]))
          limit <- min(abs(limit))

          # gaps row and order row
          orderRow <- lmdf %>%
              filter(primerid %in%
                         rownames(sccounts_sub_ct_v2)) %>%
              as.data.frame() %>%
              select(primerid, coef, dir) %>%
              arrange(coef)
          orderRow2 <- orderRow$primerid

          sccounts_sub_ct_v2 <- sccounts_sub_ct_v2[orderRow2, , drop = FALSE]
          gaprow <- as.numeric(table(orderRow$dir)[1])

          # plot
          outFile3 <- gsub(".txt", "_scHeatmap.pdf", outFile1)
          pdf(file = outFile3, width = plotWidth, height = plotHeight)
          pheatmap(mat = sccounts_sub_ct_v2[orderRow2, orderCol2],
                   breaks = c(min(sccounts_sub_ct_v2),
                              seq(from = -1*limit, to = limit, length.out = 99),
                              max(sccounts_sub_ct_v2)),
                   color = colorPalette, cellwidth = 2, cellheight = 3,
                   cluster_cols = FALSE, cluster_rows = FALSE,
                   show_rownames = TRUE, show_colnames = FALSE,
                   annotation = colannot2, annotation_colors = ann_colors,
                   gaps_row = gaprow, fontsize_col = 15, fontsize_row = 3,
                   fontsize_number = 20, border_color = NA)
          dev.off()
          }, classes = "message")

          # if number of timepoints = 2, treat as discrete
        suppressMessages(
        if(length(tps1) == 2) {
            message(paste("modeling time discrete(t=",length(tps1),"):",
                  ptid,"-", ct, sep = ""))

            # make sca object with time as discrete
            cdat <- scmetadata_sub
            cdat$Time <- as.factor(cdat$Time)

            # set reference as first timepoint
            cdat$Time <- relevel(cdat$Time,
                        ref = as.character(sort(unique(cdat$Time))[1]))

            # make SCA object
            sca <- FromMatrix(exprsArray = as.matrix(sccounts_sub_ct),
                      cData = cdat, fData = fdat)

            # Fit hurdle model
            # add adjusting factor if provided
            if(is.null(adjfac) & addCDR == FALSE) {
                form <- as.formula(paste("~Time", sep=""))
            } else if(is.null(adjfac) & addCDR == TRUE) {
                form <- as.formula(paste("~Time + CDR", sep=""))
            } else if(!is.null(adjfac) & addCDR == FALSE) {
                form <- paste(adjfac, collapse=" + ")
                form <- as.formula(paste("~Time + ", form, sep=""))
            } else if(!is.null(adjfac) & addCDR == TRUE) {
                form <- paste(adjfac, collapse=" + ")
                form <- as.formula(paste("~Time + CDR + ", form, sep=""))
            }

            zlmCond <- zlm(formula = form,
                           sca = sca,
                           method = "bayesglm",
                           ebayes = TRUE,
                           parallel = TRUE)

            # perform LRT
            lrt <- paste("Time",
                         as.character(sort(unique(cdat$Time))[2]),
                         sep = "")
            summaryGroup <- MAST::summary(object = zlmCond, doLRT = lrt)

            # get summary table
            summaryDt <- summaryGroup$datatable

            # extract p-values and coefficients
            pval <- summaryDt[summaryDt$component == "H"
                              & summaryDt$contrast == lrt,
                              c("primerid", "contrast", "Pr(>Chisq)")]
            coefDF <- summaryDt[summaryDt$component == "logFC"
                        & summaryDt$contrast == lrt,
                        c("primerid", "coef")]
            lmdf <- merge(pval, coefDF, by = "primerid")
            colnames(lmdf)[3] <- "nomP"
            lmdf$adjP <- p.adjust(lmdf$nomP, method = "BH")
            fgtp <- gsub("Time", "", lrt)
            dir <- ifelse(lmdf$coef > 0, paste("upregulated at ", fgtp, sep = ""),
                          paste("downregulated at ", fgtp, sep = ""))
            lmdf <- data.frame(lmdf, donorID = ptid, group_type = ct,
                               dir = dir, stringsAsFactors = FALSE)
            lmdf <- lmdf[order(lmdf$coef, decreasing = TRUE),]

            outFile1 <- paste(filePATH,"/",fileName,"_",ptid,"_", ct,
                              "_DEGs.txt", sep = "")
            write.table(lmdf, file = outFile1, sep = "\t",
                        quote = FALSE, row.names = FALSE)

            # take top 50 genes by direction
            topG <- lmdf %>%
                group_by(dir) %>%
                top_n(min(50, length(primerid)), abs(coef)) %>%
                as.data.frame() %>% .$primerid

            ## Plotting single cells heatmap
            ## downsampling to a random 50 or min cells per sample: for
            ## plotting purpose only
            #set.seed(seed = 456)
            colannot3 <- cdat %>%
                rownames_to_column(var = "barcodes") %>%
                select(barcodes, Sample, Time, !!donorIDcol, !!group_column) %>%
                group_by(Sample) %>%
                sample_n(min(50, length(barcodes))) %>%
                as.data.frame() %>%
                select(-Sample) %>%
                column_to_rownames(var = "barcodes")
            sccounts_sub_ct_v2 <- sccounts_sub_ct[topG, rownames(colannot3)]
            sccounts_sub_ct_v2 <- as.matrix(sccounts_sub_ct_v2)

            # order col by Time
            orderCol3 <- colannot3 %>%
                rownames_to_column() %>%
                arrange(Time) %>% .$rowname

            # scale matrix and define breaks
            sccounts_sub_ct_v2 <- sccounts_sub_ct_v2[apply(sccounts_sub_ct_v2,
                            1, function(x)  length(unique(x)) != 1), ]
            sccounts_sub_ct_v2 <- t(scale(t(sccounts_sub_ct_v2))) %>%
                as.data.frame()

            ## set limit for color scale on heatmap
            limit <- range(sccounts_sub_ct_v2)
            limit <- c(ceiling(limit[1]), floor(limit[2]))
            limit <- min(abs(limit))

            ## gaps row and order row
            orderRow <- lmdf %>%
                filter(primerid %in%
                           rownames(sccounts_sub_ct_v2)) %>%
                as.data.frame() %>%
                select(primerid, coef, dir) %>%
                arrange(coef)
            orderRow2 <- orderRow$primerid
            sccounts_sub_ct_v2 <- sccounts_sub_ct_v2[orderRow2, , drop = FALSE]
            gaprow <- as.numeric(table(orderRow$dir)[1])

            # plot
            outFile4 <- gsub(".txt", "_scHeatmap.pdf", outFile1)
            pdf(file = outFile4, width = plotWidth, height = plotHeight)
            pheatmap(mat = sccounts_sub_ct_v2[orderRow2, orderCol3],
                     breaks = c(min(sccounts_sub_ct_v2),
                        seq(from = -1*limit,
                            to = limit,
                            length.out = 99),
                        max(sccounts_sub_ct_v2)),
                     color = colorPalette,
                     cellwidth = 2,
                     cellheight = 3,
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     show_rownames = TRUE,
                     show_colnames = FALSE,
                     annotation = colannot3,
                     gaps_row = gaprow,
                     fontsize_col = 15,
                     fontsize_row = 3,
                     fontsize_number = 20,
                     border_color = NA)
            dev.off()
          } , classes = "message")

          if(length(tps1) >= 2) {
            return(value = lmdf)
          }


        }) #End celltype/group_oi

      }) #End ptid

    allres <- do.call(rbind, lapply(lmLS, function(x) do.call(rbind, x)))
    data_object@result$degs <- allres
    return(data_object)
}
