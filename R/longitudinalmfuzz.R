#' longitudinalmfuzz Function
#'
#' This function allows you to identify gene/feature trajectory over
#' longitudinal points. The function uses mfuzz package (for more information
#' refer to https://www.bioconductor.org/packages/release/bioc/html/Mfuzz.html)
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation information
#' and expression data from Bulk or single cell data.
#' @param group_column User-defined group name like 'group','celltype'
#' @param timeColumn User-defined time column name like 'Time'
#' @param donorColumn User-defined donor/participant column name like 'PTID'
#' @param timeOrder (Optional) User-defined order of time variable like
#' ('D1','D2','D3')
#' @param baseline_timepoint (Optional) If baseline donors known (like 'PTID1')
#' @param featurelist (Optional) User-defined genes/features of interest
#' @param group_oi User-defined groups to consider for example from celltypes
#' select few
#' @param mfuzz_thres \code{mfuzz:thres} threshold for excluding genes
#' @param mfuzz_min.std \code{mfuzz:min.std} threshold for minimum standard
#' deviation
#' @param max_cluster Number of clusters to explore (Default 2^n)
#' @param delta \code{mfuzz:delta} threshold for minimum standard deviation
#' @param plotsize Size of plot width and height. Default 10 (in).
#' @param cl Number of clusters. Use nCores-1 to run parallel. Default 2
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory \emph{PATH} Default, current
#' directory
#' @keywords longitudinalmfuzz
#' @return longitudinal trajectory dataframe
#' @export
#' @examples
#' \dontrun{
#' longitudinalmfuzz(data_object=palmo_obj, group_column='group',
#' timeColumn='Time', donorColumn='PTID')
#' }

longitudinalmfuzz <- function(data_object, group_column = "group",
                              timeColumn = "Time", timeOrder = NULL,
                              donorColumn = "PTID", baseline_timepoint = NULL,
                              featurelist = NULL, group_oi = NULL,
                              mfuzz_thres = 0.25, mfuzz_min.std = 0,
                              max_cluster = NULL, delta = 0.5,
                              plotsize = 10, cl = 2,
                              fileName = NULL, filePATH = NULL) {

    message(date(), ": Identifying Longitudinal Trajectories")
    ## If filename or filepath null
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    if (is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }

    ## get the data
    ann <- data_object@curated$anndata
    mat <- data_object@curated$data
    check_data <- all.equal(row.names(ann), colnames(mat))
    if (check_data == FALSE) {
        stop(date(), ": Annotation of samples (rows) and datamatrix columns
             do not match")
    }

    ## By user defined group
    ann$bygroup <- ann[, group_column]
    ann <- ann[!is.na(ann$bygroup), ]
    ann$Sample_bygroup <- paste(ann$Sample, ann$bygroup, sep = ":")

    ## Define donor column
    ann$donorGroup <- ann[, donorColumn]
    ann <- ann[!is.na(ann$donorGroup), ]
    uniDonor <- sort(unique(ann$donorGroup))

    ## Input group of interest
    ann <- ann[ann$bygroup %in% group_oi, ]
    mat <- mat[, row.names(ann)]
    ann <- ann[order(ann$donorGroup, ann$bygroup), ]
    uniGroup <- sort(unique(ann$bygroup))
    uniSample_group <- unique(ann$Sample_bygroup)
    # print(uniGroup)

    ## By user defined time
    ann$bytime <- ann[, timeColumn]
    ann <- ann[!is.na(ann$bytime), ]
    uniTime <- sort(unique(na.omit(ann$bytime)))
    if (!is.null(timeOrder)) {
        uniTime <- intersect(timeOrder, uniTime)
    }
    message(date(), ": Time order->", uniTime, ".\n>>To change time order
            use parameter timeOrder=c('t1', 't2', 't3').")

    ## User defined genelist/featurelist
    if (is.null(featurelist)) {
        geneList <- as.character(row.names(mat))
    } else {
        geneList <- intersect(featurelist, row.names(mat))
    }

    ## Define maximum clusters
    if (is.null(max_cluster)) {
        max_cluster <- 2^(length(uniTime))
    } else {
        max_cluster <- as.numeric(max_cluster)
    }

    ## by group
    res_group <- lapply(uniGroup, function(uIG) {
        ## by donor
        res_donor <- lapply(uniDonor, function(uID) {
            message(date(), "Running ::", uID, "-", uIG)
            ann_sub <- ann[ann$bygroup %in% uIG & ann$donorGroup %in% uID, ]
            mat_sub <- mat[geneList, row.names(ann_sub)]
            uniTime_sub <- intersect(uniTime, unique(ann_sub$bytime))
            op <- pboptions(type = "timer")  # default
            res <- pblapply(geneList, cl = cl, function(gL) {
                df <- data.frame(gene = as.numeric(mat_sub[gL, ]),
                                 time = ann_sub$bytime, blank = 0)
                suppressMessages(df1 <- df %>% group_by(time, blank) %>%
                                     summarise(val = median(gene)) %>%
                                     data.frame() %>%
                                     column_to_rownames(var = "time"),
                                 classes = "message")
                df1 <- df1[uniTime_sub, ]
                return(df1$val)
            })
            pboptions(op)
            res <- data.frame(do.call(rbind, res))
            colnames(res) <- uniTime_sub
            row.names(res) <- geneList
            res <- data.frame(res, check.names=FALSE, stringsAsFactors=FALSE)

            ## Define maximum cluster
            nCluster <- 2^(length(uniTime_sub))
            if (nCluster > max_cluster) {
                nCluster <- max_cluster
            }

            ## Mfuzz based clustering set.seed(1234)
            expDF <- new("ExpressionSet", exprs = as.matrix(res))
            possible_clusters <- nCluster
            labels <- uniTime_sub
            expDF <- filter.NA(expDF, thres = mfuzz_thres)
            expDF <- filter.std(expDF, min.std = mfuzz_min.std, visu = FALSE)

            ## Check baseline timepoint
            if(is.null(baseline_timepoint)) {
                expDF <- standardise(expDF)
            } else {
                print(baseline_timepoint)
                if(length(intersect(colnames(res), baseline_timepoint)) == 1) {
                  expDF <- standardise2(expDF, timepoint = baseline_timepoint)
                } else {
                  stop(date(), ": Baseline timepoint is not found")
                }
            }

            ## Estimate m
            m1 <- mestimate(expDF)
            t.cl <- mfuzz(expDF, c = possible_clusters, m = m1)
            temp <- data.frame(t.cl$centers, check.names = FALSE,
                               stringsAsFactors = FALSE)
            temp1 <- apply(temp, 1, function(x) {
                diff <- lapply(2:length(x), function(y) {
                  y1 <- x[y] - x[y - 1]
                  y1 <- sign(y1) * ifelse(abs(y1) >= delta, 1, 0)
                  return(y1)
                })
                diff <- as.numeric(unlist(diff))
                y <- paste(paste(diff, collapse = ":"), sep = "")
                return(y)
            })
            temp <- data.frame(temp, trajectory_direction = temp1,
                               cluster = row.names(temp), check.names = FALSE,
                               stringsAsFactors = FALSE)

            # Cluster information
            temp$group <- uIG
            temp$donor <- uID
            write.csv(temp, file = paste(filePATH, "/longitudinal_",
                                         fileName, "_", uID, "_", uIG, ".csv",
                                         sep = ""), quote = FALSE)

            plot_panels <- ceiling(sqrt(nCluster))
            pdf(paste(filePATH, "/longitudinal_", fileName, "_", uID, "_", uIG,
                      ".pdf", sep = ""), width = plotsize, height = plotsize)
            mfuzz.plot(expDF, cl = t.cl, time.labels = labels,
                       mfrow = c(plot_panels, plot_panels), new.window = FALSE)
            dev.off()

            mRes <- data.frame(Genes = names(t.cl[["cluster"]]),
                               cluster = t.cl[["cluster"]],
                               stringsAsFactors = FALSE)
            write.csv(mRes, file = paste(filePATH, "/longitudinal_",
                        fileName, "_", uID, "_", uIG, "_geneclustering.csv",
                        sep = ""), quote = FALSE)
            mRes <- mRes[geneList, ]
            row.names(mRes) <- geneList
            colnames(mRes)[2] <- paste(uID, uIG, sep = ":")
            return(mRes)
        })
        res_donor <- do.call(cbind, res_donor)
        res_donor <- res_donor[, setdiff(colnames(res_donor), "Genes")]
        return(res_donor)
    })
    res_group <- do.call(cbind, res_group)
    data_object@result$mclust_longitudinal <- res_group
    return(data_object)
}
