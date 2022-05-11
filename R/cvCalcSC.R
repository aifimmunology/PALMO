#' cvCalcSC Function
#'
#' This function allows to calculate Intra-donor variations in single cell data
#' over longitudinal timepoints. The coefficient of variation (CV=SD/mean) is
#' calculated in average expression data in same donor/participant and
#' corresponding user-defined group (like celltype, cluster) across longitudinal
#' timepoints.
#' @param data_object Input \emph{PALMO} S4 object. Contains annotation table and
#' expression matrix or data frame. Rows represent gene/proteins column
#' represents participant samples (same as annotation table Sample column)
#' @param meanThreshold Average expression threshold to filter lowly expressed
#' genes Default is 0.1 (log2 scale)
#' @param cvThreshold Coefficient of variation threshold to select variable and
#' stable genes. Default is 5 for bulk data. Users can use 10-20 for single cell
#' average expression data.
#' @param housekeeping_genes Optional, vector of housekeeping genes. Default is
#' c('ACTB', 'GAPDH')
#' @param cl Number of clusters. Use nCores-1 to run parallel. Default 2
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory \emph{PATH} Default, current
#' directory
#' @return PALMO object with CV list
#' @keywords cvCalcSC
#' @export
#' @examples
#' \dontrun{
#' palmo_obj=cvCalcSC(data_object=palmo_obj, meanThreshold=0.1, cvThreshold=5)
#' }

cvCalcSC <- function(data_object, meanThreshold = NULL, cvThreshold = NULL, housekeeping_genes = NULL,
    cl = 2, fileName = NULL, filePATH = NULL) {

    message(date(), ": Performing Coefficient of variance analysis\n")

    ## If filename or filepath null
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    if (is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }

    ## Assign housekeeping_genes
    if (is.null(housekeeping_genes)) {
        housekeeping_genes <- c("ACTB", "GAPDH")
        data_object@housekeeping_genes <- housekeeping_genes
    }

    ## meanThrehold
    if (is.null(meanThreshold)) {
        meanThreshold <- 0
        message(date(), ": Using mean threshold >= 0\n")
    }
    data_object@meanThreshold <- meanThreshold

    ## cvThreshold
    if (is.null(cvThreshold)) {
        cvThreshold <- 10
        message(date(), ": Using cv threshold 10\n")
    }
    data_object@cvThreshold <- cvThreshold

    ## get the data
    ann <- data_object@curated$anndata
    mat <- data_object@curated$data
    check_data <- all.equal(row.names(ann), colnames(mat))
    if (check_data == FALSE) {
        stop(date(), ": Annotation of samples (rows) and datamatrix columns do not match")
    }

    ## Define group
    ann$group_donor <- paste(ann$group, ann$PTID, sep = ":")
    sample_freq <- data.frame(table(ann$group_donor))
    sample_freq <- sample_freq[sample_freq$Freq > 1, ]
    sample_freq <- as.character(sample_freq$Var1)
    if (length(sample_freq) > 0) {
        ann <- ann[ann$group_donor %in% sample_freq, ]
        mat <- mat[, row.names(ann)]
    } else {
        stop(date(), ": Not enough group samples to perform CV analysis\n")
    }

    # Calculate CV vs Mean for all genes per celltype
    unigene <- row.names(mat)
    uniDonor <- sort(unique(ann$PTID))
    uniSamplegroup <- as.character(unique(ann$group_donor))

    # All genes CV calculations
    message(date(), ": Performing CV calculations\n")
    op <- pboptions(type = "timer")  # default
    res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
        # print(uS)
        ann_df <- ann[ann$group_donor %in% uS, ]
        df <- mat[unigene, ann_df$Sample_group]
        df <- data.frame(df,
                         na = apply(df, 1, function(x) {sum(x != 0) }),
                         mean = rowMeans(df, na.rm = TRUE),
                         sd = apply(df, 1, sd, na.rm = TRUE),
                         var = apply(df, 1, var, na.rm = TRUE),
                         stringsAsFactors = FALSE)
        df$cv <- 100 * df$sd/df$mean
        return(df$cv)
    })
    pboptions(op)
    cv_res <- do.call(cbind, res)
    row.names(cv_res) <- unigene
    colnames(cv_res) <- uniSamplegroup
    cv_res <- data.frame(cv_res, check.names = FALSE, stringsAsFactors = FALSE)
    rm(res)
    # save result
    save(cv_res, file = paste(filePATH, "/", fileName, "-CV-allgenes-raw.Rda", sep = ""))
    cv_all <- cv_res

    # Genes with minimum mean threshold CV calculations
    message(date(), ": Checking Mean Threshold \n")
    op <- pboptions(type = "timer")  # default
    res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
        # print(uS)
        ann_df <- ann[ann$group_donor %in% uS, ]
        df <- mat[unigene, ann_df$Sample_group]
        df <- data.frame(df,
                         na = apply(df, 1, function(x) {sum(x != 0)}),
                         mean = rowMeans(df,na.rm = TRUE),
                         sd = apply(df, 1, sd, na.rm = TRUE),
                         var = apply(df, 1, var, na.rm = TRUE),
                         stringsAsFactors = FALSE)
        # CV <- 100*df$sd/df$mean return(CV)
        df$cv <- 100 * df$sd/df$mean
        df$cv <- ifelse(df$mean >= meanThreshold, df$cv, NA)
        return(df$cv)
    })
    pboptions(op)
    cv_res <- do.call(cbind, res)
    row.names(cv_res) <- unigene
    colnames(cv_res) <- uniSamplegroup
    cv_res <- data.frame(cv_res, check.names = FALSE, stringsAsFactors = FALSE)
    rm(res)
    # save result
    save(cv_res, file = paste(filePATH, "/", fileName, "-CV-allgenes.Rda", sep = ""))

    # Variable genes
    message(date(), ": Performing Variable gene CV analysis\n")
    op <- pboptions(type = "timer")  # default
    res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
        ann_df <- ann[ann$group_donor %in% uS, ]
        df <- mat[unigene, ann_df$Sample_group]
        df <- data.frame(df,
                         nonZero = apply(df, 1, function(x) {sum(x != 0)}),
                         mean = rowMeans(df, na.rm = TRUE),
                         sd = apply(df, 1, sd, na.rm = TRUE),
                         var = apply(df, 1, var, na.rm = TRUE),
                         stringsAsFactors = FALSE)
        df$CV <- 100 * df$sd/df$mean
        # the CV becomes very high for data with 0
        df <- df[df$mean >= meanThreshold, ]  #minimum expression >2^0.1=1
        dp2a <- df[df$mean >= meanThreshold & df$CV > cvThreshold, c("mean", "sd", "var", "CV")]
        dp2a <- dp2a[order(dp2a$CV, dp2a$mean, decreasing = TRUE), ]
        # Find variable genes
        if (nrow(dp2a) >= 1) {
            variable_gene <- data.frame(donor = uS, gene = row.names(dp2a),
                                        dp2a, stringsAsFactors = FALSE)
            return(variable_gene)
        }
    })
    pboptions(op)
    variable_gene <- do.call(rbind, res)
    rm(res)
    # save result
    save(variable_gene, file = paste(filePATH, "/", fileName, "-CV-Variablegene.Rda", sep = ""))

    # Stable genes
    message(date(), ": Performing Stable gene CV analysis\n")
    op <- pboptions(type = "timer")  # default
    res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
        ann_df <- ann[ann$group_donor %in% uS, ]
        df <- mat[unigene, ann_df$Sample_group]
        df <- data.frame(df,
                         nonZero = apply(df, 1, function(x) { sum(x != 0) }),
                         mean = rowMeans(df, na.rm = TRUE),
                         sd = apply(df, 1, sd, na.rm = TRUE),
                         var = apply(df, 1, var, na.rm = TRUE),
                         stringsAsFactors = FALSE)
        df$CV <- 100 * df$sd/df$mean
        # the CV becomes very high for data with 0
        df <- df[df$mean >= meanThreshold, ]  #minimum expression >2^0.1=1
        dp2b <- df[df$mean >= meanThreshold & df$CV <= cvThreshold, c("mean", "sd", "var", "CV")]
        dp2b <- dp2b[order(-dp2b$mean, dp2b$CV, decreasing = FALSE), ]
        # Find stable genes
        if (nrow(dp2b) >= 1) {
            non_variable_gene <- data.frame(donor = uS, gene = row.names(dp2b),
                                            dp2b, stringsAsFactors = FALSE)
            return(non_variable_gene)
        }
    })
    pboptions(op)
    non_variable_gene <- do.call(rbind, res)
    rm(res)
    # save result
    save(non_variable_gene, file = paste(filePATH, "/",fileName,
                                         "-CV-nonVariablegene.Rda", sep = ""))

    # Housekeeping genes data
    message(date(), ": Checking Housekeeping-genes CV\n")
    op <- pboptions(type = "timer")  # default
    res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
        ann_df <- ann[ann$group_donor %in% uS, ]
        df <- mat[unigene, ann_df$Sample_group]
        df <- data.frame(df,
                         nonZero = apply(df, 1, function(x) { sum(x != 0) }),
                         mean = rowMeans(df, na.rm = TRUE),
                         sd = apply(df, 1, sd, na.rm = TRUE),
                         var = apply(df, 1, var, na.rm = TRUE),
                         stringsAsFactors = FALSE)
        df$CV <- 100 * df$sd/df$mean
        # the CV becomes very high for data with 0
        df <- df[df$mean >= meanThreshold, ]  #minimum expression >2^0.1=1
        temp <- df[housekeeping_genes, c("mean", "sd", "var", "CV")]
        thr <- round(1 + max(df$CV, na.rm = TRUE))
        # Housekeeping gene profile
        temp <- data.frame(gene = row.names(temp), donor = uS, temp, stringsAsFactors = FALSE)
        return(temp)
    })
    pboptions(op)
    hg_res <- do.call(rbind, res)
    rm(res)
    # Housekeeping genes
    hg_res <- hg_res[!is.na(hg_res$mean), ]

    # Variable genes
    rn <- data.frame(do.call(rbind, strsplit(variable_gene$donor, split = ":")), stringsAsFactors = FALSE)
    variable_gene$sample <- rn$X2
    variable_gene$group <- rn$X1
    plot2 <- ggplot(variable_gene, aes(x = mean, y = CV)) +
        geom_point() +facet_wrap(~group)

    # stable/non-variable genes
    rn <- data.frame(do.call(rbind, strsplit(non_variable_gene$donor, split = ":")), stringsAsFactors = FALSE)
    non_variable_gene$sample <- rn$X2
    non_variable_gene$group <- rn$X1
    plot3 <- ggplot(non_variable_gene, aes(x = mean, y = CV)) +
        geom_point() + facet_wrap(~group)

    if (nrow(hg_res) > 0) {
        # house-keeping genes
        rn <- data.frame(do.call(rbind, strsplit(hg_res$donor, split = ":")), stringsAsFactors = FALSE)
        hg_res$sample <- rn$X2
        hg_res$group <- rn$X1

        plot4 <- ggplot(variable_gene, aes(x = mean, y = CV)) +
            geom_point() +
            geom_point(data = non_variable_gene, aes(x = mean, y = CV), color = "red") +
            geom_point(data = hg_res, aes(x = mean, y = CV), color = "blue") +
            geom_hline(yintercept = cvThreshold) +
            facet_wrap(~group, scales = "free_y")
    } else {
        plot4 <- ggplot(variable_gene, aes(x = mean, y = CV)) +
            geom_point() +
            geom_point(data = non_variable_gene, aes(x = mean, y = CV), color = "red") +
            geom_hline(yintercept = cvThreshold) +
            facet_wrap(~group, scales = "free_y")
    }

    message(date(), ": Saving CV plots in output directory\n")
    png(paste(filePATH, "/", fileName, "-CV-distribution.png", sep = ""),
        width = 10, height = 10, res = 200, units = "in")
    print(plot4)
    dev.off()

    ## Add CV result
    data_object@result$cv_all <- cv_all
    data_object@result$cv_meanthreshold <- cv_res  #minimum mean threshold
    data_object@result$variable_gene <- variable_gene
    data_object@result$non_variable_gene <- non_variable_gene

    message(date(), ": Done. Please check output directory for results.\n")
    return(data_object)
}
