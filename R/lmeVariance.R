#' lmeVariance Function
#'
#' This function allows you to calculate inter-donor variation between
#' participants over longitudinal timepoints. It uses linear mixed model to
#' calculate variance contribution from each given feature list.
#' @param data_object Input \emph{PALMO} S4 object. It contains annotation
#' information and expression data from Bulk or single cell data.
#' @param featureSet Variance analysis carried out for the feature set provided such
#' as c('PTID', 'Time', 'Sex')
#' @param fixed_effect_var Fixed effect variables. In linear mixed model
#' fixed_effect_var included as fixed effect variables and variance contribution
#' obtained by adding them as random variables
#' @param meanThreshold Average expression threshold to filter lowly expressed
#' genes/features Default is 0
#' @param selectedFeatures User-defined gene/feature list
#' @param NA_to_zero Convert NAs to zero. Default FALSE
#' @param cl Number of clusters. Use nCores-1 to run parallel. Default 2
#' @param fileName User-defined file name, Default outputFile
#' @param filePATH User-defined output directory \emph{PATH} Default, current
#' directory
#' @return PALMO object with variance lmem_res dataframe
#' @keywords lmeVariance
#' @export
#' @examples
#' \dontrun{
#' palmo_obj=lmeVariance(data_object=palmo_obj,
#' featureSet=c('PTID','Time','Sex'))
#' }

lmeVariance <- function(data_object, featureSet, fixed_effect_var = NULL,
                        meanThreshold = NULL, selectedFeatures = NULL,
                        NA_to_zero = FALSE, cl = 2, fileName = NULL,
                        filePATH = NULL) {

    message(date(), ": Performing variance decomposition\n")
    if (is.null(fileName)) {
        fileName <- "outputFile"
    }
    if (is.null(filePATH)) {
        filePATH <- data_object@filePATH
    }

    ## meanThrehold
    if (is.null(meanThreshold)) {
        meanThreshold <- 0
        message(date(), ": Using mean threshold >= 0\n")
    }
    data_object@meanThreshold <- meanThreshold

    ## get the data
    ann <- data_object@curated$anndata
    mat <- data_object@curated$data
    check_data <- all.equal(row.names(ann), colnames(mat))
    if (check_data == FALSE) {
        stop(date(), ": Annotation of samples (rows) and datamatrix columns do not match")
    }

    ## If features selected before hand
    if (!is.null(selectedFeatures)) {
        mat <- mat[selectedFeatures, ]
    }

    ## featureset
    if (!is.null(featureSet)) {
        check_featureSet <- intersect(featureSet, colnames(ann))
        if (length(check_featureSet) != length(featureSet)) {
            stop(date(), ": Input Featureset not found")
        }
    } else {
        stop(date(), ": Please provide appropriate featureset")
    }
    data_object@featureSet <- featureSet

    ## Define formula
    form <- paste(paste("(1|", featureSet, ")", sep = ""), collapse = " + ")

    ## check fixed effect variables
    if (!is.null(fixed_effect_var)) {
        check_fixedeffect <- intersect(fixed_effect_var, colnames(ann))
        if (length(check_fixedeffect) > 0) {
            featureSet = c(featureSet, check_fixedeffect)
            form <- paste(paste("(1|", featureSet, ")", sep = ""), collapse = " + ")
            form <- paste(paste(check_fixedeffect, sep = "", collapse = " + "), " + ", form,
                sep = "", collapse = " + ")
        } else {
            stop(date(), ": Fixed effect variable (fixed_effect_var) not found")
        }
    }
    data_object@result$var_formula <- form

    ## make NAs zero if defined (avoid group error and for exploration)
    if (NA_to_zero == TRUE) {
        mat[is.na(mat)] <- 0
    }

    ## Define featureList
    featureList <- c(featureSet, "Residual")
    rowN <- row.names(mat)
    op <- pboptions(type = "timer")  # default
    suppressMessages(lmem_res <- pblapply(1:length(rowN), cl = cl, function(gn) {
        geneName <- rowN[gn]
        # print(geneName)
        df <- data.frame(exp = as.numeric(mat[geneName, ]), ann, stringsAsFactors = FALSE)

        ## remove features with not enough values
        gene_group <- lapply(1:length(featureSet), function(i) {
            df1 <- df[!is.na(df$exp) & !is.na(df[, featureSet[i]]), ]
            res <- data.frame(table(df1[, featureSet[i]]))
            res <- res[res$Freq > 1, ]
            return(nrow(res))
        })
        gene_group <- as.numeric(gene_group)
        ## Check any of attribute has only one group
        check_gene_group <- sum(gene_group == 1)
        if (check_gene_group == 0) {
            ## Define formula form <- paste(paste('(1|', featureSet,')', sep=''), collapse=' + ')
            form1 <- as.formula(paste("exp ~ ", form, sep = ""))

            # linear mixed effect model
            lmem <- lmer(formula = form1, data = df)
            lmem_re <- as.data.frame(VarCorr(lmem))
            row.names(lmem_re) <- lmem_re$grp
            lmem_re <- lmem_re[featureList, ]
            # lmem_re$vcov/sum(lmem_re$vcov)
            fix_effect <- fixef(lmem)  #get fixed effect
            lmem_re$CV <- lmem_re$sdcor/fix_effect  ##Calculate CV
            return(c(geneName, mean(df$exp, na.rm = TRUE),
                     median(df$exp, na.rm = TRUE), sd(df$exp, na.rm = TRUE),
                     max(df$exp, na.rm = TRUE),
                     (lmem_re$vcov)/sum(lmem_re$vcov)))
        }
    }), classes = "message")
    pboptions(op)
    lmem_res <- do.call(rbind, lmem_res)
    lmem_res <- data.frame(lmem_res, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(lmem_res) <- c("Gene", "mean", "median", "sd", "max", featureList)
    lmem_res$mean <- as.numeric(lmem_res$mean)
    lmem_res <- lmem_res[!is.na(lmem_res$mean), ]
    row.names(lmem_res) <- lmem_res$Gene
    # Some data converted into character hence convert into numeric
    temp <- apply(lmem_res[, -1], 1, function(x) {
        as.numeric(x)
    })
    row.names(temp) <- colnames(lmem_res)[-1]
    lmem_res <- data.frame(Gene = colnames(temp), t(temp), check.names = FALSE,
                           stringsAsFactors = FALSE)
    lmem_res <- lmem_res[order(lmem_res[, featureList[1]], decreasing = TRUE), ]
    write.csv(lmem_res, file = paste(filePATH, "/", fileName, "-Variance.csv", sep = ""))

    ## Overall gene mean result
    plot1 <- ggplot(lmem_res, aes(x = mean)) + geom_histogram() + labs(title = "Mean Expression")
    plot2 <- ggplot(lmem_res, aes(x = log10(mean + 1))) +
        geom_histogram() + labs(title = "Mean Expression (log10)")
    pdf(paste(filePATH, "/", fileName, "-Meanplot.pdf", sep = ""), width = 5, height = 5)
    print(plot1)
    print(plot2)
    dev.off()

    ## Create data matrix
    res <- lmem_res[lmem_res$max > meanThreshold, featureList]
    df <- melt(data.matrix(res))
    df$value <- df$value * 100
    df$feature <- paste(df$Var1, df$Var2, sep = "_")

    sigFeatures <- c()
    for (i in 1:length(featureSet)) {
        dfx <- res[order(res[, featureSet[i]], decreasing = TRUE), ]
        sigFeatures <- c(sigFeatures, paste(row.names(dfx[1:5, ]), featureSet[i], sep = "_"))
    }
    df1 <- df[df$feature %in% sigFeatures, ]
    plot1 <- ggplot(df, aes(x = Var2, y = value, fill = Var2)) +
        geom_violin(scale = "width") +
        geom_boxplot(width = 0.1, fill = "white") +
        # ggpubr::stat_compare_means(label.y = 140) + # Add global p-value
        geom_text_repel(data = df1, aes(x = Var2, y = value),
                        label = df1$Var1, size = 2, segment.size = 0.1,
                        segment.alpha = 0.9, max.overlaps = 20, color = "black") +
        labs(x = "FeatureList", y = "Variance Explained (%)") +
        theme_classic()
    ## Plot
    pdf(paste(filePATH, "/", fileName, "-VarianceExplained-Boxplot.pdf", sep = ""),
        width = 5, height = 5)
    print(plot1)
    dev.off()
    print(plot1)

    ## Return object
    lmem_res[, featureList] <- 100 * lmem_res[, featureList]  #in percentage
    data_object@result$variance_decomposition <- lmem_res

    return(data_object)
}
