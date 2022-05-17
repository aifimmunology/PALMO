#' outlierDetectP Function
#'
#' This function allows to identify significant abnormal event identified
#' from outlier analysis.
#' @param outlier_events Identified outlier events
#' @param z_cutoff |Z| cutoff threshold to find potential outliers (Eg.
#' z_cutoff= 2, equals to Mean/SD 2)
#' @param nGenes Number of background genes/features
#' @param group_by Column name to use for groupwise outlier analysis. Default is
#' PTID (donor or participant id).
#' @param alternative alternative hypothesis, must be one of "one.sided"
#' or "two.sided" (default)
#' @return PALMO object with outlier event p value dataframe
#' @keywords outlierDetectP
#' @return PALMO object with outlier event p value
#' @examples
#' \dontrun{
#' outlierDetectP(outlier_events=outlier_res, z_cutoff=2, nGenes=1043)
#' }

outlierDetectP <- function(outlier_events, z_cutoff = 2, nGenes,
                           group_by = "PTID", alternative="two.sided") {

    message(date(), ": Identifying significant abnormal events")
    outlier_res <- outlier_events

    res_list <- list() #define result list

    ## Calculate p values
    if(alternative == "one.sided") {
      # one-sided
      rate <- 1 - pnorm(z_cutoff)
      # events with z > z_cutoff
      upcounts <- outlier_res[outlier_res$z >= z_cutoff, ]
      if (nrow(upcounts) > 1) {
          uptable <- table(upcounts[, group_by], upcounts[, "Time"])
          df <- data.frame(uptable, stringsAsFactors = FALSE)
          colnames(df) <- c("Sample", "Time", "Freq")
          df$id <- paste(df$Sample, df$Time, sep = "")
          upevents <- as.integer(df$Freq)
          df$pvals <- p_value_for_event(events=upevents,tries=nGenes,rate=rate)
          df$signP <- -log10(df$pvals)
          plot1a <- ggplot(df, aes(x = Freq, y = signP)) + geom_point() +
            labs(x = "# Features (>z cutoff)", y = "-log10(pvalue)",
                 title = "Up (one-sided)") +
            ggrepel::geom_text_repel(data = df, aes(x = Freq, y = signP,
                              label = id), size = 2, segment.size = 0.1,
                            segment.alpha = 0.9, max.overlaps = 20) +
            theme_bw()

          plot1b <- ggplot(df, aes(x = id, y = Freq)) +
            geom_bar(stat = "identity") +
            labs(x="", y="# Features (>z cutoff)", title="Up (one-sided)") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5,
                                             vjust = 0.5))

          plot1c <- ggplot(df, aes(x = id, y = signP)) +
            geom_bar(stat = "identity") +
            labs(x="", y="-log10(pvalue)", title="Up (one-sided)") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5,
                                             vjust = 0.5))

          res_list$upcounts_p <- df
          print(plot_grid(plot1a, plot1b, plot1c, ncol = 1, align = "hv"))
          rm(df)
      } else {
        res_list$upcounts_p <- NA
      }

      # events with z < -z_cutoff
      dncounts <- outlier_res[outlier_res$z <= -z_cutoff, ]
      if (nrow(dncounts) > 1) {
          dntable <- table(dncounts[, group_by], dncounts[, "Time"])
          df <- data.frame(dntable, stringsAsFactors = FALSE)
          colnames(df) <- c("Sample", "Time", "Freq")
          df$id <- paste(df$Sample, df$Time, sep = "")
          dnevents <- as.integer(df$Freq)
          df$pvals <- p_value_for_event(events=dnevents,tries=nGenes,rate=rate)
          df$signP <- -log10(df$pvals)
          plot2a <- ggplot(df, aes(x = Freq, y = signP)) + geom_point() +
            labs(x = "# Features (>z cutoff)", y = "-log10(pvalue)",
                 title = "Down (one-sided)") +
            ggrepel::geom_text_repel(data = df, aes(x = Freq, y = signP,
                        label = id), size = 2, segment.size = 0.1,
                        segment.alpha = 0.9, max.overlaps = 20) +
            theme_bw()

          plot2b <- ggplot(df, aes(x = id, y = Freq)) +
            geom_bar(stat = "identity") +
            labs(x="", y="# Features (>z cutoff)", title="Down (one-sided)") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5,
                                             vjust = 0.5))

          plot2c <- ggplot(df, aes(x = id, y = signP)) +
            geom_bar(stat = "identity") +
            labs(x="", y="-log10(pvalue)", title="Down (one-sided)") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5,
                                             vjust = 0.5))
          res_list$dncounts_p <- df
          print(plot_grid(plot2a, plot2b, plot2c, ncol = 1, align = "hv"))
          rm(df)
      } else {
        res_list$dncounts_p <- NA
      }
      #one-sided end

    } else if(alternative == "two.sided") {

      # two-sided
      rate <- 2 * (1 - pnorm(z_cutoff))
      counts <- outlier_res[abs(outlier_res$z) >= z_cutoff, ]
      if(nrow(counts) > 1) {
        countstable <- table(counts[, group_by], counts[, "Time"])
        df <- data.frame(countstable, stringsAsFactors = FALSE)
        colnames(df) <- c("Sample", "Time", "Freq")
        df$id <- paste(df$Sample, df$Time, sep = "")
        events <- as.integer(df$Freq)
        df$pvals <- p_value_for_event(events=events, tries=nGenes, rate=rate)
        df$signP <- -log10(df$pvals)
        plot3a <- ggplot(df, aes(x = Freq, y = signP)) + geom_point() +
          labs(x = "# Features (>z cutoff)", y = "-log10(pvalue)") +
          ggrepel::geom_text_repel(data = df, aes(x = Freq, y = signP,
                      label = id), size = 2, segment.size = 0.1,
                      segment.alpha = 0.9, max.overlaps = 20) +
          theme_bw()

        plot3b <- ggplot(df, aes(x = id, y = Freq)) +
          geom_bar(stat = "identity") +
          labs(x="", y="# Features (>z cutoff)") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 0.5,
                                           vjust = 0.5))

        plot3c <- ggplot(df, aes(x = id, y = signP)) +
          geom_bar(stat = "identity") +
          labs(x="", y="-log10(pvalue)") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 0.5,
                                           vjust = 0.5))

        res_list$counts_p <- df
        rm(df)
        print(plot_grid(plot3a, plot3b, plot3c, ncol = 1, align = "hv"))
      } else {
        res_list$counts_p <- NA
      }
      #two-sided end
    }

    return(res_list)
}

#' p_value_for_event Function
#'
#' This function allows to calculate p value for identified outlier significant
#' abnormal events
#' @param events Identified outlier events
#' @param tries Number of background genes/features
#' @param rate probability distribution
#' @return outlier event p value
#' @keywords p_value_for_event
#' @examples
#' \dontrun{
#' p_value_for_event(events, tries, rate)
#' }

# identify significant abnormal event
p_value_for_event <- function(events, tries, rate) {
    pvalues <- events %>%
      lapply(function(x) binom.test(x, tries, rate,
                                    alternative = "greater")$p.value) %>%
      as.double() %>% p.adjust(method = "BH")

    return(pvalues)
}
