#' A genePlot Function
#'
#' This function allows you to perform UMAP visualization of gene of interest list.
#' @param ann Annotation table. Table must consist column Sample (Participant 
#' sample name), PTID (Participant), Time (longitudinal time points), group,
#' name of the group, group_donor (combined string using group:Sample)
#' @param data Average Expression matrix or data frame. Rows represents gene/proteins
#' column represents participant samples with group (optional).
#' @param geneName User-defined gene name
#' @param groupName User-defined group name column from annotation table
#' @keywords genePlot
#' @export
#' @examples
#' #plot <- genePlot(ann=annotation, data=ExpressionData, geneName="FOLR3", groupName="Time")

genePlot <- function(ann, data, geneName, groupName=NULL) {
  
  
  if(is.null(groupName)) {
    ann$group <- ann$Time
  } else {
    ann$group <- ann[,groupName]
  }
  
  #Get overlap
  Overlap <- intersect(row.names(ann), colnames(data))
  ann <- ann[Overlap,]
  data <- data[,Overlap]
  uniSample <- sort(unique(ann$PTID))
  
  #Plot gene/genes
  df <- data.frame(exp=as.numeric(data[geneName,]), ann, stringsAsFactors = F)
  df$PTID <- factor(df$PTID, levels = uniSample)
  plot1 <- ggplot(df, aes(x=PTID, y=exp)) +
      geom_boxplot() + geom_point(aes(color=Time)) +
      facet_grid(~group, scales="free") +
      labs(x="", y="Expression", title=geneName) +
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1))
  print(plot1)
  
  plot2 <- ggplot(df, aes(x=group, y=exp, color=PTID)) +
      geom_boxplot() +
      labs(x="", y="Expression", title=geneName) +
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1))
  print(plot2)

  #Participant
  plot3 <- ggboxplot(df, x = "PTID", y = "exp", add.params = list(shape="Time"), color="group", 
                  add = c("jitter", "boxplot"), 
                  ylab = "Expression", xlab = "PTID", 
                  title = geneName, legend = "right", outlier.shape = NA) + 
          scale_shape_manual(values = 0:length(unique(ann$Time)))
  print(plot3)
  
  plot4 <- ggline(df, x = "PTID", y = "exp", add.params = list(shape="Time",color="group"),
                  add = c("mean", "jitter"), 
                  ylab = "Expression", xlab = "PTID", 
                  title = geneName, legend = "right", outlier.shape = NA) + 
          scale_shape_manual(values = 0:length(unique(ann$Time)))
  print(plot4)
  
  plot5 <- ggplot(df, aes(x=PTID, y=exp, fill=Time)) +
    geom_boxplot(outlier.size = 0.1) +
    facet_wrap(~group) +
    labs(x="", y="Expression", title=geneName) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),
          legend.position = "right")
  print(plot5)
  
  # #Time-wise
  #plot4 <- ggline(df, x = "Time", y = "exp", add.params = list(color="PTID"), add = c("mean_se", "jitter", "boxplot"), ylab = "NPX", xlab = "PTID", title = geneName, legend = "right", outlier.shape = NA) + scale_shape_manual(values = 0:10)
  #print(plot4)
  # df$Time <- factor(df$Time, levels = unique(df$Time))
  # plot3 <- ggplot(df, aes(x=Time, y=exp, color=PTID, group=PTID)) +
  #       geom_point() +
  #       geom_line() +
  #       labs(x="", y="Expression", title=geneName) +
  #       theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1))
  # print(plot3)
      
  return(list(df=df, plot1=plot1,plot2=plot2, plot3=plot3, 
              plot4=plot4, plot5=plot5))
  
}

