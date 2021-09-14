README for longitudinalDynamics
===============

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Install package and load library](#library)
* [Quick Usage (Longitudinal data)](#usage)
    * [Plasma proteome](#olink)
    * [Flow cell-proportion](#flow)
    * [Blood count data](#cbc)
    * [Single cell RNA data](#scrna)
    * [Single cell ATAC data](#scatac)
* [Tutorials](#example-main)
    * [Tutorial-1: Plasma proteome](#example1)
    * [Tutorial-2: scRNA](#example2)
    * [Tutorial-3: scATAC](#example3)
    * [Tutorial-4: CNP data](#example4)
* [Authors](#authors)
* [Acknowledgements](#acknowledgements)
* [License](#license)

* * *

## <a name="introduction"></a> Introduction
LongitudinalDynamics `(longitudinalDynamics)` is an implementation of the longitudinal data analysis platform in R. It allows to identify intra-donor, inter-donor variations over longitudinal time points. The analysis can be done on bulk expression dataset without known celltype information or single cell with celltype information or known groups. It allows to identify stable and variable features in given donor and each celltype (or user defined group). The outlier in intra-donor can be performed to identify perturbed or change in sample currespnding to donor/participant. The stable and variable gene signature helps to identify gene signatures in dataset and perturbation specific featureset.

<br><br> ![](vignettes/img1.png) <br><br>
Fig.1 General workflow and analysis schema of **LongitudinalDynamics**. It can work with longitudinal data obtained from bulk such as clinical, bulk RNAseq, proteomic or single cell dataset from scRNAseq, scATACseq and scADT.


## <a name="library"></a> Install package and load library

To install library, simply run
   
    install.packages("longitudinalDynamics_0.1.0.tar.gz", repos = NULL, type ="source")
    library("longitudinalDynamics")

## <a name="usage"></a> Quick Usage (Longitudinal data)
### <a name="olink"></a> Plasma proteome
    #Load Plasma proteomic expression (NPX) data
    load("data/Olink_NPX_log2_Protein.Rda")
    #Clinical Metdata/annotation
    load("data/data_Annotation.Rda")
    #Run longuitudinal analysis
    long_res <- longitudinalDynamics(metadata=ann, datamatrix=data, 
                         datatype="bulk", omics="Protein",
                         fileName="olink", features=c("PTID", "Time"), 
                         meanThreshold=1, cvThreshold=5,
                         outputDirectory="output")

### <a name="flow"></a> Flow cell-proportion
    #Load Flow-proportion data
    load("data/Flow_matrix.Rda")
    #Clinical Metdata/annotation
    load("data/data_Annotation.Rda")
    #Run longuitudinal analysis
    long_res <- longitudinalDynamics(metadata=ann, datamatrix=data,
                         datatype="bulk", omics="Flow",
                         fileName="Flow", features=c("PTID", "Time"), 
                         meanThreshold=1, cvThreshold=5,
                         outputDirectory="output")

### <a name="cbc"></a> Blood count data
    #Load CBC data
    load("data/CBC_matrix.Rda")
    #Clinical Metdata/annotation
    load("data/data_Annotation.Rda")
    #Run longuitudinal analysis
    long_res <- longitudinalDynamics(metadata=ann, datamatrix=data, 
                         datatype="bulk", omics="CBC",
                         fileName="CBC", features=c("PTID", "Time"), 
                         meanThreshold=1, cvThreshold=5,
                         outputDirectory="output")

### <a name="scrna"></a> Single cell RNA data
    #Seurat object
    pbmc <- readRDS("data/scRNA-PBMC-longitudinaldata.RDS")
    #Clinical Metdata/annotation
    load("data/data_Annotation.Rda")
    metadata=ann
    
    #Celltypes observed
    cell_type <- c("CD4_Naive","CD4_TEM","CD4_TCM","CD4_CTL","CD8_Naive","CD8_TEM","CD8_TCM","Treg","MAIT","gdT","dnT", 
                   "CD4_Proliferating", "CD8_Proliferating","NK_Proliferating",
                   "NK", "NK_CD56bright",
                   "B_naive", "B_memory", "B_intermediate","Plasmablast",
                   "CD14_Mono","CD16_Mono",
                   "cDC1","cDC2","pDC","ASDC",
                   "Platelet","Eryth", "ILC","HSPC","Doublet")
    #Selected celltypes (manuscript) based on proportions >5%
    group_oi <- c("CD4_Naive","CD4_TEM","CD4_TCM","CD4_CTL","CD8_Naive","CD8_TEM","CD8_TCM","Treg","MAIT","gdT",
                  "NK", "NK_CD56bright",
                  "B_naive", "B_memory", "B_intermediate",
                  "CD14_Mono","CD16_Mono",
                  "cDC2","pDC")
    
    #Run longuitudinal analysis
    long_res <- longitudinalDynamics(metadata=ann, scrnaObj=pbmc, datatype="singlecell", omics="rna",
                         fileName="scrna", features=c("PTID", "Time"), 
                         meanThreshold=0.1, cvThreshold=10,
                         avgGroup="celltype",
                         housekeeping_genes=c("GAPDH", "ACTB"),
                         group_oi=group_oi,
                         nPC=15,
                         donorThreshold=4, groupThreshold=40, topFeatures=25,
                         method=c("intra-donor","inter-donor"),
                         outputDirectory="output")
    
### <a name="scatac"></a> Single cell ATAC data
    #Load archR based genescore matrix
    load("load("data/scATAC-genescore_archr.Rda")
    #Clinical Metdata/annotation
    load("data/data_Annotation.Rda")
    
    #Run longuitudinal analysis
    long_res <- longitudinalDynamics(metadata=ann, atacObj=atacObj, datatype="singlecell", omics="atac",
                                     fileName="scatac", features=c("PTID", "Time", "group"), 
                                     meanThreshold=0.1, cvThreshold=10, celllabel=FALSE,
                                     coding_genes=FALSE,
                                     housekeeping_genes=c("GAPDH", "ACTB"), 
                                     #group_oi=group_oi,
                                     nPC=15,
                                     donorThreshold=4, groupThreshold=40, topFeatures=25,
                                     method=c("intra-donor","inter-donor"),
                                     outputDirectory="output")

## <a name="example-main"></a> Tutorials
### <a name="example1"></a> Tutorial-1: Plasma proteome [Bulk dataset]

This tutorial allows users to explore bulk plasma proteome measured from 6 healthy donors over 10 timepoints. Plasma proteomic data available at github. 1. Olink_NPX_log2_Protein.Rda (Normalized protein expression data) 2. data_Annotation.Rda (clinical metadata). Longitudinal dataset have 6 donors (3 male and 3 females). PBMC was collected for 10 weeks. Please follow following steps.

#### 1.1: Load Library and assign parameters
   
    #Load Library
    library("longituinalDynamics")
    library("Hmisc")
    library("ggpubr")
    
#### Load data and assign paramaters

    load("Olink_NPX_log2_Protein.Rda")
    load("data_Annotation.Rda")
    #assign rownames with sample name
    row.names(ann) <- ann$Sample
    #Parameters
    metadata=ann
    datamatrix=data
    features=c("PTID", "Time")
    meanThreshold=1
    cvThreshold=5
    housekeeping_genes <- c("GAPDH", "ACTB")
    
#### Create output directory

    outputDirectory <- "output"
    filePATH <- paste(getwd(), "/",outputDirectory, sep="")
    dir.create(file.path(getwd(), outputDirectory), showWarnings = FALSE)
    
#### Sample overlap

    overlap <- intersect(metadata$Sample, colnames(datamatrix))
    metadata <- metadata[overlap,]
    datamatrix <- data.frame(datamatrix, check.names = F, stringsAsFactors = F)
    datamatrix <- datamatrix[,overlap]

#### 1.2:  Check data
#### PCA Plot

    row.has.na <- apply(datamatrix,1,function(x){any(is.na(x))})
    datamatrix_nonNA <- datamatrix[!row.has.na,]
    res.pca <- suppressMessages(prcomp(t(datamatrix_nonNA),  center=T, scale = TRUE), classes = "message")
    plot1 <- fviz_pca_ind(res.pca, col.ind = metadata$PTID, geom.ind =c("point", "text"),  labelsize = 3, addEllipses=FALSE, ellipse.level=0.95) +
             theme(axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black"), legend.position = "none") +
             theme_classic()
    print(plot1)

#### Sample variability (Correlation)

    cor_mat <- rcorr(as.matrix(datamatrix_nonNA), type="pearson")
    res <- cor_mat$r
    #Plot heatmap
    ha_col <- HeatmapAnnotation(df=data.frame(PTID=metadata$PTID))
    ht1 <- Heatmap(data.matrix(res), cluster_rows =F,  
               cluster_columns = F,
               row_split = as.character(metadata$PTID), column_split = as.character(metadata$PTID),
               na_col = "grey", col = colorRamp2(c(-1,0,0.9,1), c("black","white","pink","red")),
               row_names_max_width=unit(10, "cm"),
               column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize =6),
               top_annotation = ha_col,
               heatmap_legend_param = list(title = "Pearson",heatmap_legend_side = "right") )
    draw(ht1)

#### Remove genes with >40%NAs (optional)

    row.na <- apply(datamatrix,1,function(x) { sum(is.na(x)) })
    row.non.na <- row.na[row.na < (0.4*ncol(datamatrix))] #select features with NAs <40%
    datamatrix <- datamatrix[names(row.non.na),]
    rowN <- data.frame(row.names(datamatrix), stringsAsFactors = F)

#### 1.3:  Features contributing towards donor variations
#### Variance decomposition

    lmem_res <- lmeVariance(ann=metadata, mat=datamatrix, features=features, meanThreshold=meanThreshold)
    res <- lmem_res[,c("PTID","Time","Residual")]
    colnames(res) <- c("donor","week","Residuals")
    res <- res*100 #in percentage
    
#### Donor-specific variance contrubuting features
    
    df1 <- filter(res, donor>week & Residuals < 50)
    df1 <- df1[order(df1$donor, decreasing = T),]
    df <- melt(data.matrix(df1[1:15,]))
    df$Var2 <- factor(df$Var2, levels = rev(c("donor","week", "Residuals")))
    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
    p1 <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity", position="stack") +
        scale_fill_manual(values = c("donor"="#C77CFF", "celltype"="#00BFC4", "week"="#7CAE00", "Residuals"="grey")) +
        theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),legend.position = "right") +
        coord_flip()
    print(p1)

#### Plot the variables

    genelist <- c("FOLR3", "GH2", "MICA", "FOSB", "EIF4G1", "SRC", "DNAJB1", "SIRT2")
    uniSample <- as.character(unique(metadata$PTID))
    splots <- list()
    for(i in 1:length(genelist)) {
      geneName <- genelist[i]
      df <- data.frame(exp=as.numeric(datamatrix[geneName,]), metadata, stringsAsFactors = F)
      df$PTID <- factor(df$PTID, levels = uniSample)
      plot1 <- ggpubr::ggline(df, x = "PTID", y = "exp", 
                add.params = list(shape="Time"), add = c("mean_se", "jitter", "boxplot"), 
                ylab = "NPX", xlab = "Donor", title = geneName, legend = "none", 
                outlier.shape = NA) + scale_shape_manual(values = 0:10)
      splots[[i]] <- plot1
    }
    plot_grid(plotlist=splots, ncol= 3, align="hv")

#### 1.4:  Intra-donor variations over time
#### CV vs Mean

    cv_res <- cvCalcBulk(mat=datamatrix, ann=metadata, meanThreshold=meanThreshold, cvThreshold=cvThreshold)
    CV <- cv_res$CV
    variable_genes <- cv_res$variable_genes
    stable_genes <- cv_res$stable_genes

#### 1.5:  Outlier analysis
#### Calculate IQR

    IQR_res <- iqrBulk(ann=metadata, mat=datamatrix)

#### IQR Plot
    
    df <- melt(data.matrix(IQR_res))
    df <- df[df$value != 0,]
    df <- df[!is.na(df$Var1),]
    df <- df[order(df$value, decreasing = T),]
    plot1 <- ggplot(df, aes(x=Var2, y=value)) +
        geom_violin(scale="width") +
        #geom_boxplot(width=0.1, fill="white") +
        labs(x="", y="Z-score (>2SD)") +
        ggforce::geom_sina(size=0.5) +
        theme_classic() + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1, size=6), axis.text.y = element_text(size=6), legend.position = "right")
    print(plot1)

#### Gene plot (probable outliers)

    genelist <- c("NAA10", "IFI30", "FCAR", "TNFRSF13C", "IL15")
    genelist <- as.character(unique(df$Var1)[1:9])
    splots <- list()
    for(i in 1:length(genelist)) {
       geneName <- genelist[i]
       df <- data.frame(exp=as.numeric(datamatrix[geneName,]), metadata, stringsAsFactors = F)
       df$PTID <- factor(df$PTID, levels = uniSample)
       plot1 <- ggline(df, x = "PTID", y = "exp", add.params = list(shape="Time"), add = c("mean_se", "jitter"), ylab = "NPX", xlab = "Donor", title = geneName, legend = "right", outlier.shape = NA) +
       scale_shape_manual(values = 0:10) + theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
       splots[[i]] <- plot1
    }
    plot_grid(plotlist=splots, ncol= 3, align="hv")
    

### <a name="example2"></a> Tutorial-2: scRNA longitudinal data (n=4 and 6 weeks follow-up)

This tutorial allows users to explore single cell RNAseq data measured from 4 healthy donors over 6 timepoints (week 2-7). Single cell data available at **GEOXXX**. (1) pbmc_longitudinal_data (Normalized scRNA seurat object) (2) data_Annotation.Rda (clinical metadata). Longitudinal dataset have 4 donors (2 male and 2 females). Please follow following steps.

#### 2.1: Load Library and assign parameters
   
    #Load Library
    library("longituinalDynamics")
    library("Hmisc")
    library("ggpubr")
    
#### Load data and assign paramaters

    #scRNA seurat object
    pbmc <- readRDS("data/01-scRNA-PBMC-FinalData.RDS")
    metaData <- pbmc@meta.data
    pbmc@meta.data$Sample <- pbmc@meta.data$orig.ident

    #Clinical metadata/annotation
    load("data/data_Annotation.Rda")
    metadata=ann

    #Parameters
    dataObj <- pbmc
    features=c("PTID", "Time") 
    avgGroup="celltype"
    housekeeping_genes <- c("GAPDH", "ACTB")

    #Celltypes to be considered
    cell_type <- sort(unique(pbmc@meta.data$celltype))
    group_oi <- c("CD4_Naive","CD4_TEM","CD4_TCM","CD4_CTL","CD8_Naive","CD8_TEM","CD8_TCM","Treg","MAIT","gdT",
              "NK", "NK_CD56bright",
              "B_naive", "B_memory", "B_intermediate",
              "CD14_Mono","CD16_Mono",
              "cDC2","pDC")
    
#### Create output directory

    outputDirectory <- "output"
    filePATH <- paste(getwd(), "/",outputDirectory, sep="")
    dir.create(file.path(getwd(), outputDirectory), showWarnings = FALSE)
    
#### Sample overlap

    overlap <- intersect(metadata$Sample, dataObj@meta.data$Sample)
    metadata <- metadata[overlap,]
    #in-case subset of samples only
    dataObj <- subset(x = dataObj, subset = Sample %in% overlap)
 
#### 2.2: Aggregate data at celltypes (psuedo-bulk)
#### For single cell data merge annotation and single cell metadata

    metaData <- dataObj@meta.data
    metadata1 <- metadata[metaData$Sample,]
    metaData <- cbind(metaData, metadata1)
    dataObj@meta.data <- metaData

#### Define sample group and Calculate average expression

    dataObj@meta.data$Sample_group <- paste(dataObj@meta.data$Sample, dataObj@meta.data[,avgGroup], sep=":")
    dataObj@meta.data$Sample_group <- gsub(" ", "_", dataObj@meta.data$Sample_group)
    metaData <- dataObj@meta.data
    
#### Calculate average expression across group/celltype
    
    Idents(dataObj) <- "Sample_group"
    table(Idents(dataObj))
    scrna_avgmat <- avgExpCalc(dataObj=dataObj, group.by="Sample_group")

#### Keep genes with avgExpression > zero

    rowDF <- rowSums(scrna_avgmat)
    rowDF <- rowDF[rowDF > 0]
    mat <- scrna_avgmat[names(rowDF),]

#### Create annotation

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

#### Final matrix

    Overlap <- intersect(colnames(mat), row.names(ann))
    ann <- ann[Overlap,]
    mat <- mat[,Overlap]
    write.table(sort(unique(ann$group)), file=paste(filePATH,"/",fileName,"-group.txt", sep=""), row.names = F, col.names=F, quote=F)

#### 2.3: Check data
#### UMAP plot

    p1 <- DimPlot(object = pbmc, reduction = 'umap', group.by = "Sample", label = F)
    p2 <- DimPlot(object = pbmc, reduction = 'umap', group.by = "celltype", label = F)
    print(plot_grid(p1, p2, align="hv", ncol=2))

#### CV profile

    exp_profile <- cvSampleprofile(mat=mat, ann=ann)
    #Mean Plot
    ggplot(exp_profile, aes(x=mean)) +
      geom_histogram(binwidth=0.1, color="black", fill="white") +
      scale_x_continuous(trans='log10')

    ggplot(exp_profile, aes(x=mean)) +
      geom_histogram(binwidth=0.1, color="black", fill="white") +
      scale_x_continuous(trans='log10') +
      facet_wrap(~group)

    #CV plot
    ggplot(exp_profile, aes(x=cv)) +
      geom_histogram(binwidth=1, color="black", fill="white")
    
    #Housekeeping genes data
    exp_profile <- exp_profile[exp_profile$gene %in% c("ACTB","GAPDH"),]
    ggplot(exp_profile, aes(x=mean, y=cv)) + geom_point() +
      facet_wrap(~gene)

#### 2.4: Features contributing towards donor variations
#### Variance decomposition

    meanThreshold <- 0.1
    lmem_res <- lmeVariance(ann=ann, mat=mat, features=c(features,"group"), meanThreshold=meanThreshold, fileName=fileName, filePATH=filePATH)
    res <- lmem_res[,c("PTID","Time","group","Residual")]
    colnames(res) <- c("donor","week","celltype","Residuals")
    res <- res*100 #in percentage
    
#### Donor-specific variance

    df1 <- filter(res, donor>week & donor>celltype & Residuals < 50)
    df1 <- df1[order(df1$donor, decreasing = T),]
    df <- melt(data.matrix(df1))
    df$Var2 <- factor(df$Var2, levels = rev(c("donor","week", "celltype", "Residuals")))
    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
    p1 <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
      geom_bar(stat="identity", position="stack") +
      scale_fill_manual(values = c("donor"="#C77CFF", "celltype"="#00BFC4", "week"="#7CAE00", "Residuals"="grey")) +
      theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),legend.position = "right") +
      coord_flip()
    print(p1)  

#### Plot the variables

    genelist <- c("MTRNR2L8", "XIST", "RPS4Y1", "RPS26",
                  "JUN", "DAXX", "FOSB", "STAT1",
                  "LILRA4", "CLEC9A", "CST3", "BANK1")
    plotFunction <- function(ann, mat, geneName) {
      df <- data.frame(exp=as.numeric(mat[geneName,]), ann, stringsAsFactors = F)
      wk <- c("W2","W3","W4","W5","W6","W7")
      col <- colorRampPalette(c("cyan", "grey", "black"))(6)
      names(col) <- wk
      df$Time <- factor(df$Time, levels = wk)
      plot1 <- ggplot(df, aes(x=PTID, y=exp, fill=Time)) +
          #geom_bar(stat="identity", position="dodge") +
          #geom_point(position=position_dodge(width=0.9), aes(y=exp), size=0.2) +
          geom_boxplot(outlier.size = 0.1) +
          scale_fill_manual(values = col) +
          labs(x="", y="RNA expression", title=geneName) +
          theme_classic() + theme(axis.text.x = element_text(angle=0, hjust = 0.5, vjust = 1),legend.position = "right")
        
      plot2 <- ggplot(df, aes(x=PTID, y=exp, fill=Time)) +
          #geom_bar(stat="identity", position="dodge") +
          #geom_point(position=position_dodge(width=0.9), aes(y=exp), size=0.2) +
          geom_boxplot(outlier.size = 0.1) +
          scale_fill_manual(values = col) +
          facet_wrap(~group) +
          labs(x="", y="RNA expression", title=geneName) +
          theme_classic() + theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),legend.position = "right")
      return(list(plot1=plot1, plot2=plot2))
    }
    
    plots <- plotFunction(ann, mat, geneName="MTRNR2L8")
    print(plots$plot1)
    plots <- plotFunction(ann, mat, geneName="LILRA4")
    print(plots$plot2)

### <a name="example3"></a> Tutorial-3: scATAC
### <a name="exampl4"></a> Tutorial-4: CNP data
    
    
    
## <a name="authors"></a> Authors

[Suhas Vasaikar](https://github.com/suhasaii), [Aarthi talla](https://github.com/aarthitallaAI) and [Xiaojun Li](https://github.com/Xiaojun-Li) designed the longitudinalDynamics algorithm. [Suhas Vasaikar](https://github.com/suhasaii) implemented the longitudinalDynamics package.

## <a name="acknowledgements"></a> Acknowledgements

We thank Dr. Adam Savage, Dr. Troy Torgerson, Dr. Peter Skene and Dr. Tom Bumol for contribution.


## <a name="license"></a> License

longitudinalDynamics is licensed under the [MIT-License](https://github.com/git/git-scm.com/blob/main/MIT-LICENSE.txt).
