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
    * [Tutorial-1: Plasma proteome longitudinal data](#example1)
    * [Tutorial-2: scRNA longitudinal data](#example2)
    * [Tutorial-3: scATAC longitudinal data](#example3)
    * [Tutorial-4: CNP0001102 data longitudinal data](#example4)
    * [Tutorial-5: Differential Genes from longitudinal data](#example5)
* [Authors](#authors)
* [License](#license)

* * *

## <a name="introduction"></a> Introduction
LongitudinalDynamics `(longitudinalDynamics)` is a platform for anayzing longitudinal data from bulk as well as single cell. It allows to identify inter-, intra-donor variations in genes over longitudinal time points. The analysis can be done on bulk expression dataset without known celltype information or single cell with celltype/user-defined groups. It allows to infer stable and variable features in given donor and each celltype (or user defined group). The outlier analysis can be performed to identify techinical/biological perturbed samples in donor/participant. Further, differential analysis can be performed to deciher time-wise changes in gene expression in a celtype.

<br> ![img](vignettes/LongitudinalDynamics-workflow.png) <br>
General workflow and analysis schema of **LongitudinalDynamics**. It can work with longitudinal data obtained from bulk such as clinical, bulk RNAseq, proteomic or single cell dataset from scRNAseq, and scATACseq.


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

This tutorial allows users to explore bulk plasma proteome measured from 6 healthy donors over 10 timepoints. Plasma proteomic data available at github. 1. [Olink_NPX_log2_Protein.Rda](https://github.com/aifimmunology/longitudinalDynamics/blob/main/data/Olink_NPX_log2_Protein.Rda) (Normalized protein expression data) 2. [data_Annotation.Rda](https://github.com/aifimmunology/longitudinalDynamics/blob/main/data/data_Annotation.Rda) (clinical metadata). Longitudinal dataset includes 6 donors (3 male and 3 females). PBMC samples were collected from 6 donors over 10 weeks. To interrogate longitudinal data, please follow following steps.

#### 1.1: Load Library
   
    #Load Library and other vizualization packages
    library("longitudinalDynamics")
    library("Hmisc")
    library("ggpubr")
    
#### Load data and assign paramaters

The annotation table `metadata` must consist of column `Sample` (Participant sample name), `PTID` (Participant), `Time` (longitudinal time points). The datamatrix is an Expression data frame, where rows represents gene/proteins and column represents participant samples (same as annotation table `Sample` column).

    load("Olink_NPX_log2_Protein.Rda")
    load("data_Annotation.Rda")
    #assign rownames with sample name
    row.names(ann) <- ann$Sample
    
    #Parameters
    metadata=ann
    datamatrix=data
    featureSet=c("PTID", "Time") #variation attributed to traits
    
#### Create output directory

The output directory where user can save the result files. Default `output` directory is created.

    outputDirectory <- "output"
    filePATH <- paste(getwd(), "/",outputDirectory, sep="")
    dir.create(file.path(getwd(), outputDirectory), showWarnings = FALSE)
    
#### Sample overlap

The data matrix were overlapped with metadata for selecting available samples only.

    overlap <- intersect(metadata$Sample, colnames(datamatrix))
    metadata <- metadata[overlap,]
    datamatrix <- data.frame(datamatrix, check.names = F, stringsAsFactors = F)
    datamatrix <- datamatrix[,overlap]

#### 1.2:  Check data
#### Remove genes with >40%NAs (optional)

For downstream analysis select genes/proteins with less than 40% of missing values. Users can select cut-off for missing values as necessary.

    row.na <- apply(datamatrix,1,function(x) { sum(is.na(x)) })
    row.non.na <- row.na[row.na < (0.4*ncol(datamatrix))] #select features with NAs <40%
    datamatrix <- datamatrix[names(row.non.na),]
    rowN <- data.frame(row.names(datamatrix), stringsAsFactors = F)

#### 1.3:  Features contributing towards donor variations
#### Variance decomposition

To perform variance decomposition apply `lmeVariance` function with input metadata, and datamatrix. The `featureSet` is a list of variables to which freaction variance explained by each gene is attributed. `meanThreshold` defines the minimum average expression threshold to be used for ongitudinal dataset. Here we used normalized protein expression 1 based on mean expression profile of each gene across longitudinal samples. Residuals suggest the variance can not be explained by available feature set. The variance explained by each gene towards the featureSet of interest given in percentage.

    lmem_res <- lmeVariance(ann=metadata, mat=datamatrix, featureSet=featureSet, meanThreshold=1)

<br> <img src="vignettes/Tutorial-1-variance.png" width="50%" height="50%"> <br>

    
    res <- lmem_res[,c("PTID","Time","Residual")]
    colnames(res) <- c("donor","week","Residuals")
    res <- res*100 #in percentage
    head(res)
    
    #Features      donor      week  Residuals
    #FOLR3   99.90070 0.0000000 0.09930098
    #GH2     99.49856 0.0000000 0.50144042
    #PSPN    99.26882 0.1021076 0.62906798
    #CDHR2   99.07933 0.1157406 0.80493162
    #SSC4D   98.82794 0.0000000 1.17206000
    #XPNPEP2 98.67628 0.0000000 1.32372323
    
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
    
<br> <img src="vignettes/Tutorial-1-DonorVariance.png" width="50%" height="50%"> <br>


#### Plot the Top variables

    genelist <- c("FOLR3", "MICA", "EIF4G1", "SRC")
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

<br> <img src="vignettes/Tutorial-1-geneplot.png" width="50%" height="50%"> <br>

#### 1.4:  Intra-donor variations over time
#### CV vs Mean

    cv_res <- cvCalcBulk(mat=datamatrix, ann=metadata, meanThreshold=1, cvThreshold=5)

<br> <img src="vignettes/Tutorial-1-cvDistribution.png" width="50%" height="50%"> <br>
<br> <img src="vignettes/Tutorial-1-VariableStable-Features.png" width="80%" height="80%"> <br>

    CV <- cv_res$CV
    variable_genes <- cv_res$variable_genes
    stable_genes <- cv_res$stable_genes

#### 1.5:  Outlier analysis
#### Sample variability (Correlation)

Perform the sample correlation to find out overall correlation between longitudinal samples.

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
    
<br> <img src="vignettes/Tutorial-1-samplecorrelation.png" width="50%" height="50%"> <br>

#### Detect outliers (if any)

    outlier_res <- outlierDetect(ann=metadata, mat=datamatrix)
    
    head(outlier_res)
    #Var1     Var2    value
    #IFI30 PB1194W6 2.845471
    #DPEP2 PB1194W6 2.844629
    #FCAR PB1194W6 2.844607
    #DPEP1 PB1194W6 2.844574
    #TNFRSF13C PB1194W6 2.844410
    #KIR2DL3 PB1194W6 2.844018
    
    #Z-score Plot
    plot1 <- ggplot(outlier_res, aes(x=sample, y=Z)) +
      geom_violin(scale="width") +
      #geom_boxplot(width=0.1, fill="white") +
      labs(x="", y="Z-score (>2SD)") +
      ggforce::geom_sina(size=0.5) +
      theme_classic() + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1, size=6), axis.text.y = element_text(size=6), legend.position = "right")
    print(plot1)

<br> <img src="vignettes/Tutorial-1-IQRplot.png" width="100%" height="100%"> <br>
    
    #Stringent SD
    outlier_res <- outlierDetect(ann=metadata, mat=datamatrix, SD_threshold= 2.5)
    df <- data.frame(table(outlier_res$sample))
    df <- df[order(df$Freq, decreasing = T),]
    head(df)
    #Var1 Freq
    #PB1194W6   71
    #PB5206W9   28
    #PB2216W5   20
    #PB1051W1   17
    #PB1051W8   12
    #PB7626W1    7
    
    #### Gene plot (probable outliers)

    genelist <- c("IFI30", "DPEP2","FCAR", "TNFRSF13C", "IL15", "IL32")
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

<br> <img src="vignettes/Tutorial-1-IQRgeneplot.png" width="100%" height="100%"> <br>
    
### <a name="example2"></a> Tutorial-2: scRNA longitudinal data (n=4 and 6 weeks follow-up)

This tutorial allows users to explore single cell RNAseq data measured from 4 healthy donors over 6 time points (week 2-7). Single cell data available at **GEOXXX**. (1) pbmc_longitudinal_data (Normalized scRNA seurat object) (2) data_Annotation.Rda (clinical metadata). Longitudinal data set includes 4 donors and 24 samples. To infer iner-donor, intra-donor variations, and stable features, please follow following steps.

#### 2.1: Load Library
   
    #Load Library
    library("longitudinalDynamics")
    library("Hmisc")
    library("ggpubr")
    library("Seurat")
    
#### Load data and assign paramaters

    #scRNA seurat object
    pbmc <- readRDS("data/01-scRNA-PBMC-FinalData.RDS")
    metaData <- pbmc@meta.data
    pbmc@meta.data$Sample <- pbmc@meta.data$orig.ident
    #UMAP plot
    p1 <- DimPlot(object = pbmc, reduction = 'umap', group.by = "Sample", label = F)
    p2 <- DimPlot(object = pbmc, reduction = 'umap', group.by = "celltype", label = F)
    print(plot_grid(p1, p2, align="hv", ncol=2))
    
<br> <img src="vignettes/Tutorial-2-UMAPPlot.png" width="100%" height="100%"> <br>

#### Load annotation data

    #Clinical metadata/annotation
    load("data/data_Annotation.Rda")
    metadata=ann
    row.names(metadata) <- metadata$Sample

    #Parameters
    dataObj <- pbmc
    features=c("PTID", "Time") 
    avgGroup="celltype"
    housekeeping_genes <- c("GAPDH", "ACTB")

    #Celltypes observed in dataset
    cell_type <- sort(unique(pbmc@meta.data$celltype))
    #Celltypes selected for analysis consisting atleast >5% of cells in each celltype.
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
    print(sort(unique(ann$group)))

#### 2.3: CV profile

    #Check the mean expression and CV cross groups (celltypes)
    cv_profile <- cvprofile(mat=mat, ann=ann)

<br> <img src="vignettes/Tutorial-2-CVPlot-1.png" width="50%" height="50%"> <br>

    #Lowly expressed genes show abnormal CV, which needs to be filtered.
    cv_profile <- cvprofile(mat=mat, ann=ann, meanThreshold = 0.1)

<br> <img src="vignettes/Tutorial-2-CVPlot-2.png" width="50%" height="50%"> <br>

    #Sample Celltype Mean-CV plot (check output folder)
    cvSampleprofile(mat=mat, ann=ann, meanThreshold = 0.1, cvThreshold = 10)
    
#### 2.4: Features contributing towards donor variations
#### Variance decomposition

    meanThreshold <- 0.1
    lmem_res <- lmeVariance(ann=ann, mat=mat, featureSet=c(featureSet,"group"), meanThreshold=meanThreshold, fileName=fileName, filePATH=filePATH)
    res <- lmem_res[,c("PTID","Time","group","Residual")]
    colnames(res) <- c("PTID","Time","celltype","Residuals")
    res <- res*100 #in percentage
    
<br> <img src="vignettes/Tutorial-2-variancePlot.png" width="50%" height="50%"> <br>
    
#### Donor-specific variance

    df1 <- filter(res, donor>week & donor>celltype & Residuals < 50)
    df1 <- df1[order(df1$donor, decreasing = T),]
    df <- melt(data.matrix(df1[1:15,)) #Top15
    df$Var2 <- factor(df$Var2, levels = rev(c("donor","week", "celltype", "Residuals")))
    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
    p1 <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
      geom_bar(stat="identity", position="stack") +
      scale_fill_manual(values = c("donor"="#C77CFF", "celltype"="#00BFC4", "week"="#7CAE00", "Residuals"="grey")) +
      theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),legend.position = "right") +
      coord_flip()
    print(p1)  

<br> <img src="vignettes/Tutorial-2-Donor-variancePlot2.png" width="50%" height="50%"> <br>

    #Similar procedure applied to get Time- and celltype-attributed variance features
    
<br><img src="vignettes/Tutorial-2-Donor-variancePlot.png" width="100%" height="100%"> <br>

#### Plot the variables

    plots <- genePlot(ann, mat, geneName="LILRA4")
    print(plots$plot1)
    
<br> <img src="vignettes/Tutorial-2-celltype-LILRA4-1.png" width="50%" height="50%"> <br>
    
    print(plots$plot2)
    
<br> <img src="vignettes/Tutorial-2-celltype-LILRA4-2.png" width="50%" height="50%"> <br>
    
#### 2.5: Intra-donor variations over time
#### Calculate CV

    meanThreshold=0.1
    cvThreshold=10
    cv_res <- cvCalcSC(mat=mat, ann=ann, meanThreshold=meanThreshold, cvThreshold=cvThreshold, housekeeping_genes=housekeeping_genes, filePATH=filePATH, fileName="scrna")
    
    #Plots saved in user-defined output directory
    #Overall CV profile in celltypes (black) as well as CV for house-keeping genes (blue). Based on CV of house-keeping genes 10% CV cut-off used and genes considered stable elow 10% CV.
    
<br> <img src="vignettes/Tutorial-2-CV-distribution-HousekeepingGenes.png" width="50%" height="50%"> <br>

#### Find stable and variable features in longitudinal data

    donorThreshold <- 4
    groupThreshold <- 40 #number of donors * number of celltypes/2 (4x19/2)
    topFeatures <- 25
    var_gene <- VarFeatures(ann=ann, group_oi=group_oi, meanThreshold=meanThreshold, cvThreshold=cvThreshold, donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=topFeatures, filePATH=filePATH, fileName=fileName)

    stable_gene <- StableFeatures(ann=ann, group_oi=group_oi, meanThreshold=meanThreshold, cvThreshold=cvThreshold, donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=topFeatures, filePATH=filePATH, fileName=fileName)

    #Variable genes observed in longitidinal data (CV>10%)

<br> <img src="vignettes/Tutorial-2-Variable-Plot.png" width="100%" height="100%"> <br>
    
    #Stable genes observed in longitidinal data (CV<10%)

<br> <img src="vignettes/Tutorial-2-Stable-Plot.png" width="100%" height="100%"> <br>

#### UMAP Plot

    #Top variable and stable features used for UMAP
    rnaObj <- dimUMAPPlot(rnaObj=dataObj, nPC=15, gene_oi=var_gene, groupName=avgGroup, plotname="variable", filePATH=filePATH, fileName=fileName)
    
<br> <img src="vignettes/Tutorial-2-scRNA-UMAP-variable-Genes.png" width="100%" height="100%"> <br>
    
    rnaObj <- dimUMAPPlot(rnaObj=dataObj, nPC=15, gene_oi=stable_gene, groupName=avgGroup, plotname="stable", filePATH=filePATH, fileName=fileName)

<br> <img src="vignettes/Tutorial-2-scRNA-UMAP-stable-Genes.png" width="100%" height="100%"> <br>


#### Circular gene expression plot

    load("output/scrna-CV-allgenes-raw.Rda")
    geneList <- c("IL32","CCL5","TCF7","IL7R","LEF1") #T-cell
    res <- genecircosPlot(data=cv_res, geneList=geneList, groupBy=group_oi)
    
<br> <img src="vignettes/Tutorial-2-Tcelltype-circularPlot.png" width="50%" height="50%"> <br>

### <a name="example3"></a> Tutorial-3: scATAC Longitudinal data (n=4 and 6 weeks follow-up)
This tutorial allows users to explore single cell ATACseq genscore data measured from 4 healthy donors over 6 timepoints (week 2-7). Single cell ATAC data available at GEOXXX. (1) pbmc_scatac_archr_genescore_longitudinal_data (2) data_Annotation.Rda (clinical metadata). Longitudinal dataset have 4 donors and 18 samples. To infer the variations at single cell ATAC please follow following steps.

#### 3.1: Load Library

    #Load Library
    library("longitudinalDynamics")
    library("Hmisc")
    library("ggpubr")

#### Load data and assign paramaters

    #scATAC object
    #Load genescorematrix from archR or relevant tools (Aggregate data at celltypes (psuedo-bulk))
    load("data/AIFI-scATAC-PBMC-FinalData.Rda")
    boxplot(log2(scatac_gm[,1:50]+1), las=2)
    
    #Load annotation data
    load("data/data_Metadata.Rda")
    metadata <- ann
    row.names(metadata) <- metadata$Sample

#### Parameters

    metadata <- ann
    atacObj <- log2(scatac_gm+1)
    featureSet=c("PTID", "Time") 
    avgGroup="celltype"
    housekeeping_genes <- c("GAPDH", "ACTB")
    fileName <- "scATAC"
    #Celltypes to be considered
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

    #Get the annotation
    temp <- data.frame(do.call(rbind, strsplit(colnames(atacObj), split = ":")), stringsAsFactors = F)
    cn <- data.frame(id=colnames(atacObj), Sample=temp$X1, group=temp$X2, check.names=F, stringsAsFactors = F)
    cn <- cn[cn$Sample %in% metadata$Sample,]
    overlap <- as.character(unique(cn$Sample))
    atac_overlap <- cn$id
    
    #Sample overlap
    metadata <- metadata[overlap,]
    atacObj <- atacObj[,atac_overlap]

#### 3.2: Check data
#### Keep genes with avgExpression > zero

    rowDF <- rowSums(atacObj)
    rowDF <- rowDF[rowDF > 0]
    mat <- atacObj[names(rowDF),]

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
    print(sort(unique(ann$group)))

#### 3.3: CV profile
#### CV profile

    cv_profile <- cvprofile(mat=mat, ann=ann, housekeeping_genes=housekeeping_genes)
    
<br> <img src="vignettes/Tutorial-3-cvDistribution-1.png" width="100%" height="100%"> <br>

    cv_profile <- cvprofile(mat=mat, ann=ann, housekeeping_genes=housekeeping_genes, meanThreshold = 0.1)
    
<br> <img src="vignettes/Tutorial-3-cvDistribution-2.png" width="100%" height="100%"> <br>

    #Sample Celltype Mean-CV plot
    cv_sample_profile <- cvSampleprofile(mat=mat, ann=ann, meanThreshold = 0.1, cvThreshold = 10)
    #plots saved in output directory

#### 3.4: Features contributing towards donor variations

    meanThreshold <- 0.1
    lmem_res <- lmeVariance(ann=ann, mat=mat, featureSet=c(featureSet,"group"), meanThreshold=meanThreshold, fileName=fileName, filePATH=filePATH)
    res <- lmem_res[,c("PTID","Time","group","Residual")]
    colnames(res) <- c("PTID","Time","celltype","Residuals")
    res <- res*100 #in percentage

<br> <img src="vignettes/Tutorial-3-variancePlot.png" width="50%" height="50%"> <br> 

#### Donor-specific
    df1 <- filter(res, PTID>Time & PTID>celltype & Residuals < 50)
    df1 <- df1[order(df1$PTID, decreasing = T),]
    df <- melt(data.matrix(df1[1:15,])) #select Top15
    df$Var2 <- factor(df$Var2, levels = rev(colnames(res)))
    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
    p1 <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
      geom_bar(stat="identity", position="stack") + labs(x="Features", y="% Variance explained") +
      scale_fill_manual(values = c("PTID"="#C77CFF", "celltype"="#00BFC4", "Time"="#7CAE00", "Residuals"="grey")) +
      theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),legend.position = "right") +
      coord_flip()
    
#### Time-specific
    df2 <- filter(res, PTID<Time & Time>celltype)
    df2 <- df2[order(df2$Time, decreasing = T),]
    df <- melt(data.matrix(df2[1:15,])) #select Top15
    df$Var2 <- factor(df$Var2, levels = rev(colnames(res)))
    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
    p2 <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
      geom_bar(stat="identity", position="stack") + labs(x="Features", y="% Variance explained") +
      scale_fill_manual(values = c("PTID"="#C77CFF", "celltype"="#00BFC4", "Time"="#7CAE00", "Residuals"="grey")) +
      theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),legend.position = "right") +
      coord_flip()
    
#### celltype-specific
    df3 <- filter(res, celltype>PTID & celltype>Time & Residuals < 50)
    df3 <- df3[order(df3$celltype, decreasing = T),]
    df <- melt(data.matrix(df3[1:15,])) #select Top15
    df$Var2 <- factor(df$Var2, levels = rev(colnames(res)))
    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
    p3 <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
      geom_bar(stat="identity", position="stack") + labs(x="Features", y="% Variance explained") +
      scale_fill_manual(values = c("PTID"="#C77CFF", "celltype"="#00BFC4", "Time"="#7CAE00", "Residuals"="grey")) +
      theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),legend.position = "right") +
      coord_flip()
    
    #Plot
    plot_grid(p1,p2,p3, align="hv", ncol=3)
    
<br> <img src="vignettes/Tutorial-3-Donor-variancePlot.png" width="100%" height="100%"> <br>

    #Top genes
    plots <- genePlot(ann, mat, geneName="FIRRE", groupName="group")
    print(plots$plot1)
    print(plots$plot2)
    print(plots$plot3)
    plots <- genePlot(ann, mat, geneName="DPP6", groupName="group")
    plots <- genePlot(ann, mat, geneName="SPIB", groupName="group")

<br> <img src="vignettes/Tutorial-3-geneplot.png" width="100%" height="100%"> <br>

#### 3.4: Intra-donor variations over time

    meanThreshold=0.1
    cvThreshold=10
    cv_res <- cvCalcSC(mat=mat, ann=ann, meanThreshold=meanThreshold, 
                       cvThreshold=cvThreshold, housekeeping_genes=housekeeping_genes, 
                       filePATH=filePATH, fileName=fileName)
                       
<br> <img src="vignettes/Tutorial-3-CV-distribution.png" width="50%" height="50%"> <br>
    
#### Find stable and variable features in longitudinal data

    donorThreshold <- 4
    groupThreshold <- 28 #number of donors * number of celltypes/2 (4x14/2)
    topFeatures <- 25
    stable_gene <- StableFeatures(ann=ann, meanThreshold=meanThreshold, cvThreshold=cvThreshold, donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=topFeatures, filePATH=filePATH, fileName=fileName)
    
<br> <img src="vignettes/Tutorial-3-Stable-Plot.png" width="100%" height="100%"> <br>

    var_gene <- VarFeatures(ann=ann, meanThreshold=meanThreshold, cvThreshold=cvThreshold, donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=topFeatures, filePATH=filePATH, fileName=fileName)

<br> <img src="vignettes/Tutorial-3-Variable-Plot.png" width="100%" height="100%"> <br>

#### Circos CV plot

    load(paste(filePATH,"/",fileName,"-CV-allgenes-raw.Rda", sep=""))
    geneList <- c("HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DPA1","HLA-DRB1", "ACTB","GAPDH")
    res <- genecircosPlot(data=cv_res, geneList=geneList, colorThreshold=15)
    
    group_oi <- c("CD4_Naive","CD4_TEM","CD4_TCM","CD4_CTL",
                  "CD8_Naive","CD8_TEM","CD8_TCM","Treg","MAIT","gdT",
                  "NK", "NK_CD56bright",
                  "B_naive","B_memory", "B_intermediate",
                  "CD14_Mono","CD16_Mono",
                  "cDC2","pDC")
    res <- genecircosPlot(data=cv_res, geneList=geneList, group_oi=group_oi, colorThreshold=15)
    
<br> <img src="vignettes/Tutorial-3-celltype-circularPlot.png" width="100%" height="100%"> <br>


### <a name="example4"></a> Tutorial-4: CNP0001102 data
This tutorial allows users to explore single cell RNAseq data variability across COVID and FLU donors. PBMC from the patients were collected longitudinally. Single cell data from [Zhu et al. 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7368915/) downloaded from [here](https://db.cngb.org/search/project/CNP0001102/). Metadata is downloaded from table and can be found in the [data](https://github.com/aifimmunology/longitudinalDynamics/tree/main/data). To infer variability (inter- and Intra-) and identify stable genes, please follow following steps.

#### 4.1: Load Library
   
    #Load Library
    library("longitudinalDynamics")
    library("Hmisc")
    library("ggpubr")
    library("Seurat")
    
#### Load data and assign paramaters

    #scRNA
    pbmc <- readRDS("data/CNP0001102_Final_nCoV_0716_upload.RDS")
    pbmc@meta.data$Sample <- pbmc@meta.data$batch
    pbmc@meta.data$celltype <- gsub(" ", "_", pbmc@meta.data$cell_type)
        
    #Clinical annotations (Table S1. Clinical data of the enrolled subjects)
    metadata <- read.csv("data/CNP0001102-annotation.csv", stringsAsFactors = F)
    row.names(metadata) <- metadata$Sample

    #Exploring only COVID samples 
    metadata <- metadata[metadata$PTID %in% c("COV-1", "COV-2", "COV-3", "COV-4", "COV-5"),]
    #Exploring only FLU samples
    metadata <- metadata[metadata$PTID %in% c("IAV-1","IAV-2"),]
        
    #Parameters
    dataObj <- pbmc
    featureSet=c("PTID", "Time") 
    avgGroup="celltype"
    housekeeping_genes <- c("GAPDH", "ACTB")
    cell_type <- sort(unique(pbmc@meta.data$celltype))
    fileName="CNP0001102"
    
#### Create output directory (optional)

    outputDirectory <- "output"
    filePATH <- paste(getwd(), "/",outputDirectory, sep="")
    dir.create(file.path(getwd(), outputDirectory), showWarnings = FALSE)
    
#### Sample overlap

    overlap <- intersect(metadata$Sample, dataObj@meta.data$Sample)
    metadata <- metadata[overlap,]
    #in-case subset of samples only
    dataObj <- subset(x = dataObj, subset = Sample %in% overlap)
 
#### 4.2: Aggregate data at celltypes (psuedo-bulk)
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
    print(unique(ann$group))
    
    #[1] "Activated_CD4_T_cells" "Cycling_Plasma"        "Cycling_T_cells"       "Cytotoxic_CD8_T_cells"
    #[5] "MAIT"                  "Megakaryocytes"        "Memory_B_cells"        "Naive_B_cells"        
    #[9] "Naive_T_cells"         "NKs"                   "Plasma"                "Stem_cells"           
    #[13] "XCL+_NKs"              "DCs"                   "Monocytes" 

#### 4.3: CV profile

    cv_profile <- cvprofile(mat=mat, ann=ann)

<br> <img src="vignettes/Tutorial-4-cvDistribution.png" width="100%" height="100%"> <br>

    #Sample Celltype Mean-CV plot (output directory)
    cv_sample_profile <- cvSampleprofile(mat=mat, ann=ann, meanThreshold = 0.1, cvThreshold = 25)

#### 4.4: Features contributing towards donor variations
#### Variance decomposition

    meanThreshold <- 0.1
    lmem_res <- lmeVariance(ann=ann, mat=mat, featureSet=c(featureSet,"group"), meanThreshold=meanThreshold, fileName=fileName, filePATH=filePATH)
    res <- lmem_res[,c("PTID","Time","group","Residual")]
    colnames(res) <- c("PTID","Time","celltype","Residuals")
    res <- res*100 #in percentage
    
<br> <img src="vignettes/Tutorial-4-variancePlot.png" width="50%" height="50%"> <br>
    
    head(res[order(res$celltype, decreasing = T),])
    #        PTID         Time celltype Residuals
    #MZB1    3.634604e-08 1.242595e-09 93.91550  6.084498
    #JCHAIN  0.000000e+00 0.000000e+00 93.67179  6.328213
    #MKI67   3.013979e-01 1.236726e-01 92.19240  7.382530
    #MANF    3.197369e-01 0.000000e+00 91.82359  7.856671
    #SDF2L1  5.210284e-02 1.019559e-01 90.75454  9.091404
    #POU2AF1 2.432974e-01 3.463273e-01 90.24514  9.165233

#### Variance explained

    #Donor-specific
    df1 <- filter(res, PTID>Time & PTID>celltype & Residuals < 50)
    df1 <- df1[order(df1$PTID, decreasing = T),]
    df <- melt(data.matrix(df1))
    df$Var2 <- factor(df$Var2, levels = rev(colnames(res)))
    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
    p1 <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
      geom_bar(stat="identity", position="stack") + labs(x="Features", y="% Variance explained") +
      scale_fill_manual(values = c("PTID"="#C77CFF", "celltype"="#00BFC4", "Time"="#7CAE00", "Residuals"="grey")) +
      theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),legend.position = "right") +
      coord_flip()
    
    #Time-specific
    df2 <- filter(res, Time>PTID & Time>celltype & Residuals < 50)
    #df2 <- res[order(res$Time, decreasing = T),]
    #select Top15
    df2 <- df2[1:15,]
    df <- melt(data.matrix(df2))
    df$Var2 <- factor(df$Var2, levels = rev(colnames(res)))
    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
    p2 <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
      geom_bar(stat="identity", position="stack") + labs(x="Features", y="% Variance explained") +
      scale_fill_manual(values = c("PTID"="#C77CFF", "celltype"="#00BFC4", "Time"="#7CAE00", "Residuals"="grey")) +
      theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),legend.position = "right") +
      coord_flip()
    
    #celltype-specific
    df3 <- filter(res, celltype>PTID & celltype>Time & Residuals < 50)
    df3 <- df3[order(df3$celltype, decreasing = T),]
    #select Top15
    df3 <- df3[1:15,]
    df <- melt(data.matrix(df3))
    df$Var2 <- factor(df$Var2, levels = rev(colnames(res)))
    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
    p3 <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
      geom_bar(stat="identity", position="stack") + labs(x="Features", y="% Variance explained") +
      scale_fill_manual(values = c("PTID"="#C77CFF", "celltype"="#00BFC4", "Time"="#7CAE00", "Residuals"="grey")) +
      theme_bw() + theme(axis.text.x = element_text(angle=90, hjust = 0.5, vjust = 1),legend.position = "right") +
      coord_flip()
    
    #Plot
    plot_grid(p1,p2,p3, align="hv", ncol=3) 
    
<br> <img src="vignettes/Tutorial-4-Donor-variancePlot.png" width="100%" height="100%"> <br>


#### Plot the variables

    plots <- plotFunction(ann, mat, geneName="MKI67", groupName="group")
    print(plots$plot1)
    print(plots$plot2)
    print(plots$plot3)
    
<br> <img src="vignettes/Tutorial-4-celltype-MKI67.png" width="100%" height="100%"> <br>
    
#### 4.5: Intra-donor variations over time
#### Calculate CV

    meanThreshold=0.1
    cvThreshold=25
    cv_res <- cvCalcSC(mat=mat, ann=ann, meanThreshold=meanThreshold, cvThreshold=cvThreshold, housekeeping_genes=housekeeping_genes, filePATH=filePATH, fileName=fileName)

<br> <img src="vignettes/Tutorial-4-cvDistribution.png" width="100%" height="100%"> <br>

#### Find stable and variable features in longitudinal data

    donorThreshold <- 5 #number of donors
    groupThreshold <- 38 #number of donors * number of celltypes/2 (5x15/2)
    topFeatures <- 25
    stable_gene <- StableFeatures(ann=ann, meanThreshold=meanThreshold, cvThreshold=cvThreshold, donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=topFeatures, filePATH=filePATH, fileName=fileName)
    
<br> <img src="vignettes/Tutorial-4-Stable-Plot.png" width="100%" height="100%"> <br>

    var_gene <- VarFeatures(ann=ann, meanThreshold=meanThreshold, cvThreshold=cvThreshold, donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=topFeatures, filePATH=filePATH, fileName=fileName)
    
<br> <img src="vignettes/Tutorial-4-Variable-Plot.png" width="100%" height="100%"> <br>

#### UMAP Plot

    #Stable genes UMAP
    dimUMAPPlot(rnaObj=dataObj, nPC=15, gene_oi=unique(stable_gene$gene), groupName=avgGroup, plotname="stable", filePATH=filePATH, fileName=fileName)

<br> <img src="vignettes/Tutorial-4-scRNA-UMAP-stable-Genes.png" width="100%" height="100%"> <br>
    
    #Variable genes UMAP
    dimUMAPPlot(rnaObj=dataObj, nPC=15, gene_oi=unique(var_gene$gene), groupName=avgGroup, plotname="variable", filePATH=filePATH, fileName=fileName)

<br> <img src="vignettes/Tutorial-4-scRNA-UMAP-variable-Genes.png" width="100%" height="100%"> <br>

#### Circos CV Plot

    load(paste(filePATH,"/",fileName, "-CV-allgenes-raw.Rda", sep=""))
    geneList <- c("IRF3","MAP4K4","XPC","DNAJB6", "KLF13") #Activated CD4 T-cells
    geneList <- c("HMGN2", "IFI16", "PTGES3", "SH3KBP1", "PTBP1") #Cycling T-cells
    res <- genecircosPlot(data=cv_res, geneList=geneList, colorThreshold=cvThreshold)
    
<br> <img src="vignettes/Tutorial-4-celltype-circularPlot.png" width="100%" height="100%"> <br>
    
### <a name="example5"></a> Tutorial-5: Differential Genes from longitudinal data]

This tutorial allows users to identify differential expressed genes in direction of time-points. As an example single cell data from [Zhu et al. 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7368915/) downloaded from [here](https://db.cngb.org/search/project/CNP0001102/). Metadata is downloaded from table and can be found in the [data](https://github.com/aifimmunology/longitudinalDynamics/tree/main/data). The dataset consists of 5 Covid-19 donors, 2 Flu donors with longitudinal data and 3 controls. To explore differetial expressed gened in each celltype of each donor we used hurdle model based modeling on input data to retrive the DEGs. To infer DEGs in each celltype towards time progression (timepoints considered as continoues if more than 2), please follow following steps.

#### 5.1: load data and clinical metadata
#### Single cell object CNP0001102
    
    pbmc <- readRDS("data/CNP0001102_Final_nCoV_0716_upload.RDS")
    #Add column Sample and group as celltype
    pbmc@meta.data$Sample <- pbmc@meta.data$batch
    pbmc@meta.data$celltype <- gsub(" ", "_", pbmc@meta.data$cell_type)
    
#### Clinical annotations [Table S1. Clinical data of the enrolled subjects](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7368915/)
    
    metadata <- read.csv("data/CNP0001102-annotation.csv", stringsAsFactors = F)
    row.names(metadata) <- metadata$Sample
    
#### Load library and run

    library("longitudinalDynamics")
    #run
    DEGres <- sclongitudinalDGE(ann=metadata, dataObj=pbmc, scassay="RNA", celltypecol="celltype")
    
    >Fitting a ZLM model for donorID:  COV-1 
    >Fitting a ZLM model for donorID:  COV-2 
    >Fitting a ZLM model for donorID:  COV-3 
    >Fitting a ZLM model for donorID:  COV-4 
    >Fitting a ZLM model for donorID:  COV-5 
    >Fitting a ZLM model for donorID:  IAV-1 
    >Fitting a ZLM model for donorID:  IAV-2 

    #Plots can be seen in output directory output
    head(DEGres[order(DEGres$coef, decreasing = T),])
    
    #primerid contrast         nomP     coef         adjP donorID celltype dir
    #IGHG4   TimeD9 1.701579e-26 3.056092 1.453999e-23   IAV-2   Plasma upregulated at D9
    #JCHAIN   TimeD9 8.759407e-32 2.647757 2.245474e-28   IAV-2   Plasma upregulated at D9
    #IGHG3   TimeD9 8.470810e-22 2.485250 2.412769e-19   IAV-2   Plasma upregulated at D9
    #IGLC2   TimeD9 1.352490e-16 2.289544 8.890022e-15   IAV-2   Plasma upregulated at D9
    #IGHG1   TimeD9 1.292885e-16 2.219146 8.608601e-15   IAV-2   Plasma upregulated at D9
    #SYNE2   TimeD4 7.249010e-13 2.209541 1.215659e-09   COV-4 XCL+_NKs upregulated at D4

General analysis schema and differential results in each donor over timepoints in celltype Cytotoxic CD8 T-cells using `longitudinalDynamics` shown.

<br><br> ![](vignettes/img5a-CNP0001102-DEGs.png) <br><br>

## <a name="authors"></a> Authors

[Suhas Vasaikar](https://github.com/suhasaii), [Aarthi talla](https://github.com/aarthitallaAI) and [Xiaojun Li](https://github.com/Xiaojun-Li) designed the longitudinalDynamics algorithm. [Suhas Vasaikar](https://github.com/suhasaii) implemented the longitudinalDynamics package.

## <a name="license"></a> License

longitudinalDynamics is licensed under the [MIT-License](https://github.com/git/git-scm.com/blob/main/MIT-LICENSE.txt).
