README for longitudinalDynamics
===============

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Install Library](#library)
* [Quick Usage (Longitudinal data)](#usage)
    * [Plasma proteome](#olnk)
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
longitudinalDynamics is an implementation of the longitudinal data analysis platform in R. It allows to identify intra-donor, inter-donor variations over longitudinal time points. The analysis can be done on bulk expression dataset without known celltype information or single cell with celltype information or known groups. It allows to identify stable and variable features in given donor and each celltype (or user defined group). The outlier in intra-donor can be performed to identify perturbed or change in sample currespnding to donor/participant. The stable and variable gene signature helps to identify gene signatures in dataset and perturbation specific featureset.

## <a name="library"></a> Install Library

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
### <a name="ex1"></a> Tutorial-1: Plasma proteome
### <a name="ex1"></a> Tutorial-2: scRNA
### <a name="ex1"></a> Tutorial-3: scATAC
### <a name="ex1"></a> Tutorial-4: CNP data
    
    
    
## <a name="authors"></a> Authors

[Suhas Vasaikar](https://github.com/suhasaii), [Aarthi talla](https://github.com/aarthitallaAI) and [Xiaojun Li](https://github.com/Xiaojun-Li) designed the longitudinalDynamics algorithm. [Suhas Vasaikar](https://github.com/suhasaii) implemented the longitudinalDynamics package.

## <a name="acknowledgements"></a> Acknowledgements

We thank Dr. Adam Savage, Dr. Troy Torgerson, Dr. Peter Skene and Dr. Tom Bumol for contribution.


## <a name="license"></a> License

longitudinalDynamics is licensed under the [MIT-License](https://github.com/git/git-scm.com/blob/main/MIT-LICENSE.txt).
