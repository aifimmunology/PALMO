## Libraries

import rpy2.robjects.packages as rpackages
import pandas as pd
import numpy as np
import os

# working with dataframes
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.vectors import StrVector

# TODO: suppress system warnings from imports
from rpy2.robjects.packages import importr

utils = importr('utils')
# getting ready to download these packages...
importr('PALMO')
importr('Hmisc')
rpackages.importr('ggpubr')


def rpy_2py(df):
    ''' convenient function that will convert rpy object to a pandas data.frame 
    '''
    with localconverter(ro.default_converter + pandas2ri.converter):
        return (ro.conversion.rpy2py(df))


def py2_rpy(df):
    ''' convenient function that will convert a pandas data.frame to rpy object that can then be used 
        for R functions 
    '''
    with localconverter(ro.default_converter + pandas2ri.converter):
        return (ro.conversion.py2rpy(df))


def load_data(data_fp, meta_fp, datatype):
    ''' 
    '''
    if datatype == 'bulk':
        palmo_obj = ro.r('''
            library("PALMO")
            load("{}")
            load("{}")
            palmo_obj <- createPALMOobject(anndata=ann, data=data)
        '''.format(data_fp, meta_fp))
    elif datatype == 'singlecell':
        palmo_obj = ro.r('''
            pbmc <- readRDS("{}")
            metaData <- pbmc@meta.data 
            pbmc@meta.data$Sample <- pbmc@meta.data$orig.ident
            pbmc@meta.data$celltype <- gsub(" ", "_", pbmc@meta.data$celltype)

            # load annotation 
            load("{}")
            palmo_obj <- createPALMOobject(anndata=ann, data=pbmc)              
        '''.format(data_fp, meta_fp))

    return palmo_obj


def prep_data(palmo_obj, datatype, na_cutoff):
    """
    """
    if datatype == "bulk":
        palmo_obj = ro.r('''
            palmo_obj <- annotateMetadata(data_object=palmo_obj, sample_column="Sample", donor_column="PTID", time_column="Time")
            palmo_obj <- mergePALMOdata(data_object=palmo_obj, datatype="{dt}")
            palmo_obj <- checkReplicates(data_object=palmo_obj, mergeReplicates=T)
            palmo_obj <- naFilter(data_object=palmo_obj, na_cutoff={na})
            '''.format(dt=datatype, na=float(na_cutoff)))
    if datatype == "singlecell":
        palmo_obj = ro.r('''
            palmo_obj <- annotateMetadata(data_object=palmo_obj, sample_column="Sample", donor_column="PTID", time_column="Time")
            palmo_obj <- mergePALMOdata(data_object=palmo_obj, datatype="{dt}")
            palmo_obj <- avgExpCalc(data_object=palmo_obj, assay="RNA", group_column="{ct}")
            palmo_obj <- checkReplicates(data_object=palmo_obj, mergeReplicates=T)
        '''.format(dt=datatype, ct='celltype'))
    return palmo_obj


def run_lmeVariance(palmo_obj, mean_threshold, feature_set, datatype):
    '''
    '''
    if datatype == 'bulk':
        # load lmeVariance object function and run
        palmo_obj = ro.r('''
            featureSet <- c("{f1}", "{f2}")
            palmo_obj <- lmeVariance(data_object=palmo_obj, featureSet=featureSet, meanThreshold={mean}, fileName="olink")
            '''.format(f1=feature_set[0],
                       f2=feature_set[1],
                       mean=mean_threshold))
    elif datatype == 'singlecell':
        palmo_obj = ro.r('''
            featureSet <- c("{f1}", "{f2}", "{f3}")
            palmo_obj <- lmeVariance(data_object=palmo_obj, cl=4, featureSet=featureSet, meanThreshold={mean}, fileName="scrna")
            '''.format(f1=feature_set[0],
                       f2=feature_set[1],
                       f3=feature_set[2],
                       mean=mean_threshold))
    return palmo_obj


def run_outlier_detection(palmo_obj, z_cutoff):
    '''
    '''
    # robject of outlierDetect function
    palmo_obj = ro.r('''
        palmo_obj <- outlierDetect(data_object=palmo_obj, z_cutoff={})
        
        # assign cutoffs 
        '''.format(z_cutoff))
    return palmo_obj


def run_cvCalc(palmo_obj, meanThreshold, cvThreshold, datatype):
    """ Intro-donor variations over time 
        TODO: label top 10 genes, and most stable ones 
    """
    if datatype == 'bulk':
        palmo_obj = ro.r('''
            palmo_obj <- cvCalcBulk(data_object=palmo_obj, meanThreshold={m}, cvThreshold={c})
            '''.format(m=meanThreshold, c=cvThreshold))
    elif datatype == 'singlecell':
        palmo_obj = ro.r('''
            palmo_obj <- cvCalcSC(data_object=palmo_obj, meanThreshold={m}, cvThreshold={c})
            '''.format(m=meanThreshold, c=cvThreshold))
    return palmo_obj


def run_cvcalc_scprofile(palmo_obj, mean_threshold, housekeeping_genes):
    palmo_obj = ro.r('''
    palmo_obj <- cvCalcSCProfile(data_object=palmo_obj, meanThreshold={avg}, housekeeping_genes=c("GAPDH", "ACTB"), fileName="scrna")
    '''.format(avg=mean_threshold))
    return palmo_obj


def find_variable_features(palmo_obj, cv_threshold):
    palmo_obj = ro.r('''
        donorThreshold <- 4 
        groupThreshold <- 40 
        celltype_oi <- c("CD4_Naive","CD4_TEM","CD4_TCM","CD4_CTL","CD8_Naive",
                  "CD8_TEM","CD8_TCM","Treg","MAIT","gdT",
                  "NK", "NK_CD56bright",
                  "B_naive", "B_memory", "B_intermediate",
                  "CD14_Mono","CD16_Mono",
                  "cDC2","pDC")
        palmo_obj <- VarFeatures(data_object=palmo_obj, group_oi=celltype_oi,
            cvThreshold={cv}, groupThreshold=groupThreshold, topFeatures=25, 
            fileName="scrna")
    '''.format(cv=cv_threshold))
    return palmo_obj


def find_stable_features(palmo_obj, cv_threshold):
    '''
    '''
    palmo_obj = ro.r('''
    donorThreshold <- 4 
    groupThreshold <- 40 
    celltype_oi <- c("CD4_Naive","CD4_TEM","CD4_TCM","CD4_CTL","CD8_Naive",
                  "CD8_TEM","CD8_TCM","Treg","MAIT","gdT",
                  "NK", "NK_CD56bright",
                  "B_naive", "B_memory", "B_intermediate",
                  "CD14_Mono","CD16_Mono",
                  "cDC2","pDC")
    StableFeatures(data_object=palmo_obj, group_oi=celltype_oi, cvThreshold={}, donorThreshold=donorThreshold, groupThreshold=groupThreshold, topFeatures=25, fileName="scrna")
    '''.format(cv_threshold))
    return palmo_obj


def run(data_filepath, metadata_filepath, datatype, z_cutoff, mean_threshold,
        cv_threshold, na_threshold, housekeeping_genes, feature_set):
    ''' main function 
        1. load datasets 
        2. run variance custom function 
        3. run outlier detection custom function  
    '''
    print('datatype is..{}'.format(datatype))
    palmo_obj = load_data(data_filepath, metadata_filepath, datatype)
    palmo_obj = prep_data(palmo_obj, datatype, na_threshold)
    if datatype == 'bulk':
        palmo_obj = run_lmeVariance(palmo_obj, mean_threshold, feature_set,
                                    datatype)
        palmo_obj = run_cvCalc(palmo_obj,
                               meanThreshold=mean_threshold,
                               cvThreshold=cv_threshold,
                               datatype=datatype)
        palmo_obj = run_outlier_detection(palmo_obj, z_cutoff)

        # extract curated data.frames and results

        datamatrix = rpy_2py(ro.r("palmo_obj@curated$data"))
        metadata = rpy_2py(ro.r("palmo_obj@curated$ann"))
        decomp_var_df = rpy_2py(
            ro.r("palmo_obj@result$variance_decomposition"))
        cv_all = rpy_2py(ro.r("palmo_obj@result$cv_all"))
        outlier_res = rpy_2py(ro.r("palmo_obj@result[['outlier_res']]"))
        """
        workdir = '/Users/james.harvey/workplace/bulk'
        datamatrix.to_csv('{}/data.csv'.format(workdir))
        decomp_var_df.to_csv('{}/var_decomp.csv'.format(workdir))
        cv_all.to_csv('{}/cv_res.csv'.format(workdir))
        outlier_res.to_csv('{}/outlier.csv'.format(workdir))
        """
        return (datamatrix, metadata, decomp_var_df, cv_all, outlier_res)
    if datatype == 'singlecell':
        """
        work_dir = '{}/data'.format(os.getcwd())
        data = pd.read_csv('{}/data.csv'.format(work_dir))
        metadata = pd.read_csv('{}/metadata.csv'.format(work_dir))
        var_decomp = pd.read_csv('{}/var_decomp.csv'.format(work_dir))
        var_genes = pd.read_csv('{}/var_genes.csv'.format(work_dir))
        var_genes['is_var_gene'] = 1
        stable_genes = pd.read_csv('{}/non_var_genes.csv'.format(work_dir))
        stable_genes['is_var_gene'] = 0
        cv_res = pd.concat([var_genes, stable_genes])
        outlier_res = pd.DataFrame()
        """
        palmo_obj = run_cvcalc_scprofile(palmo_obj, mean_threshold,
                                         housekeeping_genes)
        palmo_obj = run_lmeVariance(palmo_obj, mean_threshold, feature_set,
                                    datatype)
        palmo_obj = run_cvCalc(palmo_obj,
                               datatype=datatype,
                               meanThreshold=mean_threshold,
                               cvThreshold=cv_threshold)
        palmo_obj = find_variable_features(palmo_obj, cv_threshold)
        palmo_obj = find_stable_features(palmo_obj, cv_threshold)
        # extract cureated data.frames and results
        data = rpy_2py(ro.r("palmo_obj@curated$data"))
        metadata = rpy_2py(ro.r("palmo_obj@curated$anndata"))
        cv_res = rpy_2py(ro.r("palmo_obj@result$cv_meanthreshold"))
        var_decomp = rpy_2py(ro.r("palmo_obj@result$variance_decomposition"))
        var_genes = rpy_2py(ro.r("palmo_obj@result$variable_gene"))
        stable_genes = rpy_2py(ro.r("palmo_obj@result$non_variable_gene"))
        cv_res = pd.concat([var_genes, stable_genes])
        outlier_res = pd.DataFrame()

        return (data, metadata, var_decomp, cv_res, outlier_res)


if __name__ == "__main__":
    import sys
    datamatrix_filepath = sys.argv[1]
    metadata_filepath = sys.argv[2]
    datatype = sys.argv[3]
    z_cutoff = float(sys.argv[4])
    mean_threshold = float(sys.argv[5])
    cv_threshold = float(sys.argv[6])
    na_threshold = float(sys.argv[7])
    housekeeping_genes = sys.argv[8]
    feature_set = sys.argv[9]
    run(datamatrix_filepath=datamatrix_filepath,
        metadata_filepath=metadata_filepath,
        datatype=datatype,
        z_cutoff=z_cutoff,
        mean_threshold=mean_threshold,
        cv_threshold=cv_threshold,
        na_threshold=na_threshold,
        housekeeping_genes=housekeeping_genes,
        feature_set=feature_set)
