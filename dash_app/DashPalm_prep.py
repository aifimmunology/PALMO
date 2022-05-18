## Libraries

import rpy2.robjects.packages as rpackages
import pandas as pd
import numpy as np

# working with dataframes
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.vectors import StrVector

# TODO: suppress system warnings from imports
from rpy2.robjects.packages import importr

utils = importr('utils')
# getting ready to download these packages...
importr('PALM')
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
    if datatype == 'bulk'
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
            palmo_obj <- avgExpCalc(data_object=palmo_obj, assay="RNA", group_column="{dt}")
            palmo_obj <- checkReplicates(data_object=palmo_obj, mergeReplicates=T)
        '''.format(dt=datatype))
    return palmo_obj


def run_lmeVariance(palmo_obj, mean_threshold, feature_set, output_dir):
    '''
    '''
    # load lmeVariance object function and run
    palmo_obj = ro.r('''
        featureSet <- c("{f1}", "{f2}")
        palmo_obj <- lmeVariance(data_object=palmo_obj, featureSet=featureSet, meanThreshold={mean}, fileName="olink")
        '''.format(f1=feature_set[0], f2=feature_set[1], mean=mean_threshold))
    return palmo_obj


def run_outlier_detection(palmo_obj, z_cutoff, output_dir):
    '''
    '''
    # robject of outlierDetect function
    palmo_obj = ro.r('''
        palmo_obj <- outlierDetect(data_object=palmo_obj, z_cutoff={})
        
        # assign cutoffs 
        '''.format(z_cutoff))
    return palmo_obj


def run_cvCalcBulk(palmo_obj, meanThreshold, cvThreshold):
    """ Intro-donor variations over time 
        TODO: label top 10 genes, and most stable ones 
    """

    palmo_obj = ro.r('''
        palmo_obj <- cvCalcBulk(data_object=palmo_obj, meanThreshold={m}, cvThreshold={c})
        '''.format(m=meanThreshold, c=cvThreshold))
    return palmo_obj


def run(data_filepath,
        metadata_filepath,
        datatype,
        do_outlier,
        z_cutoff,
        mean_threshold,
        cv_threshold,
        na_threshold,
        housekeeping_genes,
        feature_set,
        output_dir=None):
    ''' main function 
        1. load datasets 
        2. run variance custom function 
        3. run outlier detection custom function  
    '''
    palmo_obj = load_data(data_filepath, metadata_filepath, datatype)
    palmo_obj = prep_data(palmo_obj, datatype, na_threshold)
    if datatype == 'bulk':
        palmo_obj = run_lmeVariance(palmo_obj, mean_threshold, feature_set,
                                    output_dir)
        palmo_obj = run_cvCalcBulk(palmo_obj,
                                   meanThreshold=mean_threshold,
                                   cvThreshold=cv_threshold)
        palmo_obj = run_outlier_detection(palmo_obj, z_cutoff, output_dir)

        # extract curated data.frames and results
        datamatrix = rpy_2py(ro.r("palmo_obj@curated$data"))
        metadata = rpy_2py(ro.r("palmo_obj@curated$ann"))
        decomp_var_df = rpy_2py(
            ro.r("palmo_obj@result$variance_decomposition"))
        cv_all = rpy_2py(ro.r("palmo_obj@result$cv_all"))
        outlier_res = rpy_2py(ro.r("palmo_obj@result[['outlier_res']]"))
        #decomp_var_df.to_csv(
        #    '/Users/james.harvey/workplace/variance_input.csv')
        return (datamatrix, metadata, decomp_var_df, cv_all, outlier_res)
    elif datatype == 'singlecell':
        palmo_obj = run_cvcalc_scprofile(palmo_obj, housekeeping_genes)

if __name__ == "__main__":
    import sys
    datamatrix_filepath = sys.argv[1]
    metadata_filepath = sys.argv[2]
    datatype = sys.argv[3]
    do_outlier = bool(sys.argv[4])
    z_cutoff = float(sys.argv[5])
    mean_threshold = float(sys.argv[6])
    cv_threshold = float(sys.argv[7])
    na_threshold = float(sys.argv[8])
    housekeeping_genes = list(sys.argv[9])
    output_dir = sys.argv[10]
    feature_set = sys.argv[11]
    import pdb; pdb.set_trace()
    run(datamatrix_filepath=datamatrix_filepath,
        metadata_filepath=metadata_filepath,
        datatype=datatype,
        do_outlier=do_outlier,
        z_cutoff=z_cutoff,
        mean_threshold=mean_threshold,
        cv_threshold=cv_threshold,
        na_threshold=na_threshold,
        output_dir=output_dir,
        feature_set=feature_set)
