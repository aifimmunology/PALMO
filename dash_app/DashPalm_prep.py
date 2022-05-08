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


def sample_overlap(mat, meta):
    ''' overlap the datamtrix with emtadata for selecting available samples only 
    '''
    # create bridging robject
    r = ro.r

    # load robject functions - intersect & colnames
    r.intersect = ro.r['intersect']
    r.colnames = ro.r['colnames']
    overlap = r.intersect(meta[0], r.colnames(mat))

    # convert rpy object to py object and overlap data
    mat_py = rpy_2py(mat)
    meta_py = rpy_2py(meta)
    mat_py = mat_py.loc[:, list(overlap)]
    meta_py = meta_py.loc[meta_py['Sample'].isin(list(overlap)), ]
    meta_py.index = ro.r['ann'][0]
    return (mat_py, meta_py)


def apply_na_threshold(data, threshold=0.4):
    '''
    '''
    # calculate how many NA's in each row
    data['number_NA'] = data.isna().sum(axis=1)
    filtered_data = data.loc[data['number_NA'] <= len(data.columns) *
                             threshold, ]
    filtered_data.drop(columns=['number_NA'], inplace=True)
    return filtered_data


def load_data(matrix_fp, meta_fp, na_threshold):
    ''' 
    '''
    # create bridging robject
    r = ro.r

    # loads robjects named "data"  and "ann"
    datamatrix_r = ro.r['load'](matrix_fp)
    datamatrix_r = ro.r['data']  # users can choose data here
    metadata_r = ro.r['load'](meta_fp)
    metadata_r = ro.r['ann']  # users could choose data here

    datamatrix, metadata = sample_overlap(datamatrix_r, metadata_r)

    # TODO: should this cut-off be added to dash app UI?
    # DEFAULT - 0.4 of number of columns
    datamatrix = apply_na_threshold(datamatrix, na_threshold)

    return (datamatrix, metadata)


def run_lmeVariance(matrix, metadata, mean_threshold, feature_set, output_dir):
    '''
    '''
    # load lmeVariance object function and run
    lmeVariance = ro.r['lmeVariance']

    # convert to R StrVector
    featureSet = ro.StrVector(feature_set)
    lmem_res = rpy_2py(
        lmeVariance(ann=py2_rpy(metadata),
                    mat=py2_rpy(matrix),
                    featureSet=featureSet,
                    meanThreshold=mean_threshold,
                    filePATH=output_dir))
    mean_vals = lmem_res.loc[:, 'mean']
    lmem_res['mean_log10'] = np.log10(lmem_res['mean'] + 1)
    return lmem_res


def run_outlier_detection(metadata, datamatrix, z_cutoff, output_dir):
    '''
    '''
    # robject of outlierDetect function
    outlierDetect = ro.r['outlierDetect']
    outlier_res = outlierDetect(ann=py2_rpy(metadata),
                                mat=py2_rpy(datamatrix),
                                z_cutoff=z_cutoff,
                                filePATH=output_dir)
    outlier_respy = rpy_2py(outlier_res)

    # add categorical variables
    outlier_respy.loc[outlier_respy['z'] > 0, 'zgroup'] = "> Z"
    outlier_respy.loc[outlier_respy['z'] < 0, 'zgroup'] = "< -Z"

    outlier_respy['z_cutoff'] = z_cutoff
    return outlier_respy


def gen_dash_app():
    '''
    '''
    pass


def run_cvCalcBulk(mat, ann, meanThreshold, cvThreshold):
    """ Intro-donor variations over time 
        TODO: label top 10 genes, and most stable ones 
    """

    unigene = mat.index.tolist()
    uniSample = ann['PTID'].unique().tolist()

    variable_gene = pd.DataFrame()
    stable_gene = pd.DataFrame()

    for i in list(range(0, len(uniSample))):  #len(uniSample))):
        uS = uniSample[i]
        meta_df = ann.loc[ann['PTID'].eq(uS), ]
        if len(meta_df) > 1:
            samples = meta_df['Sample'].unique().tolist()

            df = mat.loc[mat.index.isin(unigene), samples]

            # column of how many NAs exist
            df['NAs'] = df.isna().sum(axis=1)

            # mean column (rowwise calc)
            df['mean'] = df[samples].mean(axis=1)

            # row-wise sd
            df['sd'] = df[samples].std(axis=1)

            # row-wise var
            df['var'] = df[samples].var(axis=1)

            df['CV'] = (100 * df['sd']) / df['mean']

            cutoff = 0.5 * len(meta_df)
            dp1 = df.copy(deep=True)
            dp1 = dp1.loc[dp1['NAs'] <= cutoff, ]
            dp1 = dp1.loc[np.abs(dp1['mean']) >= meanThreshold, ]

            dp2a = dp1.loc[np.abs(dp1['CV']) > cvThreshold, ]
            if len(dp2a) > 0:
                # subset to just the first 10 labels
                dp2a.sort_values(by=['CV', 'mean'],
                                 key=abs,
                                 inplace=True,
                                 ascending=False)
                if len(dp2a) > 10:
                    dp2a_label = dp2a[0:10]
                variable_gene = pd.concat(
                    [variable_gene, dp2a[['mean', 'sd', 'var', 'CV']]], axis=1)
                variable_gene['donor'] = uS
                variable_gene['gene'] = variable_gene.index

            dp2b = dp1.loc[np.abs(dp1['CV']) <= cvThreshold, ]
            if len(dp2b) > 0:
                dp2b.sort_values(by=['mean', 'CV'],
                                 key=abs,
                                 inplace=True,
                                 ascending=False)
                stable_gene = pd.concat(
                    [stable_gene, dp2b[['mean', 'sd', 'var', 'CV']]], axis=1)
                stable_gene['donor'] = uS
                stable_gene['gene'] = stable_gene.index

            temp = pd.DataFrame(data={'var': df['CV']})
            temp['gene'] = temp.index
            temp.rename(columns={'var': '{}_var'.format(uS)}, inplace=True)
            if (i == 0):
                res = temp.copy(deep=True)
            else:
                res = pd.merge(res, temp, on='gene')

    # remove blanks
    variable_gene = variable_gene.loc[~variable_gene['donor'].isna(), ]
    stable_gene = stable_gene.loc[~stable_gene['donor'].isna(), ]

    res_var = res.copy(deep=True)
    res_var.index = res_var['gene'].values
    res_var.columns = res_var.columns.str.rstrip('_var')

    uniq_sample = list(
        set(uniSample).intersection(set(res_var.columns.tolist())))

    res_var['Mean'] = res_var[uniq_sample].mean(axis=1)
    res_var['Median'] = res_var[uniq_sample].median(axis=1)
    return res_var


def run(datamatrix_filepath,
        metadata_filepath,
        datatype,
        do_outlier,
        z_cutoff,
        mean_threshold,
        cv_threshold,
        na_threshold,
        feature_set,
        output_dir=None):
    ''' main function 
        1. load datasets 
        2. run variance custom function 
        3. run outlier detection custom function  
    '''
    if datatype == 'bulk':
        datamatrix, metadata = load_data(datamatrix_filepath,
                                         metadata_filepath, na_threshold)
        lmem_df = run_lmeVariance(datamatrix, metadata, mean_threshold,
                                  feature_set, output_dir)
        cv_res = run_cvCalcBulk(mat=datamatrix,
                                ann=metadata,
                                meanThreshold=mean_threshold,
                                cvThreshold=cv_threshold)
        outlier_res_py = run_outlier_detection(metadata, datamatrix, z_cutoff,
                                               output_dir)
        return (datamatrix, metadata, lmem_df, cv_res, outlier_res_py)
    elif datatype == 'singlecell' and omics == 'rna':
        pass


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
    output_dir = sys.argv[9]
    feature_set = sys.argv[10]
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
