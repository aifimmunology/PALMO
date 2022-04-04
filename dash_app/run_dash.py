''' run_dash.py

'''

# load libraries 
from jupyter_dash import JupyterDash
import hisepy as hp
import os
import dash
from dash import html, dcc
from dash import dash_table
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from dash.dependencies import Input, Output, State
import pandas as pd
import socket
import diskcache
from dash.long_callback import DiskcacheLongCallbackManager

# additional libraries that wasn't handled properly in prep code 
from plotly.subplots import make_subplots 
import rpy2.robjects as robjects 
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import numpy as np 

# import prep module 
import DashPalm_prep as dpp 

# argument handler 
import argparse 
import pathlib

cache = diskcache.Cache("./cache")
long_callback_manager = DiskcacheLongCallbackManager(cache)

# bootstrap - style sheet 
external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css']
app = JupyterDash(__name__, long_callback_manager=long_callback_manager, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)
server = app.server

# default aesthetics for all plots 
pio.templates.default = "plotly_white"
colors = {
    'aiwhite': 'rgb(64, 85, 100)',
    'aiblue': 'rgb(0, 48, 87)',
    'aigreen': 'rgb(120, 156, 74)',
    'aiorange': 'rgb(238, 118, 35)',
    'aidarkorange': 'rgb(176, 101, 51)',
    'aigrey': 'rgb(124, 125, 127)',
    'ailightgrey': 'rgb(231, 230, 230)',
}


cvThreshold = 5 
meanThreshold = 1
######################
# --- Arg Parser --- # 
######################
""" other arguments 
parser.add_argument('--omics', '-o', nargs='?', default="plasmaproteome"
                    help='User defined name like RNA, ATAC, Proteomics, FLOW') 
parser.add_argument('--filename', '-fn', nargs='?', 
                    help='User defined filename') 
parser.add_argument('--column_sep', '-cs', nargs='?',
                   help='Separator of "PTID" and "Time" in "Sample" column of Annotation table like column_sep="W" for PTID1W1, column_sep=":" for PTID1W1:Tcell)'
parser.add_argument('--cluster_by', '-cb', nargs='?',
                   help='for sample correlation cluster columns by ("donor", "group")')
parser.add_argument('--coding_genes', '-cg', nargs='?',
                   help='Selecting protein coding/user-defined gene list only') 
parser.add_argument('--avg_group', '-ag', nargs='?',
                   help='Group label to be used to calculate average gene expression by group label') 
parser.add_argument('--housekeeping_genes', '-hg', nargs='?',
                   help='Optional list of housekeeping genes to focus on Default is NULL)') 
parser.add_argument('--group_oi', '-go', nargs='?',
                   help='Group of interest to focus on, Default is NULL') 
parser.add_argument('--
"""
def parse_args(): 
    
    parser = argparse.ArgumentParser() 

    # create arguments 
    parser.add_argument('--metadata', '-m',  required=True, 
                        help='Annotation table. Table must consist column Sample (Participant sample name), PTID (Participant), Time (longitudinal timepoints)') 
    parser.add_argument('--data', '-d', required=True, 
                        help='Expression matrix or data frame. Rows represents gene/proteins column represents participant samples (same as annotation table Sample column). For single cell, Single cell RNA Seurat object, if datatype is single cell RNA and Single cell ATAC genescore matrix or data frame')
    parser.add_argument('--datatype', '-dt', nargs='?', default='bulk',
                        choices=['bulk', 'singlecell'],
                        help='Data input can be bulk or singlecell')
    parser.add_argument('--do_outlier', '-out', type=bool, 
                        help='Specify whether to perform outlier analysis')
    parser.add_argument('--z_cutoff', '-zo', type=float, 
                        help='|z| cutoff threshold to find potential outliers')
    parser.add_argument('--mean_threshold', '-mt', nargs='?', type=float,
                        help='Average expression threshold to filter lowly expressed') 
    parser.add_argument('--cv_threshold', '-cvt', nargs='?',  type=float,
                        help='Coefficient of variation threshold to select variable and stable genes Default is 10 for single cell RNA (100*SD/mean)')
    parser.add_argument('--na_threshold', '-nat', nargs='?', type=float, default=0.4,
                        help='Number of NAs in data (numeric value or NULL). Default, 40% * number of columns.')
    parser.add_argument('--output_dir', '-o', 
                        help='user-defined output directory') 
    parser.add_argument('--feature_set', '-fs', nargs='*', type=str, default=["PTID", "Time"],
                        help='Variance analysis carried out on the featureSet provided such as c("PTID", "Time", "Sex")')
    parser.add_argument('--method', type=str, default='spearman', 
                        help='Sample correlation analysis') 
    return parser.parse_args() 


########################## 
# --- Variance Plots --- #
########################## 

def prep_var_contribute_df(tbl : pd.DataFrame): 
    ''' Takes table of residuals and does some formatting and subsetting
    '''
    rmelt = robjects.r['melt']
    rdata_matrix = robjects.r['data.matrix']
    featureList = ['PTID', 'Time'] + ['Residual'] # args.here
    # meanThreshold = meanThreshold 
    tbl = tbl.loc[tbl['max'] > meanThreshold, featureList]
    res_tbl = tbl.rename(columns={'PTID' : 'donor', 'Time':'week', 'Residual':'Residuals'})
    res_tbl = res_tbl * 100 
    df1 = res_tbl.loc[(res_tbl['donor'] > res_tbl['week']) & (res_tbl['Residuals'] < 50), ]
    df1.sort_values(by='donor', ascending=False, inplace=True) 

    # convert object and perform R's melt, and data.matrix functions 
    with localconverter(ro.default_converter + pandas2ri.converter):
        df1_r = ro.conversion.py2rpy(df1)
    var_df_r = rmelt(rdata_matrix(df1_r))

    # convert back to pandas data.frame 
    with localconverter(ro.default_converter + pandas2ri.converter):
        return(ro.conversion.rpy2py(var_df_r)) 

@app.callback(
    Output('variance-out', 'data'),
    Input('feature-button', 'n_clicks'),
    [State('dropdown-select', 'value')]
)
def subset_var_df(n, var_chosen): 
    '''
    '''
    var_df = prep_var_contribute_df(lmem_py) 
    var_df.sort_values(by='value', ascending=False, inplace=True) 
    if n == 0: 
        top_vars = var_df[0:15]['Var1'].unique().tolist()
        dff = var_df[var_df['Var1'].isin(top_vars)] 
    elif n > 0: 
        dff = var_df[var_df['Var1'].isin(var_chosen)]
    return dff.to_json(date_format='iso', orient='split') 


@app.callback(
    Output('variance-contribution', 'figure'),
    Input('variance-out', 'data'),
    suppress_callback_exceptions=True,
    prevent_initial_callback=False
)
def var_contribute_plot(dff): 
    ''' separate function that renames and subsets data 
    '''
    dff = pd.read_json(dff, orient='split') 
    
    # sort by donor # 
    donors_sorted = dff.loc[dff['Var2'].eq('donor'), ]
    donors_sorted.sort_values(by='value', ascending=False, inplace=True) 
    # dff.sort_values(by='value', ascending=True, inplace=True) 
    var_fig = px.bar(dff, x='value', y='Var1', color='Var2', orientation='h', 
                     category_orders = {'Var1': donors_sorted['Var1'].tolist()},
                     labels= {'Var1' : 'Features',
                            'value' : '% variance explained'}
                    ) 
    var_fig.update_traces().update_layout(title_x=0.5, legend_title_text='FeatureList') 
    return (var_fig)


def violin_box_plot(): 
    ''' creates violin with box plot overlayed 
    TODO: plot sigFeature datapoints 
    '''
    featureList = ['PTID', 'Time'] + ['Residual'] 
    # meanThreshold = meanThreshold 
    rmelt = robjects.r['melt']
    rdata_matrix = robjects.r['data.matrix']
    res = lmem_py.loc[lmem_py['max'] > meanThreshold, featureList]

    # convert object to r 
    with localconverter(ro.default_converter + pandas2ri.converter):
        res_r = ro.conversion.py2rpy(res)
    df = rmelt(rdata_matrix(res_r))

    with localconverter(ro.default_converter + pandas2ri.converter):
        df_py = ro.conversion.rpy2py(df)

    df_py['value'] = df_py['value'] * 100 

    df_py['feature'] = df_py[['Var1', 'Var2']].agg('_'.join, axis=1)
    #df_py.to_csv('/home/jupyter/dash_var.csv') 
    
    # list features 
    feature_list = df_py['Var2'].unique().tolist() 
    fig= go.Figure()
    for f in feature_list: 
        fig.add_trace(go.Violin(x=df_py.loc[df_py['Var2'].eq(f), 'Var2'],
                                y = df_py.loc[df_py['Var2'].eq(f), 'value'],
                                name = f,
                                box_visible=True, 
                                opacity=0.6,
                                spanmode='hard',
                                showlegend=True
                               )) 
    fig.update_traces(width=0.5) 
    fig.update_layout(legend_title_text='featureList')
    
    return fig 


@app.callback(
    Output('download-var-df','data'),
    Input('var_btn', 'n_clicks'), 
    State('variance-out', 'data')
) 
def download_var_test(n_clicks, dff):
    '''
    '''    
    if n_clicks != None: 
        dff = pd.read_json(dff, orient='split') 
        # create json blob that will hold these values
        return dcc.send_data_frame(dff.to_csv, 'variance_viz_data.csv')
    
    
######################### 
# --- Outlier Plots --- #
######################### 

@app.callback(
    Output('download-outlier-df', 'data'), 
    Input('outlier_dl_button', 'n_clicks'), 
    State('outlier-input', 'data')
)
def download_outlier_data(n_clicks, dff): 
    '''
    '''
    if n_clicks != None: 
        dff = pd.read_json(dff, orient='split') 
        return dcc.send_data_frame(dff.to_csv, 'outlier_viz_data.csv') 
    
    
def prep_outlierP_data(outlier_df, z_cutoff, z_score_subset, nGenes, groupby): 
    '''
    '''
    
    # load custom function to calcualte p-values 
    p_value_for_event = robjects.r['p_value_for_event']
    rpnorm = robjects.r['pnorm']
    rate = 1 - rpnorm(z_cutoff)[0]
    if z_score_subset == 'above':
        data = outlier_df.loc[outlier_df['z'] >= z_cutoff, ] 
    elif z_score_subset == 'below': 
        data = outlier_df.loc[outlier_df['z'] <= -z_cutoff, ]
    else: 
        data = outlier_df.loc[np.abs(outlier_df['z']) >= z_cutoff, ]
        rate = 2*(1-rpnorm(z_cutoff)[0])
        
    # create a data.frame with freq #'s 
    # first, create a table and convert to a data.frame
    tab = data.groupby([groupby, 'Time']).size()
    df = tab.unstack().reset_index()
    df = pd.wide_to_long(df, stubnames="W", i='PTID', j='Time').reset_index() 

    # now do some cleaning; and fill NAs
    df['tmp'] = 'W'
    df['Time'] = df['Time'].astype(str) 
    df['Time'] = df[['tmp','Time']].agg(''.join, axis=1) 
    df.fillna(value=0, inplace=True)

    # rename and drop columns 
    df.rename(columns={groupby:'Sample', 'W':'Freq'}, inplace=True) 
    df.drop(columns='tmp', inplace=True) 
    df['pvals'] = p_value_for_event(events=df['Freq'].tolist(), tries=nGenes, rate=rate) 
    df['signP'] = -np.log10(df['pvals'])

    # create an id column 
    df['id'] = df[['Sample','Time']].agg(''.join, axis=1) 
    df.sort_values(by=['Time','Sample'], inplace=True) 
    
    return df 


@app.callback( 
    Output('outlier-input', 'data'),
    Input('z-submit','n_clicks')
)
def run_and_cache_outlier_input(n_clicks): 
    '''
        TODO: add some signaling to all methods that utilize this as input 
        TODO: cache results by just appending to this df 
    '''
    if n_clicks == 0: 
        return outlier_res_py.to_json(date_format='iso', orient='split') 
    else: 
        return outlier_res_py.to_json(date_format='iso', orient='split') 




## TODO: parameter: featureSet (mp) 
## TODO: handle different ann and matrices (lp)
## TODO: better color aesthetics 
## TODO: create a "select plot" dropdown
@app.callback(
    Output('number-features-bar', 'figure'), 
    Input('outlier-input', 'data'),
    [State('sample-selector', 'value')]
) 
def plot_no_features(df, select_samples): 
    '''
    '''
    # create column for number of features 
    nofeatures_df = pd.read_json(df, orient='split') 
    nofeatures_df['no_features'] = 1
    nofeatures_df = nofeatures_df.groupby(['Sample', 'zgroup'])['no_features'].sum().reset_index() 
    
    if select_samples != None: 
        nofeatures_df = nofeatures_df.loc[nofeatures_df['Sample'].isin(select_samples), ]
    nofeat = px.bar(nofeatures_df, x='Sample', y='no_features', color='zgroup', width=500, barmode='group', text='no_features')
    nofeat.update_traces(texttemplate='%{text:.2s}', textposition='outside')
    return nofeat



@app.callback(
    Output('p-scatter', 'figure'), 
    Input('outlier-input', 'data'),
    Input('z-score-cutoff','value'),
    [State('sample-selector', 'value')]
)
def outlier_p_scatter_plot(data_input, z_subset, select_samples): 
    '''
    '''
    d_input = pd.read_json(data_input, orient='split') 
    
    # if none, we can assume the df is from the initial callback. which should only have a single value under z cut-off column
    input_df = prep_outlierP_data(d_input, z_cutoff=2, z_score_subset=z_subset, nGenes=1042, groupby="PTID")
    if select_samples != None: 
        input_df = input_df.loc[input_df['id'].isin(select_samples), ]
    fig = px.scatter(input_df, x='Freq', y='signP', hover_name='id',
                 labels={'Freq': '# Features', 
                         'signP' : '-log10(p-value)'
                        }) 
    return fig


@app.callback(
    Output('p-bar', 'figure'),
    Input('outlier-input', 'data'),
    Input('z-score-cutoff', 'value'),
    [State('sample-selector', 'value')]
)
def outlier_p_bar_plot(data_input, z_subset,  select_samples): 
    '''
    '''
    d_input = pd.read_json(data_input, orient='split') 
    # if none, we can assume the df is from the initial callback. which should only have a single value under z cut-off column
 
    input_df = prep_outlierP_data(d_input, z_cutoff=2, z_score_subset=z_subset, nGenes=1042, groupby="PTID")
    input_df['sample_time'] =  input_df[['Sample','Time']].agg(''.join, axis=1) 
    input_df.sort_values(by=['sample_time'],inplace=True) 
    
    if select_samples != None: 
        input_df = input_df.loc[input_df['sample_time'].isin(select_samples), ]
    fig = px.bar(input_df, x='sample_time', y='signP',
                labels={'sample_time' : 'Sample', 
                        'signP': '-log10(pvalue)'
                       }
                )
    return fig 


@app.callback(
    Output('outlier-plot', 'figure'), 
    Input('center-measure', 'value'),
    Input('outlier-input', 'data'),
    [State('sample-selector', 'value')]
)
def outlier_detect_plot(measure_col, data_input, these_samples): 
    ''' 
    TODO: order the samples correctly 
    '''
    # sort values 
    d_input = pd.read_json(data_input, orient='split') 
    d_input.sort_values(by='Sample', inplace=True)
    if these_samples != None: 
         d_input = d_input.loc[d_input['Sample'].isin(these_samples), ]
    
    out_fig = go.Figure() 
    z_group_list = list(d_input['zgroup'].unique())
    z_group_list.sort() 
    for this_group in z_group_list: 
        violin_input = d_input.loc[d_input['zgroup'].eq(this_group), ].sort_values(by='Sample') 
        
        out_fig.add_trace(go.Violin(x=violin_input['Sample'],
                                    y=violin_input[measure_col],
                                    legendgroup=this_group,
                                    name=this_group,
                                    box_visible=True,
                                    spanmode='hard',
                                    opacity=0.6,

                 ))
    
    #out_fig = px.violin(d_input, x='Sample', y=measure_col, color='zgroup', box=True) 
    out_fig.update_traces(width=0.5) 
    sample_order = d_input['Sample'].unique().tolist() 
    sample_order.sort()
    out_fig.update_xaxes(categoryorder='array', categoryarray=sample_order) 
    # out_fig.update_layout(legend_traceorder="reversed")
    return out_fig 



##################################
# --- Expression Level Plots --- # 
################################## 

@app.callback(
    Output('download-expression-df', 'data'), 
    Input('expression_dl_button', 'n_clicks'), 
    State('expression-out', 'data')
)
def download_outlier_data(n_clicks, dff): 
    '''
    '''
    if n_clicks != None: 
        dff = pd.read_json(dff, orient='split') 
        return dcc.send_data_frame(dff.to_csv, 'expression_lvl_viz_data.csv') 
    
@app.callback(
    Output('expression-out', 'data'), 
    Input('this-gene', 'value')
)
def store_gene_data(this_gene): 
    ''' Caches and stores data for callbacks 
    '''
    this_gene_vals = datamatrix[datamatrix.index == this_gene].values
    ann = metadata.copy(deep=True) 
    ann['exp'] = list(this_gene_vals[0]) 
    return ann.to_json(date_format='iso', orient='split') 


@app.callback(
    Output('var-gene-plot', 'figure'),
    Input('expression-out', 'data')
)
def var_gene_plot(df): 
    '''
    '''
    gene_df = pd.read_json(df, orient='split')     
    gfig = px.scatter(gene_df, x='PTID', y='exp', symbol='Time', color='Time')
    boxf = px.box(gene_df, x='PTID', y='exp')
    combined_plot = go.Figure(data=boxf.data + gfig.data)
    return combined_plot 


@app.callback( 
    Output('geneplot', 'figure'),
    Input('expression-out', 'data') 
)
def geneplot(df): 
    '''
    '''
    ann = pd.read_json(df, orient='split') 
    fig = px.scatter(ann, x='PTID', y='exp', color='Time',  facet_col='Time')
    fig.for_each_annotation(lambda a: a.update(text="")) # remove annotation from each subplot 
    return fig 


##############################
# --- correlation matrix --- #
##############################
@app.callback(
    Output('corr-matrix', 'data'), 
    Output('corr-fig', 'figure'),
    Input('corr-dropdown', 'value')
)
def correlation_matrix_plot(correlation_method='spearman'): 
    cols = datamatrix.columns.tolist() 
    cor_res = datamatrix.corr(method=correlation_method) 
    data_input = cor_res[cols].to_numpy().transpose()
    fig = go.Figure(go.Heatmap(
            z=data_input, 
            colorscale=[[0, 'white'], [0.5, 'rgb(250, 138, 130)'], [1.0, 'rgb(255, 17, 0)']],
            x=cols,
            y=cols,
            dx=1,
            dy=1
        ))
    return (pd.DataFrame(data_input).to_json(orient='split'), fig)


@app.callback(
    Output('download-corr-df','data'),
    Input('corr_btn', 'n_clicks'),
    State('corr-matrix', 'data')
) 
def download_correlation(n_clicks, dff):
    '''
    '''    
    if n_clicks != None: 
        dff = pd.read_json(dff, orient='split') 
        # create json blob that will hold these values
        return dcc.send_data_frame(dff.to_csv, 'correlation_matrix_viz_data.csv')

##################################
# --- Intra-donor variations --- #
##################################

# NOTE: hard-coded vars here 
def prep_cvcalc_bulk(meanThreshold, cvThreshold): 
    ''' Intro-donor variations over time 
    
    TODO: label top 10 genes, and most stable ones 
    '''
    mat=datamatrix
    ann=metadata

    unigene = mat.index.tolist()
    uniSample = ann['PTID'].unique().tolist()

    variable_gene = pd.DataFrame()
    stable_gene = pd.DataFrame() 

    for i in list(range(0,len(uniSample))): #len(uniSample))): 
        uS = uniSample[i]
        meta_df = ann.loc[ann['PTID'].eq(uS),] 
        if len(meta_df) > 1:
            samples = meta_df['Sample'].unique().tolist()

            df = mat.loc[mat.index.isin(unigene), samples] 

            # column of how many NAs exist 
            df['NAs'] =  df.isna().sum(axis=1) 

            # mean column (rowwise calc) 
            df['mean'] = df[samples].mean(axis=1)

            # row-wise sd 
            df['sd'] = df[samples].std(axis=1)

            # row-wise var 
            df['var'] = df[samples].var(axis=1) 

            df['CV'] = (100* df['sd']) / df['mean'] 

            cutoff = 0.5 * len(meta_df) 
            dp1 = df.copy(deep=True) 
            dp1 = dp1.loc[dp1['NAs'] <= cutoff,]
            dp1 = dp1.loc[np.abs(dp1['mean']) >= meanThreshold, ]


            dp2a = dp1.loc[np.abs(dp1['CV']) > cvThreshold, ]
            if len(dp2a) > 0: 
                # subset to just the first 10 labels 
                dp2a.sort_values(by=['CV', 'mean'], key=abs, inplace=True, ascending=False)
                if len(dp2a) > 10: 
                    dp2a_label = dp2a[0:10]
                variable_gene = pd.concat([variable_gene, dp2a[['mean','sd','var','CV']]], axis=1)
                variable_gene['donor'] = uS
                variable_gene['gene'] = variable_gene.index

            dp2b = dp1.loc[np.abs(dp1['CV']) <= cvThreshold,]
            if len(dp2b) > 0: 
                dp2b.sort_values(by=['mean','CV'], key=abs, inplace=True, ascending=False)
                stable_gene = pd.concat([stable_gene, dp2b[['mean','sd','var','CV']]], axis=1)
                stable_gene['donor'] = uS
                stable_gene['gene'] = stable_gene.index

            temp = pd.DataFrame(data={'var':df['CV']})
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

    uniq_sample = list(set(uniSample).intersection(set(res_var.columns.tolist())))

    res_var['Mean'] = res_var[uniq_sample].mean(axis=1)
    res_var['Median'] = res_var[uniq_sample].median(axis=1) 
    return res_var

@app.callback(
    Output('cv-density', 'figure'), 
    Input('gene-df', 'data') 
)
def cv_density_plot(gene_df): 
    '''
    '''
    dff = pd.read_json(gene_df, orient='split') 
    uniSample = metadata['PTID'].unique().tolist()
    #input_data = prep_cvcalc_bulk(meanThreshold, cvThreshold)
    # uniq_sample = list(set(uniSample).intersection(set(gene_df.columns.tolist())))
    bar_plot_input = pd.wide_to_long(dff, stubnames='PTID', i='gene', j= 'var').reset_index()
    
    fig = px.histogram(bar_plot_input, x='PTID',labels={'PTID':'CV'}).update_layout(yaxis_title='density')
    return fig  


@app.callback(
    Output('cvplot', 'figure'),
    Output('gene-df', 'data'),
    Input('gene-type', 'value')
)
def cvplot(gene_type): 
    ''' here we can either pick the genes of interest for heatmap
        or we could explore by making the mean & cv threshold interactive. (think a slider) 
    '''
    mat_input = prep_cvcalc_bulk(meanThreshold, cvThreshold) # for now, just use the default 
    
    if (gene_type == 'variable'): 
        mat = mat_input.sort_values(by='Median', ascending=False)
        mat = mat.loc[np.abs(mat['Median']) > cvThreshold, ]
        col_scale = [[0, 'rgb(210, 201, 242)'], [0.1, 'rgb(233, 237, 14)'], [0.2, 'rgb(224, 36, 29)'], [1, 'rgb(82, 22, 20)']]
    elif (gene_type == 'stable'): 
        mat = mat_input.sort_values(by='Median', ascending=True) 
        mat = mat.loc[np.abs(mat['Median']) <= cvThreshold, ]
        col_scale = [[0, 'rgb(210, 201, 242)'], [0.3,'rgb(100, 79, 179)'],  [0.7, 'rgb(82, 68, 242)'], [0.9, 'rgb(233, 237, 14)'], [1, 'rgb(82, 22, 20)']]
    
    # subset to first 50 entries then create heat map 
    mat = mat[0:50] 
    heat_fig = go.Heatmap(
        z=mat[['PTID1','PTID2','PTID3','PTID4','PTID5','PTID6']].values,
        x=['PTID1','PTID2','PTID3','PTID4','PTID5','PTID6'],
        y=mat['gene'].unique().tolist(),
        colorscale=col_scale
    )
    d1 = mat['PTID1'].values 
    d2 = mat['PTID2'].values
    d3 = mat['PTID3'].values
    d4 = mat['PTID4'].values 
    d5 = mat['PTID5'].values 
    d6 = mat['PTID6'].values 

    #b1 = annotate_box = go.Figure(go.Box(y=d1))
    b1= go.Box(y=d1)
    b2= go.Box(y=d2) 
    b3= go.Box(y=d3)
    b4= go.Box(y=d4)
    b5= go.Box(y=d5)
    b6= go.Box(y=d6)

    master_fig = make_subplots(rows=2, cols=1,
                               row_heights=[0.2, 0.8],
                               vertical_spacing=0,
                               specs=[[{'type':'box'}],
                                    [{'type':'xy'}]])
    master_fig.add_traces([b1,b2,b3,b4,b5,b6], rows=1, cols=1).update_layout(showlegend=False)
    master_fig.update_xaxes(showticklabels=False, row=1) 
    master_fig.add_traces([heat_fig], rows=2, cols=1) 
    return master_fig, mat.to_json(orient='split')


@app.callback(
    Output('download-cv-gene-df','data'),
    Input('cv_gene_btn', 'n_clicks'),
    State('gene-df', 'data')
) 
def download_cv_gene(n_clicks, dff):
    '''
    '''    
    if n_clicks != None: 
        dff = pd.read_json(dff, orient='split') 
        # create json blob that will hold these values
        return dcc.send_data_frame(dff.to_csv, 'correlation_matrix_viz_data.csv')


############################   
# --- Run Prep & Store --- # 
############################
import DashPalm_prep as dpp 


"""running=[
        (Output("button_id", "disabled"), True, False),
        (Output("cancel_button_id", "disabled"), False, True),
        (
            Output("paragraph_id", "style"),
            {"visibility": "hidden"},
            {"visibility": "visible"},
        ),
        (
            Output("progress_bar_graph", "style"),
            {"visibility": "visible"},
            {"visibility": "hidden"},
        ),
    ],'
"""
@app.long_callback(
    output= [
        Output('input-datamatrix', 'data'), 
        Output('input-metadata', 'data'), 
        Output('input-var', 'data'), 
        Output('input-outlier', 'data')
    ],
    inputs=[
        Input('run_app_btn' ,'n_clicks'),
        State('metadata-entry', 'value'), 
        State('datamatrix-entry', 'value'), 
        State('params-dtype', 'value'), 
        State('params-outlier', 'value'), 
        State('params-z-score', 'value'), 
        State('params-mean', 'value'), 
        State('params-cv', 'value'), 
        State('params-na', 'value'), 
        State('params-features', 'value'), 
        State('params-output', 'value')
    ],
    running=[
        (Output("run_app_btn", "disabled"), True, False)
    ],
    prevent_initial_call=True
)
def parse_params(run_n_clicks, 
                 metadata_fpath, 
                 datamatrix_fpath, 
                 datatype, 
                 run_outlier, 
                 z_score_cutoff, 
                 mean_cutoff, 
                 cv_cutoff, 
                 na_cutoff, 
                 feature_list, 
                 output_dir):
    if run_n_clicks != None: 
        print(run_n_clicks) 
        # TODO: function to check valid parameters 
        print('about to run Dash prep...') 
        global datamatrix, metadata, lmem_py, outlier_res_py 
        (datamatrix, metadata, lmem_py, outlier_res_py) = dpp.run(datamatrix_filepath= datamatrix_fpath, 
                                                                  metadata_filepath= metadata_fpath,
                                                                  datatype= datatype,
                                                                  do_outlier= run_outlier, 
                                                                  z_cutoff= z_score_cutoff,
                                                                  mean_threshold= mean_cutoff, 
                                                                  cv_threshold= cv_cutoff, 
                                                                  na_threshold= na_cutoff, 
                                                                  output_dir= output_dir, 
                                                                  feature_set= ['PTID','Time'])

        # convert objects to json and return 
        return (datamatrix.to_json(orient='split'), metadata.to_json(orient='split'), lmem_py.to_json(orient='split'), outlier_res_py.to_json(orient='split'))
    else: 
        return pd.DataFrame().to_json(), pd.DataFrame().to_json(), pd.DataFrame().to_json(), pd.DataFrame().to_json()

    
####################
# --- Misc/Dev --- #
####################

@app.callback(
    Output('gene-list-options', 'options'),
    Input('input-datamatrix', 'data')
)
def gen_exp_lvl_vals(dm_df): 
    '''
    '''

    try: 
        dff = pd.read_json(dm_df, orient='split') 
        return df.index.unique().tolist() 
    except: 
        print('nope')
        return []

    
    
######################    
# --- Tab handler ---#
######################


def parameters_layout(): 
    '''
    '''
    return html.Div([
        html.H4("Metadata filepath entry", style=header_button_text),     
        dcc.Input(
            value='/home/jupyter/PALM_DASHv2/PALM/data/data_Metadata.Rda',
            id='metadata-entry',
            type='text', 
            style=center_component_style
        ),

        html.H4("data matrix filepath entry", style=header_button_text),
        dcc.Input(
            value='/home/jupyter/PALM_DASHv2/PALM/data/Olink_NPX_log2_Protein.Rda',
            id='datamatrix-entry', 
            type='text',
            style=center_component_style
        ),

        html.H4("Choose datatype", style=header_button_text), 
        dcc.RadioItems(options=['bulk'],
                       value='bulk',
                       id='params-dtype',
                       style=center_component_style),

        html.H4("Choose whether to run outlier analysis", style=header_button_text), 
        dcc.RadioItems(options=['True', 'False'],
                       id='params-outlier', 
                       value='True',
                       style=center_component_style),

        html.H4("Choose z-score cutoff", style=header_button_text),
        dcc.Input( 
            id='params-z-score',
            type='number',
            value=2,
            style=center_component_style
        ),

        html.H4("Choose mean threshold", style=header_button_text), 
        dcc.Input( 
            id='params-mean',
            type='number', 
            value=1,
            style=center_component_style),

        html.H4("Choose CV threshold", style=header_button_text),
        dcc.Input(
            id='params-cv',
            type='number',
            value=5,
            style=center_component_style), 

        html.H4("Choose NA Threshold", style=header_button_text), 
        dcc.Input(
            id='params-na',
            type='number', 
            value=0.4,
            style=center_component_style), 

        html.H4("Choose name of feature set", style=header_button_text), 
        dcc.Input(
            id='params-features', 
            type='text',
            value='PTID Time',
            style=center_component_style), 

        html.H4("Choose where to save outputs", style=header_button_text), 
        dcc.Input(
            id='params-output',
            type='text',
            value='/home/jupyter',
            style=center_component_style),

        html.Button("Run App", id='run_app_btn', style={'align-items' : 'center',
                                                   'width': '200px',
                                                   'display':'grid',
                                                   'justify-content' : 'center',
                                                   'margin' : '0 auto',
                                                   'background' : colors['aigreen']})
    ])



def need_to_run_layout(): 
    '''
    '''
    
    return html.Div([
        html.H4("Please go to the submit parameters tab first and run the app to generate data") 
    ]) 


def correlation_layout(): 
    '''
    '''
    return html.Div([
            html.H4("Choose the correlation coefficient method", style={'textAlign' : 'Center', 
                                                                        'color' : colors['aiblue'], 
                                                                        'padding-top' : '20px',
                                                                        'display':'grid',
                                                                        'justify-content' : 'center'}), 
            dcc.Dropdown(id='corr-dropdown', 
                         options = ['pearson','spearman'],
                         value='spearman',
                         style={'width':'200px',
                                'margin' : '0 auto',
                                'align-items' : 'center'}
                        ),
            html.H4("Download dataset", style={'textAlign' : 'Center', 
                                                     'color' : colors['aiblue'], 
                                                     'padding-top' : '40px',
                                                     'display':'grid',
                                                     'justify-content' : 'center'}), 
            html.Button("Download CSV", id='corr_btn', style={'align-items' : 'center',
                                                               'width': '200px',
                                                               'display':'grid',
                                                               'justify-content' : 'center',
                                                               'margin' : '0 auto',
                                                               'background' : colors['aigreen']}), 
            dcc.Download(id='download-corr-df'), 


            dcc.Graph(id='corr-fig', style={'height': '800px',
                                        'width': '800px',
                                        'margin' : '0 auto'}),
            ])

    
def exp_layout(): 
    '''
    ''' 
    return html.Div([
        html.H4('Pick gene of interest', style={'textAlign' : 'Center', 
                                                'color' : colors['aiblue'],
                                                'padding-top' : '20px'}), 
        dcc.Dropdown(datamatrix.index.unique(),
                     value='ABL1', 
                     id='this-gene',
                     style={'width': '200px', 'margin' : '0 auto'}),  

        html.H4("Download dataset", style={'textAlign' : 'Center', 
                                                 'color' : colors['aiblue'], 
                                                 'padding-top' : '20px',
                                                 'display':'grid',
                                                 'justify-content' : 'center'}
               ),
        html.Button('Download CSV', id='expression_dl_button',
                   style={'align-items' : 'center',
                         'width': '200px',
                         'display':'grid',
                         'justify-content' : 'center',
                         'margin' : '0 auto',
                         'background' : colors['aigreen']}
                   ), 
        dcc.Download(id='download-expression-df'), 

        dcc.Graph(id='var-gene-plot', style={'display': 'inline-block', 'width' : '50%'}),     
        dcc.Graph(id='geneplot', style={'display': 'inline-block', 'width':'50%'})
    ])
  
def var_layout(): 
    '''
    '''
    return html.Div([
        html.H4("Select feature(s) of interest", style={'textAlign' : 'Center', 
                                                               'color' : colors['aiblue'], 
                                                               'padding-top' : '20px',
                                                               'display':'grid',
                                                               'justify-content' : 'center'}), 

        dcc.Dropdown(id = 'dropdown-select',
                     options= prep_var_contribute_df(lmem_py)['Var1'].unique().tolist(), 
                     multi=True),
        html.Button(id='feature-button', n_clicks=0, children='Update feature plot',
                    style={'align-items' : 'center',
                           'width': '200px',
                           'display':'grid',
                           'justify-content' : 'center',
                           'margin' : '0 auto',
                           'background' : colors['aigreen']}),

        html.H4("Download dataset", style={'textAlign' : 'Center', 
                                                 'color' : colors['aiblue'], 
                                                 'padding-top' : '40px',
                                                 'display':'grid',
                                                 'justify-content' : 'center'}), 
        html.Button("Download CSV", id='var_btn', style={'align-items' : 'center',
                                                       'width': '200px',
                                                       'display':'grid',
                                                       'justify-content' : 'center',
                                                       'margin' : '0 auto',
                                                       'background' : colors['aigreen']}), 
        dcc.Download(id='download-var-df'), 


        dcc.Graph(id='variance-contribution'),
        dcc.Graph(figure=violin_box_plot()),

    ])


def intra_var_layout(): 
    '''
    '''
    return html.Div([
        html.H4("Choose to plot variable or stable genes", style={'textAlign' : 'Center', 
                                                                 'color' : colors['aiblue'], 
                                                                 'padding-top' : '40px',
                                                                 'display':'grid',
                                                                 'justify-content' : 'center'}), 
        dcc.RadioItems(['stable','variable'],
                'stable',
                id='gene-type',
                style={'width': '200px', 
                       'margin' : '0 auto',
                       'align-items' : 'center'}),

        html.H4("Download dataset", style={'textAlign' : 'Center', 
                                                 'color' : colors['aiblue'], 
                                                 'padding-top' : '40px',
                                                 'display':'grid',
                                                 'justify-content' : 'center'}), 
        html.Button("Download CSV", id='cv_gene_btn', style={'align-items' : 'center',
                                                           'width': '200px',
                                                           'display':'grid',
                                                           'justify-content' : 'center',
                                                           'margin' : '0 auto',
                                                           'background' : colors['aigreen']}), 
        dcc.Download(id='download-cv-gene-df'), 
        dcc.Graph(id='cvplot', style={'height': '150vh'}),
        dcc.Graph(id='cv-density'),
    ])


def outlier_layout(): 
    '''
    '''
    return html.Div([
        html.H4('Select a sample of interest', style= {'textAlign' : 'Center', 
                                                       'color' : colors['aiblue'], 
                                                       'padding-top' : '20px',
                                                       'display':'grid',
                                                       'justify-content' : 'center'}), 
        dcc.Dropdown(
            id='sample-selector',
            options=outlier_res_py['Sample'].unique().tolist(),
            multi=True, 
            placeholder='Select samples of interest',
            style={'width':'200px',
                   'margin' : '0 auto',
                   'align-items' : 'center'
                  }
        ),


        # insert button
        html.Button('Update Plots', 
                    id='z-submit', 
                    n_clicks=0,
                    style={'width': '200px', 
                           'margin' : '0 auto',
                           'align-items' : 'center',
                           'display':'grid',
                           'background' : colors['aigreen']}
        ),


        html.H4("Download dataset", style={'textAlign' : 'Center', 
                                         'color' : colors['aiblue'], 
                                         'padding-top' : '20px',
                                         'display':'grid',
                                         'justify-content' : 'center'}
        ),
        html.Button('Download CSV', id='outlier_dl_button', 
                    style={'align-items' : 'center',
                           'width': '200px',
                           'display':'grid',
                           'justify-content' : 'center',
                           'margin' : '0 auto',
                           'background' : colors['aigreen']}), 
        dcc.Download(id='download-outlier-df'),
        dcc.RadioItems(
            ['meanDev', 'z'],
            'z',
            id='center-measure',
        ),

        dcc.Store(id='outlier-input'),   
            # TODO: z-score cutoff text box 
            # TODO: update plot button 

        html.Div(children=[

            dcc.Graph(id='outlier-plot', style={'display': 'inline-block', 'width':'50%'}),
            # TODO: option to just graph all groups together 
            dcc.Graph(id='number-features-bar', style={'display': 'inline-block', 'width': '50%'}),
           # dcc.Graph(id='number-features-pbar'),
        ]),

        html.H5("Choose to plot where |Z| is above or below z-cutoff", style={'textAlign' : 'Center', 
                                                   'color' : colors['aiblue'], 
                                                   'padding-top' : '20px',
                                                   'display':'grid',
                                                   'justify-content' : 'center'}),
        dcc.RadioItems(
                ['below', 'above', 'all'],
                'all',
                id='z-score-cutoff',
                style={'width': '200px', 
                       'margin' : '0 auto',
                       'align-items' : 'center'}
        ), 
        html.Div(children=[    
            dcc.Graph(id='p-scatter', style={'display': 'inline-block', 'width': '50%'}),
            dcc.Graph(id='p-bar', style={'display':'inline-block', 'width': '50%' }) 
        ])
    ])


@app.callback(Output('tabs-content', 'children'),
              Input('tabs-graph', 'value'),
              Input('run_app_btn', 'n_clicks')
)
def render_content(tab, run_n_clicks):
    print(run_n_clicks) 
    print('within render func') 
    if run_n_clicks == None: 
        if tab == 'correlation':
            return need_to_run_layout()
        elif tab == 'expression level': 
            return need_to_run_layout()
        elif tab == 'variance': 
            return need_to_run_layout()
        elif tab == 'intra-donor-variation': 
            return need_to_run_layout()
        elif tab == 'outliers': 
            return need_to_run_layout()
    else: 
        if tab == 'correlation':
            return correlation_layout()
        elif tab == 'expression level': 
            return exp_layout() 
        elif tab == 'variance': 
            return var_layout() 
        elif tab == 'intra-donor-variation': 
            return intra_var_layout()
        elif tab == 'outliers': 
            return outlier_layout() 



######################    
# --- App Layout --- #  
###################### 


# style tabs 
tabs_styles = {
    'height': '44px'
}
tab_style = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontWeight': 'bold',
    'color': 'white'
}

tab_selected_style = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': colors['aigreen'],
    'color': 'white',
    'padding': '6px'
}

header_button_text = {'textAlign' : 'Center', 
                       'color' : colors['aiblue'], 
                       'padding-top' : '20px',
                       'display':'grid',
                       'justify-content' : 'center'}

center_component_style = {'width':'200px',
                          'margin' : '0 auto',
                          'align-items' : 'center',
                          'display':'grid',}


def run_app(): 
    app.layout = html.Div(children=[
        html.H2('Platform for Analyzing Longitudinal Multi-omics data',
                style={
                    'textAlign': 'Center',
                    'color': colors['aigreen']
                }),
        dcc.Tabs(id='tabs-graph', value='parameters', children=[

            # landing page 
            # we need to render this page at initial loading to gain access to the run_var_btn id 
            dcc.Tab(label='Submit Parameters', value='parameters', style=tab_style, selected_style=tab_selected_style, children=[
                 html.Div([
                    html.H4("Metadata filepath entry", style=header_button_text),     
                    dcc.Input(
                        value='/home/jupyter/PALM_DASHv2/PALM/data/data_Metadata.Rda',
                        id='metadata-entry',
                        type='text', 
                        style=center_component_style
                    ),

                    html.H4("data matrix filepath entry", style=header_button_text),
                    dcc.Input(
                        value='/home/jupyter/PALM_DASHv2/PALM/data/Olink_NPX_log2_Protein.Rda',
                        id='datamatrix-entry', 
                        type='text',
                        style=center_component_style
                    ),

                    html.H4("Choose datatype", style=header_button_text), 
                    dcc.RadioItems(options=['bulk'],
                                   value='bulk',
                                   id='params-dtype',
                                   style=center_component_style),

                    html.H4("Choose whether to run outlier analysis", style=header_button_text), 
                    dcc.RadioItems(options=['True', 'False'],
                                   id='params-outlier', 
                                   value='True',
                                   style=center_component_style),

                    html.H4("Choose z-score cutoff", style=header_button_text),
                    dcc.Input( 
                        id='params-z-score',
                        type='number',
                        value=2,
                        style=center_component_style
                    ),

                    html.H4("Choose mean threshold", style=header_button_text), 
                    dcc.Input( 
                        id='params-mean',
                        type='number', 
                        value=1,
                        style=center_component_style),

                    html.H4("Choose CV threshold", style=header_button_text),
                    dcc.Input(
                        id='params-cv',
                        type='number',
                        value=5,
                        style=center_component_style), 

                    html.H4("Choose NA Threshold", style=header_button_text), 
                    dcc.Input(
                        id='params-na',
                        type='number', 
                        value=0.4,
                        style=center_component_style), 

                    html.H4("Choose name of feature set", style=header_button_text), 
                    dcc.Input(
                        id='params-features', 
                        type='text',
                        value='PTID Time',
                        style=center_component_style), 

                    html.H4("Choose where to save outputs", style=header_button_text), 
                    dcc.Input(
                        id='params-output',
                        type='text',
                        value='/home/jupyter',
                        style=center_component_style),

                    html.Button("Run App", id='run_app_btn', style={'align-items' : 'center',
                                                               'width': '200px',
                                                               'display':'grid',
                                                               'justify-content' : 'center',
                                                               'margin' : '0 auto',
                                                               'background' : colors['aigreen']})
                ])
            ]),


            ## correlation 
            dcc.Tab(label='Correlation', value='correlation', style=tab_style, selected_style=tab_selected_style),


            ## expression levels 
            dcc.Tab(label='Expression Level', value='expression level', style=tab_style, selected_style=tab_selected_style),

            ## variance tab 
            dcc.Tab(label='Variance', value='variance', style=tab_style, selected_style=tab_selected_style),

            # intra-donor-variation 
            dcc.Tab(label='Intra-donor Variation', value='intra-donor-variation', style=tab_style, selected_style=tab_selected_style),

            # outliers tab 
            dcc.Tab(label='Outliers', value='outliers', style=tab_style, selected_style=tab_selected_style)
        ],
                 colors={
                        'border':colors['aiwhite'],
                        'primary' : colors['aiwhite'], 
                        'background' : colors['aiblue']}),

        # id that holds all the tab content 
        html.Div(id='tabs-content'),

        # storing data 
        dcc.Store(id='input-datamatrix'),
        dcc.Store(id='input-metadata'), 
        dcc.Store(id='input-var'),
        dcc.Store(id='input-outlier'),
        dcc.Store(id='expression-out'),
        dcc.Store(id='corr-matrix'),
        dcc.Store(id='variance-out'),
        dcc.Store(id='var-dl-counter'),
        dcc.Store(id='gene-df')
    ])


    del app.config._read_only["requests_pathname_prefix"]
    app.run_server(mode="inline", debug=True)
    return

if __name__ == "__main__": 
    run_app()








