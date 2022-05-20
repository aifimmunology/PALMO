''' run_dash.py

'''
#from jupyter_dash.comms import _send_jupyter_config_comm_request # jupyter-only
#_send_jupyter_config_comm_request() # jupyter-only

#from jupyter_dash import JupyterDash # jupyter-only
#JupyterDash.infer_jupyter_proxy_config() # jupyter-only

# load libraries
import os
import dash
from dash import html, dcc
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from dash.dependencies import Input, Output, State
import pandas as pd
import diskcache
from dash.long_callback import DiskcacheLongCallbackManager
import dash_bootstrap_components as dbc

# additional libraries that wasn't handled properly in prep code
from plotly.subplots import make_subplots
import rpy2.robjects as robjects
import numpy as np

# import prep module
import DashPalm_prep as dpp

# long callback manager set up
cache = diskcache.Cache("./cache")
long_callback_manager = DiskcacheLongCallbackManager(cache)

# bootstrap - style sheet
external_stylesheets = [
    'https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css'
]
#app = JupyterDash(__name__, external_stylesheets=external_stylesheets,long_callback_manager=long_callback_manager) # jupyter-only
app = dash.Dash(__name__,
                external_stylesheets=external_stylesheets,
                long_callback_manager=long_callback_manager)

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

##########################
# --- Variance Plots --- #
##########################


def prep_var_contribute_df(tbl: pd.DataFrame, mean_threshold):
    ''' Takes table of residuals and does some formatting and subsetting
    '''
    featureList = ['PTID', 'Time'] + ['Residual']
    res_tbl = tbl.loc[tbl['max'] > mean_threshold, featureList]
    # res_tbl.sort_values(by='PTID', ascending=False, inplace=True)

    # reshape
    res_tbl['genes'] = res_tbl.index
    var_df = pd.melt(res_tbl, id_vars='genes', value_vars=featureList)
    return var_df


@app.callback(Output('scrna_variance_decomp', 'figure'),
              Input('input-var', 'data'))
def scrna_variance_decomp_plot(var_input):
    var_decomp = pd.read_json(var_input, orient='split')
    var_decomp.index = var_decomp['Gene']
    var_decomp = pd.melt(var_decomp,
                         id_vars='Gene',
                         value_vars=['PTID', 'Time', 'Residual', 'celltype'])
    feature_list = var_decomp['variable'].unique().tolist()
    fig = go.Figure()
    for f in feature_list:
        fig.add_trace(
            go.Violin(x=var_decomp.loc[var_decomp['variable'].eq(f),
                                       'variable'],
                      y=var_decomp.loc[var_decomp['variable'].eq(f), 'value'],
                      name=f,
                      box_visible=True,
                      opacity=0.6,
                      spanmode='hard',
                      showlegend=True))
    fig.update_traces(width=0.5)
    fig.update_layout(legend_title_text='featureList')
    return fig


@app.callback(Output('scrna_features_decomp', 'figure'),
              Input('input-var', 'data'))
def scrna_features_decomp_plot(var_input):
    var_decomp = pd.read_json(var_input, orient='split')
    var_decomp.index = var_decomp['Gene']
    var_decomp = pd.melt(var_decomp,
                         id_vars='Gene',
                         value_vars=['PTID', 'Time', 'Residual', 'celltype'])
    donors_sorted = var_decomp.loc[var_decomp['variable'].eq('donor'), ]
    var_decomp = var_decomp.loc[
        var_decomp['Gene'].isin(['XIST', 'TSIX', 'JUN']), ]
    donors_sorted.sort_values(by='value', ascending=False, inplace=True)
    var_fig = px.bar(var_decomp,
                     x='value',
                     y='Gene',
                     color='variable',
                     orientation='h',
                     category_orders={
                         'Gene': donors_sorted['Gene'].tolist(),
                         'variable': ['PTID', 'Time', 'celltype', 'Residual']
                     },
                     labels={
                         'Gene': 'Features',
                         'value': '% variance explained'
                     })
    var_fig.update_traces().update_layout(title_x=0.5,
                                          legend_title_text='FeatureList')
    return var_fig


@app.callback(Output('variance-out', 'data'),
              Input('feature-button', 'n_clicks'), Input('input-var', 'data'),
              Input('params-mean', 'value'),
              [State('dropdown-select', 'value')])
def subset_var_df(n, input_var, mean_threshold, var_chosen):
    '''
    '''
    lmem_py = pd.read_json(input_var, orient='split')
    var_df = prep_var_contribute_df(lmem_py, mean_threshold)
    if (n == 0):
        ptid_entries = var_df.loc[var_df['variable'].eq('PTID'), ].copy(
            deep=True)
        top_vars = ptid_entries.sort_values(
            by='value', ascending=False)[0:15]['genes'].unique().tolist()
        dff = var_df.loc[var_df['genes'].isin(top_vars), ]
    elif n > 0:
        dff = var_df[var_df['genes'].isin(var_chosen)]
    return dff.to_json(orient='split')


@app.callback(Output('variance-contribution', 'figure'),
              Input('variance-out', 'data'),
              suppress_callback_exceptions=True,
              prevent_initial_callback=False)
def var_contribute_plot(dff):
    ''' separate function that renames and subsets data 
    '''
    dff = pd.read_json(dff, orient='split')
    # sort by donor
    donors_sorted = dff.loc[dff['variable'].eq('PTID'), ]
    donors_sorted.sort_values(by='value', ascending=False, inplace=True)
    var_fig = px.bar(
        dff,
        x='value',
        y='genes',
        color='variable',
        orientation='h',
        category_orders={'genes': donors_sorted['genes'].tolist()},
        labels={
            'genes': 'Features',
            'value': '% variance explained'
        })
    var_fig.update_traces().update_layout(title_x=0.5,
                                          legend_title_text='FeatureList')
    return (var_fig)


@app.callback(Output('violin-box-var-plot', 'figure'),
              Input('input-var', 'data'), Input('params-mean', 'value'))
def violin_box_plot(input_var, mean_threshold):
    ''' creates violin with box plot overlayed
    TODO: plot sigFeature datapoints
    '''
    lmem_py = pd.read_json(input_var, orient='split')
    featureList = ['PTID', 'Time'] + ['Residual']
    res = lmem_py.loc[lmem_py['max'] > mean_threshold, featureList]
    res['genes'] = res.index
    df_py = pd.melt(res,
                    id_vars='genes',
                    value_vars=['PTID', 'Time', 'Residual'])
    df_py['value'] = df_py['value'] * 100
    df_py['feature'] = df_py[['genes', 'variable']].agg('_'.join, axis=1)

    # list features
    feature_list = df_py['variable'].unique().tolist()
    fig = go.Figure()
    for f in feature_list:
        fig.add_trace(
            go.Violin(x=df_py.loc[df_py['variable'].eq(f), 'variable'],
                      y=df_py.loc[df_py['variable'].eq(f), 'value'],
                      name=f,
                      box_visible=True,
                      opacity=0.6,
                      spanmode='hard',
                      showlegend=True))
    fig.update_traces(width=0.5)
    fig.update_layout(legend_title_text='featureList')

    return fig


@app.callback(Output('download-var-df', 'data'), Input('var_btn', 'n_clicks'),
              State('variance-out', 'data'))
def download_var_test(n_clicks, dff):
    '''
    '''
    if n_clicks is not None:
        dff = pd.read_json(dff, orient='split')
        # create json blob that will hold these values
        return dcc.send_data_frame(dff.to_csv, 'variance_viz_data.csv')


#########################
# --- Outlier Plots --- #
#########################


@app.callback(Output('download-outlier-df', 'data'),
              Input('outlier_dl_button', 'n_clicks'),
              State('outlier-input-state', 'data'))
def download_outlier_data(n_clicks, dff):
    '''
    '''
    if n_clicks is not None:
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
        rate = 2 * (1 - rpnorm(z_cutoff)[0])

    # create a data.frame with freq #'s
    # first, create a table and convert to a data.frame
    tab = data.groupby([groupby, 'Time']).size()
    df = tab.unstack().reset_index()
    df = pd.wide_to_long(df, stubnames="W", i='PTID', j='Time').reset_index()

    # now do some cleaning; and fill NAs
    df['tmp'] = 'W'
    df['Time'] = df['Time'].astype(str)
    df['Time'] = df[['tmp', 'Time']].agg(''.join, axis=1)
    df.fillna(value=0, inplace=True)

    # rename and drop columns
    df.rename(columns={groupby: 'Sample', 'W': 'Freq'}, inplace=True)
    df.drop(columns='tmp', inplace=True)
    df['pvals'] = p_value_for_event(events=df['Freq'].tolist(),
                                    tries=nGenes,
                                    rate=rate)
    df['signP'] = -np.log10(df['pvals'])

    # create an id column
    df['id'] = df[['Sample', 'Time']].agg(''.join, axis=1)
    df.sort_values(by=['Time', 'Sample'], inplace=True)

    return df


@app.callback(Output('outlier-input-state', 'data'),
              Input('z-submit', 'n_clicks'), Input('input-outlier', 'data'),
              [State('sample-selector', 'value')])
def run_and_cache_outlier_input(n_clicks, data_input, select_samples):
    outlier_res_py = pd.read_json(data_input, orient='split')
    if n_clicks == 0:
        return outlier_res_py.to_json(date_format='iso', orient='split')
    else:
        outlier_res_py = outlier_res_py.loc[
            outlier_res_py['Sample'].isin(select_samples), ]
        return outlier_res_py.to_json(date_format='iso', orient='split')


@app.callback(Output('number-features-bar', 'figure'),
              Input('outlier-input-state', 'data'),
              Input('params-z-score', 'value'),
              [State('sample-selector', 'value')])
def plot_no_features(df, z_cutoff, select_samples):
    '''
    '''
    # create column for number of features
    nofeatures_df = pd.read_json(df, orient='split')
    nofeatures_df['no_features'] = 1
    nofeatures_df.loc[nofeatures_df['z'] > 0, 'zgroup'] = '> Z'
    nofeatures_df.loc[nofeatures_df['z'] < 0, 'zgroup'] = '< -Z'

    nofeatures_df = nofeatures_df.groupby(
        ['Sample', 'zgroup'])['no_features'].sum().reset_index()

    if select_samples is not None:
        nofeatures_df = nofeatures_df.loc[
            nofeatures_df['Sample'].isin(select_samples), ]
    nofeat = px.bar(nofeatures_df,
                    x='Sample',
                    y='no_features',
                    color='zgroup',
                    width=500,
                    barmode='group',
                    text='no_features')
    nofeat.update_traces(texttemplate='%{text:.2s}', textposition='outside')
    return nofeat


@app.callback(Output('p-scatter', 'figure'),
              Input('outlier-input-state', 'data'),
              Input('z-score-cutoff', 'value'), Input('params-z-score',
                                                      'value'),
              [State('sample-selector', 'value')])
def outlier_p_scatter_plot(data_input, z_subset, z_score_cutoff,
                           select_samples):
    '''
    '''
    d_input = pd.read_json(data_input, orient='split')
    input_df = prep_outlierP_data(d_input,
                                  z_cutoff=z_score_cutoff,
                                  z_score_subset=z_subset,
                                  nGenes=1042,
                                  groupby="PTID")
    if select_samples is not None:
        input_df = input_df.loc[input_df['id'].isin(select_samples), ]
    fig = px.scatter(input_df,
                     x='Freq',
                     y='signP',
                     hover_name='id',
                     labels={
                         'Freq': '# Features',
                         'signP': '-log10(p-value)'
                     })
    return fig


@app.callback(Output('p-bar', 'figure'), Input('outlier-input-state', 'data'),
              Input('z-score-cutoff', 'value'), Input('params-z-score',
                                                      'value'),
              [State('sample-selector', 'value')])
def outlier_p_bar_plot(data_input, z_subset, z_score_cutoff, select_samples):
    '''
    '''
    d_input = pd.read_json(data_input, orient='split')
    input_df = prep_outlierP_data(d_input,
                                  z_cutoff=z_score_cutoff,
                                  z_score_subset=z_subset,
                                  nGenes=1042,
                                  groupby="PTID")
    input_df['sample_time'] = input_df[['Sample', 'Time']].agg(''.join, axis=1)
    input_df.sort_values(by=['sample_time'], inplace=True)

    if select_samples is not None:
        input_df = input_df.loc[input_df['sample_time'].isin(select_samples), ]
    fig = px.bar(input_df,
                 x='sample_time',
                 y='signP',
                 labels={
                     'sample_time': 'Sample',
                     'signP': '-log10(pvalue)'
                 })
    return fig


@app.callback(Output('outlier-plot', 'figure'), Input('center-measure',
                                                      'value'),
              Input('outlier-input-state', 'data'),
              Input('params-z-score', 'value'),
              [State('sample-selector', 'value')])
def outlier_detect_plot(measure_col, data_input, z_cut_off, these_samples):
    ''' 
    TODO: order the samples correctly 
    '''
    # sort values
    d_input = pd.read_json(data_input, orient='split')
    d_input.sort_values(by='Sample', inplace=True)
    if these_samples is not None:
        d_input = d_input.loc[d_input['Sample'].isin(these_samples), ]

    out_fig = go.Figure()
    d_input.loc[d_input['z'] > 0, 'zgroup'] = '> Z'
    d_input.loc[d_input['z'] < 0, 'zgroup'] = '< -Z'
    z_group_list = list(d_input['zgroup'].unique())
    z_group_list.sort()
    for this_group in z_group_list:
        violin_input = d_input.loc[
            d_input['zgroup'].eq(this_group), ].sort_values(by='Sample')
        out_fig.add_trace(
            go.Violin(
                x=violin_input['Sample'],
                y=violin_input[measure_col],
                legendgroup=this_group,
                name=this_group,
                box_visible=True,
                spanmode='hard',
                opacity=0.6,
            ))

    out_fig.update_traces(width=0.5)
    sample_order = d_input['Sample'].unique().tolist()
    sample_order.sort()
    out_fig.update_xaxes(categoryorder='array', categoryarray=sample_order)
    return out_fig


##################################
# --- Expression Level Plots --- #
##################################


@app.callback(Output('scrna_exp_gene_plot', 'figure'),
              Input('input-datamatrix', 'data'), Input('input-metadata',
                                                       'data'))
def plot_scrna_gene_exp(input_data, input_metadata):
    data = pd.read_json(input_data, orient='split')
    metadata = pd.read_json(input_metadata, orient='split')
    data.index = data['Unnamed: 0'].values
    del data['Unnamed: 0']
    this_gene = 'LILRA4'
    this_gene_vals = data.loc[data.index == this_gene, ].values
    metadata['exp'] = list(this_gene_vals[0])

    # plot
    gfig = px.strip(metadata, x='PTID', y='exp', color='Time')
    boxf = px.box(metadata, x='PTID', y='exp')
    combined_plot = go.Figure(data=boxf.data + gfig.data)
    return combined_plot


@app.callback(Output('scrna_exp_gene_time_plot', 'figure'),
              Input('input-datamatrix', 'data'), Input('input-metadata',
                                                       'data'))
def plot_scrna_gene_bytime(input_data, input_metadata):
    data = pd.read_json(input_data, orient='split')
    metadata = pd.read_json(input_metadata, orient='split')
    data.index = data['Unnamed: 0'].values
    del data['Unnamed: 0']
    this_gene = 'LILRA4'
    this_gene_vals = data.loc[data.index == this_gene, ].values
    metadata['exp'] = list(this_gene_vals[0])
    fig = px.strip(metadata, x='PTID', y='exp', color='Time', facet_col='Time')
    fig.for_each_annotation(
        lambda a: a.update(text=""))  # remove annotation from each subplot
    return fig


@app.callback(Output('download-expression-df', 'data'),
              Input('expression_dl_button', 'n_clicks'),
              State('expression-out', 'data'))
def download_outlier_data(n_clicks, dff):
    '''
    '''
    if n_clicks is not None:
        dff = pd.read_json(dff, orient='split')
        return dcc.send_data_frame(dff.to_csv, 'expression_lvl_viz_data.csv')


@app.callback(Output('expression-out', 'data'), Input('this-gene', 'value'),
              Input('input-datamatrix', 'data'), Input('input-metadata',
                                                       'data'))
def store_gene_data(this_gene, input_data, input_metadata):
    ''' Caches and stores data for callbacks 
    '''
    datamatrix = pd.read_json(input_data, orient='split')
    metadata = pd.read_json(input_metadata, orient='split')
    this_gene_vals = datamatrix[datamatrix.index == this_gene].values
    ann = metadata.copy(deep=True)
    ann['exp'] = list(this_gene_vals[0])
    return ann.to_json(orient='split')


@app.callback(Output('var-gene-plot', 'figure'), Input('expression-out',
                                                       'data'))
def var_gene_plot(df):
    '''
    '''
    gene_df = pd.read_json(df, orient='split')
    gfig = px.scatter(gene_df, x='PTID', y='exp', symbol='Time', color='Time')
    boxf = px.box(gene_df, x='PTID', y='exp')
    combined_plot = go.Figure(data=boxf.data + gfig.data)
    return combined_plot


@app.callback(Output('geneplot', 'figure'), Input('expression-out', 'data'))
def geneplot(df):
    '''
    '''
    ann = pd.read_json(df, orient='split')
    fig = px.scatter(ann, x='PTID', y='exp', color='Time', facet_col='Time')
    fig.for_each_annotation(
        lambda a: a.update(text=""))  # remove annotation from each subplot
    return fig


##############################
# --- correlation matrix --- #
##############################


@app.callback(Output('corr-matrix', 'data'), Output('corr-fig', 'figure'),
              Input('input-datamatrix', 'data'), Input('corr-dropdown',
                                                       'value'))
def correlation_matrix_plot(data_input, correlation_method='spearman'):
    datamatrix = pd.read_json(data_input, orient='split')
    cols = datamatrix.columns.tolist()
    cor_res = datamatrix.corr(method=correlation_method)
    data_input = cor_res[cols].to_numpy().transpose()
    fig = go.Figure(
        go.Heatmap(z=data_input,
                   colorscale=[[0, 'white'], [0.5, 'rgb(250, 138, 130)'],
                               [1.0, 'rgb(255, 17, 0)']],
                   x=cols,
                   y=cols,
                   dx=1,
                   dy=1))
    return (pd.DataFrame(data_input).to_json(orient='split'), fig)


@app.callback(Output('download-corr-df', 'data'),
              Input('corr_btn', 'n_clicks'), State('corr-matrix', 'data'))
def download_correlation(n_clicks, dff):
    '''
    '''
    if n_clicks is not None:
        dff = pd.read_json(dff, orient='split')
        # create json blob that will hold these values
        return dcc.send_data_frame(dff.to_csv,
                                   'correlation_matrix_viz_data.csv')


##################################
# --- Intra-donor variations --- #
##################################
"""
@app.callback(Output('scrna_cv_profile', 'figure'),
              Input('input-cvprof', 'data'), Input('params-cv', 'value'))
def scrna_cv_profile_plot(input_cv, cv_threshold):
    print('trying to create this scrna plot man')
    cv_input = pd.read_json(input_cv, orient='split')
    # housekeeping_genes = ['GAPDH', 'ACTB']
    this_gene = 'ASDC'
    var_genes = cv_input.loc[cv_input['is_var_gene'].eq(1), ]
    stable_genes = cv_input.loc[cv_input['is_var_gene'].eq(0), ]
    fig = px.scatter(var_genes.loc[var_genes['group'].eq(this_gene), ],
                     x='mean',
                     y='CV',
                     color_discrete_sequence=['black'])
    fig.add_trace(
        go.Scatter(mode='markers',
                   name='non_variable_genes',
                   x=stable_genes.loc[stable_genes['group'].eq(this_gene),
                                      'mean'],
                   y=stable_genes.loc[stable_genes['group'].eq(this_gene),
                                      'CV']))
    # add horizontal line
    fig.add_hline(y=cv_threshold)

    fig.add_trace(
       go.Scatter(mode='markers',
                  name='housekeeping_genes',
                   x=cv_input.loc[cv_input['group'].isin(housekeeping_genes),
                                  'mean'],
                   y=cv_input.loc[cv_input['group'].isin(housekeeping_genes),
                                  'CV'],
                   marker={'color': 'blue'}))

    return fig
"""


# NOTE: hard-coded vars here
def prep_cvcalc_bulk(meanThreshold, cvThreshold, datamatrix, metadata):
    ''' Intro-donor variations over time 
    
    TODO: label top 10 genes, and most stable ones 
    '''
    mat = datamatrix
    ann = metadata

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


@app.callback(Output('cv-density', 'figure'), Input('gene-df', 'data'),
              Input('input-metadata', 'data'))
def cv_density_plot(gene_df, metadata_input):
    '''
    '''
    metadata = pd.read_json(metadata_input, orient='split')
    dff = pd.read_json(gene_df, orient='split')
    uniSample = metadata['PTID'].unique().tolist()
    bar_plot_input = pd.wide_to_long(dff, stubnames='PTID', i='gene',
                                     j='var').reset_index()

    fig = px.histogram(bar_plot_input, x='PTID', labels={
        'PTID': 'CV'
    }).update_layout(yaxis_title='density')
    return fig


@app.callback(Output('cvplot', 'figure'), Output('gene-df', 'data'),
              Input('gene-type', 'value'), Input('input-metadata', 'data'),
              Input('input-datamatrix', 'data'), Input('params-mean', 'value'),
              Input('params-cv', 'value'))
def cvplot(gene_type, metadata_input, matrix_input, mean_threshold,
           cv_threshold):
    ''' here we can either pick the genes of interest for heatmap
        or we could explore by making the mean & cv threshold interactive. (think a slider) 
    '''
    datamatrix = pd.read_json(matrix_input, orient='split')
    metadata = pd.read_json(metadata_input, orient='split')
    mat_input = prep_cvcalc_bulk(mean_threshold, cv_threshold, datamatrix,
                                 metadata)  # for now, just use the default

    if (gene_type == 'variable'):
        mat = mat_input.sort_values(by='Median', ascending=False)
        mat = mat.loc[np.abs(mat['Median']) > cv_threshold, ]
        col_scale = [[0, 'rgb(210, 201, 242)'], [0.1, 'rgb(233, 237, 14)'],
                     [0.2, 'rgb(224, 36, 29)'], [1, 'rgb(82, 22, 20)']]
    elif (gene_type == 'stable'):
        mat = mat_input.sort_values(by='Median', ascending=True)
        mat = mat.loc[np.abs(mat['Median']) <= cv_threshold, ]
        col_scale = [[0, 'rgb(210, 201, 242)'], [0.3, 'rgb(100, 79, 179)'],
                     [0.7, 'rgb(82, 68, 242)'], [0.9, 'rgb(233, 237, 14)'],
                     [1, 'rgb(82, 22, 20)']]

    # subset to first 50 entries then create heat map
    mat = mat[0:50]
    heat_fig = go.Heatmap(
        z=mat[['PTID1', 'PTID2', 'PTID3', 'PTID4', 'PTID5', 'PTID6']].values,
        x=['PTID1', 'PTID2', 'PTID3', 'PTID4', 'PTID5', 'PTID6'],
        y=mat['gene'].unique().tolist(),
        colorscale=col_scale)
    d1 = mat['PTID1'].values
    d2 = mat['PTID2'].values
    d3 = mat['PTID3'].values
    d4 = mat['PTID4'].values
    d5 = mat['PTID5'].values
    d6 = mat['PTID6'].values

    b1 = go.Box(y=d1)
    b2 = go.Box(y=d2)
    b3 = go.Box(y=d3)
    b4 = go.Box(y=d4)
    b5 = go.Box(y=d5)
    b6 = go.Box(y=d6)

    master_fig = make_subplots(rows=2,
                               cols=1,
                               row_heights=[0.2, 0.8],
                               vertical_spacing=0,
                               specs=[[{
                                   'type': 'box'
                               }], [{
                                   'type': 'xy'
                               }]])
    master_fig.add_traces([b1, b2, b3, b4, b5, b6], rows=1,
                          cols=1).update_layout(showlegend=False)
    master_fig.update_xaxes(showticklabels=False, row=1)
    master_fig.add_traces([heat_fig], rows=2, cols=1)
    return master_fig, mat.to_json(orient='split')


@app.callback(Output('download-cv-gene-df', 'data'),
              Input('cv_gene_btn', 'n_clicks'), State('gene-df', 'data'))
def download_cv_gene(n_clicks, dff):
    '''
    '''
    if n_clicks is not None:
        dff = pd.read_json(dff, orient='split')
        # create json blob that will hold these values
        return dcc.send_data_frame(dff.to_csv,
                                   'correlation_matrix_viz_data.csv')


############################
# --- Run Prep & Store --- #
############################
def get_datapath(data_entry):
    '''
    '''
    print('getting filepath...')
    if data_entry == 'Bulk Plasma':
        return (
            '/Users/james.harvey/workplace/repos/PALM/dash_app/data/Olink_NPX_log2_Protein.Rda'
        )
    if data_entry == "scRNA":
        return '/Users/james.harvey/workplace/repos/PALM/dash_app/data/Olink_NPX_log2_Protein.Rda'


def get_metapath():
    return '/Users/james.harvey/workplace/repos/PALM/dash_app/data/data_Metadata.Rda'


"""
@app.callback(
    Output('input-datamatrix', 'data'),
    Output('input-metadata', 'data'),
    Output('input-var', 'data'),
    Output('input-outlier', 'data'),
    Output('params-mean', 'value'),
    Output('params-cv', 'value'),
    Output('params-z-score', 'value'),
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
)
def generate_bulk_results(metadata_fpath, data_fpath, datatype,
                 run_outlier, z_score_cutoff, mean_cutoff, cv_cutoff,
                 na_cutoff, feature_list, output_dir):
    (datamatrix, metadata, lmem_py, cv_res,
         outlier_res_py) = dpp.run(data_filepath=get_datapath(data_fpath),
                                   metadata_filepath=metadata_fpath,
                                   datatype=datatype,
                                   do_outlier=run_outlier,
                                   z_cutoff=z_score_cutoff,
                                   mean_threshold=mean_cutoff,
                                   cv_threshold=cv_cutoff,
                                   na_threshold=na_cutoff,
                                   housekeeping_genes=['GAPDH', 'ACTB'],
                                   output_dir=output_dir,
                                   feature_set=['PTID', 'Time'])
        # convert objects to json and return
        print('returning json blobs...')
        return (datamatrix.to_json(orient='split'),
                metadata.to_json(orient='split'),
                lmem_py.to_json(orient='split'),
                outlier_res_py.to_json(orient='split'), mean_cutoff, cv_cutoff,
                z_score_cutoff)
"""


@app.long_callback(
    output=[
        Output('input-datamatrix', 'data'),
        Output('input-metadata', 'data'),
        Output('input-var', 'data'),
        Output('input-outlier', 'data'),
        # Output('input-cvprof', 'data'),
        Output('params-mean', 'value'),
        Output('params-cv', 'value'),
        Output('params-z-score', 'value')
    ],
    inputs=[
        Input('run_app_btn', 'n_clicks'),
        # State('metadata-entry', 'value'),
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
    running=[(Output("run_app_btn", "disabled"), True, False)],
    prevent_initial_call=True)
def parse_params(run_n_clicks, data_fpath, datatype, run_outlier,
                 z_score_cutoff, mean_cutoff, cv_cutoff, na_cutoff,
                 feature_list, output_dir):
    print('hm...{}'.format(run_n_clicks))
    if (run_n_clicks is not None):
        (datamatrix, metadata, lmem_py, cv_res,
         outlier_res_py) = dpp.run(data_filepath=get_datapath(data_fpath),
                                   metadata_filepath=get_metapath(),
                                   datatype=datatype,
                                   do_outlier=run_outlier,
                                   z_cutoff=z_score_cutoff,
                                   mean_threshold=mean_cutoff,
                                   cv_threshold=cv_cutoff,
                                   na_threshold=na_cutoff,
                                   housekeeping_genes=['GAPDH', 'ACTB'],
                                   output_dir=output_dir,
                                   feature_set=['PTID', 'Time'])
        # convert objects to json and return
        print(run_n_clicks)
        print('returning json blobs...')
        return (datamatrix.to_json(orient='split'),
                metadata.to_json(orient='split'),
                lmem_py.to_json(orient='split'),
                outlier_res_py.to_json(orient='split'), mean_cutoff, cv_cutoff,
                z_score_cutoff)

    else:
        return pd.DataFrame().to_json(), pd.DataFrame().to_json(
        ), pd.DataFrame().to_json(), pd.DataFrame().to_json(), pd.DataFrame(
        ).to_json(), 0, 0, 0


######################
# --- Tab handler ---#
######################


def need_to_run_layout():
    return html.Div([
        html.H4("""Please go to the submit parameters tab first and run the app
                to generate visualizations""",
                className='need-to-run-message')
    ])


def need_to_run_layout2():
    return html.Div([
        html.H4(
            """PPPPPlease go to the submit parameters tab first and run the app
                to generate visualizations""",
            className='need-to-run-message')
    ])


@app.callback(Output('tabs-content', 'children'), Input('tabs-graph', 'value'),
              Input('run_app_btn', 'n_clicks'), Input('input-outlier', 'data'),
              Input('input-datamatrix', 'data'), Input('input-var', 'data'),
              Input('params-dtype', 'value'))
def render_content(tab, run_n_clicks, outlier_input, datamatrix_input,
                   input_var, dtype):
    if run_n_clicks is None:
        if tab in [
                'correlation', 'expression_level', 'variance',
                'intra-donor-variation', 'outliers'
        ]:
            return need_to_run_layout()
    else:
        if tab == 'correlation':
            if dtype == 'bulk':
                return html.Div([
                    html.H4("Choose the correlation coefficient method",
                            className='header-button-text'),
                    dcc.Dropdown(id='corr-dropdown',
                                 options=['pearson', 'spearman'],
                                 value='spearman',
                                 className='graph-component-align'),
                    html.H4("Download dataset",
                            className='header-button-text'),
                    html.Button("Download CSV",
                                id='corr_btn',
                                className='button-default-style'),
                    dcc.Download(id='download-corr-df'),
                    dcc.Graph(id='corr-fig', className='matrix-graph-style'),
                ])
            elif dtype == 'singlecell':
                return need_to_run_layout2()
        if tab == 'expression_level':
            datamatrix = pd.read_json(datamatrix_input, orient='split')
            if dtype == 'bulk':
                return html.Div([
                    html.H4('Pick gene of interest',
                            className='header-button-text'),
                    dcc.Dropdown(datamatrix.index.unique(),
                                 value='ABL1',
                                 id='this-gene',
                                 className='custom-dropdown'),
                    html.H4("Download dataset",
                            className='header-button-text'),
                    html.Button('Download CSV',
                                id='expression_dl_button',
                                className='button-default-style'),
                    dcc.Download(id='download-expression-df'),
                    dcc.Graph(id='var-gene-plot',
                              className='side-by-side-graph'),
                    dcc.Graph(id='geneplot', className='side-by-side-graph')
                ])
            elif dtype == 'singlecell':
                return html.Div([
                    dcc.Graph(id='scrna_exp_gene_plot'),
                    dcc.Graph(id='scrna_exp_gene_time_plot')
                ])
        if tab == 'variance':
            if dtype == 'bulk':
                lmem_py = pd.read_json(input_var, orient='split')
                return html.Div([
                    html.H4("Select feature(s) of interest",
                            className='header-button-text'),
                    dcc.Dropdown(id='dropdown-select',
                                 options=prep_var_contribute_df(
                                     lmem_py, 0)['genes'].unique().tolist(),
                                 multi=True),
                    html.Button(id='feature-button',
                                n_clicks=0,
                                children='Update feature plot',
                                className='button-default-style'),
                    html.H4("Download dataset",
                            className='header-button-text'),
                    html.Button("Download CSV",
                                id='var_btn',
                                className='button-default-style'),
                    dcc.Download(id='download-var-df'),
                    dcc.Graph(id='variance-contribution'),
                    dcc.Graph(id='violin-box-var-plot'),
                ])
            elif dtype == 'singlecell':
                return html.Div([
                    dcc.Graph(id='scrna_features_decomp'),
                    dcc.Graph(id='scrna_variance_decomp')
                ])
        if tab == 'intra-donor-variation':
            if dtype == 'bulk':
                return html.Div([
                    html.H4("Choose to plot variable or stable genes",
                            className='header-button-text'),
                    dcc.RadioItems(['stable', 'variable'],
                                   'stable',
                                   id='gene-type',
                                   className='graph-component-align'),
                    html.H4("Download dataset",
                            className='header-button-text'),
                    html.Button("Download CSV",
                                id='cv_gene_btn',
                                className='button-default-style'),
                    dcc.Download(id='download-cv-gene-df'),
                    dcc.Graph(id='cvplot', className='custom-cv-graph'),
                    dcc.Graph(id='cv-density'),
                ])
            elif dtype == 'singlecell':
                return need_to_run_layout()
        if tab == 'outliers':
            if dtype == 'bulk':
                outlier_res_py = pd.read_json(outlier_input, orient='split')
                sample_list = outlier_res_py['Sample'].unique().tolist()
                sample_list.sort()
                return html.Div([
                    html.H4('Select a sample of interest',
                            className='header-button-text'),
                    dcc.Dropdown(id='sample-selector',
                                 options=sample_list,
                                 multi=True,
                                 placeholder='Select samples of interest',
                                 className='graph-component-align'),

                    # insert button
                    html.Button('Update Plots',
                                id='z-submit',
                                n_clicks=0,
                                className='button-default-style'),
                    html.H4("Download dataset",
                            className='header-button-text'),
                    html.Button('Download CSV',
                                id='outlier_dl_button',
                                className='button-default-style'),
                    dcc.Download(id='download-outlier-df'),
                    dcc.RadioItems(
                        ['meanDev', 'z'],
                        'z',
                        id='center-measure',
                    ),
                    html.Div(children=[
                        dcc.Graph(id='outlier-plot',
                                  className='side-by-side-graph'),
                        dcc.Graph(id='number-features-bar',
                                  className='side-by-side-graph'),
                    ]),
                    html.H5(
                        "Choose to plot where |Z| is above or below z-cutoff",
                        className='header-button-text'),
                    dcc.RadioItems(['below', 'above', 'all'],
                                   'all',
                                   id='z-score-cutoff',
                                   className='graph-component-align'),
                    html.Div(children=[
                        dcc.Graph(id='p-scatter',
                                  className='side-by-side-graph'),
                        dcc.Graph(id='p-bar', className='side-by-side-graph')
                    ])
                ])
            elif dtype == 'singlecell':
                return need_to_run_layout()


######################
# --- App Layout --- #
######################


def run_app():
    app.layout = html.Div(children=[
        html.H2('Platform for Analyzing Longitudinal Multi-omics data',
                className='custom-h2-header'),
        dcc.Tabs(
            id='tabs-graph',
            value='parameters',
            children=[
                # we need to render this page at initial loading to gain access to the run_var_btn id
                dcc.Tab(
                    label='Submit Parameters',
                    value='parameters',
                    className='default-tab-style',
                    selected_className='tab-selected-style',
                    children=[
                        html.Div([

                            #dbc.Col([
                            #    html.H4("Metadata filepath entry",
                            #            className='header-button-text'),
                            #    dcc.Input(value='{}/data/data_Metadata.Rda'.
                            #              format(os.getcwd()),
                            #              id='metadata-entry',
                            #              type='text',
                            #              className='center-component-style')
                            #]),
                            dbc.Col([
                                html.H4("Choose datatype",
                                        className='header-button-text'),
                                dcc.RadioItems(
                                    options=['bulk', 'singlecell'],
                                    value='singlecell',
                                    id='params-dtype',
                                    className='center-component-style')
                            ]),
                            dbc.Row([
                                dbc.Col([
                                    html.H4("Choose dataset entry",
                                            className='header-button-text'),
                                    dcc.Dropdown(
                                        options=['Bulk Plasma', 'scRNA'],
                                        id='datamatrix-entry')
                                ]),
                                dbc.Col([
                                    html.
                                    H4("Choose whether to run outlier analysis",
                                       className='header-button-text'),
                                    dcc.RadioItems(
                                        options=['True', 'False'],
                                        id='params-outlier',
                                        value='True',
                                        className='center-component-style')
                                ]),
                                dbc.Col([
                                    html.H4("Choose where to save outputs",
                                            className='header-button-text'),
                                    dcc.Input(
                                        id='params-output',
                                        type='text',
                                        value=os.getcwd(),
                                        className='center-component-style')
                                ])
                            ]),
                            html.H4("Choose z-score cutoff",
                                    className='header-button-text'),
                            dcc.Input(id='params-z-score',
                                      type='number',
                                      value=2,
                                      className='center-component-style'),
                            html.H4("Choose mean threshold",
                                    className='header-button-text'),
                            dcc.Input(id='params-mean',
                                      type='number',
                                      value=1,
                                      className='center-component-style'),
                            html.H4("Choose CV threshold",
                                    className='header-button-text'),
                            dcc.Input(id='params-cv',
                                      type='number',
                                      value=5,
                                      className='center-component-style'),
                            html.H4("Choose NA Threshold",
                                    className='header-button-text'),
                            dcc.Input(id='params-na',
                                      type='number',
                                      value=0.4,
                                      className='center-component-style'),
                            html.H4("Choose name of feature set",
                                    className='header-button-text'),
                            dcc.Input(id='params-features',
                                      type='text',
                                      value='PTID Time',
                                      className='center-component-style'),
                            html.Button("Run App",
                                        id='run_app_btn',
                                        className='button-default-style'),
                        ])
                    ]),

                ## correlation
                dcc.Tab(label='Correlation',
                        value='correlation',
                        className='default-tab-style',
                        selected_className='tab-selected-style'),

                ## expression levels
                dcc.Tab(label='Expression Level',
                        value='expression_level',
                        className='default-tab-style',
                        selected_className='tab-selected-style'),

                ## variance tab
                dcc.Tab(label='Variance',
                        value='variance',
                        className='default-tab-style',
                        selected_className='tab-selected-style'),

                # intra-donor-variation
                dcc.Tab(label='Intra-donor Variation',
                        value='intra-donor-variation',
                        className='default-tab-style',
                        selected_className='tab-selected-style'),

                # outliers tab
                dcc.Tab(label='Outliers',
                        value='outliers',
                        className='default-tab-style',
                        selected_className='tab-selected-style')
            ],
            colors={
                'border': colors['aiwhite'],
                'primary': colors['aiwhite'],
                'background': colors['aiblue']
            }),

        # id that holds all the tab content
        html.Div(id='tabs-content'),

        # storing data
        dcc.Store(id='input-datamatrix'),
        dcc.Store(id='input-metadata'),
        dcc.Store(id='outlier-input-state'),
        dcc.Store(id='input-var'),
        # dcc.Store(id='input-cvprof'),
        dcc.Store(id='input-outlier'),
        dcc.Store(id='expression-out'),
        dcc.Store(id='corr-matrix'),
        dcc.Store(id='variance-out'),
        dcc.Store(id='var-dl-counter'),
        dcc.Store(id='gene-df')
    ])

    # del app.config._read_only["requests_pathname_prefix"]
    app.run_server(debug=True,
                   host=os.getenv('IP', '0.0.0.0'),
                   port=int(os.getenv('PORT', 4442)))
    return


if __name__ == "__main__":
    run_app()
