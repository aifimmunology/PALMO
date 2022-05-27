# dash stuff
from dash import Dash, html, dcc
import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

# plotting
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

# data manipulation
import pandas as pd
import numpy as np
import rpy2.robjects as robjects

# custom prep module
import DashPalmPrep as dpp

# debugging
import pdb
import os
import time

# bootstrap - style sheet
external_stylesheets = [
    'https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css'
]
#app = JupyterDash(__name__, external_stylesheets=external_stylesheets,long_callback_manager=long_callback_manager) # jupyter-only
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.config['suppress_callback_exceptions'] = True

THISCOLORSCALE = px.colors.qualitative.Dark2
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
##############################
# --- correlation matrix --- #
##############################


@app.callback(Output('corr-matrix', 'data'),
              Output('corr-fig', 'figure'),
              Input('input-data', 'data'),
              Input('corr-dropdown', 'value'),
              prevent_initial_call=True)
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
              Input('input-outlier', 'data'),
              State('z-submit', 'n_clicks'),
              [State('sample-selector', 'value')],
              prevent_initial_call=True)
def run_and_cache_outlier_input(data_input, n_clicks, select_samples):
    if n_clicks is None:
        raise PreventUpdate
    else:
        outlier_res_py = pd.read_json(data_input, orient='split')
        if n_clicks == 0:
            return outlier_res_py.to_json(date_format='iso', orient='split')
        else:
            outlier_res_py = outlier_res_py.loc[
                outlier_res_py['Sample'].isin(select_samples), ]
            return outlier_res_py.to_json(date_format='iso', orient='split')


@app.callback(Output('number-features-bar', 'figure'),
              Input('outlier-input-state', 'data'),
              State('params-z-score', 'value'),
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
                    text='no_features',
                    color_discrete_sequence=THISCOLORSCALE)
    nofeat.update_traces(texttemplate='%{text:.2s}', textposition='outside')
    return nofeat


@app.callback(Output('p-scatter', 'figure'),
              Input('outlier-input-state', 'data'),
              Input('z-score-cutoff', 'value'), State('params-z-score',
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
                     },
                     color_discrete_sequence=THISCOLORSCALE)
    return fig


@app.callback(Output('p-bar', 'figure'), Input('outlier-input-state', 'data'),
              Input('z-score-cutoff', 'value'), State('params-z-score',
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
                 },
                 color_discrete_sequence=THISCOLORSCALE)
    return fig


@app.callback(Output('outlier-plot', 'figure'), Input('center-measure',
                                                      'value'),
              Input('outlier-input-state', 'data'),
              State('params-z-score', 'value'),
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


@app.callback(Output('sample-selector', 'options'),
              Input('input-outlier', 'data'))
def define_sample_options(outlier_data):
    outlier = pd.read_json(outlier_data, orient='split')
    sample_list = outlier['Sample'].unique().tolist()
    sample_list.sort()
    return sample_list


##################################
# --- Expression Level Plots --- #
##################################


@app.callback(Output('this-gene-scrna', 'options'),
              Input('this-gene-scrna', 'value'), State('var-decomp', 'data'))
def define_gene_options(this_gene, var_decomp):
    var_decomp = pd.read_json(var_decomp, orient='split')
    gene_list = var_decomp['Gene'].unique().tolist()
    return gene_list


@app.callback(Output('download-scrna-exp-df', 'data'),
              Input('scrna_exp_dl_btn', 'n_clicks'),
              State('scrna-subset-data', 'data'))
def download_scrna_exp(n_clicks, dff):
    '''
    '''
    if n_clicks is not None:
        dff = pd.read_json(dff, orient='split')
        # create json blob that will hold these values
        return dcc.send_data_frame(dff.to_csv, 'scrna_exp_data.csv')


@app.callback(Output('scrna_exp_gene_plot', 'figure'),
              Output('scrna-subset-data', 'data'),
              Input('this-gene-scrna', 'value'), State('input-data', 'data'),
              State('input-metadata', 'data'))
def plot_scrna_gene_exp(this_gene, input_data, input_metadata):
    data = pd.read_json(input_data, orient='split')
    metadata = pd.read_json(input_metadata, orient='split')
    data.index = data['Unnamed: 0'].values
    del data['Unnamed: 0']
    this_gene_vals = data.loc[data.index == this_gene, ].values
    metadata['exp'] = list(this_gene_vals[0])

    # plot
    gfig = px.strip(metadata,
                    x='PTID',
                    y='exp',
                    color='Time',
                    color_discrete_sequence=THISCOLORSCALE)
    boxf = px.box(metadata,
                  x='PTID',
                  y='exp',
                  color_discrete_sequence=THISCOLORSCALE)
    combined_plot = go.Figure(data=boxf.data + gfig.data)
    return combined_plot, metadata.to_json(orient='split')


@app.callback(
    Output('scrna_exp_gene_time_plot', 'figure'),
    Input('this-gene-scrna', 'value'),
    State('input-data', 'data'),
    State('input-metadata', 'data'),
)
def plot_scrna_gene_bytime(this_gene, input_data, input_metadata):
    data = pd.read_json(input_data, orient='split')
    metadata = pd.read_json(input_metadata, orient='split')
    data.index = data['Unnamed: 0'].values
    del data['Unnamed: 0']
    this_gene_vals = data.loc[data.index == this_gene, ].values
    metadata['exp'] = list(this_gene_vals[0])
    fig = px.strip(metadata,
                   x='PTID',
                   y='exp',
                   color='Time',
                   facet_col='Time',
                   color_discrete_sequence=THISCOLORSCALE)
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


@app.callback(Output('expression-out', 'data'),
              Input('this-gene', 'value'),
              Input('input-data', 'data'),
              Input('input-metadata', 'data'),
              prevent_initial_call=True)
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
    gfig = px.scatter(gene_df,
                      x='PTID',
                      y='exp',
                      symbol='Time',
                      color='Time',
                      color_discrete_sequence=THISCOLORSCALE)
    boxf = px.box(gene_df,
                  x='PTID',
                  y='exp',
                  color_discrete_sequence=THISCOLORSCALE)
    combined_plot = go.Figure(data=boxf.data + gfig.data)
    return combined_plot


@app.callback(Output('geneplot', 'figure'), Input('expression-out', 'data'))
def geneplot(df):
    '''
    '''
    ann = pd.read_json(df, orient='split')
    fig = px.scatter(ann,
                     x='PTID',
                     y='exp',
                     color='Time',
                     facet_col='Time',
                     color_discrete_sequence=THISCOLORSCALE)
    fig.for_each_annotation(
        lambda a: a.update(text=""))  # remove annotation from each subplot
    return fig


@app.callback(Output('this-gene', 'options'), Input('var-decomp', 'data'))
def define_gene_options(var_decomp):
    var_decomp = pd.read_json(var_decomp, orient='split')
    gene_list = var_decomp['Gene'].unique().tolist()
    return gene_list


##########################
# --- Variance Plots --- #
##########################


@app.callback(Output('download-scrna-variance', 'data'),
              Input('scrna_variance_dl_btn', 'n_clicks'),
              State('scrna_gene_subset_output', 'data'))
def download_scrna_variance(n_clicks, dff):
    '''
    '''
    if n_clicks is not None:
        dff = pd.read_json(dff, orient='split')
        # create json blob that will hold these values
        return dcc.send_data_frame(dff.to_csv, 'scrna_variance_data.csv')


@app.callback(Output('scrna_variance_decomp', 'figure'),
              Input('var-decomp', 'data'))
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
              Output('scrna_gene_subset_output', 'data'),
              Input('var-decomp', 'data'),
              Input('var-select-btn', 'n_clicks'),
              State('gene-dropdown-select', 'value'),
              prevent_initial_call=True)
def scrna_features_decomp_plot(var_input, clicks, var_chosen):
    if clicks is None:
        raise PreventUpdate
    else:
        var_decomp = pd.read_json(var_input, orient='split')
        var_decomp.index = var_decomp['Gene']
        var_decomp = pd.melt(
            var_decomp,
            id_vars='Gene',
            value_vars=['PTID', 'Time', 'Residual', 'celltype'])
        donors_sorted = var_decomp.loc[var_decomp['variable'].eq('donor'), ]
        if clicks > 0:
            var_decomp = var_decomp.loc[var_decomp['Gene'].isin(var_chosen), ]
        else:
            # choose some default genes - top 15
            top_vars = var_decomp.sort_values(
                by='value', ascending=False)[0:15]['Gene'].unique().tolist()
            var_decomp = var_decomp.loc[var_decomp['Gene'].isin(top_vars), ]
        donors_sorted.sort_values(by='value', ascending=False, inplace=True)
        var_fig = px.bar(var_decomp,
                         x='value',
                         y='Gene',
                         color='variable',
                         orientation='h',
                         category_orders={
                             'Gene': donors_sorted['Gene'].tolist(),
                             'variable':
                             ['PTID', 'Time', 'celltype', 'Residual']
                         },
                         labels={
                             'Gene': 'Features',
                             'value': '% variance explained'
                         },
                         color_discrete_sequence=THISCOLORSCALE)
        var_fig.update_traces().update_layout(title_x=0.5,
                                              legend_title_text='FeatureList')
        return var_fig, var_decomp.to_json(orient='split')


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


@app.callback(Output('variance-out', 'data'),
              Input('feature-button', 'n_clicks'),
              Input('var-decomp', 'data'),
              State('params-mean', 'value'),
              [State('gene-dropdown-select', 'value')],
              prevent_initial_call=True)
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
        },
        color_discrete_sequence=THISCOLORSCALE)
    var_fig.update_traces().update_layout(title_x=0.5,
                                          legend_title_text='FeatureList')
    return (var_fig)


@app.callback(Output('violin-box-var-plot', 'figure'),
              Input('var-decomp', 'data'), State('params-mean', 'value'))
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


##################################
# --- Intra-donor variations --- #
##################################


@app.callback(Output('download-scrna-cv-gene-df', 'data'),
              Input('scrna_cv_gene_btn', 'n_clicks'),
              State('scrna_stablevar_output', 'data'))
def download_scrna_cv_gene(n_clicks, dff):
    '''
    '''
    if n_clicks is not None:
        dff = pd.read_json(dff, orient='split')
        # create json blob that will hold these values
        return dcc.send_data_frame(dff.to_csv, 'scrna_cv_gene_data.csv')


@app.callback(Output('cvprof-hg-fig', 'figure'), Input('input-cvprof', 'data'))
def scrna_hg_cvprofile(input_data):
    cv_input = pd.read_json(input_data, orient='split')
    housekeeping_genes = ['GAPDH', 'ACTB']
    housekeep_input = cv_input.loc[cv_input['gene'].isin(housekeeping_genes), ]
    fig = px.scatter(housekeep_input,
                     x='mean',
                     y='CV',
                     facet_col='gene',
                     color_discrete_sequence=THISCOLORSCALE)
    return fig


@app.callback(Output('this-gene-intra', 'options'),
              Input('input-cvprof', 'data'))
def define_group_options(cvprof):
    cvprof = pd.read_json(cvprof, orient='split')
    gene_list = cvprof['group'].unique().tolist()
    return gene_list


@app.callback(Output('scrna_cv_profile', 'figure'),
              Input('input-cvprof', 'data'),
              Input('this-gene-intra', 'value'),
              State('params-cv', 'value'),
              prevent_initial_call=True)
def scrna_cv_profile_plot(input_cv, this_gene, cv_threshold):
    cv_input = pd.read_json(input_cv, orient='split')
    housekeeping_genes = ['GAPDH', 'ACTB']
    var_genes = cv_input.loc[cv_input['is_var_gene'].eq(1), ]
    stable_genes = cv_input.loc[cv_input['is_var_gene'].eq(0), ]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(mode='markers',
                   name='variable_genes',
                   x=var_genes.loc[var_genes['group'].eq(this_gene), 'mean'],
                   y=var_genes.loc[var_genes['group'].eq(this_gene), 'CV'],
                   marker={'color': 'black'}))

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

    fig = px.histogram(bar_plot_input,
                       x='PTID',
                       labels={
                           'PTID': 'CV'
                       },
                       color_discrete_sequence=THISCOLORSCALE).update_layout(
                           yaxis_title='density')
    return fig


@app.callback(Output('cvplot', 'figure'),
              Output('gene-df', 'data'),
              Input('gene-type', 'value'),
              Input('input-metadata', 'data'),
              Input('input-data', 'data'),
              State('params-mean', 'value'),
              State('params-cv', 'value'),
              prevent_initial_call=True)
def cvplot(gene_type, metadata_input, matrix_input, mean_threshold,
           cv_threshold):
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
        return dcc.send_data_frame(dff.to_csv, 'cv_gene_data.csv')


#################################
# --- Stable/Var gene plots --- #
#################################


@app.callback(Output('input_stablevar', 'data'),
              Output('stablevar-selector', 'options'),
              Output('stablevar-selector', 'value'),
              Input('scrna_stablevar_select', 'value'),
              Input('subset-stablevar', 'n_clicks'),
              State('input-cvprof', 'data'),
              [State('stablevar-selector', 'value')])
def subset_scrna_stablevar(var_or_stable, clicks, data, genes_chosen):
    input_df = pd.read_json(data, orient='split')
    ctx = dash.callback_context
    action = ctx.triggered[0]['prop_id'].split('.')[0]
    if var_or_stable == "Stable":
        stable_df = input_df.loc[input_df['is_var_gene'].eq(0), ]
        selectable_genes = stable_df['gene'].unique().tolist()
        selectable_genes.sort()
        if (clicks == 0) or (action == 'scrna_stablevar_select'):
            return stable_df.to_json(orient='split'), selectable_genes, []
        else:
            stable_df = stable_df.loc[stable_df['gene'].isin(genes_chosen), ]
            return stable_df.to_json(
                orient='split'), selectable_genes, genes_chosen
    elif var_or_stable == "Variable":
        var_df = input_df.loc[input_df['is_var_gene'].eq(1), ]
        selectable_genes = var_df['gene'].unique().tolist()
        selectable_genes.sort()
        if (clicks == 0) or (action == 'scrna_stablevar_select'):
            return var_df.to_json(orient='split'), selectable_genes, []
        else:
            var_df = var_df.loc[var_df['gene'].isin(genes_chosen), ]
            return var_df.to_json(
                orient='split'), selectable_genes, genes_chosen


@app.callback(Output('scrna_stable_plot', 'figure'),
              Output('scrna_stablevar_output', 'data'),
              Input('input_stablevar', 'data'),
              State('scrna_stablevar_select', 'value'))
def scrna_variable_gene_plot(dff, is_stable):

    def _calc_top25(gene_df):
        # remove duplicates
        top_genes = gene_df.drop_duplicates(subset='gene')

        # then sort by top 25 freqs
        top_genes.sort_values(by='freq_x', inplace=True, ascending=False)
        # return list of top genes
        return top_genes['gene'][0:25].unique().tolist()

    input_df = pd.read_json(dff, orient='split')
    if is_stable == 'Stable':
        var_genes = input_df.loc[input_df['is_var_gene'].eq(0), ]
        # calculate freq
        var_genes['freq'] = 1
        freqs = var_genes.groupby(
            ['gene'])['freq'].sum().reset_index().sort_values(by='freq',
                                                              ascending=False)
        var_genes = pd.merge(freqs, var_genes, on='gene')
        this_col_scale = [[0, 'rgb(210, 201, 242)'], [1, 'rgb(82, 68, 242)']]
        top_genes = _calc_top25(var_genes)
        var_genes = var_genes.loc[var_genes['gene'].isin(top_genes), ]
        var_genes.sort_values(by='donor', inplace=True, ascending=True)
        fig = go.Figure(data=go.Heatmap(z=var_genes['CV'],
                                        x=var_genes['donor'],
                                        y=var_genes['gene'],
                                        colorscale=this_col_scale))
        # fig.update_layout(zaxis_title={'title': 'CV'})
        return fig, var_genes.to_json(orient='split')
    elif is_stable == "Variable":
        var_genes = input_df.loc[input_df['is_var_gene'].eq(1), ]
        var_genes['freq'] = 1
        freqs = var_genes.groupby(
            ['gene'])['freq'].sum().reset_index().sort_values(by='freq',
                                                              ascending=False)
        var_genes = pd.merge(freqs, var_genes, on='gene')
        this_col_scale = [[0,
                           'rgb(210, 201, 242)'], [0.1, 'rgb(233, 237, 14)'],
                          [0.2, 'rgb(224, 36, 29)'], [1, 'rgb(82, 22, 20)']]
        top_genes = _calc_top25(var_genes)
        var_genes = var_genes.loc[var_genes['gene'].isin(top_genes), ]
        var_genes.sort_values(by='donor', inplace=True, ascending=True)
        fig = go.Figure(data=go.Heatmap(z=var_genes['CV'],
                                        x=var_genes['donor'],
                                        y=var_genes['gene'],
                                        colorscale=this_col_scale))
        #fig.update_layout(z='CV')
        return fig, var_genes.to_json(orient='split')


######################
# --- scRNA UMAP --- #
######################
@app.callback(Output('celltype_umap', 'figure'), Output('pc_fig', 'figure'),
              Input('umap_input', 'data'))
def plot_umap(dff):
    umap_input = pd.read_json(dff, orient='split')
    umap_input.sort_values(by='predicted.celltype.l2.y', inplace=True)
    # plot umap
    umap_celltype = px.scatter(umap_input,
                               x='UMAP_1',
                               y='UMAP_2',
                               color='predicted.celltype.l2.y',
                               height=600,
                               labels={'predicted.celltype.l2.y': 'celltype'})
    umap_celltype.update_layout(xaxis=dict(showgrid=False),
                                yaxis=dict(showgrid=False))
    umap_celltype.update_traces(marker=dict(size=1.5))

    # plot PC1 PC2
    pc_fig = px.scatter(umap_input,
                        x='PC_1',
                        y='PC_2',
                        color='predicted.celltype.l2.y',
                        height=600,
                        labels={'predicted.celltype.l2.y': 'celltype'})
    pc_fig.update_layout(xaxis=dict(showgrid=False),
                         yaxis=dict(showgrid=False))
    pc_fig.update_traces(marker=dict(size=1.5))
    return umap_celltype, pc_fig


##################################
# --- Prep and data processing --- #
##################################
@app.callback(Output('gene-dropdown-select', 'options'),
              Input('var-decomp', 'data'))
def define_gene_options(var_decomp):
    var_decomp = pd.read_json(var_decomp, orient='split')
    gene_list = var_decomp['Gene'].unique().tolist()
    return gene_list


@app.callback(Output('output-content', 'children'),
              Output('input-data', 'data'),
              Output('input-metadata', 'data'),
              Output('var-decomp', 'data'),
              Output('input-cvprof', 'data'),
              Output('input-outlier', 'data'),
              Output('umap_input', 'data'), [Input('run_app_btn', 'n_clicks')],
              State('params-dtype', 'value'),
              State('params-z-score', 'value'),
              State('params-mean', 'value'),
              State('params-cv', 'value'),
              State('params-na', 'value'),
              prevent_initial_call=True)
def render_tabs(click1, dtype, zscore_cutoff, mean_cutoff, cv_cutoff,
                na_cutoff):
    """ handles loading in the correct tabs """
    ctx = dash.callback_context
    action = ctx.triggered[0]['prop_id'].split('.')[0]
    if (action == 'run_app_btn') and (click1 is not None):
        (data, metadata, var_decomp, cv_res, outlier,
         umap) = run_prep(dtype, zscore_cutoff, mean_cutoff, cv_cutoff,
                          na_cutoff)
        # if user doesn't specify z-score: return empty data.frame

        # get gene options
        tabs = []
        for num in dtype_dict[dtype]:
            content = get_tab_content(dtype, num)
            if (zscore_cutoff is None) and (num == 'Outliers'):
                outlier = pd.DataFrame()
                continue
            tabs.append(
                dcc.Tab(label=num,
                        value=num,
                        children=[content],
                        className='default-tab-style',
                        selected_className='tab-selected-style'))
        # return to some default tab
        return (dcc.Tabs(id='tab',
                         value='tab1',
                         children=tabs,
                         colors={
                             'border': colors['aiwhite'],
                             'primary': colors['aiwhite'],
                             'background': colors['aiblue']
                         }), data.to_json(orient='split'),
                metadata.to_json(orient='split'),
                var_decomp.to_json(orient='split'),
                cv_res.to_json(orient='split'),
                outlier.to_json(orient='split'), umap.to_json(orient='split'))


def run_prep(datatype, zscore_cutoff, mean_cutoff, cv_cutoff, na_cutoff):
    meta_fpath = get_metapath()
    (data, metadata, var_decomp, cv_res, outlier,
     umap) = dpp.run(data_filepath=get_datapath(datatype),
                     metadata_filepath=meta_fpath,
                     datatype=determine_datatype(datatype),
                     z_cutoff=zscore_cutoff,
                     mean_threshold=mean_cutoff,
                     cv_threshold=cv_cutoff,
                     na_threshold=na_cutoff,
                     housekeeping_genes=['GAPDH', 'ACTB'],
                     feature_set=['PTID', 'Time'])
    return (data, metadata, var_decomp, cv_res, outlier, umap)


# structure to determine number/names of tabs
dtype_dict = {
    'Bulk Plasma': [
        'Correlation', 'Expression Level', 'Variance', 'Intra-donor Variation',
        'Outliers'
    ],
    'Single-Cell':
    ['Expression Level', 'Variance', 'Intra-donor Variation', 'UMAP']
}

# some sort of tab handler maps tab-dtype to content
bulk_tab_content = {
    'Correlation':
    html.Div([
        html.Br(),
        html.Button("Download CSV",
                    id='corr_btn',
                    className='button-default-style'),
        html.H4("Choose the correlation coefficient method",
                className='header-button-text'),
        dcc.Dropdown(id='corr-dropdown',
                     options=['pearson', 'spearman'],
                     value='spearman',
                     className='graph-component-align'),
        dcc.Download(id='download-corr-df'),
        dcc.Graph(id='corr-fig', className='matrix-graph-style'),
    ]),
    'Expression Level':
    html.Div([
        html.Br(),
        html.Button('Download CSV',
                    id='expression_dl_button',
                    className='button-default-style'),
        html.H4('Pick gene of interest', className='header-button-text'),
        dcc.Dropdown(value='ABL1', id='this-gene',
                     className='custom-dropdown'),
        dcc.Download(id='download-expression-df'),
        dcc.Graph(id='var-gene-plot', className='side-by-side-graph'),
        dcc.Graph(id='geneplot', className='side-by-side-graph')
    ]),
    'Variance':
    html.Div([
        html.Br(),
        html.Button("Download CSV",
                    id='var_btn',
                    className='button-default-style'),
        html.H4("Select feature(s) of interest",
                className='header-button-text'),
        dcc.Dropdown(id='gene-dropdown-select',
                     multi=True,
                     className='custom-dropdown'),
        html.Button(id='feature-button',
                    n_clicks=0,
                    children='Update feature plot',
                    className='button-default-style'),
        dcc.Download(id='download-var-df'),
        dcc.Graph(id='variance-contribution'),
        dcc.Graph(id='violin-box-var-plot'),
    ]),
    'Intra-donor Variation':
    html.Div([
        html.Br(),
        html.Button("Download CSV",
                    id='cv_gene_btn',
                    className='button-default-style'),
        html.H4("Choose to plot variable or stable genes",
                className='header-button-text'),
        dcc.RadioItems(['stable', 'variable'],
                       'stable',
                       id='gene-type',
                       className='graph-component-align'),
        dcc.Download(id='download-cv-gene-df'),
        dcc.Graph(id='cvplot', className='custom-cv-graph'),
        dcc.Graph(id='cv-density'),
    ]),
    'Outliers':
    html.Div([
        html.Br(),
        html.Button('Download CSV',
                    id='outlier_dl_button',
                    className='button-default-style'),
        html.H4('Select a sample of interest', className='header-button-text'),
        dcc.Dropdown(id='sample-selector',
                     multi=True,
                     placeholder='Select samples of interest',
                     className='graph-component-align'),

        # insert button
        html.Button('Update Plots',
                    id='z-submit',
                    className='button-default-style'),
        dcc.Download(id='download-outlier-df'),
        dcc.RadioItems(
            ['meanDev', 'z'],
            'z',
            id='center-measure',
        ),
        html.Div(children=[
            dcc.Graph(id='outlier-plot', className='side-by-side-graph'),
            dcc.Graph(id='number-features-bar',
                      className='side-by-side-graph'),
        ]),
        html.H5("Choose to plot where |Z| is above or below z-cutoff",
                className='header-button-text'),
        dcc.RadioItems(['below', 'above', 'all'],
                       'all',
                       id='z-score-cutoff',
                       className='graph-component-align'),
        html.Div(children=[
            dcc.Graph(id='p-scatter', className='side-by-side-graph'),
            dcc.Graph(id='p-bar', className='side-by-side-graph')
        ])
    ])
}

scrna_tab_content = {
    'Expression Level':
    html.Div([
        html.Button("Download CSV",
                    id='scrna_exp_dl_btn',
                    className='button-default-style'),
        dcc.Download(id='download-scrna-exp-df'),
        html.H4('Pick gene of interest', className='header-button-text'),
        dcc.Dropdown(value='ABL1',
                     id='this-gene-scrna',
                     className='custom-dropdown'),
        dcc.Graph(id='scrna_exp_gene_plot', className='side-by-side-graph'),
        dcc.Graph(id='scrna_exp_gene_time_plot',
                  className='side-by-side-graph')
    ]),
    'Variance':
    html.Div([
        html.Br(),
        html.Button("Download CSV",
                    id='scrna_variance_dl_btn',
                    className='button-default-style'),
        dcc.Download(id='download-scrna-variance'),
        html.H4("Select feature(s) of interest",
                className='header-button-text'),
        dcc.Dropdown(id='gene-dropdown-select',
                     multi=True,
                     className='custom-dropdown'),
        html.Button(id='var-select-btn',
                    children='Update feature plot',
                    className='button-default-style'),
        dcc.Graph(id='scrna_features_decomp'),
        dcc.Graph(id='scrna_variance_decomp')
    ]),
    'Intra-donor Variation':
    html.Div([
        html.Br(),
        html.Button("Download CSV",
                    id='scrna_cv_gene_btn',
                    className='button-default-style'),
        html.H4("Select between variale or stable",
                className='header-button-text'),
        dcc.RadioItems(['Stable', 'Variable'],
                       'Variable',
                       id='scrna_stablevar_select',
                       className='graph-component-align'),
        dcc.Dropdown(id='stablevar-selector',
                     multi=True,
                     placeholder='Select genes of interest',
                     className='custom-dropdown'),
        html.Button(id='subset-stablevar',
                    n_clicks=0,
                    children='Update Plot',
                    className='button-default-style'),
        dcc.Download(id='download-scrna-cv-gene-df'),
        dcc.Graph(id='scrna_stable_plot'),
        html.H4("Select celltype of interest", className='header-button-text'),
        dcc.Dropdown(id='this-gene-intra',
                     value='ASDC',
                     className='custom-dropdown'),
        dcc.Graph(id='scrna_cv_profile'),
        # dcc.Graph(id='cvprof-hg-fig'),
    ]),
    'UMAP':
    html.Div([
        dcc.Graph(id='pc_fig', className='side-by-side-graph'),
        dcc.Graph(id='celltype_umap', className='side-by-side-graph'),
    ])
}


def get_tab_content(datatype, this_tab):
    if datatype == 'Bulk Plasma':
        return bulk_tab_content[this_tab]
    elif datatype == "Single-Cell":
        return scrna_tab_content[this_tab]


def get_datapath(data_entry):
    '''
    '''
    if data_entry == 'Bulk Plasma':
        return ('{}/data/Olink_NPX_log2_Protein.Rda'.format(os.getcwd()))
    if data_entry == "Single-Cell":
        return ('{}/data/AIFI-scRNA-PBMC-FinalData.RDS'.format(os.getcwd()))


def get_metapath():
    return '{}/data/data_Metadata.Rda'.format(os.getcwd())


def determine_datatype(user_input):
    """ based on users' submission, return proper datatype"""
    if user_input == 'Bulk Plasma':
        return 'bulk'
    elif user_input == 'Single-Cell':
        return 'singlecell'


@app.callback(Output('params-content', 'children'),
              Input('params-dtype', 'value'))
def render_params(dtype):
    """ returns div object that holds submission for parameters """
    if dtype == 'Bulk Plasma':
        return html.Div(children=[
            dbc.Row([]),
            dbc.Row([
                dbc.Col([
                    html.H4("Choose CV threshold", className='params-header'),
                    html.P('(Required)', className='optional-text'),
                    dcc.Input(id='params-cv',
                              type='number',
                              value=5,
                              className='center-component-style'),
                ]),
                dbc.Col([
                    html.H4("Choose mean threshold",
                            className='params-header'),
                    html.P('(Required)', className='optional-text'),
                    dcc.Input(id='params-mean',
                              type='number',
                              value=1,
                              className='center-component-style'),
                ]),
                dbc.Col([
                    html.H4("Choose z-score cutoff",
                            className='params-header'),
                    html.P('(Optional)', className='optional-text'),
                    dcc.Input(id='params-z-score',
                              type='number',
                              className='center-component-style'),
                ]),
                dbc.Col([
                    html.H4("Choose NA Threshold", className='params-header'),
                    html.P('(Optional)', className='optional-text'),
                    dcc.Input(id='params-na',
                              type='number',
                              className='center-component-style')
                ])
            ],
                    style={"background-color": colors['aiblue']})
        ])
    elif dtype == "Single-Cell":
        return html.Div(children=[
            dbc.Row([
                dbc.Col([
                    dcc.Input(id='params-z-score',
                              type='number',
                              className='center-component-style',
                              style={'display': 'none'})
                ]),
                dbc.Col([
                    dcc.Input(id='params-na',
                              type='number',
                              className='center-component-style',
                              style={'display': 'none'})
                ])
            ]),
            dbc.Row([
                dbc.Col([
                    html.H4("Choose CV threshold", className='params-header'),
                    html.P('(Required)', className='optional-text'),
                    dcc.Input(id='params-cv',
                              type='number',
                              value=10,
                              className='center-component-style')
                ]),
                dbc.Col([
                    html.H4("Choose mean threshold",
                            className='params-header'),
                    html.P('(Required)', className='optional-text'),
                    dcc.Input(id='params-mean',
                              type='number',
                              value=0.1,
                              className='center-component-style'),
                ])
            ]),
        ],
                        style={"background-color": colors['aiblue']})


submit_params_page = html.Div(children=[
    html.Div(children=[
        html.H2('Platform for Analyzing Longitudinal Multi-omics data',
                className='custom-h2-header'),
        html.H4("Choose datatype", className='params-header'),
        dcc.Dropdown(id='params-dtype',
                     options=list(dtype_dict.keys()),
                     className='custom-dropdown'),
        html.Div(id='params-content'),
        html.Button(
            "Run App", id='run_app_btn', className='button-default-style'),
    ],
             style={'background-color': colors['aiblue']}),
    dbc.Spinner(id='loading-content',
                children=[html.Div(id='output-content')],
                delay_hide=10,
                fullscreen=True)
])

data_store = html.Div(children=[
    dcc.Store('input-data'),
    dcc.Store('scrna_stablevar_output'),
    dcc.Store('scrna_gene_subset_output'),
    dcc.Store('input-metadata'),
    dcc.Store('var-decomp'),
    dcc.Store('scrna-subset-data'),
    dcc.Store('umap_input'),
    dcc.Store('input-cvprof'),
    dcc.Store('input-outlier'),
    dcc.Store('corr-matrix'),
    dcc.Store('download-var-df'),
    dcc.Store('variance-out'),
    dcc.Store('expression-out'),
    dcc.Store('outlier-input-state'),
    dcc.Store('gene-df'),
    dcc.Store('input_stablevar')
])
app.layout = html.Div([submit_params_page, data_store])

if __name__ == '__main__':
    app.run_server(debug=True,
                   host=os.getenv('IP', '0.0.0.0'),
                   port=int(os.getenv('PORT', 4444)))