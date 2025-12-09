#for importing, formatting and data manipulation
import pandas as pd
import numpy as np
import glob
import datetime
#from time import time
#from datetime import datetime
#from datetime import timedelta
import tempfile
from qiime2 import Artifact
import zipfile
import yaml
import re
import os


#for plotting
import matplotlib.pyplot as plt
import seaborn as sns
#sns.set(style="whitegrid")
import plotly.express as px
from IPython.display import display
from upsetplot import plot
#import pyupset as pyu
pd.set_option('display.max_rows', 15)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from pandas.plotting import register_matplotlib_converters
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
from statsmodels.tsa.stattools import acf, pacf
from statsmodels.tsa.statespace.sarimax import SARIMAX
register_matplotlib_converters()

#for statistical analyses
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest
from skbio.diversity import alpha_diversity
from skbio.stats.distance import permanova
from skbio import DistanceMatrix
from scipy.spatial.distance import cdist
from skbio.stats.composition import clr
from skbio.stats.composition import alr
from skbio.stats.composition import ilr
from skbio.diversity.alpha import chao1
from scipy.stats import zscore
from scipy.stats import entropy

from skbio.stats.composition import multiplicative_replacement
from skbio.stats.distance import DistanceMatrix
from scipy.spatial.distance import pdist, squareform
from skbio.diversity import beta_diversity
from collections import OrderedDict
from matplotlib.patches import Patch
from matplotlib import gridspec

from skbio.diversity.alpha import shannon
import networkx as nx
from matplotlib.patches import Rectangle

# Special thanks to Alex Manuele https://github.com/alexmanuele
def consolidate_tables(year=all):
    if year == all:
        table_list = glob.glob('{0}/*/02-PROKs/DADA2/table.qza'.format('/Users/Diana/Documents/escuela/phd/ch2/splitted_data'))
        print("Found all yearly tables")

    else:
        table_list = glob.glob('{0}/02-PROKs/DADA2/table.qza'.format('/Users/Diana/Documents/escuela/phd/ch2/splitted_data/'+year))
        print("Found year "+year+" tables")

    dataframes = []
    for table_path in table_list:
        with tempfile.TemporaryDirectory() as tempdir:
            #load table, dump contents to tempdir
            table = Artifact.load(table_path)
            #Make sure the tables are all FeatureFrequency type
            assert str(table.type) == 'FeatureTable[Frequency]', "{0}: Expected FeatureTable[Frequency], got {1}".format(table_path, table.type)
            Artifact.extract(table_path, tempdir)
            #get the provenance form the tempdir and format it for DF
            prov = '{0}/{1}/provenance/'.format(tempdir, table.uuid)
            action = yaml.load(open("{0}action/action.yaml".format(prov), 'r'), Loader=yaml.BaseLoader)
            paramlist = action['action']['parameters']
            paramlist.append({'table_uuid': "{}".format(table.uuid)})
            paramdict = {}
            for record in paramlist:
                paramdict.update(record)

            # Get the data into a dataframe
              #Biom data
            df = table.view(pd.DataFrame).unstack().reset_index()
            df.columns = ['feature_id', 'sample_name', 'feature_frequency']
            df['table_uuid'] = ["{}".format(table.uuid)] * df.shape[0]
              #param data
            pdf = pd.DataFrame.from_records([paramdict])
              #merge params into main df
            df = df.merge(pdf, on='table_uuid')


            #I like having these columns as the last three. Makes it more readable
            cols = df.columns.tolist()
            reorder = ['sample_name', 'feature_id', 'feature_frequency']
            for val in reorder:
                cols.append(cols.pop(cols.index(val)))
            df = df[cols]
            df['table_path'] = [table_path] * df.shape[0]
            df['sample_name'] = df['sample_name'].str.replace('-', '.')
            df['sample_name'] = df['sample_name'].str.replace(r'^V4V5\.', '', regex=True)

            sample_names = df['sample_name'].astype(str)
            parts = sample_names.str.extract(r'^([^\.]+)\.([^\.]+)')
            df['sampleid'] = parts[0] + '.' + parts[1]
            dataframes.append(df)

    #Stick all the dataframes together
    #outputfile="merged_all_tables.tsv"
    df = pd.concat(dataframes)

    #we have some problematic sample names
    df.sampleid = df.sampleid.replace({'BB19.35CSb': 'BB19.35CS', 'BB19.35BSb': 'BB19.35BS'})
    df['sampleid'] = df['sampleid'].str.replace(r'\.a(\d{1,2}[A-Za-z]{1,2})', r'.\1', regex=True)
    #df.to_csv(comm+'/merged_all_tables.tsv', sep='\t', index=False)
    print("Successfully saved all tables.")
    return df

def decorticate_sampleid (df, samplecol='sample_name'):
    #define dict for depth_code
    depth_num = {
        "A": 1,
        "B": 5,
        "C": 10,
        "D": 60,
        "E": 30
    }


    df = df[['nouveau', 'feature_id', 'feature_frequency']].copy()
    df.rename(columns={'nouveau':'sample_name'}, inplace=True)
    df['sample_name'] = df['sample_name'].str.replace(r'\.0(\d)', r'.\1', regex=True) #remove the leading zero

    df = df[df.feature_frequency != 0] #remove null rows

    # This regex first looks for 'SL', then 'S', then 'L'. If none is found, we fill with 'W'.
    df["size_code"] = df[samplecol].str.extract(r'(SL|S|L)')
    df["size_code"] = df["size_code"].fillna('W')

    df["depth_code"] = df[samplecol].str.extract(r'[1-9][0-9]?([A-E])')
    df['depth']= df['depth_code'].map(depth_num)
    df["weekn"] = df[samplecol].str.extract(r'\.([1-9][0-9]?)[A-E]')
    df['weekn'] = pd.to_numeric(df['weekn'])
    df['depth'] = pd.to_numeric(df['depth'])

    df['Total'] = df['feature_frequency'].groupby(df[samplecol]).transform('sum')
    df['ratio'] = df['feature_frequency']/df['Total']
    df['nASVs'] = df['feature_id'].groupby(df[samplecol]).transform('count')
    df['weekdepth'] = df["weekn"].astype(str) + df["depth"].astype(str)
    df['avg'] = df['nASVs'].groupby(df['weekdepth']).transform('mean')
    df['diff'] = df['nASVs'] - df['avg']

    df["year"] = '20' + df[samplecol].str.extract(r'BB(\d{2})(?=\.)')
    df['datetime'] = pd.to_datetime(
    df['year'].astype(str) + '-W' + df['weekn'].astype(int).astype(str) + '-1',
    format='%Y-W%W-%w')

    df=df.rename(columns={samplecol: "sampleid"})

    return df


def merge_metadata(df, all_md):
    #df = pd.read_csv('02-PROKs/'+'/merged_all_tables.tsv', sep='\t')

    depth_num = {
        "A": 1,
        "B": 5,
        "C": 10,
        "D": 60,
        "E": 30
    }
    tables = df[['sample_name', 'feature_id', 'feature_frequency']].copy()
    tables.rename(columns={'sample_name':'sampleid'}, inplace=True)

    all_md['sampleid'] = all_md['sampleid'].str.replace('_', '.')
    merged = pd.merge(tables,all_md, on='sampleid', how='outer') #all_md is the metadata file
    merged = merged[merged.feature_frequency != 0]

    merged["size_code"] = merged["sampleid"].str.extract(r'[1-9][0-9]?[A-E]([L-S])')
    merged["size_code"] = merged["size_code"].fillna('W')
    merged["depth_code"] = merged["sampleid"].str.extract(r'[1-9][0-9]?([A-E])')
    merged['depth']= merged['depth_code'].map(depth_num)
    merged["weekn"] = merged["sampleid"].str.extract(r'\.([1-9][0-9]?)[A-E]')
    merged['weekn'] = pd.to_numeric(merged['weekn'])
    merged['depth'] = pd.to_numeric(merged['depth'])
    group_dates = merged.groupby('weekn', as_index=False)['date'].first()
    merged = merged.drop(columns='date').merge(group_dates, on='weekn', how='left')

    merged['Total'] = merged['feature_frequency'].groupby(merged['sampleid']).transform('sum')
    merged['ratio'] = merged['feature_frequency']/merged['Total']
    merged['nASVs'] = merged['feature_id'].groupby(merged['sampleid']).transform('count')
    merged['weekdepth'] = merged["weekn"].astype(str) + merged["depth"].astype(str)
    merged['avg'] = merged['nASVs'].groupby(merged['weekdepth']).transform('mean')
    merged['diff'] = merged['nASVs'] - merged['avg']

    print('Set up metadata ...')

    #merged.to_csv(comm+'/merged_asvs_metadata.tsv', sep = '\t')
    print('Saved merged_asvs_metadata.tsv')

    return merged

#use this function if keeping chloro and bacteria together
def pick_metadata_taxo(comm, merged, depth='all', size_fraction='both', year='all', R='all', F='all', txsubset = 'all'):
#make df of features/composition+run+comm

    depth = depth
    year = year
    size_fraction = size_fraction
    txsubset = txsubset

    files = glob.glob('{0}/20**/taxa*/classification/*/data/taxonomy.tsv'.format('/Users/Diana/Documents/escuela/phd/ch2/bb_data'))
    taxos = []
#    if not os.path.exists(path+composition):
#        os.mkdir(path+composition)
    for filename in files:
        tax = pd.read_csv(filename, sep='\t')
        taxos.append(tax)

    print('Appended all taxonomies to taxos')
    taxos = pd.concat(taxos)
    taxos = taxos.rename(columns={"Feature ID": "feature_id"}, errors="raise")
    taxos = taxos.drop_duplicates()

    separated = merged.merge(taxos, how='left', on='feature_id') #merged excludes features of frequency = 0
    separated = separated.drop_duplicates()

    if depth != 'all':
        separated = separated[separated["depth"] == depth]
    if size_fraction != 'both':
        separated = separated[separated["size_fraction"] == size_fraction]

    separated[['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = separated['Taxon'].str.split('; ', expand=True)
    cols = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    for col in cols:
        separated[col] = separated[col].fillna('Unassigned')

    print('Saved separated by metadata dataframe.')
    print('Community is '+comm)

    contaminants = None
    if comm in ['02-PROKs', '02-EUKS']:
        searchfor = ["Cyanobacteria", "Chloroplast"]
        contaminants = separated[separated.Taxon.str.contains('|'.join(searchfor))]
        separated = separated[~separated.Taxon.str.contains('|'.join(searchfor))]
        separated = separated.reset_index(drop=True)
        print('Removed cyanobacteria and chloroplast from '+comm)

        #re-calculate ratios when removing chloroplast for 16S or 18S
        separated['Total'] = separated['feature_frequency'].groupby(separated['sampleid']).transform('sum')
        separated['ratio'] = separated['feature_frequency']/separated['Total']
        separated['nASVs'] = separated['feature_id'].groupby(separated['sampleid']).transform('count')
        separated['weekdepth'] = separated["weekn"].astype(str) + separated["depth"].astype(str)
        separated['avg'] = separated['nASVs'].groupby(separated['weekdepth']).transform('mean')
        separated['diff'] = separated['nASVs'] - separated['avg']

    elif comm == 'chloroplast':
        #run these lines to switch to chloroplast comm
        searchfor = ["Cyanobacteria", "Chloroplast"]
        contaminants = separated[separated.Taxon.str.contains('|'.join(searchfor))]
        separated = contaminants.copy()
        separated = separated.reset_index(drop=True)
        comm = 'chloroplast'
        print('Switched to cyanobacteria and chloroplast')

        #add phytorep taxonomy
        cp_tax = pd.read_csv('/Users/Diana/Documents/escuela/phd/ch2/taxonomy.tsv', sep='\t')
        cp_tax = cp_tax.rename(columns={"Feature ID": "feature_id", "Taxon": "PRTaxon", "Confidence":"PRConfidence"})
        separated = pd.merge(separated, cp_tax, on="feature_id", how="left")
        separated['PRSpecies'] = separated['PRTaxon'].str.split('|').str[-1]

        #re-calculate ratios when removing chloroplast for 16S
        separated['Total'] = separated['feature_frequency'].groupby(separated['sampleid']).transform('sum')
        separated['ratio'] = separated['feature_frequency']/separated['Total']
        separated['nASVs'] = separated['feature_id'].groupby(separated['sampleid']).transform('count')
        separated['weekdepth'] = separated["weekn"].astype(str) + separated["depth"].astype(str)
        separated['avg'] = separated['nASVs'].groupby(separated['weekdepth']).transform('mean')
        separated['diff'] = separated['nASVs'] - separated['avg']

    return separated, contaminants


#use this one if picking chloroplast, or removing chloroplast from 16S
def pick_metadata(
    comm,
    merged,
    depth='all',
    size_fraction='both',
    year='all',                 # reserved (not used below)
    txsubset='all',             # reserved (not used below)
    base_dir='/Users/Diana/Documents/escuela/phd/ch2/bb_data',
    chloroplast_pattern=r'Chloroplast|Cyanobacteria'  # edit if needed
):

    # --- collect all taxonomy.tsv files ---
    files = glob.glob(f'{base_dir}/20**/taxa*/classification/*/data/taxonomy.tsv')
    if not files:
        raise FileNotFoundError(f'No taxonomy.tsv found under: {base_dir}')

    taxos = []
    for fn in files:
        # Minimal, robust read; coerce Feature ID to str
        t = pd.read_csv(fn, sep='\t', dtype={'Feature ID': 'string'}, engine='python')
        taxos.append(t)

    taxos = pd.concat(taxos, ignore_index=True)
    if 'Feature ID' not in taxos.columns:
        raise KeyError("Expected 'Feature ID' column in taxonomy.tsv files.")
    if 'Taxon' not in taxos.columns:
        raise KeyError("Expected 'Taxon' column in taxonomy.tsv files.")

    taxos = taxos.rename(columns={"Feature ID": "feature_id"})
    taxos['feature_id'] = taxos['feature_id'].astype('string')
    taxos = taxos.drop_duplicates()

    # --- merge with composition table ---
    if 'feature_id' not in merged.columns:
        raise KeyError("merged DataFrame must contain 'feature_id' column.")
    separated = merged.merge(taxos, how='left', on='feature_id').drop_duplicates()

    # --- basic filters ---
    if depth != 'all':
        separated = separated.loc[separated['depth'] == depth].copy()
    if size_fraction != 'both':
        separated = separated.loc[separated['size_fraction'] == size_fraction].copy()
    # (You can add `year` or `txsubset` logic here if desired.)

    # --- split taxonomy into ranks (handles missing gracefully) ---
    if 'Taxon' in separated.columns:
        parts = separated['Taxon'].fillna('Unassigned').str.split(r';\s*', expand=True)
        # Ensure we have 7 columns
        while parts.shape[1] < 7:
            parts[parts.shape[1]] = 'Unassigned'
        parts.columns = ['Domain','Phylum','Class','Order','Family','Genus','Species']
        for col in parts.columns:
            parts[col] = parts[col].fillna('Unassigned')
        separated = pd.concat([separated, parts], axis=1)

    print('Saved separated-by-metadata dataframe.')
    print('Community is', comm)

    contaminants = None
    # Case A: remove cyano/chloro from 16S/18S
    if comm in ['02-PROKs', '02-EUKS']:
        mask_contam = separated['Taxon'].fillna('').str.contains(chloroplast_pattern, case=False, na=False, regex=True)
        contaminants = separated.loc[mask_contam].copy()
        separated = separated.loc[~mask_contam].copy()
        print(f'Removed Cyanobacteria/Chloroplast from {comm}')

        if {'feature_frequency','sampleid','feature_id','weekn','depth'}.issubset(separated.columns):
            separated['Total'] = separated.groupby('sampleid')['feature_frequency'].transform('sum')
            separated['ratio'] = separated['feature_frequency'] / separated['Total'].where(separated['Total'] != 0, pd.NA)
            separated['nASVs'] = separated.groupby('sampleid')['feature_id'].transform('nunique')
            separated['weekdepth'] = separated['weekn'].astype(str) + separated['depth'].astype(str)
            separated['avg'] = separated.groupby('weekdepth')['nASVs'].transform('mean')
            separated['diff'] = separated['nASVs'] - separated['avg']

    elif comm == 'chloroplast':
        mask_cp = separated['Taxon'].fillna('').str.contains(chloroplast_pattern, case=False, na=False, regex=True)
        contaminants = separated.loc[mask_cp].copy()   # original matching rows
        separated = contaminants.copy()
        print('Switched to Cyanobacteria/Chloroplast community')

        pr_path = '/Users/Diana/Documents/escuela/phd/ch2/taxonomy.tsv'
        if os.path.exists(pr_path):
            cp_tax = pd.read_csv(pr_path, sep='\t', dtype={'Feature ID': 'string'}, engine='python')
            cp_tax = cp_tax.rename(columns={"Feature ID": "feature_id", "Taxon": "PRTaxon", "Confidence": "PRConfidence"})
            separated = pd.merge(separated, cp_tax, on="feature_id", how="left")
            separated['PRSpecies'] = separated['PRTaxon'].fillna('').str.split('|').str[-1]

        if {'feature_frequency','sampleid','feature_id','weekn','depth'}.issubset(separated.columns):
            separated['Total'] = separated.groupby('sampleid')['feature_frequency'].transform('sum')
            separated['ratio'] = separated['feature_frequency'] / separated['Total'].where(separated['Total'] != 0, pd.NA)
            separated['nASVs'] = separated.groupby('sampleid')['feature_id'].transform('nunique')
            separated['weekdepth'] = separated['weekn'].astype(str) + separated['depth'].astype(str)
            separated['avg'] = separated.groupby('weekdepth')['nASVs'].transform('mean')
            separated['diff'] = separated['nASVs'] - separated['avg']

    return separated, contaminants




def load_df():

    filenames = glob.glob('/Users/Diana/Documents/escuela/phd/ch2/bb_data/20**/METADATA.txt')
    #load all metadata and concatenate them into one dataframe
    md = []
    for filename in filenames:
        df = pd.read_csv(filename, sep='\t')
        md.append(df)
        print (filename)

    md = pd.concat(md)

    #drop empty columns and rows
    md.dropna(how='all', axis=1, inplace=True) #empty cols
    md.dropna(how='all', inplace=True) #empty rows

    return md

def make_defract(all_md, separated):
    #make sure all size codes are indicated
    all_md["size_code"] = all_md["sampleid"].str.extract(r'[1-9][0-9]?[A-E]([L-S])')
    all_md["size_code"] = all_md["size_code"].fillna('W')

    #only keep values from weeks 1 to 16
    sep_SL = all_md[
        (all_md.size_code != "W") &
        (all_md.size_code != "P")
    ]
    sep_SL = sep_SL.drop(sep_SL[sep_SL.weekn > 16].index)

    #sum [DNA] of small and large size fractions
    sep_SL['[DNAt]'] = sep_SL.groupby(['weekn', 'depth'])['[DNA]ng/uL'].transform('sum')

    #separate small and large size fraction
    sep_S = sep_SL[sep_SL.size_code == 'S']
    sep_L = sep_SL[sep_SL.size_code == 'L']

    #calculate DNA proportion per size fraction
    sep_SL['DNApr'] = sep_SL['[DNA]ng/uL']/sep_SL['[DNAt]']

    #merge with separated on common columns to get corresponding rel. abundances
    sep_SL = sep_SL[['sampleid', 'DNApr', '[DNAt]']].copy()
    sepSLRA = pd.merge(separated, sep_SL, on=['sampleid'], how='left') #all_md is the metadata file

    #exclude ASVs from the whole water
    sep_SLRA = sepSLRA[
        (sepSLRA.size_code != 'W') &
        (sepSLRA.size_code != 'P')
    ]

    #calculate corrected per sample ratio, and corrected feature frequency of de-fractionated samples
    sep_SLRA['Newfeature_frequency'] = sep_SLRA['feature_frequency'] * sep_SLRA['DNApr']
    sep_SLRA['Newff'] = sep_SLRA.groupby(['feature_id', 'weekn', 'depth'])['Newfeature_frequency'].transform('sum')


    #sep_SLRA = sep_SLRA.drop(['sampleid', 'size_code'], axis=1)
    sep_SLRA['sampleid'] = "BB22." + sep_SLRA['weekn'].astype(str) + sep_SLRA['depth_code'] + "SL"

    #uncomment the line below if keeping small and large original sample
    #sep_SLRA['size_code'] = sep_SLRA['size_code'] + '-DFr'

    #uncomment the line above if merging smallandlarge
    sep_SLRA['size_code'] = 'SL'

    #drop unecessary columns which might rise merging conflicts
    sep_SLRA = sep_SLRA.drop(['feature_frequency', 'Total', 'ratio', 'nASVs', 'weekdepth', 'avg',
                              'diff', '[DNA]ng/uL',
                              'Newfeature_frequency'], axis=1)
    sep_SLRA.rename(columns={'Newff':'feature_frequency'}, inplace=True)
    sep_SLRA = sep_SLRA.drop_duplicates()

    #recalculate ratios
    sep_SLRA['Total'] = sep_SLRA['feature_frequency'].groupby(sep_SLRA['sampleid']).transform('sum')
    sep_SLRA['ratio'] = sep_SLRA['feature_frequency']/sep_SLRA['Total']
    sep_SLRA['nASVs'] = sep_SLRA['feature_id'].groupby(sep_SLRA['sampleid']).transform('nunique')

    sep_SLRA = sep_SLRA.drop_duplicates()

    #make new df dependingg on plotting needs
    sep_WO = separated[separated.size_code == "W"]
    sep_WO = sep_WO.drop_duplicates()

    sep_PO = separated[separated.size_code == "P"]
    sep_PO = sep_PO.drop_duplicates()

    sep_S = separated[separated.size_code == "S"]
    sep_L = separated[separated.size_code == "L"]


    sep_WO.reset_index(inplace=True, drop=True)
    sep_SLRA.reset_index(inplace=True, drop=True)

    #newseparated = pd.concat([sep_SLRA.reset_index(drop=True), sep_WO.reset_index(drop=True)], axis=0).reset_index(drop=True)
    newseparated = pd.concat([sep_SLRA, sep_WO, sep_PO, sep_L, sep_S], ignore_index=True)

    newseparated['weekdepth'] = newseparated["weekn"].astype(str) + newseparated["depth"].astype(str)
    newseparated['avg'] = newseparated['nASVs'].groupby(newseparated['weekdepth']).transform('mean')
    newseparated['diff'] = newseparated['nASVs'] - newseparated['avg']

    newseparated["rank"] = newseparated.groupby("sampleid")["ratio"].rank(method="average", ascending=False)
    newseparated["ranktot"] = newseparated['rank'] / newseparated['nASVs']

    #calculate shannon diversity index
    grouped = newseparated.groupby('sampleid')['feature_frequency'].apply(list)
    diversity = grouped.apply(shannon)
    newseparated['shannon_diversity'] = newseparated['sampleid'].map(diversity)

    #calculate dnaconc of SL (debugged by chatGPT)
    cleaned_data = newseparated.drop_duplicates(subset=["weekn", "depth", "date", "size_code", "[DNA]ng/uL"])

    sl_fill_values = (
        cleaned_data[cleaned_data["size_code"].isin(["S", "L"])]
        .groupby(["weekn", "depth", "date"], as_index=False)["[DNA]ng/uL"]
        .sum()
    )
    sl_fill_values["size_code"] = "SL"
    newseparated.loc[
        (newseparated["size_code"] == "SL") & (newseparated["[DNA]ng/uL"].isna()),
        "[DNA]ng/uL",] = newseparated.merge(
        sl_fill_values,
        on=["weekn", "depth", "date", "size_code"],
        how="left")["[DNA]ng/uL_y"]

    return newseparated

import pandas as pd

def calculate_normal_range(series, method='std', multiplier=1.0, lower_quantile=0.05, upper_quantile=0.95):
    if not pd.api.types.is_numeric_dtype(series):
        raise ValueError("The series must be numeric")

    if method == 'std':
        mean_val = series.mean()
        std_val = series.std()
        lower = mean_val - multiplier * std_val
        upper = mean_val + multiplier * std_val
    elif method == 'quantile':
        lower = series.quantile(lower_quantile)
        upper = series.quantile(upper_quantile)

    return lower, upper


def make_int_plot(col2plot, depth='all'):
    lab_md = pd.read_csv("/Users/Diana/Downloads/allmetadata_edit_Oct17_2024(allmetadata).csv")
    if depth != 'all':
        lab_md = lab_md[lab_md.depth == depth]

    lab_md['time_string'] = pd.to_datetime(lab_md['time_string'])
    lab_md.set_index('time_string', inplace=True)
    lab_md_agg = lab_md.groupby(lab_md.index).mean() #for duplicate samples take the avg

    ts = lab_md_agg[col2plot]

    # Calculate the normal range for the selected column.
    # Here we use the standard deviation method with a multiplier of 1.0.
    lower_threshold, upper_threshold = calculate_normal_range(ts, method='std', multiplier=1.0)

    week_numbers = ts.index.isocalendar().week

    # --------------------------
    # Create an interactive Plotly figure
    # --------------------------
    fig = go.Figure()

    # Add the time series trace (markers connected by a line)
    fig.add_trace(go.Scatter(
        x=ts.index,
        y=ts,
        mode='lines+markers',
        marker=dict(size=6),  # Adjust marker size as desired
        line=dict(width=2),
        name=col2plot,
        customdata=week_numbers,
        hovertemplate='Date: %{x|%Y-%m-%d}<br>'
    ))

    # Add horizontal dashed line for the lower normal threshold
    fig.add_shape(
        type='line',
        x0=ts.index.min(),
        x1=ts.index.max(),
        y0=lower_threshold,
        y1=lower_threshold,
        line=dict(dash='dash'),
    )

    # Add horizontal dashed line for the upper normal threshold
    fig.add_shape(
        type='line',
        x0=ts.index.min(),
        x1=ts.index.max(),
        y0=upper_threshold,
        y1=upper_threshold,
        line=dict(dash='dash'),
    )

    # Update layout with titles and axis labels
    fig.update_layout(
        xaxis_title='Date',
        yaxis_title=col2plot,
        hovermode='closest'
    )

    # Display the interactive plot
    fig.show()

def all_depths(column_to_plot, mdpath, hidden=None):

    colors = px.colors.qualitative.Plotly
    lab_md = pd.read_csv(mdpath)
    if hidden is not None:
            lab_md = lab_md[~lab_md['depth'].isin(hidden)]

    lab_md['time_string'] = pd.to_datetime(lab_md['time_string'])
    lab_md.set_index('time_string', inplace=True)
    # Create a new figure
    fig = go.Figure()

    # Loop through each unique depth and add a trace
    unique_depths = sorted(lab_md['depth'].unique())
    for i, depth in enumerate(unique_depths):
        # Filter data for this depth
        group = lab_md[lab_md['depth'] == depth].copy()
        group.sort_index(inplace=True)

        # Compute week numbers and day names for hover information
        week_numbers = group.index.isocalendar().week
        customdata = list(zip(week_numbers))

        # Calculate the normal range for this depth's series
        lower_threshold, upper_threshold = calculate_normal_range(group[column_to_plot], method='std', multiplier=1.0)

    # Assign a color for this depth trace
        color = colors[i % len(colors)]

        fig.add_trace(go.Scatter(
            x=group.index,
            y=group[column_to_plot],
            mode='lines+markers',
            marker=dict(size=6, color=color),  # Set marker color
            line=dict(width=2, color=color),     # Set line color
            name=f"Depth {depth}",
            customdata=customdata,
            hovertemplate='Date: %{x|%Y-%m-%d}<br>' +
                          'Week: %{customdata}<br>' +
                          'Depth: ' + str(depth) + '<br>' +
                          'Value: %{y:.2f}<extra></extra>'
        ))

    # Update layout (adjust title and axis labels as needed)
    fig.update_layout(
        xaxis_title='Date',
        yaxis_title=column_to_plot,
        hovermode='closest'
    )

    fig.show()

def reformat_md(md):
    #create a dictionary for months
    month_dic = {
        "Jan": 1,
        "Feb": 2,
        "Mar": 3,
        "Apr": 4,
        "May": 5,
        "Jun": 6,
        "Jul": 7,
        "Aug": 8,
        "Sep": 9,
        "Oct": 10,
        "Nov": 11,
        "Dec": 12
    }
    month_season = {
        "Jan": "Winter",
        "Feb": "Winter",
        "Mar": "Spring",
        "Apr": "Spring",
        "May": "Spring",
        "Jun": "Summer",
        "Jul": "Summer",
        "Aug": "Summer",
        "Sep": "Autumn",
        "Oct": "Autumn",
        "Nov": "Autumn",
        "Dec": "Winter"
    }

    #add month to a new column
    md['month_name'] = md['date'].str.split('-').str[1]

    #add month number
    md['month']= md['month_name'].map(month_dic)

    #add day number
    md['day'] = md['date'].str.split('-').str[0]
    md[["year", "month", "day"]] = md[["year", "month", "day"]].apply(pd.to_numeric)

    #remove symbol for better handling of data
    #md.rename(columns={"Week#": "Weekn"}, inplace=True)
    #md.rename(columns={"Depth": "depth"}, inplace=True) #to match dfo

    #change to int to remove decimals from date columns
    md.year = md.year.apply(int)
    md.depth = md.depth.apply(int)
    md.weekn = md.weekn.apply(int)

    #change to str to aggregate them into time_string to match dfos formatting of the date
    md.year = md.year.apply(str)
    md.month = md.month.apply(str)
    md.day = md.day.apply(str)
    #add leading zero to match date format in dfo metadata
    md['month'] = md['month'].str.zfill(2)
    md['day'] = md['day'].str.zfill(2)

    #add leading zero to match date format in dfo metadata
    md['month'] = md['month'].str.zfill(2)
    md['day'] = md['day'].str.zfill(2)

    md['time_string'] = md[['year', 'month', 'day']].agg('-'.join, axis=1)

    dfo_md = pd.read_csv("/Users/Diana/Documents/escuela/phd/ch2/bb_data/bbmp_aggregated_profiles.csv")
    bio_niskin = pd.read_csv("/Users/Diana/Documents/escuela/phd/ch2/bb_data/BBMP_Data_2022.csv")#
    #dfo_metadata_y14 = pd.read_csv("/Users/Diana/Documents/escuela/phd/bb_data/2019/data_export/trim-analysis/dfo_metadata_y14.tsv", sep='\t')

    #change to str to aggregate them into time_string
    bio_niskin.year = bio_niskin.year.apply(str)
    bio_niskin.month = bio_niskin.month.apply(str)
    bio_niskin.day = bio_niskin.day.apply(str)
    #add leading zero to match date format in dfo metadata
    bio_niskin['month'] = bio_niskin['month'].str.zfill(2)
    bio_niskin['day'] = bio_niskin['day'].str.zfill(2)

    bio_niskin['time_string'] = bio_niskin[['year', 'month', 'day']].agg('-'.join, axis=1)

    #make a new column for time_string without the time
    dfo_md['time_string_time'] = dfo_md['time_string']
    dfo_md['time_string'] = dfo_md['time_string'].str.split(' ').str[0]

    #renaming columns to ensure correct merging
    dfo_md.rename(columns={"depth":"bbmpdepth","pressure": "depth", "year_time": "year", "month_time": "month", "day_time": "day"}, inplace=True)

    #change to int to remove decimals from date columns
    cols = ['year', 'depth', 'month', 'day']
    md[cols] = md[cols].apply(pd.to_numeric, errors='ignore', axis=1)
    dfo_md[cols] = dfo_md[cols].apply(pd.to_numeric, errors='ignore', axis=1)
    bio_niskin[cols] = bio_niskin[cols].apply(pd.to_numeric, errors='ignore', axis=1)

    #make a season column
    md['season'] = ''

    for month, season in month_season.items():
        md.loc[md['month_name'] == month, 'season'] = season

    #merging party
    merged = pd.merge(md, dfo_md, on=["year", "month", "depth", "day"], how="left")
    allyears = pd.merge(md, dfo_md, on=["year", "month", "depth", "day"], how="outer")

    #add nutrient data
    preall_md= pd.merge(allyears, bio_niskin, on=["day", "month", "year", 'depth'], how="outer")
    all_md = pd.merge(merged, bio_niskin, on=["day", "month", "year", 'depth'], how="left")

    #split dfs by depth
    shallow_depths = [1, 5, 10]
    shallow = all_md[all_md["depth"] < 30]
    #shallow = shallow.groupby(['year', 'month', "day"]).mean().reset_index()
    deep = all_md[all_md.depth == 60]

    #split dfs by season
    year_season = preall_md.groupby(by = ['year','season']).mean().reset_index()

    Winter = year_season.loc[year_season['season'] == 'Winter',:]
    Spring = year_season.loc[year_season['season'] == 'Spring',:]
    Summer = year_season.loc[year_season['season'] == 'Summer',:]
    Autumn = year_season.loc[year_season['season'] == 'Autumn',:]

    #save output as csv
    all_md.to_csv('allmetadata.csv')

    return all_md

def plot_adiv(df, depth=None):
    #select only cols of interest for plotting alpha diversity
    copy_of = df[['sampleid', 'feature_id', 'feature_frequency', 'year',
             'weekn', 'depth', 'depth_code', '[DNA]ng/uL', 'size_code', 'month', 'day',
            'season', 'time_string_x', 'time_string', 'temperature', 'Chlorophyll A', 'date', 'Total',
            'nASVs', 'avg', 'diff', 'Taxon', 'Confidence']].copy()
    #drop any null rows
    df_clean = copy_of.dropna(subset=['feature_id'])
    #make sure the feature asv count per sample is accurate
    df_clean['unique_feature_count'] = df_clean.groupby('sampleid')['feature_id'].transform('nunique')
    # for plotting make datetime:
    df_clean['time_string_x'] = pd.to_datetime(df_clean['time_string_x'])

    if depth != 'all':
        df_clean = df_clean[df_clean.depth == depth]

    # Group by date, depth, and size_code, then compute unique feature_id count per group
    ts_data = df_clean.groupby(['time_string_x', 'depth_code', 'size_code'])['feature_id'].nunique().reset_index(name='unique_feature_count')
    plt.figure(figsize=(12, 6))
    sns.lineplot(
        data=ts_data,
        x='time_string_x',
        y='unique_feature_count',
        hue='size_code',
        markers=True,
        dashes=False
    )

    plt.xlabel('Date')
    plt.ylabel('nASVs')
    plt.tight_layout()
    plt.show()


def pca_and_correlation_analysis(df, selected_columns, n_components=2, plot=True, depth=1):
    # Filter to the selected columns and drop rows with missing data
    data = df[df.depth == depth]
    data = data[selected_columns].dropna()

    # Standardize the data so that each variable has mean=0 and std=1
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)

    # Perform PCA
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(data_scaled)

    # Create a DataFrame for the principal components
    pc_columns = [f'PC{i+1}' for i in range(n_components)]
    pc_df = pd.DataFrame(principal_components, columns=pc_columns, index=data.index)

    # Create a DataFrame for the loadings (coefficients of the linear combination)
    loadings = pca.components_.T  # shape: (n_features, n_components)
    loadings_df = pd.DataFrame(loadings, index=selected_columns, columns=pc_columns)

    # Compute the correlation matrix for the selected columns
    corr_matrix = data.corr()

    if plot:
        # Plot the correlation matrix as a heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', fmt=".2f")
        plt.show()

    pca_result = {
        'pca_object': pca,
        'explained_variance': pca.explained_variance_ratio_,
        'loadings': loadings_df,
        'principal_components': pc_df
    }

    return pca_result, corr_matrix


def correlation_pairplot(df, selected_columns, depth):
    # Filter the DataFrame for the selected columns and drop rows with missing data
    data = df[df.depth == depth]
    data = data[selected_columns].dropna()

    # Create the pairplot with histograms on the diagonal
    pairplot = sns.pairplot(data, diag_kind='hist')
    # Adjust the title and layout
    plt.show()


import pandas as pd
import plotly.express as px
import itertools

def taxbarplot_by_year(comm, table, level, depth, topn, include_other=True):
    """
    Create one stacked bar plot per year for the top `topn` taxa at a given depth,
    using only size_codes 'SL' and 'W'. In each plot, the x-axis is week number (1-53)
    and the y-axis shows relative abundance ('ratio'). If include_other is True (default),
    taxa not in the global top list are grouped into an "Other" category. If False, only the top
    taxa are included.

    Color mapping is computed globally so that the same taxon always has the same color across years.

    Parameters:
      comm (str): Community identifier (e.g., 'chloroplast' to trigger a taxonomic level change).
      table (pd.DataFrame): Input DataFrame containing at least the columns:
                            'time_string', 'sampleid', 'size_code', 'depth', 'weekn', level,
                            'feature_frequency', and 'ratio'.
      level (str): Name of the column to use for taxonomy (e.g. 'Taxon').
      depth (int): The depth (as an integer) to be plotted.
      topn (int): Number of top taxa to display.
      include_other (bool): If True, taxa not in the global top list are recoded as "Other".
                            If False, those rows are dropped.

    Returns:
      figures (dict): A dictionary with keys = year and values = Plotly figure objects.
    """
    # (1) If working with chloroplast, switch the taxonomic level
    if comm == 'chloroplast':
        level = 'PRTaxon'

    # (2) Replace unassigned taxa with the feature_id if needed
    table.loc[table[level].isin(['Unassigned', 'g__uncultured', 's__uncultured']), level] = table['feature_id']

    # (3) Update the size_code column if needed (assumed done before calling this function)
    # Filter to only use size_codes 'SL' and 'W'
    df = table[table['size_code'].isin(['SL', 'W'])].copy()

    # (4) Filter to the selected depth
    df = df[df['depth'] == depth].copy()

    # (6) Compute a global top-n taxa list (across all years) based on overall abundance.
    # Here we use 'feature_frequency' to rank taxa.
    agg_all = df.groupby(level).agg({'feature_frequency': 'sum'}).reset_index()
    global_top = agg_all.sort_values('feature_frequency', ascending=False).head(topn)[level].tolist()

    # (7) Build a global palette mapping dictionary using a chosen color sequence.
    # Use Plotly's sequential Plasma palette.
    palette_source = px.colors.sequential.Plasma
    if topn > len(palette_source):
        # If topn is larger than available colors, cycle through the palette.
        palette_colors = list(itertools.islice(itertools.cycle(palette_source), topn))
    else:
        palette_colors = palette_source[:topn]
    palette_dict = {taxon: color for taxon, color in zip(global_top, palette_colors)}
    # Define a color for "Other"
    palette_dict["Other"] = "#A9A9A9"  # DarkGray

    # (8) Get a sorted list of years (excluding any NaN)
    years = sorted(df['year'].dropna().unique())

    figures = {}
    # Loop over each year to generate a plot
    for yr in years:
        df_year = df[df['year'] == yr].copy()
        if df_year.empty:
            continue

        # (a) Depending on include_other, recode or filter taxa:
        if include_other:
            df_year.loc[~df_year[level].isin(global_top), level] = "Other"
        else:
            df_year = df_year[df_year[level].isin(global_top)]

        # (b) Group by week number and taxon, summing the relative abundance ('ratio')
        grouped = df_year.groupby(['weekn', level])['ratio'].sum().reset_index()

        # (c) Ensure every week (1-53) appears: pivot and reindex to get full weeks.
        pivot = grouped.pivot(index='weekn', columns=level, values='ratio').fillna(0)
        pivot = pivot.reindex(range(1, 54), fill_value=0).reset_index()
        melted = pivot.melt(id_vars='weekn', var_name=level, value_name='ratio')

        # (d) Create the stacked bar plot with x = week number.
        fig = px.bar(melted,
                     x='weekn',
                     y='ratio',
                     color=level,
                     labels={'weekn': 'Week', 'ratio': 'Relative Abundance'},
                     color_discrete_map=palette_dict)
        fig.update_xaxes(dtick=1)

        figures[yr] = fig
        fig.show()

    return figures



def taxo_barplot(df, tax_level, topn=10, order='most',
                 abundance_col='feature_frequency', week_col='weekn'):

    df = df.copy()

    df[week_col] = pd.to_numeric(df[week_col], errors='coerce')

    tax_sum = df.groupby(tax_level)[abundance_col].sum().reset_index()
    if order == 'most':
        tax_sum = tax_sum.sort_values(abundance_col, ascending=False)
    elif order == 'least':
        tax_sum = tax_sum.sort_values(abundance_col, ascending=True)
    else:
        raise ValueError("order must be either 'most' or 'least'")

    top_taxa = tax_sum[tax_level].head(topn).tolist()

    df['tax_plot'] = df[tax_level].where(df[tax_level].isin(top_taxa), 'Other')

    total_per_week = df.groupby(week_col)[abundance_col].sum().reset_index().rename(columns={abundance_col: 'total'})

    # Aggregate data by week and tax_plot (summing abundance)
    df_weekly = df.groupby([week_col, 'tax_plot'])[abundance_col].sum().reset_index()
    df_weekly = df_weekly.merge(total_per_week, on=week_col)

    # Calculate relative abundance (this will be between 0 and 1)
    df_weekly['relative_abundance'] = df_weekly[abundance_col] / df_weekly['total']

    # Remove the "Other" category from plotting
    df_weekly_filtered = df_weekly[df_weekly['tax_plot'] != 'Other'].copy()

    df_weekly_filtered.sort_values(week_col, inplace=True)

    fig = px.bar(
        df_weekly_filtered,
        x=week_col,
        y='relative_abundance',
        color='tax_plot',
        title=f'Taxonomic Bar Plot at {tax_level} Level (Top {topn} {order} abundant) by Week',
        labels={'relative_abundance': 'Relative Abundance', week_col: 'Week'}
    )

    fig.update_layout(xaxis=dict(type='linear'))

    return fig


def make_news_memory_efficient(separated):

    # Convert key columns to categorical to save memory
    for col in ['size_code', 'weekn', 'depth', 'year']:
        if separated[col].dtype != 'category':
            separated[col] = separated[col].astype('category')

    # Filter rows that are used in subsequent operations and drop unwanted columns early
    mask = ~separated['size_code'].isin(["W", "P"])
    sep_SL = separated.loc[mask].copy()

    sample_dna = sep_SL[['sampleid', 'weekn', 'depth', 'year', '[DNA]ng/uL']].drop_duplicates('sampleid')

    group_sum = sample_dna.groupby(['year', 'weekn', 'depth'], observed=True)['[DNA]ng/uL'].sum().reset_index()

    group_sum.rename(columns={'[DNA]ng/uL': '[DNAt]'}, inplace=True)

    sep_SL = sep_SL.merge(group_sum, on=['year', 'weekn', 'depth'], how='left')
    sep_SL['DNApr'] = sep_SL['[DNA]ng/uL'] / sep_SL['[DNAt]']

    # Prepare a minimal lookup DataFrame for merging; drop duplicates immediately
    merge_df = sep_SL[['sampleid', 'DNApr', '[DNAt]']].drop_duplicates()

    # Free memory from sep_SL if not needed further
    del sep_SL

    # Merge back with the original DataFrame using a lightweight lookup table
    sepSLRA = separated.merge(merge_df, on='sampleid', how='left')

    # Drop rows not needed and filter further in place
    sepSLRA = sepSLRA.loc[~sepSLRA['size_code'].isin(["W", "P"])].copy()

    # Compute corrected feature frequency with groupby, then drop the intermediate column
    sepSLRA['Newfeature_frequency'] = sepSLRA['feature_frequency'] * sepSLRA['DNApr']
    sepSLRA['Newff'] = sepSLRA.groupby(['feature_id', 'weekn', 'depth', 'year'], observed=True)['Newfeature_frequency'].transform('sum')

    # Update sampleid and set size_code in one step
    sepSLRA['sampleid'] = (
    "BB" + sepSLRA['year'].astype(str).str[-2:] + "." +
    sepSLRA['weekn'].astype(str) +
    sepSLRA['depth_code'] + "SL")
    sepSLRA['size_code'] = 'SL'

    # Drop unneeded columns early
    cols_to_drop = ['feature_frequency', 'Total', 'ratio', 'nASVs', 'weekdepth', 'avg', 'year',
                    'diff', '[DNA]ng/uL', 'Newfeature_frequency']
    sepSLRA.drop(columns=cols_to_drop, inplace=True, errors='ignore')
    sepSLRA.rename(columns={'Newff': 'feature_frequency'}, inplace=True)
    sepSLRA.drop_duplicates(inplace=True)

    # Calculate totals and ratios in place
    sepSLRA['Total'] = sepSLRA.groupby('sampleid', observed=True)['feature_frequency'].transform('sum')
    sepSLRA['ratio'] = sepSLRA['feature_frequency'] / sepSLRA['Total']
    sepSLRA['nASVs'] = sepSLRA.groupby('sampleid', observed=True)['feature_id'].transform('nunique')
    sepSLRA.drop_duplicates(inplace=True)

    # Process other groups; filtering only when needed
    sep_WO = separated.loc[separated['size_code'] == "W"].drop_duplicates().copy()
    sep_PO = separated.loc[separated['size_code'] == "P"].drop_duplicates().copy()
    sep_S  = separated.loc[separated['size_code'] == "S"].copy()
    sep_L  = separated.loc[separated['size_code'] == "L"].copy()

    # Reset indices without creating extra copies
    sep_WO.reset_index(drop=True, inplace=True)
    sepSLRA.reset_index(drop=True, inplace=True)

    # Concatenate all parts at once
    newseparated = pd.concat([sepSLRA, sep_WO, sep_PO, sep_L, sep_S], ignore_index=True)

    # Compute combined columns and group statistics
    newseparated['weekdepth'] = newseparated["weekn"].astype(str) + newseparated["depth"].astype(str)
    newseparated['avg'] = newseparated.groupby('weekdepth', observed=True)['nASVs'].transform('mean')
    newseparated['diff'] = newseparated['nASVs'] - newseparated['avg']
    newseparated["rank"] = newseparated.groupby("sampleid", observed=True)["ratio"].rank(method="average", ascending=False)
    newseparated["ranktot"] = newseparated['rank'] / newseparated['nASVs']

    # Calculate shannon diversity index per sample (assuming shannon is defined)
    diversity = newseparated.groupby('sampleid', observed=True)['feature_frequency'].apply(lambda freqs: shannon(freqs))
    newseparated['shannon_diversity'] = newseparated['sampleid'].map(diversity)

    # For DNA concentration in SL, compute and merge as needed
    cleaned_data = newseparated.drop_duplicates(subset=["weekn", "depth", 'year', "size_code", "[DNA]ng/uL"])
    sl_fill_values = (
        cleaned_data.loc[cleaned_data["size_code"].isin(["S", "L"])]
        .groupby(["weekn", "depth", 'year'], observed=True)["[DNA]ng/uL"]
        .sum()
        .reset_index()
    )
    sl_fill_values["size_code"] = "SL"

    newseparated = newseparated.merge(
        sl_fill_values,
        on=["weekn", "depth", 'year', "size_code"],
        how="left",
        suffixes=('', '_fill')
    )
    newseparated.loc[newseparated["[DNA]ng/uL"].isna(), "[DNA]ng/uL"] = newseparated.loc[newseparated["[DNA]ng/uL"].isna(), "[DNA]ng/uL_fill"]
    newseparated.drop(columns=["[DNA]ng/uL_fill"], inplace=True)

    return newseparated

##### for metadata analysis;
def find_non_numeric_columns(df, columns):
    problem_cols = []
    for col in columns:
        # Select values that are NOT numeric
        non_numeric_vals = df[~df[col].apply(lambda x: isinstance(x, (int, float, np.number)) or pd.isna(x))][col]
        if not non_numeric_vals.empty:
            problem_cols.append(col)
            print(f"Column '{col}' has non-numeric values. Examples:")
            print(non_numeric_vals.head(5))
            print()
    if not problem_cols:
        print("All columns are numeric or contain only NaNs.")
    return problem_cols

def list_outliers_per_season_year(heatmap_data_scaled, z_threshold=2):
    """
    Given a DataFrame of z-scored anomalies (variables x season_year),
    returns a dict mapping each season_year to a list of variable names
    that are outliers (absolute z-score > threshold).
    """
    outliers = {}
    for season_year in heatmap_data_scaled.columns:
        # Select variables where anomaly magnitude exceeds threshold
        vars_outlier = heatmap_data_scaled.index[heatmap_data_scaled[season_year].abs() > z_threshold].tolist()
        outliers[season_year] = vars_outlier
    return outliers

def identify_outlier_season_years(heatmap_data_scaled, z_threshold=2):
    """
    Given a DataFrame of z-scored anomalies (variables x season_year),
    returns a Series counting how many variables exceed the threshold
    (positive or negative) per season_year.
    """
    # Boolean mask of outliers per variable-season_year
    outliers_mask = heatmap_data_scaled.abs() > z_threshold

    # Count how many variables are outliers per season_year (column)
    outlier_counts = outliers_mask.sum(axis=0)

    # Sort descending
    outlier_counts = outlier_counts.sort_values(ascending=False)

    return outlier_counts

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def zscore_by_season(heatmap_df: pd.DataFrame) -> pd.DataFrame:
    """
    Season-specific z-scoring at a given depth.
    rows = variables, cols = 'Season-YYYY'
    For each season, compute mean/std across all years and z-score columns of that season.
    """
    seasons = pd.Series(heatmap_df.columns).str.split('-', n=1, expand=True)[0]
    z = pd.DataFrame(index=heatmap_df.index, columns=heatmap_df.columns, dtype=float)

    for season in seasons.unique():
        cols = heatmap_df.columns[seasons.values == season]
        block = heatmap_df[cols]  # variables x years for this season

        mu = block.mean(axis=1, skipna=True)
        sd = block.std(axis=1, ddof=1, skipna=True).replace(0, np.nan)

        z[cols] = (block.sub(mu, axis=0)).div(sd, axis=0)

    return z


def plot_variable_with_seasonal_norm(df, variable, depth, years='2014-2022'):
    df = df.copy()
    df['date'] = pd.to_datetime(df['date'])
    df = df.sort_values('date')

    # Filter by depth
    df = df[df['depth_code'] == depth]

    # Parse years input and filter
    if isinstance(years, str) and '-' in years:
        start_year, end_year = map(int, years.split('-'))
        df = df[(df['date'].dt.year >= start_year) & (df['date'].dt.year <= end_year)]
    elif isinstance(years, list):
        df = df[df['date'].dt.year.isin(years)]

    # Convert variable to numeric, coerce errors to NaN
    df[variable] = pd.to_numeric(df[variable], errors='coerce')

    # Extract month for seasonal norm
    df['month'] = df['date'].dt.month

    # Calculate seasonal norm = average per month across all selected years at this depth
    seasonal_norm = df.groupby('month')[variable].mean()

    # Resample to monthly frequency to explicitly show gaps (NaNs for missing months)
    df = df.set_index('date')
    df_var = df[variable].resample('M').mean()

    # Map seasonal norm back to dates in resampled df
    df_var = df_var.to_frame()
    df_var['month'] = df_var.index.month
    df_var['seasonal_norm'] = df_var['month'].map(seasonal_norm)

    # Plot with markers and lines, matplotlib will break lines on NaNs
    plt.figure(figsize=(14,6))
    plt.plot(df_var.index, df_var[variable], marker='o', linestyle='-', label='Observed', alpha=0.7)
    plt.plot(df_var.index, df_var['seasonal_norm'], color='red', linewidth=2, label='Seasonal Norm (Monthly Avg)')
    plt.xlabel('Date')
    plt.ylabel(variable)
    plt.title(f'{variable} Time Series depth {depth} ({years})')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'figures/{variable}_plot_depth{depth}.png', dpi=300)
    plt.show()

def shannon_diversity(counts):
    counts = counts[counts > 0]
    return entropy(counts, base=np.e)

def evenness(counts):
    counts = counts[counts > 0]
    if len(counts) <= 1:
        return 0.0
    return shannon_diversity(counts) / np.log(len(counts))

def rarefy_and_calc_metrics(counts, depth):
    if counts.sum() < depth:
        return np.nan, np.nan, np.nan  # Cannot rarefy
    # Subsample without replacement
    rarefied = np.random.choice(counts.index.repeat(counts.values), size=depth, replace=False)
    rare_counts = pd.Series(rarefied).value_counts()

    richness = rare_counts.count()
    shannon = shannon_diversity(rare_counts.values)
    even = evenness(rare_counts.values)

    return richness, shannon, even

def repeat_rarefy(table, depth, n_iter=100):
    results = []

    for i in range(n_iter):
        metrics = table.apply(lambda x: rarefy_and_calc_metrics(x, depth), axis=1)
        metrics_df = pd.DataFrame(metrics.tolist(), index=table.index, columns=['richness', 'shannon', 'evenness'])
        results.append(metrics_df)

    combined = pd.concat(results, axis=0, keys=range(n_iter))
    summary = combined.groupby(level=1).agg(['mean', 'std'])
    return summary

def build_comparison_networks(df,
                              excluded_years=[2019],
                              weeks=[36, 37, 38, 39, 40, 41],
                              max_asvs=50,
                              tax_level='Genus',
                              size_codes=['W', 'SL'],
                              depth_value=1,
                              corr_threshold=0.6):

    # 1. Filter base data
    df['datetime'] = pd.to_datetime(df['datetime'])
    df['weekn'] = df['datetime'].dt.isocalendar().week
    df['year'] = df['datetime'].dt.year

    base_mask = (
        (df['depth'] == depth_value) &
        (df['size_code'].isin(size_codes)) &
        (df['weekn'].isin(weeks))
    )

    df_excluded = df[base_mask & (df['year'].isin(excluded_years))]
    df_rest = df[base_mask & (~df['year'].isin(excluded_years))]

    def make_asv_table(sub_df):
        pivot = sub_df.groupby(['sampleid', 'feature_id'])['[DNAt]'].sum().unstack(fill_value=0)
        return pivot

    def preprocess_and_filter(asv_table):
        # Keep top ASVs by total abundance
        top_asvs = asv_table.sum(axis=0).sort_values(ascending=False).head(max_asvs).index
        filtered = asv_table[top_asvs]
        return filtered

    def build_network(asv_matrix):
        corr = asv_matrix.corr(method='spearman')
        np.fill_diagonal(corr.values, 0)
        G = nx.Graph()
        for i in corr.index:
            for j in corr.columns:
                if abs(corr.loc[i, j]) >= corr_threshold:
                    G.add_edge(i, j, weight=corr.loc[i, j])
        return G

    def assign_taxonomy_labels(G, mapping):
        for node in G.nodes:
            taxon = mapping.get(node, "Unknown")
            G.nodes[node]['taxon'] = taxon

    def plot_network(G, title, filename):
        pos = nx.spring_layout(G, seed=42)
        taxa = [G.nodes[n]['taxon'] for n in G.nodes]
        unique_taxa = list(set(taxa))
        color_map = dict(zip(unique_taxa, sns.color_palette("hsv", len(unique_taxa))))
        node_colors = [color_map[G.nodes[n]['taxon']] for n in G.nodes]

        plt.figure(figsize=(10, 10))
        nx.draw_networkx_nodes(G, pos, node_size=100, node_color=node_colors, alpha=0.85)
        nx.draw_networkx_edges(G, pos, alpha=0.4)
        nx.draw_networkx_labels(G, pos, labels={n: G.nodes[n]['taxon'] for n in G.nodes}, font_size=6)
        plt.title(title)
        plt.axis('off')

        # Save the figure
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.show()

    # Build ASV tables
    asv_excluded = preprocess_and_filter(make_asv_table(df_excluded))
    asv_rest = preprocess_and_filter(make_asv_table(df_rest))

    # Correlation-based networks
    G_excluded = build_network(asv_excluded)
    G_rest = build_network(asv_rest)

    # Map feature_id to taxon level
    feature_to_tax = df.drop_duplicates('feature_id').set_index('feature_id')[tax_level].to_dict()

    assign_taxonomy_labels(G_excluded, feature_to_tax)
    assign_taxonomy_labels(G_rest, feature_to_tax)

    # Plot
    plot_network(G_excluded, f"Network: Year(s) {excluded_years}", f"network_years_{'_'.join(map(str, excluded_years))}.png")
    plot_network(G_rest, f"Network: All Years Except {excluded_years}", f"network_all_except_{'_'.join(map(str, excluded_years))}.png")


    return G_excluded, G_rest

def build_comparison_networks_flexible(df,
                                       excluded_years=[2019],
                                       weeks=[36, 37, 38, 39, 40, 41],
                                       max_asvs=50,
                                       tax_level='Genus',
                                       size_codes=['W', 'SL'],
                                       depth_value=1,
                                       corr_threshold=0.6,
                                       strict_size_code=True,
                                       output_prefix="network"):

    df['datetime'] = pd.to_datetime(df['datetime'])
    df['weekn'] = df['datetime'].dt.isocalendar().week
    df['year'] = df['datetime'].dt.year


    if strict_size_code:
        size_mask = df['size_code'].isin(size_codes)
    else:
        size_mask = df['size_code'].notna()  # allow all available codes

    base_mask = (
        (df['depth'] == depth_value) &
        size_mask &
        (df['weekn'].isin(weeks))
    )

    df_excluded = df[base_mask & (df['year'].isin(excluded_years))]
    df_rest = df[base_mask & (~df['year'].isin(excluded_years))]

    def make_asv_table(sub_df):
        if sub_df.empty:
            return pd.DataFrame()
        pivot = sub_df.groupby(['sampleid', 'feature_id'])['[DNAt]'].sum().unstack(fill_value=0)
        return pivot

    def preprocess_and_filter(asv_table):
        if asv_table.empty:
            return pd.DataFrame()
        top_asvs = asv_table.sum(axis=0).sort_values(ascending=False).head(max_asvs).index
        return asv_table[top_asvs]

    def build_network(asv_matrix):
        if asv_matrix.empty or asv_matrix.shape[1] < 2:
            return nx.Graph()

        # Apply log1p transformation
        log_asv = np.log1p(asv_matrix)

        # Pearson correlation on log-transformed data
        corr = log_asv.corr(method='pearson')
        np.fill_diagonal(corr.values, 0)

        G = nx.Graph()
        for i in corr.index:
            for j in corr.columns:
                if abs(corr.loc[i, j]) >= corr_threshold:
                    G.add_edge(i, j, weight=corr.loc[i, j])
        return G

    def assign_taxonomy_labels(G, mapping):
        for node in G.nodes:
            taxon = mapping.get(node, "Unknown")
            G.nodes[node]['taxon'] = taxon

    def plot_network(G, title, filename):
        if G.number_of_nodes() == 0:
            print(f"[!] Skipping plot for {title} — no valid nodes.")
            return
        pos = nx.spring_layout(G, seed=42)
        taxa = [G.nodes[n]['taxon'] for n in G.nodes]
        unique_taxa = list(set(taxa))
        color_map = dict(zip(unique_taxa, sns.color_palette("hsv", len(unique_taxa))))
        node_colors = [color_map[G.nodes[n]['taxon']] for n in G.nodes]

        plt.figure(figsize=(10, 10))
        nx.draw_networkx_nodes(G, pos, node_size=100, node_color=node_colors, alpha=0.85)
        nx.draw_networkx_edges(G, pos, alpha=0.4)
        nx.draw_networkx_labels(G, pos, labels={n: G.nodes[n]['taxon'] for n in G.nodes}, font_size=6)
        plt.title(title)
        plt.axis('off')
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.show()

    # Process
    asv_excluded = preprocess_and_filter(make_asv_table(df_excluded))
    asv_rest = preprocess_and_filter(make_asv_table(df_rest))

    G_excluded = build_network(asv_excluded)
    G_rest = build_network(asv_rest)

    # Get taxonomy map
    feature_to_tax = df.drop_duplicates('feature_id').set_index('feature_id')[tax_level].to_dict()
    assign_taxonomy_labels(G_excluded, feature_to_tax)
    assign_taxonomy_labels(G_rest, feature_to_tax)

    # Plot and save
    year_label = '_'.join(map(str, excluded_years))
    size_label = '_'.join(size_codes) if strict_size_code else 'auto'

    plot_network(G_excluded, f"Excluded Year(s): {year_label}", f"{output_prefix}_excluded_{year_label}_{size_label}.png")
    plot_network(G_rest, f"All Except: {year_label}", f"{output_prefix}_rest_not_{year_label}_{size_label}.png")

    return G_excluded, G_rest

def build_taxonomic_networks(df,
                              excluded_years=[2019],
                              weeks=[36, 37, 38, 39, 40, 41],
                              max_nodes=50,
                              tax_level='Genus',
                              size_codes=['W', 'SL'],
                              depth_value=1,
                              corr_threshold=0.6,
                              strict_size_code=True,
                              output_prefix="network"):

    # Prepare datetime and filtering columns
    df['datetime'] = pd.to_datetime(df['datetime'])
    df['weekn'] = df['datetime'].dt.isocalendar().week
    df['year'] = df['datetime'].dt.year

    # Filter size code
    size_mask = df['size_code'].isin(size_codes) if strict_size_code else df['size_code'].notna()

    # Base filter for all data
    base_mask = (
        (df['depth'] == depth_value) &
        size_mask &
        (df['weekn'].isin(weeks))
    )

    df_excluded = df[base_mask & (df['year'].isin(excluded_years))]
    df_rest = df[base_mask & (~df['year'].isin(excluded_years))]

    # Function to create a taxon-level abundance table
    def make_tax_table(sub_df):
        if sub_df.empty:
            return pd.DataFrame()
        pivot = sub_df.groupby(['sampleid', tax_level])['ratio'].sum().unstack(fill_value=0)
        return pivot

    # Keep only top N most abundant taxonomic groups
    def filter_top_taxa(table):
        if table.empty:
            return pd.DataFrame()
        top_taxa = table.sum(axis=0).sort_values(ascending=False).head(max_nodes).index
        return table[top_taxa]

    # Build correlation-based network
    def build_network(table):
        if table.empty or table.shape[1] < 2:
            return nx.Graph()
        log_table = np.log1p(table)
        corr = log_table.corr(method='pearson')
        np.fill_diagonal(corr.values, 0)
        G = nx.Graph()
        for i in corr.index:
            for j in corr.columns:
                if abs(corr.loc[i, j]) >= corr_threshold:
                    G.add_edge(i, j, weight=corr.loc[i, j])
        return G

    # Plot and save network
    def plot_network(G, title, filename):
        if G.number_of_nodes() == 0:
            print(f"[!] Skipping plot for {title} — no valid nodes.")
            return
        pos = nx.spring_layout(G, seed=42)
        nodes = list(G.nodes)
        color_map = dict(zip(nodes, sns.color_palette("hsv", len(nodes))))
        node_colors = [color_map[n] for n in nodes]

        plt.figure(figsize=(6, 6))
        nx.draw_networkx_nodes(G, pos, node_size=100, node_color=node_colors, alpha=0.85)
        nx.draw_networkx_edges(G, pos, alpha=0.4)
        nx.draw_networkx_labels(G, pos, font_size=7)
        plt.title(title)
        plt.axis('off')
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.show()

    # Process both groups
    tax_table_excluded = filter_top_taxa(make_tax_table(df_excluded))
    tax_table_rest = filter_top_taxa(make_tax_table(df_rest))

    G_excluded = build_network(tax_table_excluded)
    G_rest = build_network(tax_table_rest)

    # Plot and save
    year_label = '_'.join(map(str, excluded_years))
    size_label = '_'.join(size_codes) if strict_size_code else 'auto'

    plot_network(G_excluded,
                 f"{tax_level} network – year(s): {year_label}",
                 f"{output_prefix}_{tax_level.lower()}_excluded_{year_label}_{size_label}.png")

    plot_network(G_rest,
                 f"{tax_level} network – all but: {year_label}",
                 f"{output_prefix}_{tax_level.lower()}_not_{year_label}_{size_label}.png")

    return G_excluded, G_rest

def repeat_rarefaction(df, sample_col='sampleid', feature_col='feature_id', count_col='feature_frequency', depth=20000, repeats=100, random_state=None):
    """
    Perform repeat-rarefaction on a long-format DataFrame.

    Parameters:
      df: long-format DataFrame.
      sample_col: column name for sample IDs.
      feature_col: column name for feature IDs.
      count_col: column name for counts.
      depth: rarefaction depth.
      repeats: number of rarefaction iterations.
      random_state: optional random seed.

    Returns:
      A DataFrame of the average rarefied counts (samples x features) across repeats.
    """
    df_wide = df.pivot_table(index=sample_col, columns=feature_col, values=count_col, fill_value=0)

    all_reps = []
    for i in range(repeats):
        # Optionally vary the random seed for each repeat.
        rep = df_wide.apply(lambda row: rarefy_sample(row.values, depth, random_state=(random_state+i if random_state is not None else None)), axis=1)
        all_reps.append(rep)

    # Convert each replicate (a DataFrame) into a list and then average across replicates.
    rep_dfs = [pd.DataFrame(list(rep), index=df_wide.index, columns=df_wide.columns) for rep in all_reps]
    mean_rarefied = sum(rep_dfs) / repeats
    return mean_rarefied

def analyze_depth_distances(df, depth_code='all', size_codes=None, metric='euclidean',
                            std_multiplier=2, output_dir='figures'):
    '''
    Analyze normalized distances between consecutive weekly samples for specified depth codes.

    Parameters:
    -----------
    df: pd.DataFrame
        DataFrame containing columns ['sampleid','feature_id','ratio','size_code','depth_code','year','weekn'].
    depth_code: str
        The depth code to filter by, or 'all' to process all depth codes.
    size_codes: list of str, optional
        List of size codes to include (default ['W','SL']).
    metric: str
        Distance metric ('euclidean' or 'braycurtis').
    std_multiplier: float
        Multiplier for standard deviation threshold.
    output_dir: str
        Directory to save figures.

    Returns:
    --------
    dict or list
        If depth_code is 'all', returns dict mapping each depth code to list of sample IDs above threshold.
        If depth_code is specific, returns list of sample IDs above threshold.
    '''
    os.makedirs(output_dir, exist_ok=True)
    df = df.copy()
    # Ensure year is numeric and fill missing years based on sampleid
    df['year'] = pd.to_numeric(df['year'], errors='coerce')
    missing = df['year'].isna()
    if missing.any():
        df.loc[missing, 'year'] = (
            df.loc[missing, 'sampleid'].str.extract(r'BB(\\d{2})')[0]
            .astype(float) + 2000
        )
    # Default size codes
    if size_codes is None:
        size_codes = ['W', 'SL']
    # Determine depths to process
    depths = df['depth_code'].unique() if depth_code == 'all' else [depth_code]
    results = {}
    for d in depths:
        # Filter for size and depth
        df_f = df.query('size_code in @size_codes and depth_code == @d')
        if df_f.empty:
            continue
        # Pivot to samples x features
        df_wide = df_f.pivot_table(index='sampleid', columns='feature_id', values='ratio', fill_value=0)
        # CLR transformation
        data_repl = multiplicative_replacement(df_wide.values)
        data_clr = clr(data_repl)
        df_clr = pd.DataFrame(data_clr, index=df_wide.index, columns=df_wide.columns)
        # Compute distances
        if metric == 'euclidean':
            dists = squareform(pdist(df_clr.values, metric='euclidean'))
            ids = df_clr.index
        else:
            bc_dm = beta_diversity('braycurtis', df_wide.values, ids=df_wide.index)
            dists = bc_dm.data
            ids = bc_dm.ids
        dist_matrix = DistanceMatrix(dists, ids)
        dist_df = dist_matrix.to_data_frame()
        # Convert to long format
        dist_long = (
            dist_df.reset_index()
            .melt(id_vars='index', var_name='sample2', value_name='distance')
        )
        dist_long = dist_long[dist_long['index'] != dist_long['sample2']]
        dist_long.rename(columns={'index': 'sample1'}, inplace=True)
        # Metadata for timepoints and dates
        meta = df_f[['sampleid', 'year', 'weekn']].drop_duplicates().copy()
        meta['year'] = pd.to_numeric(meta['year'], errors='coerce')
        meta['weekn'] = pd.to_numeric(meta['weekn'], errors='coerce')
        meta['timepoint'] = meta['year'] * 52 + meta['weekn']
        meta['date'] = pd.to_datetime(
            meta['year'].astype(int).astype(str) + '-W' +
            meta['weekn'].astype(int).astype(str).str.zfill(2) + '-1',
            format='%Y-W%W-%w'
        )
        # Merge metadata
        merged = (
            dist_long
            .merge(
                meta[['sampleid', 'timepoint', 'date']]
                .rename(columns={'sampleid': 'sample1', 'timepoint': 'tp1', 'date': 'date1'}),
                on='sample1'
            )
            .merge(
                meta[['sampleid', 'timepoint', 'date']]
                .rename(columns={'sampleid': 'sample2', 'timepoint': 'tp2', 'date': 'date2'}),
                on='sample2'
            )
        )
        # Filter consecutive weeks
        df_consec = merged[merged['tp2'] == merged['tp1'] + 1].copy()
        # Normalize distances
        min_d = df_consec['distance'].min()
        max_d = df_consec['distance'].max()
        df_consec['norm_distance'] = 2 * (df_consec['distance'] - min_d) / (max_d - min_d) - 1
        mean = df_consec['norm_distance'].mean()
        std = df_consec['norm_distance'].std()
        # Plot
        sns.set_theme(context='talk', style='whitegrid', font_scale=1.3)
        plt.figure(figsize=(12, 8))
        sns.lineplot(data=df_consec, x='date2', y='norm_distance', marker='o', linewidth=2)
        plt.axhline(mean, color='red', linewidth=2)
        x_min = df_consec['date2'].min()
        x_max = df_consec['date2'].max()
        plt.fill_between([x_min, x_max], [mean - std_multiplier * std] * 2, [mean + std_multiplier * std] * 2, color='red', alpha=0.3)
        # Extract above-threshold sample IDs
        threshold = mean + std_multiplier * std
        above = df_consec[df_consec['norm_distance'] > threshold]
        above_ids = above['sample2'].unique()
        # Save figure
        plt.xlabel('Date')
        plt.ylabel('Normalized Distance')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'depth_{d}_normalized_distances.png'), bbox_inches='tight')
        plt.close()
        # Store results
        results[d] = list(above_ids)
    # Return results
    return results if depth_code == 'all' else results.get(depth_code, [])

SEASON_ORDER = ['Winter', 'Spring', 'Summer', 'Fall']
SEASON_REGEX = r'^(Winter|Spring|Summer|Fall)-(\d{4})$'

def _assign_season(ts: pd.Timestamp, mode: str = "astronomical") -> str:
    if mode not in {"astronomical", "meteorological"}:
        raise ValueError("mode must be 'astronomical' or 'meteorological'")

    y = ts.year
    if mode == "astronomical":
        if (ts >= pd.Timestamp(year=y-1, month=12, day=21)) and (ts < pd.Timestamp(year=y, month=3, day=20)):
            return 'Winter'
        elif (ts >= pd.Timestamp(year=y, month=3, day=20)) and (ts < pd.Timestamp(year=y, month=6, day=21)):
            return 'Spring'
        elif (ts >= pd.Timestamp(year=y, month=6, day=21)) and (ts < pd.Timestamp(year=y, month=9, day=22)):
            return 'Summer'
        elif (ts >= pd.Timestamp(year=y, month=9, day=22)) and (ts < pd.Timestamp(year=y, month=12, day=21)):
            return 'Fall'
        else:
            return 'Winter'
    else:  # meteorological, i.e. split  by month instead of solstice
        m = ts.month
        if m in (12, 1, 2):
            return 'Winter'
        elif m in (3, 4, 5):
            return 'Spring'
        elif m in (6, 7, 8):
            return 'Summer'
        else:  # 9,10,11
            return 'Fall'

def _make_season_year(row, mode: str = "astronomical") -> str:
    #create a new column of season + year identification for each sample
    #december is the following year's winter i.e. dec 2012 = winter 2013
    y = row['date'].year
    season = row['season']
    if season == 'Winter' and row['date'].month == 12:
        y += 1
    return f"{season}-{y}"

def _season_from_label(col: str): #this will extract the SEASON and the YEAR is distinct columns from the new 'season_year' column created previously
    s = str(col).strip()
    m = re.match(SEASON_REGEX, s)
    if m:
        return m.group(1), int(m.group(2))
    return (np.nan, np.nan)


def plot_seasonal_outlier_heatmap(df, variables, years='all', season_mode: str = "astronomical"):
    df = df.copy()
    df['date'] = pd.to_datetime(df['date'])

    # Restrict to 2014–2022 unless overridden
    df = df[(df['date'].dt.year >= 2014) & (df['date'].dt.year <= 2022)]

    # Ensure variables exist & numeric
    for v in variables:
        if v not in df.columns:
            df[v] = np.nan
    df[variables] = df[variables].apply(pd.to_numeric, errors='coerce')

    # Optional year filter
    if years != 'all':
        df = df[df['date'].dt.year.isin(years)]

    # === Add nutrient & biomass ratios ===
    def _pick_first_existing(*candidates):
        for c in candidates:
            if c in df.columns:
                return c
        return None

    # Pick columns across naming schemes (CERC vs BIO). All coerced numeric above.
    col_NO3 = _pick_first_existing('Nitrate_BIO_bbmp', 'NO3_BIO_bbmp')
    col_NO2 = _pick_first_existing('Nitrite_BIO_bbmp', 'NO2_BIO_bbmp')
    col_NH4 = _pick_first_existing('Ammonia_BIO_bbmp', 'NH4_BIO_bbmp', 'NH3_BIO_bbmp')
    col_PO4 = _pick_first_existing('Phosphate_BIO_bbmp', 'PO4_BIO_bbmp')
    col_Si  = _pick_first_existing('Silicate_BIO_bbmp', 'Si_BIO_bbmp')
    col_Chl = _pick_first_existing('Chlorophyll_A_BIO_bbmp', 'HPLC_CHLA_BIO_bbmp')

    # Build DIN & DIP (µM assumed; if units differ, ratios still dimensionless but interpret cautiously)
    DIN = 0
    n_found = 0
    for c in [col_NO3, col_NO2, col_NH4]:
        if c:
            DIN = DIN + df[c].astype(float)
            n_found += 1
    DIN = DIN if n_found else np.nan

    DIP = df[col_PO4].astype(float) if col_PO4 else np.nan
    Si  = df[col_Si].astype(float)  if col_Si  else np.nan
    NO3 = df[col_NO3].astype(float) if col_NO3 else np.nan
    Chl = df[col_Chl].astype(float) if col_Chl else np.nan

    # Safe division helper
    def _safe_div(a, b):
        with np.errstate(divide='ignore', invalid='ignore'):
            out = a / b
        # Replace inf, -inf, and invalid values (like NaN) with np.nan
        if isinstance(out, pd.Series):
            out = out.replace([np.inf, -np.inf], np.nan)
        else:
            out = np.where(np.isfinite(out), out, np.nan)
        return out
    # Ratios
    df['ratio_DIN_DIP'] = _safe_div(DIN, DIP)
    df['ratio_Si_DIN']  = _safe_div(Si, DIN)
    df['ratio_NO3_Chl'] = _safe_div(NO3, Chl)

    # If you also have Total N / Total P columns, use them for N:P; otherwise fall back to DIN:DIP
    col_TN = _pick_first_existing('TN_mean_CERC_nuts', 'TotalN_BIO_bbmp', 'TN_BIO_bbmp')
    col_TP = _pick_first_existing('TP_mean_CERC_nuts', 'TotalP_BIO_bbmp', 'TP_BIO_bbmp')
    if col_TN and col_TP:
        df['ratio_N_P'] = _safe_div(df[col_TN].astype(float), df[col_TP].astype(float))
    else:
        df['ratio_N_P'] = df['ratio_DIN_DIP']  # pragmatic fallback

    # Add new ratios to variables so they appear in the heatmap (keep original order first)
    ratio_vars = ['ratio_N_P', 'ratio_DIN_DIP', 'ratio_Si_DIN', 'ratio_NO3_Chl']
    for rv in ratio_vars:
        if rv not in variables:
            variables.append(rv)

    # Season and Season-Year (now using chosen season_mode)
    df['season'] = df['date'].apply(lambda d: _assign_season(d, mode=season_mode))
    df['season_year'] = df.apply(lambda r: _make_season_year(r, mode=season_mode), axis=1)

    heatmaps = {}

    # Helper: chronological sort key for 'Season-YYYY'
    def _sort_key(season_year_label: str):
        season, year = _season_from_label(season_year_label)
        return (year if not pd.isna(year) else 0,
                SEASON_ORDER.index(season) if season in SEASON_ORDER else -1)

    for depth in df['depth_code'].unique():
        df_depth = df[df['depth_code'] == depth].copy()

        # 1) Seasonal averages (Season-Year) for each variable
        seasonal_avgs = df_depth.groupby('season_year')[variables].mean()

        # 2) Clean/normalize the season-year index
        seasonal_avgs.index = seasonal_avgs.index.astype(str).str.strip()

        # 3) Order columns by chronological Season-Year
        ordered_cols = sorted(seasonal_avgs.index, key=_sort_key)
        seasonal_avgs = seasonal_avgs.reindex(ordered_cols)

        # 4) Variables as rows, Season-Year as columns
        heatmap_data = seasonal_avgs.T.copy()

        # --- Debug: show any columns that don't look like Season-YYYY ---
        bad_cols = [c for c in heatmap_data.columns if re.match(SEASON_REGEX, str(c).strip()) is None]
        if bad_cols:
            print(f"[depth {depth}] Columns not matching 'Season-YYYY' and kept as-is (will be skipped in season z-score):")
            for b in bad_cols:
                print("   ", b)

        # 5) Season-specific z-scoring (Fall vs Falls, etc.)
        z = pd.DataFrame(index=heatmap_data.index, columns=heatmap_data.columns, dtype=float)
        col_to_season = {c: _season_from_label(c)[0] for c in heatmap_data.columns}

        for season in SEASON_ORDER:
            season_cols = [c for c, s in col_to_season.items() if s == season]
            if not season_cols:
                continue
            block = heatmap_data[season_cols]
            mu = block.mean(axis=1, skipna=True)
            sd = block.std(axis=1, ddof=1, skipna=True).replace(0, np.nan)
            z[season_cols] = (block.sub(mu, axis=0)).div(sd, axis=0)
            z = z.dropna(how='all')

        # 6) Plot
        if z.empty or z.dropna(how='all').empty:
            print(f"No data to plot for depth {depth}. Skipping heatmap.")
            continue

        plt.figure(figsize=(max(10, len(z.columns)*0.8), len(variables)*0.4 + 2))
        sns.heatmap(z, center=0, cmap='coolwarm',
                    cbar_kws={'label': 'Season-specific z-score', 'shrink': 0.5})
        plt.title(f'Depth {depth} ({season_mode.title()} seasons)')
        plt.xlabel('Season-Year')
        plt.ylabel('Variable')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        outname = f"figures/heatmap_depth_{depth}_{season_mode}.png"
        plt.savefig(outname, dpi=300, bbox_inches='tight')
        plt.close()

        heatmaps[depth] = z

    return heatmaps

def plot_alpha_richness_by_depth(
    df,
    depths=None,
    size_codes=('SL','W'),
    rarefy_depth=10000,
    metric='mean',                  # 'mean' or 'median' to match your precomputed col
    week_col='weekn',
    date_col='datetime',
    depth_col='depth',
    size_col='size_code',
    year_col='year',                # will be created from date_col if missing
    output_dir='figures/alpha_trends',
    filename_stub='alpha_depth{depth}_rarefy{r}_metric{metric}.png'
):
    """
    Creates one PNG per depth:
      - x = weekn (1..53)
      - y = richness (repeat-rarefied), pulled from richness_{metric}_{rarefy_depth}
      - hue/style = year; yearly lines are dotted
      - thick solid black line = across-years weekly mean (the 'normal' trend)

    Assumes df already contains repeat-rarefied richness columns like:
        richness_mean_{rarefy_depth} or richness_median_{rarefy_depth}
    """
    os.makedirs(output_dir, exist_ok=True)

    # Which richness column to use
    raref_col = f"richness_{metric}_{int(rarefy_depth)}"
    if raref_col not in df.columns:
        raise KeyError(
            f"Column '{raref_col}' not found. "
            f"Expected a precomputed repeat-rarefied richness column like "
            f"'richness_mean_{rarefy_depth}' or 'richness_median_{rarefy_depth}'."
        )

    # Ensure required columns
    needed = {depth_col, size_col, week_col, raref_col}
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise KeyError(f"Dataframe is missing required columns: {missing}")

    # Add year if needed
    dfx = df.copy()
    if year_col not in dfx.columns:
        if date_col not in dfx.columns:
            raise KeyError(
                f"'{year_col}' not in df and cannot create it because '{date_col}' is missing."
            )
        dfx[date_col] = pd.to_datetime(dfx[date_col])
        dfx[year_col] = dfx[date_col].dt.year

    # Default depths
    if depths is None:
        depths = sorted(pd.unique(dfx[depth_col].dropna()))

    # Nice seaborn defaults
    sns.set(style="whitegrid", context="talk")

    for depth in depths:
        # Filter depth + size
        dfd = dfx[(dfx[depth_col] == depth) & (dfx[size_col].isin(size_codes))].copy()
        if dfd.empty:
            continue

        # Keep only rows with a valid week number
        dfd = dfd[dfd[week_col].notna()].copy()

        # Compute per-(year, week) average richness (averaging over replicates/size-codes if both present)
        yearly_week = (
            dfd.groupby([year_col, week_col], observed=True, as_index=False)[raref_col]
               .mean()
               .sort_values([year_col, week_col])
        )

        if yearly_week.empty:
            continue

        # Across-years weekly mean (the "normal trend")
        weekly_mean = (
            yearly_week.groupby(week_col, observed=True, as_index=False)[raref_col]
                       .mean()
                       .sort_values(week_col)
                       .rename(columns={raref_col: 'weekly_mean'})
        )

        # Build dashed mapping so ALL yearly lines are dotted
        years = sorted(yearly_week[year_col].unique())
        dash_map = {y: (2, 2) for y in years}  # dotted pattern
        # Palette for years
        palette = sns.color_palette("tab10", n_colors=len(years))
        pal_map = {y: palette[i % len(palette)] for i, y in enumerate(years)}

        # Nice seaborn defaults: white background, no grid
        sns.set(style="white", context="talk")

        # Plot with square figure size
        plt.figure(figsize=(11, 7))

        # All yearly dotted lines
        sns.lineplot(
            data=yearly_week,
            x=week_col,
            y=raref_col,
            hue=year_col,
            style=year_col,
            dashes=dash_map,
            palette=pal_map,
            linewidth=1.7,
            alpha=0.9,
            legend=True,
            errorbar=None,
        )

        # Overlay thick solid mean line
        plt.plot(
            weekly_mean[week_col],
            weekly_mean['weekly_mean'],
            linewidth=3.5,
            color='black',
            label='Across-years mean'
        )

        # Labels & cosmetics
        plt.title(f"Alpha diversity (richness) • Depth {depth} m • Rarefied @ {rarefy_depth}")
        plt.xlabel("Week number")
        plt.ylabel(f"Richness (repeat-rarefied {metric} @ {rarefy_depth})")
        plt.xlim(int(np.nanmin(yearly_week[week_col])), int(np.nanmax(yearly_week[week_col])))

        # Legend outside the plot
        handles, labels = plt.gca().get_legend_handles_labels()
        if 'Across-years mean' not in labels:
            handles.append(plt.Line2D([], [], color='black', linewidth=3.5))
            labels.append('Across-years mean')
        plt.legend(
            handles, labels,
            title='Year',
            bbox_to_anchor=(1.02, 1),
            loc='upper left',
            frameon=False
        )

        plt.tight_layout()

        # Make sure the mean appears in legend
        handles, labels = plt.gca().get_legend_handles_labels()
        if 'Across-years mean' not in labels:
            handles.append(plt.Line2D([], [], color='black', linewidth=3.5))
            labels.append('Across-years mean')
        plt.legend(handles, labels, title='Year', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True)

        # Save
        outpath = os.path.join(
            output_dir,
            filename_stub.format(depth=str(depth), r=int(rarefy_depth), metric=metric)
        )
        plt.savefig(outpath, dpi=300, bbox_inches='tight')
        plt.close()

    return



def pcaplot(separated, depth, comm, columnperm, spc, colrow):

    if depth == 'all':
        df = separated.copy()
    else:
        df=separated[separated.depth==depth]


    if 'SL' in separated['size_code'].unique():
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W', 'SL']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}
        dicsc = pd.Series(df.size_code.values,index=df.sampleid).to_dict()
        color_rows_sc = {k: palette_dict[v] for k, v in dicsc.items()}
        seriescr = pd.Series(color_rows_sc)

    else:
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}
        dicsc = pd.Series(df.size_code.values,index=df.sampleid).to_dict()
        color_rows_sc = {k: palette_dict[v] for k, v in dicsc.items()}
        seriescr = pd.Series(color_rows_sc)

    #month palette code
    if colrow == 'Month':
        df['Month'] = df['date'].str.split('-').str[1]
        months = ['Jan', 'Feb', 'Mar', 'May', 'Apr']
        palette_colors = sns.color_palette("flare")
        palette_dict_month = {monthname: color for monthname, color in zip(months, palette_colors)}
        dic = pd.Series(df.Month.values,index=df.sampleid).to_dict()
        color_rows_month = {k: palette_dict_month[v] for k, v in dic.items()}
        seriesmonthcr = pd.Series(color_rows_month)
        dfcolors = pd.DataFrame({'Month': seriesmonthcr,'Size code':seriescr})

    else:
        df['weekn2'] = df['weekn'].astype(str)
        weeks = list(df['weekn'].unique())
        weeks.sort()
        weeks = [str(x) for x in weeks]
        palette_colors = sns.color_palette("flare", 16)
        palette_dict_weekn = {weekname: color for weekname, color in zip(weeks, palette_colors)}
        dic = pd.Series(df.weekn2.values,index=df.sampleid).to_dict()
        color_rows_weekn = {k: palette_dict_weekn[v] for k, v in dic.items()}
        seriesweekn = pd.Series(color_rows_weekn)
        dfcolors = pd.DataFrame({'Weekn': seriesweekn, 'Size code': seriescr})

    #select data to pivot
    topiv = df[['feature_id', 'feature_frequency', 'sampleid']].copy()
    topiv = topiv.drop_duplicates()

    sfdpiv= topiv.pivot(index='sampleid', columns='feature_id', values='feature_frequency')
    sfdpiv=sfdpiv.fillna(0)
    sfdclr=sfdpiv.mask(sfdpiv==0).fillna(0.1)
    clr_transformed_array = clr(sfdclr)
    samples = sfdpiv.index
    asvs = sfdpiv.columns

    #Creating the dataframe with the clr transformed data, and assigning the sample names
    clr_transformed = pd.DataFrame(clr_transformed_array, columns=asvs)
    #Assigning the asv names
    clr_transformed['samples'] = samples
    clr_transformed = clr_transformed.set_index('samples')
    clr_transformed.head()

    #calculate distance matrix
    dist = cdist(clr_transformed, clr_transformed, 'euclid')
    distance_matrix = pd.DataFrame(dist, columns=samples)
    distance_matrix['samples'] = samples
    distance_matrix = distance_matrix.set_index('samples')

    #format for pca
    dm = DistanceMatrix(distance_matrix)

    pca = PCA(n_components=2)
    pca_features = pca.fit_transform(distance_matrix)

    ####
    sns.set(rc={"figure.figsize":(4, 3)})
    sns.set_style("whitegrid", {'axes.grid' : False})
    plot_df = pd.DataFrame(data = pca_features, columns = ['dim1', 'dim2'], index = sfdpiv.index)
    plot_df['dim1'] = plot_df['dim1']/1000
    plot_df['dim2'] = plot_df['dim2']/1000
    if depth =='all':
        plot_df2 = pd.merge(plot_df,df[['sampleid','size_code','depth']],on='sampleid', how='left')
    else:
        plot_df2 = pd.merge(plot_df,df[['sampleid','size_code','weekn']],on='sampleid', how='left')


    ##divide into pre-post bloom
    def get_stage(weekNb):
        if weekNb < 8:
            return 'Pre-bloom'
        elif weekNb >= 8:
            return 'Bloom'

    if depth != 'all':
        plot_df2['Time'] = plot_df2['weekn'].apply(get_stage)

    plot_df2 = plot_df2.rename(columns={'size_code': 'Size code'})

    pc1v = round(pca.explained_variance_ratio_[0]*100)
    pc2v = round(pca.explained_variance_ratio_[1]*100)

    #plot_df2 = plot_df2.drop_duplicates()
    #dfperm = plot_df2.set_index('sampleid')

    #permanova2 = permanova(dm, dfperm, columnperm)
    #results = permanova2(999)

    #plot

    if depth == 'all':
        var2 = 'depth'
    else:
        var2 = 'Time'

    sns.set_style("white")
    ax=sns.scatterplot(x = 'dim1', y = 'dim2', hue= 'Size code', style=var2, data = plot_df2,
                       palette=palette_dict) #, size = 'Week_Group')#,palette=sns.color_palette("dark:salmon_r", as_cmap=True))
    plt.ylabel('PCo2') #(' + str(pc2v) + '% variance explained)')
    plt.xlabel('PCo1') #(' + str(pc1v) +'% variance explained)')
    ax.set_title('Depth ' + str(depth) + 'm', loc='left', weight='bold')
    plt.legend(frameon=False)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    sns.despine()
    plt.savefig('outputs/'+comm+'/D'+str(depth)+spc+'_PCAplot.png', dpi=200, bbox_inches="tight")
    plt.clf()
    plt.cla()
    plt.close()

    print ( "Components = ", pca.n_components_ , ";\nTotal explained variance = ",
      round(pca.explained_variance_ratio_.sum(),5)  )

    print ("Components 1 and 2 are", pca.explained_variance_ratio_)

    # Retrieve Loadings
    loadings = pca.components_

    # Summarize Loadings by Metadata Category
    metadata_groups = plot_df2[var2].unique()
    metadata_contributions = {}

    for group in metadata_groups:
        group_variables = plot_df2.loc[plot_df2[var2] == group, 'sampleid']
        group_loadings = np.abs(loadings[:, [list(distance_matrix.columns).index(var) for var in group_variables]]).mean(axis=1)
        metadata_contributions[group] = group_loadings

    ##clustermap
    ax = sns.clustermap(distance_matrix, method="complete", cmap='RdBu', annot=True,
               yticklabels=True, row_colors = dfcolors,
               annot_kws={"size": 7}, figsize=(15,12));
    if colrow == 'Month':
        handles1 = [Patch(facecolor=palette_dict_month[key]) for key in palette_dict_month]
        plt.legend(handles1, palette_dict_month, title='Month',
                   bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper left')

    else:
        handles1 = [Patch(facecolor=palette_dict_weekn[key]) for key in palette_dict_weekn]
        plt.legend(handles1, palette_dict_weekn, title='Week',
                bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper left')


    plt.savefig('outputs/'+comm+'/D'+str(depth)+spc+'_clustermap.png', dpi=200, bbox_inches="tight")
    plt.clf()
    plt.cla()
    plt.close()


    plot_df2['PCo 1'] = pc1v
    plot_df2['PCo 2'] = pc2v
    plot_df2.rename(columns={'Size code': 'size_code'}, inplace=True)
    plot_df2.to_csv('R_results/R_testing_vis/'+str(depth)+'_'+comm+'for_R.csv')

    return pca, pca_features, sfdclr, distance_matrix

def plot_prevalence_abundance_faceted(
    df, taxonomy_df, *,
    abundance_col="ratio",
    feature_col="feature_id",
    sample_col="sampleid",
    phylum_col="Phylum",
    datetime_col="datetime",
    depth=None,
    season=None,
    phyla=None,
    min_presence=0.0,
    color_by=None,                # None → one color per facet (by phylum)
    logx=True,
    point_alpha=0.35,
    point_size=18,
    save_fig=None,
    save_legend=None,
    legend_title=None,
    legend_ncol=1,
    random_state=42
):
    # --- prepare & compute ---
    df = df.copy()
    df[datetime_col] = pd.to_datetime(df[datetime_col], errors="coerce")
    if "year" not in df.columns:
        df["year"] = df[datetime_col].dt.year
    if depth is not None and "depth" in df.columns:
        df = df[df["depth"] == depth]
    if season is not None and "season" in df.columns:
        df = df[df["season"].astype(str).str.lower() == str(season).lower()]

    total_samples = df[sample_col].nunique()
    if total_samples == 0:
        raise ValueError("No samples after filtering.")

    stats = (df.groupby(feature_col, observed=True)
               .agg(total_abundance=(abundance_col, "sum"),
                    n_samples_present=(abundance_col, lambda x: (x > min_presence).sum()))
               .reset_index())
    stats["prevalence"] = stats["n_samples_present"] / total_samples

    tax = taxonomy_df[[feature_col, phylum_col]].drop_duplicates()
    stats = stats.merge(tax, on=feature_col, how="left").dropna(subset=[phylum_col])

    if phyla:
        stats = stats[stats[phylum_col].isin(phyla)]
    stats = stats.sample(frac=1.0, random_state=random_state)

    # --- palettes & legend setup (respect input order) ---
    if phyla:
        # keep only those phyla present and preserve the user-specified order
        phylum_list = [p for p in phyla if p in stats[phylum_col].unique()]
    else:
        # if user didn’t specify, fall back to sorted order
        phylum_list = list(stats[phylum_col].unique())

    # assign colors in that order
    phylum_palette = sns.color_palette("tab20", n_colors=len(phylum_list))
    phylum_to_color = dict(zip(phylum_list, phylum_palette))


    if (color_by is not None) and (color_by not in stats.columns):
        # attach from df if needed
        if color_by in df.columns:
            stats = stats.merge(df[[feature_col, color_by]].drop_duplicates(),
                                on=feature_col, how="left")
        else:
            raise ValueError(f"color_by='{color_by}' not found.")

    # common x limits
    xvals = stats["total_abundance"].replace(0, np.nan).dropna()
    xmin = xvals.min() if not xvals.empty else 1
    xmax = xvals.max() if not xvals.empty else 1

    # --- build grid (5 per row, shared axes) ---
    if color_by is None:
        g = sns.FacetGrid(
            stats, col=phylum_col, col_wrap=5, height=3.5,
            sharex=True, sharey=True
        )

        def facet_scatter(data, color=None, **kwargs):
            ph = data[phylum_col].iloc[0]
            sns.scatterplot(
                data=data, x="total_abundance", y="prevalence",
                color=phylum_to_color.get(ph, "gray"),
                alpha=point_alpha, s=point_size, edgecolor="none", **kwargs
            )
    else:
        # hue by season/year/etc. (legend saved separately)
        levels = sorted(pd.Series(stats[color_by].dropna().unique()).tolist())
        hue_palette = sns.color_palette("tab20", n_colors=len(levels))
        level_to_color = {lvl: col for lvl, col in zip(levels, hue_palette)}

        g = sns.FacetGrid(
            stats, col=phylum_col, col_wrap=5, height=3.5,
            sharex=True, sharey=True, hue=color_by,
            hue_order=levels, palette=level_to_color
        )

        def facet_scatter(data, color=None, **kwargs):
            sns.scatterplot(
                data=data, x="total_abundance", y="prevalence",
                alpha=point_alpha, s=point_size, edgecolor="none", **kwargs
            )

    g.map_dataframe(facet_scatter)

    # --- axes styling (common) ---
    n = len(phylum_list)
    ncols = 5
    nrows = int(np.ceil(n / ncols))
    for i, ax in enumerate(g.axes.flatten()[:n]):
        # limits & scales
        if logx:
            ax.set_xscale("log")
        ax.set_xlim(left=max(xmin, 1e-9))  # safe lower bound for log
        ax.set_ylim(0, 1.05)

        # strip default titles
        ax.set_title("")

        # custom strip: full-width rectangle + text with phylum name (no 'p__')
        label = str(phylum_list[i]).replace("p__", "")
        ax.add_patch(Rectangle(
            (0, 1.02), 1.0, 0.14, transform=ax.transAxes, clip_on=False,
            facecolor="#f0f0f0", edgecolor="none"
        ))
        ax.text(0.5, 1.09, label, transform=ax.transAxes,
                ha="center", va="center", fontsize=11, weight="bold")

        # show y only on leftmost column
        col_idx = i % ncols
        if col_idx != 0:
            ax.set_ylabel("")
            ax.tick_params(labelleft=False)

        # show x only on bottom row
        row_idx = i // ncols
        if row_idx != (nrows - 1):
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)

    # remove any extra empty axes (when panels not multiple of 5)
    for ax in g.axes.flatten()[n:]:
        ax.remove()

    # remove legend from the grid
    if g._legend is not None:
        g._legend.remove()

    # figure title & layout
    g.fig.suptitle("Abundance vs. Prevalence by Phylum", y=1.02)
    g.fig.tight_layout()

    # --- save/show main figure ---
    if save_fig:
        g.fig.savefig(save_fig, dpi=300, bbox_inches="tight")
        plt.close(g.fig)
    else:
        plt.show()

    # --- separate legend figure ---
    if save_legend:
        if color_by is None:
            handles = [
                plt.Line2D([], [], marker='o', linestyle='', color='w',
                           markerfacecolor=phylum_to_color[p], markersize=6, label=p.replace("p__",""))
                for p in phylum_list
            ]
            labels = [p.replace("p__","") for p in phylum_list]
            title = legend_title or "Phylum"
        else:
            handles = [
                plt.Line2D([], [], marker='o', linestyle='', color='w',
                           markerfacecolor=g._legend_data[color_by].get(lvl, "#666")
                           if hasattr(g, "_legend_data") else level_to_color[lvl],
                           markersize=6, label=str(lvl))
                for lvl in levels
            ]
            labels = [str(l) for l in levels]
            title = legend_title or str(color_by).capitalize()

        fig_leg = plt.figure(figsize=(max(2.0, 0.35 * max(1, len(labels))), 1.2))
        fig_leg.legend(handles=handles, labels=labels, loc="center",
                       frameon=False, ncol=legend_ncol, title=title)
        fig_leg.tight_layout()
        fig_leg.savefig(save_legend, dpi=300, bbox_inches="tight")
        plt.close(fig_leg)

    return stats


SEASON_ORDER = ['Winter', 'Spring', 'Summer', 'Fall']
_SEASON_TO_RANK = {s: i for i, s in enumerate(SEASON_ORDER)}

def add_season_columns(df: pd.DataFrame,
                       datetime_col: str = 'datetime',
                       mode: str = 'astronomical',
                       season_col: str = 'season',
                       season_year_col: str = 'season_year',
                       copy: bool = False) -> pd.DataFrame:
    """
    Add 'season' and 'season_year' columns using an existing `_assign_season(ts, mode)`.

    - `season` uses astronomical or meteorological logic (per your _assign_season).
    - `season_year` rolls December to the *next* year's Winter (e.g., 2012-12-28 -> Winter-2013).
    - Categorical ordering is applied to `season` for sensible sorting/plotting.

    Returns the modified DataFrame (or a copy if copy=True).
    """
    out = df.copy() if copy else df

    # ensure datetime dtype
    out[datetime_col] = pd.to_datetime(out[datetime_col], errors='coerce')
    if out[datetime_col].isna().any():
        bad = out[out[datetime_col].isna()].index[:5].tolist()
        raise ValueError(f"Some {datetime_col} values could not be parsed to datetime (example bad indices: {bad})")

    # season via your helper
    out[season_col] = out[datetime_col].apply(lambda ts: _assign_season(ts, mode=mode))

    # year, with December→next-year for Winter
    year = out[datetime_col].dt.year
    month = out[datetime_col].dt.month
    is_dec_winter = (out[season_col] == 'Winter') & (month == 12)
    adj_year = year.where(~is_dec_winter, year + 1)

    out[season_year_col] = out[season_col] + '-' + adj_year.astype(str)

    # nice categorical for plotting/sorting
    out[season_col] = pd.Categorical(out[season_col], categories=SEASON_ORDER, ordered=True)

    # (optional) add a sortable helper for season_year if you ever need stable ordering
    out['_season_rank'] = out[season_col].map(_SEASON_TO_RANK).astype('Int64')
    out['_season_year_sortkey'] = list(zip(adj_year, out['_season_rank']))

    return out


# Optional: if you already have _PREFIX, re-use it; otherwise define lightweight prefixes here
_PREFIX = {"Phylum":"p__", "Class":"c__", "Order":"o__", "Family":"f__", "Genus":"g__"}

_UNASSIGNED_RE = re.compile(r"^\s*(unassigned|unknown|unclassified|na|none|null)?\s*$", re.IGNORECASE)

def _is_unassigned(value: str) -> bool:
    if value is None:
        return True
    v = str(value).strip()
    if not v:
        return True
    # treat prefix-only like 'p__' as unassigned
    if re.match(r"^[kpcfoags]__?$", v, flags=re.IGNORECASE):
        return True
    # treat strings that are literally "Unassigned"/"Unknown"/etc. as unassigned
    if _UNASSIGNED_RE.match(v):
        return True
    # treat 'x__Unassigned' or 'x__' with blank tail as unassigned
    if re.match(r"^[a-z]__\s*(unassigned|unknown|unclassified)?\s*$", v, re.IGNORECASE):
        return True
    return False

def _clean_rank_label(text: str) -> str:
    """Remove common rank prefixes and prettify underscores/spaces."""
    if text is None:
        return ""
    s = str(text).strip()
    s = re.sub(r"^[a-z]__\s*", "", s, flags=re.IGNORECASE)  # drop p__/c__/o__/...
    s = re.sub(r"_+", " ", s).strip()
    return s

def _build_fallback_label_map(taxonomy_csv: str, *, feature_col: str, rank_preference=("Order","Class","Phylum")):
    """
    Returns:
      label_map: {feature_id -> "<label> (<short-id>)"}
      rank_used_map: {feature_id -> rank name actually used, or 'Unassigned'}
    """
    want_cols = {feature_col, *rank_preference}
    try:
        tdf = pd.read_csv(taxonomy_csv, usecols=lambda c: c in want_cols, dtype={feature_col: "string"})
    except Exception:
        tdf = pd.read_csv(taxonomy_csv, dtype={feature_col: "string"})

    have_ranks = [r for r in rank_preference if r in tdf.columns]
    if feature_col not in tdf.columns:
        raise KeyError(f"[taxonomy] Expected feature column '{feature_col}' in {taxonomy_csv}")

    tdf = tdf.drop_duplicates(subset=feature_col).copy()

    label_map = {}
    rank_used_map = {}

    for _, row in tdf.iterrows():
        fid = row[feature_col]
        if pd.isna(fid):
            continue

        short_id = str(fid)[:6]
        chosen_label = None
        chosen_rank = "Unassigned"

        for rank in have_ranks:
            raw = row.get(rank)
            if _is_unassigned(raw):
                continue
            cleaned = _clean_rank_label(raw)
            if cleaned and not _is_unassigned(cleaned):
                chosen_label = cleaned
                chosen_rank = rank
                break

        if not chosen_label:
            chosen_label = "Unassigned"

        # Always append the short feature ID
        full_label = f"{chosen_label} ({short_id})"

        label_map[fid] = full_label
        rank_used_map[fid] = chosen_rank

    return label_map, rank_used_map


def _union_topk_depth(df, tax_key, depth, seasons, years, top_k):
    sub = df[(df["depth"]==depth) & (df["season"].isin(seasons)) & (df["year"].isin(years))]
    taxa = []
    for (y, s), g in sub.groupby(["year","season"]):
        top = (g.sort_values("median_rel_abund", ascending=False)
                 .head(top_k)[tax_key].tolist())
        taxa.extend(top)
    return list(OrderedDict.fromkeys(taxa))

def _global_color_map(taxa_union):
    base = plt.rcParams['axes.prop_cycle'].by_key().get('color', [])
    if len(taxa_union) > len(base):
        base = list(base) + list(plt.cm.get_cmap('tab20').colors)
    # to hex for portability
    import matplotlib.colors as mcolors
    return {t: mcolors.to_hex(base[i % len(base)]) for i, t in enumerate(taxa_union)}

def _build_legend_label_map(taxonomy_csv, feature_col="feature_id", genus_col="Genus", family_col="Family"):
    """
    Returns dict: feature_id -> legend label like 'g__Genus_abc12' or 'f__Family_abc12'.
    If both genus/family missing, uses 'u__Unassigned_abc12'.
    """
    t = pd.read_csv(taxonomy_csv, usecols=lambda c: c in {feature_col, genus_col, family_col})
    t = t.drop_duplicates(subset=feature_col).copy()
    def _label(row):
        asv = str(row.get(feature_col, ""))[:5] if pd.notna(row.get(feature_col)) else "na"
        gen = str(row.get(genus_col)) if pd.notna(row.get(genus_col)) and str(row.get(genus_col)).strip() else ""
        fam = str(row.get(family_col)) if pd.notna(row.get(family_col)) and str(row.get(family_col)).strip() else ""
        if gen:
            base = f"{gen}"
        elif fam:
            base = f"{fam}"
        else:
            base = "u__Unassigned"
        return f"{base}_{asv}"
    lbl = {row[feature_col]: _label(row) for _, row in t.iterrows()}
    return lbl

# ---------- main ----------
def grouped_by_season_stacked_by_years_per_depth_w_legcsv(
    baseline_long_csv,
    taxonomy_csv="notebooks/newseparated_sep23.csv",
    tax_key="feature_id",
    season_order=("winter","spring","summer","autumn"),
    top_k=5,
    as_percent=True,
    save_path=None,
    legend_csv_path=None,
    year_label_fontsize=8,
    season_label_fontsize=10,
    group_gap=0.8,
    within_gap=0.08,
    exclude_depths=(30,),
    year_label_step=2,
    # ---------- NEW: tell the legend which rank(s) to try, in fallback order ----------
    legend_rank_preference=("Order","Class","Phylum")
):
    df = pd.read_csv(baseline_long_csv)
    req = {"year","season","depth", tax_key, "median_rel_abund"}
    miss = req - set(df.columns)
    if miss: raise KeyError("Missing columns: %s" % miss)

    # order and filter (unchanged) ...
    df = df[df["season"].isin(season_order)].copy()
    df["season"] = pd.Categorical(df["season"], categories=list(season_order), ordered=True)
    years  = sorted(df["year"].dropna().unique())
    depths = sorted([d for d in df["depth"].dropna().unique() if d not in exclude_depths],
                    key=lambda x: (pd.isna(x), x))

    # per-depth taxa + colors (unchanged) ...
    per_depth_taxa, all_taxa = {}, []
    for dp in depths:
        per_depth_taxa[dp] = _union_topk_depth(df, tax_key, dp, season_order, years, top_k)
        all_taxa += per_depth_taxa[dp]
    global_union = list(OrderedDict.fromkeys(all_taxa))
    cmap = _global_color_map(global_union)

    # ---------- CHANGED: taxonomy-based labels with fallback ----------
    label_map, rank_used_map = _build_fallback_label_map(
        taxonomy_csv,
        feature_col=tax_key,
        rank_preference=legend_rank_preference
    )

    # figure & geometry (unchanged) ...
    nrows = len(depths)
    fig = plt.figure(figsize=(12, 3.8*nrows), facecolor="white")
    gs  = gridspec.GridSpec(nrows, 2, width_ratios=[4, 1], wspace=0.25, hspace=0.45)
    fig.subplots_adjust(bottom=0.12)
    n_years   = len(years)
    n_seasons = len(season_order)
    target_group_width = 2.2
    bar_w = max(0.1, min(0.9, (target_group_width - (n_years-1)*within_gap) / max(n_years,1)))
    group_span    = target_group_width + group_gap
    group_centers = np.arange(n_seasons) * group_span
    offsets       = (np.arange(n_years) - (n_years-1)/2.0) * (bar_w + within_gap)

    legend_rows = []
    for r, dp in enumerate(depths):
        ax = fig.add_subplot(gs[r, 0])
        ax_legend = fig.add_subplot(gs[r, 1]); ax_legend.axis("off")
        ax.set_facecolor("white")

        sub_d = df[df["depth"]==dp]
        taxa = per_depth_taxa.get(dp, [])
        if sub_d.empty or not taxa:
            ax.axis("off"); continue

        # draw grouped stacks (unchanged) ...
        for s_idx, s in enumerate(season_order):
            s_data = (sub_d[sub_d["season"]==s]
                      .pivot_table(index="year", columns=tax_key,
                                   values="median_rel_abund", aggfunc="mean")
                      .reindex(index=years, columns=taxa).fillna(0.0))
            if as_percent: s_data *= 100.0

            for j, yr in enumerate(years):
                vec = s_data.loc[yr].values if yr in s_data.index else np.zeros(len(taxa))
                x = group_centers[s_idx] + offsets[j]
                bottom = 0.0
                for t_idx, t in enumerate(taxa):
                    val = vec[t_idx]
                    if val == 0: continue
                    ax.bar(x, val, width=bar_w, bottom=bottom, color=cmap.get(t, "#999999"), edgecolor="none")
                    bottom += val

                if s == season_order[0] and (j % year_label_step == 0):
                    yy = str(int(yr))[-2:]
                    ax.text(
                        x, -0.12, yy,
                        transform=ax.get_xaxis_transform(),  # y in axes coords
                        ha="center", va="top",
                        fontsize=year_label_fontsize,
                        rotation=0,
                        clip_on=False
                    )

        ax.set_title(f"Depth {dp}", fontsize=12, pad=8)
        ax.set_ylabel("Relative abundance" + (" (%)" if as_percent else ""))
        ax.set_xticks(group_centers)
        ax.set_xticklabels([s.capitalize() for s in season_order], fontsize=season_label_fontsize)
        ax.grid(False)

        # ---------- CHANGED: legend title reflects fallback chain ----------
        handles = [Patch(facecolor=cmap[t], edgecolor='none') for t in taxa]
        row_labels = [label_map.get(t, f"u__Unassigned_{str(t)[:5]}") for t in taxa]
        pref_str = " → ".join(legend_rank_preference)
        ncols_leg = 1 if len(taxa) <= 20 else 2
        ax_legend.legend(
            handles, row_labels,
            title=f"Legend (pref: {pref_str})\nDepth {dp}",
            fontsize=7, title_fontsize=8,
            loc="upper left", frameon=False, ncol=ncols_leg
        )

        # record legend rows (now with rank_used)
        for t, lbl in zip(taxa, row_labels):
            legend_rows.append({
                "depth": dp,
                tax_key: t,
                "legend_label": lbl,
                "color_hex": cmap.get(t, "#999999"),
                "rank_used": rank_used_map.get(t, "Unassigned")
            })

    plt.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")

    if legend_csv_path:
        leg_df = pd.DataFrame(legend_rows)
        if not leg_df.empty:
            agg = (leg_df.groupby(tax_key, as_index=False)
                   .agg(legend_label=("legend_label","first"),
                        color_hex=("color_hex","first"),
                        rank_used=("rank_used","first"),
                        depths=("depth", lambda x: ",".join(map(str, sorted(set(x))))) )
                   .sort_values([tax_key]))
            os.makedirs(os.path.dirname(legend_csv_path) or ".", exist_ok=True)
            agg.to_csv(legend_csv_path, index=False)
            print(f"[legend] Saved legend table: {legend_csv_path}")
        else:
            print("[legend] No taxa to write.")
    return fig
