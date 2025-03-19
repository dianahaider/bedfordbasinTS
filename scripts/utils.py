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

from skbio.diversity.alpha import shannon

# Special thanks to Alex Manuele https://github.com/alexmanuele
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
            dataframes.append(df)

            # Adding table_id, forward and reverse trim columns
            #df['table_id'] = str(table_path.split('/')[-3]) #add a table_id column
            #df['forward_trim'], df['reverse_trim'] = df['table_id'].str.split('R', 1).str
            #df['forward_trim'] = df['forward_trim'].map(lambda x: x.lstrip('F'))
            #df["forward_trim"] = pd.to_numeric(df["forward_trim"])
            #df["reverse_trim"] = pd.to_numeric(df["reverse_trim"])

    #Stick all the dataframes together
    #outputfile="merged_all_tables.tsv"
    df = pd.concat(dataframes)
    #df.to_csv(comm+'/merged_all_tables.tsv', sep='\t', index=False)
    print("Successfully saved all tables.")
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

def pick_metadata(merged, depth='all', size_fraction='both', year='all', R='all', F='all', txsubset = 'all'):
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



    #separated['total'] = separated.groupby(['table_id','sample-id'])['feature_frequency'].transform('sum')
    #separated['ratio'] = separated['feature_frequency']/(separated['total'])
    #separated_taxonomies = separated.copy()

    #make a dictionary with keys for id-ing the taxon belonging to this sub-community
    #separated_dic = pd.Series(separated.Taxon.values,separated.feature_id.values).to_dict()
    print('Saved separated by metadata dataframe.')

    return separated

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
                              'diff', 'extraction_date', '[DNA]ng/uL', 'A260/280',
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

def all_depths(column_to_plot, hidden=None):

    colors = px.colors.qualitative.Plotly
    lab_md = pd.read_csv("/Users/Diana/Downloads/allmetadata_edit_Oct17_2024(allmetadata).csv")
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
