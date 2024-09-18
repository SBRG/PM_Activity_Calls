import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join,isdir
from scipy import stats
from scipy.signal import savgol_filter
from get_plates_from_study import get_plates_from_study
import scipy.stats as st
import seaborn as sns
import re

# All wells in the biolog plate
wells =     ['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12',
             'B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12',
             'C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12',
             'D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12',
             'E01','E02','E03','E04','E05','E06','E07','E08','E09','E10','E11','E12',
             'F01','F02','F03','F04','F05','F06','F07','F08','F09','F10','F11','F12',
             'G01','G02','G03','G04','G05','G06','G07','G08','G09','G10','G11','G12',
             'H01','H02','H03','H04','H05','H06','H07','H08','H09','H10','H11','H12']

def get_files_dir(path):
    """
    Takes in path and returns all files and directories in a path
    """
    files = [f for f in listdir(path) if isfile(join(path, f))]
    directories = [f for f in listdir(path) if isdir(join(path, f))]
    
    return files,directories
    
    
def curate_files(onlyfiles,path,strain_names_to_change):
    """
    Corrects column names from "field xyz" to strain, media, etc
    """
    for file in onlyfiles:
        temp = pd.read_csv(path+'/'+file)
        #temp = temp.drop(['Unnamed: 0.1','Unnamed: 0'])
        columns = temp.columns.tolist()
        temp = temp.rename(columns={'Field 1':'Media','Field 2':
                            'Strain','Field 3':'Specie'})
        if('Field 4' in temp.columns):
            temp = temp.drop(['Field 4'],axis=1)

        if(temp['Strain'].unique().tolist()[0] in strain_names_to_change.keys()):
            if(temp['Strain'].unique().tolist()[0]!='e.coli'):
                temp['Strain'] = strain_names_to_change[temp['Strain'].unique().tolist()[0]]
            if(temp['Strain'].unique().tolist()[0]!='e.coli'):
                temp['Strain'] = temp['Media']
                plate = temp['Plate type'].unique().tolist()
                if('PM01' in plate or 'PM02' in plate):
                    temp['Media'] = 'M9 minimal media no C source'
                elif('PM03' in plate or 'PM04' in plate):
                    temp['Media'] = 'IF0-Na-Succinate-Fe-Citrate'

        temp['Strain'] = temp['Strain'].apply(clean_strain)   
        temp['Media'] = temp['Media'].apply(clean_media)

        temp.to_csv(path+'/'+file)


def clean_strain(strain):
    # Remove special characters and spaces
    cleaned_strain = re.sub(r'[^A-Za-z0-9]+', '', strain)
    # Remove specified strings
    cleaned_strain = re.sub(r'R[1-9]', '', cleaned_strain)
    return cleaned_strain

def clean_media(media):
    if('m9' in media.lower()):
        clean_media = 'M9 minimal media no C source'
    elif('if0' in media.lower()):
        clean_media = 'IF0-Na-Succinate-Fe-Citrate'
    elif('mhb' in media.lower()):
        clean_media = 'MHB-CA'
    elif('rpmi' in media.lower()):
        clean_media = 'RPMI'
    return clean_media

def clean_plate(plate):
    temp_plate = re.sub(r'[^A-Za-z0-9]', '', plate)
    temp_plate = temp_plate.lower()
    if('pm01' in temp_plate or temp_plate=='pm1'):
        clean_plate = 'PM01'
    elif('pm02' in temp_plate or temp_plate=='pm2'):
        clean_plate = 'PM02A'
    elif('pm03' in temp_plate or temp_plate=='pm3'):
        clean_plate = 'PM03B'
    elif('pm04' in temp_plate or temp_plate=='pm4'):
        clean_plate = 'PM04A'
    elif('pm05' in temp_plate or temp_plate=='pm5'):
        clean_plate = 'PM05'
    elif('pm05' in temp_plate or temp_plate=='pm5'):
        clean_plate = 'PM05'
    elif('pm06' in temp_plate or temp_plate=='pm6'):
        clean_plate = 'PM06'
    elif('pm07' in temp_plate or temp_plate=='pm7'):
        clean_plate = 'PM07'
    elif('pm08' in temp_plate or temp_plate=='pm8'):
        clean_plate = 'PM08'
    elif('pm09' in temp_plate or temp_plate=='pm9'):
        clean_plate = 'PM09'
    elif('pm10' in temp_plate or temp_plate=='pm10'):
        clean_plate = 'PM10'
    elif('pm11' in temp_plate or temp_plate=='pm11c'):
        clean_plate = 'PM11C'
    elif('pm12' in temp_plate or temp_plate=='pm12b'):
        clean_plate = 'PM12B'
    return clean_plate
    

def fix_plate_entries(plate_dataframe):

    plate_dataframe['Strain'] = plate_dataframe['Strain'].apply(clean_strain)
    #plate_dataframe['Media'] = plate_dataframe['Media'].apply(clean_media)

    return plate_dataframe

def get_plate_files(onlyfiles,plate,path):
    """
    Collect all file names for a particular plate
    """
    plate_files = []
    for file in onlyfiles:
        temp = pd.read_csv(path+'/'+file)
        temp_plate = temp['Plate Type'].unique().tolist().pop()
        if(temp_plate in plate):
            plate_files.append(file)
    return plate_files


def get_strain_media(column_name):
    """
    Gets the strain and media type for a column name in the plate datatype
    """
    strain = ''
    media = ''
    flag = 0
    for char in column_name:
        if(char=='_'):
            flag = flag+1
        if(flag>1):
            break
        if(flag<1):
            strain = strain+char
        if(flag>0):
            media = media+char
    
    media = media[1:]
    
    return strain,media



def get_plates_from_study(path):
    """
    get list of all plates in a study/sample
    """

    onlyfiles = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    plate = []

    for file in onlyfiles:
        temp = pd.read_csv(path+'/'+file)

        for column in temp.columns:
            if(column.lower()=='plate type'):
                temp = temp.rename(columns={column:'plate'})
        
        if(not len(temp['plate'].unique().tolist())):
           print(file)
        plate.append(temp['plate'].unique().tolist().pop())

    plate = list(set(plate))
    
    return plate


def label_replicates(group):
    group['Replicates'] = 'R' + (group.groupby(['Strain', 'Well', 'Media']).cumcount() + 1).astype(str)
    return group

def get_plate_layout(plate):
    plate_layout = pd.read_table('Plate_layouts/'+plate+'.tsv',index_col='Well')
    return plate_layout

def make_plate_datatype(file_list,path,plate):
    # Create an empty DataFrame
    
    time_interval = np.linspace(0,48,193)
    combined_df = pd.DataFrame(columns=['Strain','Plate','Well', 'Media','Compound','KEGG ID','CAS ID'] + [f'{t}hrs' for t in time_interval])

    for file_name in file_list:
        # Read the data from the file
        data = pd.read_csv(path+'/'+file_name)
        strain = str(data['Strain'].unique().tolist().pop())
        media = str(data['Media'].unique().tolist().pop())


        #plate_layout = pd.read_table('Plate_layouts/'+plate+'.tsv',index_col='Well')
        plate_layout = get_plate_layout(plate)
        for well in wells:
            # Create a dictionary to store the data for the current file
            file_data = {
                'Strain': strain,
                'Well': well,
                'Media': media,
                'Plate': plate,
                'Compound': plate_layout.loc[well]['Compound'],
                'KEGG ID':plate_layout.loc[well]['KEGG ID'],
                'CAS ID':plate_layout.loc[well]['CAS ID']
            }

            # Iterate through the columns, which are the time points
            for col in range(0,len(time_interval)):
                # Extract the time value (e.g., '1hr', '1.5hr') and store it in the dictionary
                file_data[f'{time_interval[col]}hrs'] = data[well].values[col]

            # Append the data for the current file to the combined DataFrame
            combined_df = pd.concat([combined_df,pd.DataFrame([file_data])],ignore_index=True)

   
    combined_df['Replicates'] = ''
    # Sort the DataFrame by 'Strain', 'Well', and 'Media' to ensure consistent labeling
    combined_df = combined_df.sort_values(['Strain', 'Well','Media'])

    # Group the DataFrame by 'Strain', 'Well', and 'Media' and apply the labeling function
    combined_df = combined_df.groupby(['Strain', 'Well','Media'], group_keys=False).apply(label_replicates)

    # Set the columns you want at the beginning
    column_order = ['Strain','Plate','Well', 'Media', 'Replicates','Compound','KEGG ID','CAS ID']

    # Extract the remaining columns
    other_columns = [col for col in combined_df.columns if col not in column_order]

    # Combine the columns in the desired order
    new_column_order = column_order + other_columns

    # Reorder the DataFrame
    combined_df = combined_df[new_column_order]
    combined_df = combined_df.reset_index(drop=True)

    return combined_df



def remove_no_replicates(plate_dataframe):
    # Group by the unique combinations of Strain, Media, and Well
    groups = plate_dataframe.groupby(['Strain', 'Media', 'Well'])

    # Initialize an empty list to store the Strain names to drop
    strains_to_drop = []

    # Iterate through the groups
    for group_name, group in groups:
        if len(group) == 1:
            strains_to_drop.append(group_name[0])
            plate_dataframe = plate_dataframe.drop(group.index)

    return plate_dataframe, strains_to_drop



## Remove wells with low replicate correlations (<90%)

"""
This function is not used in the current version, but can be invoked to reject samples with low
replicate correlation
"""
    
def check_correlation(dataframe):
    """
    Colect all replicates for a strain + well and select ones with best correlation
    else reject the strain+well
    """
    from scipy.signal import savgol_filter
    
    filtered_dataframe = pd.DataFrame(index=dataframe.index)
    indices_to_drop = []
    for column in dataframe.columns:
        if('A01' in column):
            continue
        indices_to_keep = []
        df = dataframe[column]
        df2 = df.copy(deep=True)
        ## If no replicates are available, drop the column
        if(df2.size==df2.shape[0]):
            indices_to_drop.append(column)
            continue
        else:
            ## Apply a savgol filter to smoothen signal to remove noise for correlation
            for col in range(0,df2.shape[1]):
                df2.iloc[:,col] = savgol_filter(df2.iloc[:,col], 50, 3)   
            corr = df2.corr()
            corr = corr[corr>0.7]
            corr = corr[corr<1]
            for i in range(0,corr.shape[0]):
                for j in range(0,corr.shape[1]):
                    if(corr.iloc[i,j]==np.max(corr)[1] and ~pd.isna(np.max(corr)[-1])):
                        indices_to_keep.append(i)
                        indices_to_keep.append(j)
                        indices_to_keep = list(set(indices_to_keep))
                        signal = (df.iloc[:,indices_to_keep[0]] + df.iloc[:,indices_to_keep[1]])/2
                        filtered_dataframe[column] = signal
                    elif(pd.isna(np.max(corr)[-1])):
                        indices_to_drop.append(column)
        indices_to_keep = []
    indices_to_drop = list(set(indices_to_drop))
    
    return filtered_dataframe,indices_to_drop



"""
Note : These functions are deprecated and are only here for book keeping
These functions generate the final, human-readable csv sheets from the processed results
"""

def curate_kinetic_frame(kinetic_frame,plate,specie,plate_layout):
    
    """
    Generate the kinetic summary csv file
    """
    
    columns = ['Plate','Well','Strain','Specie','Compound','Max Resp',
              'Max Resp Rate','Time till max resp rate','AUC','Z-score',"Growth(1)/No Growth(0)",
               'KEGG ID','CAS ID']
    
    kinetic_summary = pd.DataFrame(columns=columns)
    
    for col in range(0,kinetic_frame.shape[1]):
        
        kinetic_summary.loc[col,'Plate'] = plate
        kinetic_summary.loc[col,'Specie'] = specie
        kinetic_summary.loc[col,'Strain'] = kinetic_frame.columns[col][:-4]
        kinetic_summary.loc[col,'Well'] = kinetic_frame.columns[col][-3:]
        kinetic_summary.loc[col,'Compound'] = plate_layout.loc[kinetic_frame.columns[col][-3:]]['Compound']
        kinetic_summary.loc[col,'Max Resp'] = kinetic_frame.iloc[0,col]
        kinetic_summary.loc[col,'Max Resp Rate'] = kinetic_frame.iloc[1,col]
        kinetic_summary.loc[col,'Time till max resp rate'] = kinetic_frame.iloc[2,col]
        kinetic_summary.loc[col,'AUC'] = kinetic_frame.iloc[3,col]
        kinetic_summary.loc[col,'Z-score'] = kinetic_frame.iloc[4,col]
        kinetic_summary.loc[col,'Growth(1)/No Growth(0)'] = kinetic_frame.iloc[5,col]
        kinetic_summary.loc[col,'KEGG ID'] = plate_layout.loc[kinetic_frame.columns[col][-3:]]['KEGG ID']
        kinetic_summary.loc[col,'CAS ID'] = plate_layout.loc[kinetic_frame.columns[col][-3:]]['CAS ID']
    
    return kinetic_summary


def curate_plate_frame(plate_frame,plate,specie,plate_layout):
    
    """
    Generate the plate summary csv file
    """
    
    columns = ['Plate','Strain','Specie','Well','Compound','KEGG ID','CAS ID']
    
    i = 0
    while(i<=48):
        columns.append('Time_'+str(i)+'_hrs')
        i = i+0.25
    
    plate_summary = pd.DataFrame(columns = columns)
    smooth_plate_summary = pd.DataFrame(columns = columns)
    
    for col in range(0,plate_frame.shape[1]):
        plate_summary.loc[col,'Plate'] = plate
        smooth_plate_summary.loc[col,'Plate'] = plate
        plate_summary.loc[col,'Strain'] = plate_frame.columns[col][:-4]
        smooth_plate_summary.loc[col,'Strain'] = plate_frame.columns[col][:-4]
        plate_summary.loc[col,'Well'] = plate_frame.columns[col][-3:]
        smooth_plate_summary.loc[col,'Well'] = plate_frame.columns[col][-3:]
        plate_summary.loc[col,'Specie'] = specie
        smooth_plate_summary.loc[col,'Specie'] = specie
        plate_summary.loc[col,'Compound'] = plate_layout.loc[plate_frame.columns[col][-3:]]['Compound']
        smooth_plate_summary.loc[col,'Compound'] = plate_layout.loc[plate_frame.columns[col][-3:]]['Compound']
        plate_summary.loc[col,'KEGG ID'] = plate_layout.loc[plate_frame.columns[col][-3:]]['KEGG ID']
        smooth_plate_summary.loc[col,'KEGG ID'] = plate_layout.loc[plate_frame.columns[col][-3:]]['KEGG ID']
        plate_summary.loc[col,'CAS ID'] = plate_layout.loc[plate_frame.columns[col][-3:]]['CAS ID']
        smooth_plate_summary.loc[col,'CAS ID'] = plate_layout.loc[plate_frame.columns[col][-3:]]['CAS ID']
        plate_summary.iloc[col,7:] = plate_frame.iloc[:,col]
        smooth_plate_summary.iloc[col,7:] = plate_frame.iloc[:,col]
    
    return plate_summary,smooth_plate_summary


def curate_growth_summary(growth_matrix,plate,specie,plate_layout):
    
    """
    Generate the growth summary csv file
    """
    
    strains = []
    plates = []
    species = []
    wells = growth_matrix.index.tolist()
    compounds = plate_layout.loc[wells]['Compound']
    kegg_ids = plate_layout.loc[wells]['KEGG ID'].unique()
    cas_ids = plate_layout.loc[wells]['CAS ID']
    
    columns = ['Plate','Strain','Specie','Well','Compound','KEGG ID','CAS ID','Growth(1)/No Growth(0)/NA(0.5)']
    growth_summary = pd.DataFrame(columns=columns)
    
    j = 0
    for strain in growth_matrix.columns.tolist():
        k = 0
        for i in range(j,j+len(wells)):
            growth_summary.loc[i,'Plate'] = plate
            growth_summary.loc[i,'Strain'] = strain
            growth_summary.loc[i,'Specie'] = specie
            growth_summary.loc[i,'Well'] = wells[k]
            growth_summary.loc[i,'Growth(1)/No Growth(0)/NA(0.5)']= growth_matrix.loc[wells[k],strain]
            growth_summary.loc[i,'Compound'] = plate_layout.loc[wells[k]]['Compound']
            growth_summary.loc[i,'KEGG ID'] = plate_layout.loc[wells[k]]['KEGG ID']
            growth_summary.loc[i,'CAS ID'] = plate_layout.loc[wells[k]]['CAS ID']
            k = k+1
        j = j+96
    
    return growth_summary



def create_summary_table(plate,specie,strains):
    
    """
    Generate the summary sheet - all high quality growth call strains and plates
    """
    summary_table = pd.DataFrame(columns=['Plate','Strain','Specie'])
    for strain in range(0,len(strains)):
        summary_table.loc[strain,'Plate']=plate
        summary_table.loc[strain,'Strain']=strains[strain]
        summary_table.loc[strain,'Specie']=specie
    return summary_table


def remove_bad_controls(df,bad_controls_df):

    result = pd.merge(df, bad_controls_df, on=['Plate', 'Strain'], how='left', indicator=True)

    # Filter the rows where the indicator column is 'left_only' (only in df1)
    df1_filtered = result[result['_merge'] == 'left_only']

    # Drop the indicator column
    df1_filtered = df1_filtered.drop(columns=['_merge'])
    
    bad_control_frame = df.merge(bad_controls_df, on=['Plate', 'Strain'], how='inner')
    return df1_filtered, bad_control_frame


def filter_bad_controls(bad_kinetics,threshold=200):
    # Initialize an empty dictionary to store the results
    result_dict = {}
    
    # Filter rows where 'Max Resp' > 200
    filtered_df = bad_kinetics[(bad_kinetics['Well']=='A01') & (bad_kinetics['Max Resp'] > 200)]
    
    # Iterate through the filtered DataFrame and populate the dictionary
    for index, row in filtered_df.iterrows():
        plate = row['Plate']
        strain = row['Strain']
        
        if plate in result_dict:
            if strain not in result_dict[plate]:
                result_dict[plate].append(strain)
        else:
            result_dict[plate] = [strain]
    
    return result_dict
