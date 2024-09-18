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
from scipy.stats import zscore
from multiprocessing import Pool
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from utils import get_plate_layout
from scipy.stats import ttest_ind


def perform_z_test(group):
    z_scores = zscore(group['Max Resp'], ddof=1)
    p_values = 2 * (1 - norm.cdf(abs(z_scores)))
    return p_values

def get_kinetic_parameters(signal,time):
    """
    Extract the max resp, rate, time and auc for a respiration signal
    """
    smoothened_signal = savgol_filter(signal, 50, 3)
    dt = []
    for i in range(0,np.shape(smoothened_signal)[0]-1):
        dt.append((smoothened_signal[i+1]-smoothened_signal[i])/(time[i+1]-time[i]))
    max_resp_rate = np.max(dt)
    time_till = time[np.argmax(smoothened_signal)]
    max_val = np.max(smoothened_signal)
    auc = np.trapz(smoothened_signal,time,dx=0.25)   
    return max_val,max_resp_rate,time_till,auc

def get_kinetic_dataframe(plate_dataframe):
    
    kinetic_dataframe = pd.DataFrame(index=plate_dataframe.index,columns=['Strain','Well','Plate','Media','Replicates','Compound','KEGG ID','CAS ID','Max Resp','Max Resp Rate','Time till max resp rate','AUC'])

    for i in range(0,plate_dataframe.shape[0]):
        kinetic_dataframe.iloc[i,0] = plate_dataframe['Strain'][i]
        kinetic_dataframe.iloc[i,1] = plate_dataframe['Well'][i]
        kinetic_dataframe.iloc[i,2] = plate_dataframe['Plate'][i]
        kinetic_dataframe.iloc[i,3] = plate_dataframe['Media'][i]
        kinetic_dataframe.iloc[i,4] = plate_dataframe['Replicates'][i]
        kinetic_dataframe.iloc[i,5] = plate_dataframe['Compound'][i]
        kinetic_dataframe.iloc[i,6] = plate_dataframe['KEGG ID'][i]
        kinetic_dataframe.iloc[i,7] = plate_dataframe['CAS ID'][i]
        max_val,max_resp_rate,time_till,auc = get_kinetic_parameters(plate_dataframe.iloc[i,8:],np.linspace(0,48,193))
        kinetic_dataframe.iloc[i,8] = max_val
        kinetic_dataframe.iloc[i,9] = max_resp_rate
        kinetic_dataframe.iloc[i,10] = time_till
        kinetic_dataframe.iloc[i,11] = auc

    return kinetic_dataframe

def make_growth_calls(kinetic_dataframe,max_resp_threshold=120,alpha=0.05,negative_control=True):

    if(negative_control):
        # Define a significance level

        # Create an empty 'Growth' column and initialize with zeros
        kinetic_dataframe['Growth'] = 0

        # Filter for control wells (Well is 'A01')
        control_group = kinetic_dataframe[(kinetic_dataframe['Well'] == 'A01')&(kinetic_dataframe['Max Resp']<max_resp_threshold)]

        # Create an empty list to store p-values
        p_values = []

        # Iterate through each row and perform a one-sided z-test
        for index, row in kinetic_dataframe.iterrows():
            # Perform a one-sided z-test using control group statistics
            z_score = (row['Max Resp'] - control_group['Max Resp'].mean()) / (control_group['Max Resp'].std() / (len(control_group) ** 0.5))
            p_value = 1 - norm.cdf(z_score)
            p_values.append(p_value)

        # Apply Benjamini-Hochberg correction to p-values
        p_adjusted = multipletests(p_values, method='fdr_bh')[1]

        for i, p_adj in enumerate(p_adjusted):
            if (p_adj < alpha and kinetic_dataframe['Max Resp'].iloc[i] > control_group['Max Resp'].mean() and kinetic_dataframe['Max Resp'].iloc[i]>max_resp_threshold):
                kinetic_dataframe.at[i, 'Growth'] = 1


    else:
        # Create an empty 'Growth' column and initialize with zeros
        kinetic_dataframe['Growth'] = 0

        # Filter for control wells (Well is 'A01')
        control_group = kinetic_dataframe[kinetic_dataframe['Max Resp']<max_resp_threshold]

        # Create an empty list to store p-values
        p_values = []

        # Iterate through each row and perform a one-sided z-test
        for index, row in kinetic_dataframe.iterrows():
            # Perform a one-sided z-test using control group statistics
            z_score = (row['Max Resp'] - control_group['Max Resp'].mean()) / (control_group['Max Resp'].std() / (len(control_group) ** 0.5))
            p_value = 1 - norm.cdf(z_score)
            p_values.append(p_value)

        # Apply Benjamini-Hochberg correction to p-values
        p_adjusted = multipletests(p_values, method='fdr_bh')[1]

        for i, p_adj in enumerate(p_adjusted):
            if (p_adj < alpha and kinetic_dataframe['Max Resp'].iloc[i] > control_group['Max Resp'].mean() and kinetic_dataframe['Max Resp'].iloc[i]>max_resp_threshold):
                kinetic_dataframe.at[i, 'Growth'] = 1

        # Display the resulting DataFrame
    return kinetic_dataframe

def check_significance(group,control_group):
    control_group_max_resp = pd.to_numeric(control_group['Max Resp'], errors='coerce')
    group_max_resp = pd.to_numeric(group['Max Resp'], errors='coerce')
    
    # Perform Welch's t-test with a greater alternative hypothesis
    t_stat, p_value = ttest_ind(group_max_resp, control_group_max_resp, equal_var=False, alternative='greater')
    
    growth = 1 if p_value < 0.05 else 0
    
    return growth

def make_growth_calls_strong_controls(kinetic_dataframe,alpha=0.05):

    control_group = kinetic_dataframe[kinetic_dataframe['Well'] == 'A01']
    non_controls = kinetic_dataframe[kinetic_dataframe['Well'] != 'A01']

    # Create a new DataFrame to store the results
    result_df = pd.DataFrame(columns=['Strain', 'Media', 'Plate', 'Well', 'Growth'])

    # Iterate over unique combinations of Strain, Media, and Plate
    for key, group in non_controls.groupby(['Strain', 'Media', 'Plate', 'Well']):
        growth = check_significance(group,control_group)
        result_df = result_df.append({'Strain': group['Strain'].tolist()[0], 'Media': group['Media'].tolist()[0], 'Plate': group['Plate'].tolist()[0], 'Well': group['Well'].tolist()[0], 'Growth': growth}, ignore_index=True)

    # Check for cases with growth = 0 and average max_resp > 200, and reassign them
    for key, group in result_df[result_df['Growth'] == 0].groupby(['Strain', 'Media', 'Plate', 'Well']):
        original_group = kinetic_dataframe[(kinetic_dataframe['Strain'] == group['Strain'].tolist()[0]) & (kinetic_dataframe['Media'] == group['Media'].tolist()[0]) & (kinetic_dataframe['Plate'] == group['Plate'].tolist()[0]) & (kinetic_dataframe['Well'] == group['Well'].tolist()[0])]
        
        if original_group['Max Resp'].mean() > 200:
            result_df.loc[group.index, 'Growth'] = 1


    # Create a new column "Growth" in the 'kinetics' DataFrame
    kinetic_dataframe['Growth'] = 0  # Initialize all values to 0

    # Iterate over unique combinations of Strain, Media, Plate, and Well
    for key, group in result_df.groupby(['Strain', 'Media', 'Plate', 'Well']):
        # Find the corresponding rows in the 'kinetics' DataFrame
        condition = (
            (kinetic_dataframe['Strain'] == group['Strain'].tolist()[0]) &
            (kinetic_dataframe['Media'] == group['Media'].tolist()[0]) &
            (kinetic_dataframe['Plate'] == group['Plate'].tolist()[0]) &
            (kinetic_dataframe['Well'] == group['Well'].tolist()[0])
        )

        # Update the "Growth" column in the 'kinetics' DataFrame
        kinetic_dataframe.loc[condition, 'Growth'] = group['Growth'].values[0]

    return result_df,kinetic_dataframe

def majority_growth_func(group):
    total_entries = len(group)
    growth_sum = group['Growth'].sum()

    if growth_sum / total_entries >= 0.5:
        group['Growth'] = 1
    elif (total_entries - growth_sum) / total_entries >= 0.5:
        group['Growth'] = 0
    else:
        group['Growth'] = 0.5  # In cases where there's an equal number of 1s and 0s

    return group

def make_growth_datatype(kinetic_dataframe,plate,negative_control=True):

    # Group by columns 'Strain', 'Well', 'Plate', and 'Media'
    grouped = kinetic_dataframe.groupby(['Strain', 'Well', 'Plate', 'Media'])

    # Count the number of 1s and 0s in each group
    count_1s = grouped['Growth'].sum()
    count_0s = grouped['Growth'].count() - count_1s

    # Calculate majority growth with a tie case
    majority_growth = (count_1s > count_0s).astype(float)
    majority_growth[count_1s == count_0s] = 0.5

    # Reset the index
    majority_growth = majority_growth.reset_index()
    #majority_growth.rename(columns={'Growth': 'Majority_Growth'}, inplace=True)

    bad_controls = []

    if(negative_control):
        controls = majority_growth.loc[majority_growth['Well']=='A01']
        bad_controls_df = controls.loc[controls['Growth']==1]
        bad_controls = bad_controls_df['Strain'].tolist()

    plate_layout = get_plate_layout(plate)
    majority_growth = majority_growth.merge(plate_layout, on='Well', how='left')
    return majority_growth,bad_controls

def logistic(x, A, u, d, v):
    """
    Fits a growth curve parameterized logistic function
    """
    y = (A / (1 + np.exp( ( ((4 * u)/A) * (d - x) ) + 2 ))) 
    return y

def gompertz(x, A, u, d, v):
    """
    Fits a growth curve parameterized gompertz function
    """
    y = (A * np.exp( -np.exp( (((u * np.e)/A) * (d - x)) + 1 ) ) )
    return y

def richards(x, A, u, d, v):
    """
    Fits a growth curve parameterized richards function
    """
    y = (A * pow(1 + (v + (np.exp(1 + v) * np.exp( (u/A) * (1 + v) * (1 + (1/v)) * (d - x) ) ) ),-(1/v)))
    return y

def smooth_plate(dataframe):
    smooth_frame = pd.DataFrame(index=dataframe.index,columns=dataframe.columns)
    for col in range(0,dataframe.shape[1]):
        smooth_frame.iloc[:,col] = savgol_filter(dataframe.iloc[:,col], 50, 3)
    return smooth_frame

def fit_curve(data,time):
    """
    Fits the best function out of the three above to each well
    in the plate object, if no function fits all the parameters are assigned nan.unique()
    """
    from scipy.optimize import curve_fit
    from scipy.signal import savgol_filter
    import math
    rms = []
    parameters = []
    fit_data = np.zeros([np.size(data),3])
    
    data = savgol_filter(data, 50, 3)#Smoothen data to remove noise
        
    try:
        params, pop = curve_fit(logistic,time,data,maxfev=10000,bounds=(0,[1.1*np.max(data),10*np.max(data),np.max(time),1]))
        parameters.append(params)
        fit_data[:,0] = logistic(time,*params)
        rms.append(np.sqrt(np.average((fit_data[:,0] - data)**2)))
    except:
        rms.append(math.inf)
        parameters.append(['na','na','na','na'])
            
    try:
        params, pop = curve_fit(gompertz,time,data,maxfev=10000,bounds=(0,[1.1*np.max(data),10*np.max(data),np.max(time),1]))
        fit_data[:,1] = logistic(time,*params)
        rms.append(np.sqrt(np.average((fit_data[:,1] - data)**2)))
        parameters.append(params)
    except:
        rms.append(math.inf)
        parameters.append(['na','na','na','na'])
    try:
        params, pop = curve_fit(richards,time,data,maxfev=10000,bounds=(0,[1.1*np.max(data),10*np.max(data),np.max(time),1]))
        fit_data[:,2] = logistic(time,*params)
        rms.append(np.sqrt(np.average((fit_data[:,2] - data)**2)))
        parameters.append(params)
    except:
        rms.append(math.inf)
        parameters.append(['na','na','na','na'])
            
    min_rms = np.argmin(rms)

    return fit_data[:,min_rms],parameters[min_rms]



