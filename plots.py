import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
import os

wells =     ['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12',
             'B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12',
             'C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12',
             'D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12',
             'E01','E02','E03','E04','E05','E06','E07','E08','E09','E10','E11','E12',
             'F01','F02','F03','F04','F05','F06','F07','F08','F09','F10','F11','F12',
             'G01','G02','G03','G04','G05','G06','G07','G08','G09','G10','G11','G12',
             'H01','H02','H03','H04','H05','H06','H07','H08','H09','H10','H11','H12']


def create_donut_plot(AUC,size,rad):
    """
    Generate a donut plot to compare kinetic parameters
    """
    import matplotlib as mpl
    
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.patches as mpatches
    from matplotlib.collections import PatchCollection

    colors = ['#000052','#0c44ac','#faf0ca','#ed0101','#970005'] 

    fig,ax = plt.subplots(figsize=(20,5))
    
    j = 0
    for strain in range(0,AUC.shape[0]):
        outer_colors = []
        cm = LinearSegmentedColormap.from_list('custom', colors,N=2*int(max(AUC.iloc[strain,:]-min(AUC.iloc[strain,:]))))
        cm.set_bad(color='white')
        for i in AUC.iloc[strain,:]:
            outer_colors.append(cm(int(i)))
        

        
        if(strain==0):
            ax.pie(np.ones([len(AUC.iloc[strain,:]),1])[:,0],radius=rad-strain*size,labels=AUC.columns,colors=outer_colors,
               wedgeprops=dict(width=size, edgecolor='w'),textprops={'fontsize': 5})
            ax.text(1.5, 0.5+j, AUC.index[strain])

        else:
            ax.pie(np.ones([len(AUC.iloc[strain,:]),1])[:,0],radius=rad-strain*size,colors=outer_colors,
               wedgeprops=dict(width=size, edgecolor='w'),textprops={'fontsize': 5})
            ax.text(1.5, 0.5+j, AUC.index[strain])

            
        j = j+0.2

            
            
def plate_charts(growth_matrix,kinetic_parameters,plate_layout,plate):
    """
    Make kinetic charts for each plate in a project
    """
    kinetic_matrix = pd.DataFrame(columns=growth_matrix.columns,index=growth_matrix.index)
    asymptote_matrix = pd.DataFrame(columns=growth_matrix.columns,index=growth_matrix.index)
    temp_growth_matrix = pd.DataFrame(columns=growth_matrix.columns,index=growth_matrix.index)
    for well in kinetic_matrix.index:
        for strain in kinetic_matrix.columns:
            if(growth_matrix.loc[well,strain]!=0.5):
                kinetic_matrix.loc[well,strain]=kinetic_parameters.loc['Max Resp Rate',strain+'_'+well][0]
                asymptote_matrix.loc[well,strain]=kinetic_parameters.loc['Max Resp',strain+'_'+well][0]
                temp_growth_matrix.loc[well,strain]=growth_matrix.loc[well,strain]
                
    indices = []
    for index in kinetic_matrix.index:
        indices.append(plate_layout.loc[index]['Compound']+'_'+index)
    kinetic_matrix.index = indices
    asymptote_matrix.index = indices
                
    high_var_strains = asymptote_matrix.var().sort_values(ascending=False).index[0:5]
    high_var_wells = asymptote_matrix.T.var().sort_values(ascending=False).index[0:10]
    low_var_wells = asymptote_matrix.T.var().sort_values(ascending=True).index[0:10]
    
    sns.set(rc={'figure.figsize':(20,10)})
    sns.set_style('white')
    sns.boxplot(data=kinetic_matrix.T)
    plt.xticks(rotation=90)
    plt.title('Max Respiration Rates')
    plt.ylabel('strains')
    plt.xlabel('wells')
    plt.savefig('../Stats_Data/'+plate+'/'+"resp_rates.pdf",dpi=200)
    plt.clf()

    plt.figure(figsize=(20,10))
    df = asymptote_matrix.loc[high_var_wells.append(low_var_wells),high_var_strains].T
    create_donut_plot(df.fillna(100),0.1,1)
    plt.title('Max Respiration Observed')
    plt.savefig('../Stats_Data/'+plate+'/'+"well_variance.pdf",dpi=200)



def make_max_resp_histogram(kinetic_frame,ax,negative_control=True):

    if(negative_control):
        control_wells = kinetic_frame.loc[kinetic_frame['Well']=='A01']
        other_wells = kinetic_frame.loc[kinetic_frame['Well']!='A01']
        ax.hist(control_wells['Max Resp'],color='thistle',alpha=0.5,density = True)
        ax.hist(other_wells['Max Resp'],color='darkviolet',alpha=0.5,density = True)
        ax.set_ylabel('Max Resp prob',fontweight='bold')
        ax.set_xlabel('Wells',fontweight='bold')
        ax.legend(['Controls','Others'])
    else:
        other_wells = kinetic_frame.loc[kinetic_frame['Well']!='A01']
        ax.hist(other_wells['Max Resp'],color='darkviolet',alpha=0.5,density = True)
        ax.set_ylabel('Max Resp prob',fontweight='bold')
        ax.set_xlabel('Wells',fontweight='bold')



# def make_plate_heatmap(growth_dataframe, ax, cbar=True):
#     #growth_df = growth_dataframe[['Compound', 'Growth']].set_index('Compound')

#     cmap = sns.color_palette(['red', 'lightgray', 'green'])
#     ax = sns.heatmap(growth_dataframe, ax=ax, cbar=cbar, cmap=cmap, xticklabels=True, yticklabels=True, annot=True, fmt='.1g')

def make_plate_heatmap(growth_dataframe, ax, cbar=True):
    # Define custom colors for 0, 0.5, and 1
    colors = ['white', 'red', 'violet']
    custom_cmap = ListedColormap(colors)

    # Map the values to the range [0, 1, 2]
    mapped_data = growth_dataframe.replace({0: 0, 0.5: 1, 1: 2})

    # Create a heatmap without annotations, using the custom colormap
    ax = sns.heatmap(mapped_data, ax=ax, cbar=cbar, cmap=custom_cmap, xticklabels=True, yticklabels=True, annot=False
                    ,linewidths=1, linecolor='black')
    # Set the aspect ratio to be equal to make cells square
    #ax.set_aspect(1.2)
    # Clear x and y axis labels
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='both', which='major', labelsize=40)
    ax.set_xticklabels(ax.get_xticklabels(), weight='bold',size=10,rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), weight='bold',size=10)


def plot_bad_controls(controls_df,controls_kinetics_df,ax):

    no_growth_wells = controls_kinetics_df.loc[(controls_kinetics_df['Growth']==0)&(controls_kinetics_df['Well']!='A01')]['Well'].tolist()
    growth_wells = controls_kinetics_df.loc[(controls_kinetics_df['Growth']==1)&(controls_kinetics_df['Well']!='A01')]['Well'].tolist()
    control_wells = ['A01']
    no_growth_signals = controls_df.loc[controls_df['Well'].isin(no_growth_wells)]
    growth_signals = controls_df.loc[controls_df['Well'].isin(growth_wells)]
    control_signals = controls_df.loc[controls_df['Well'].isin(control_wells)]
    time = np.linspace(0,48,193)
    signal_columns = [str(w)+'hrs' for w in time]

    for sig in range(0,no_growth_signals.shape[0]):
        ax.plot(time,no_growth_signals.iloc[sig,:][signal_columns],label='no g')

    for sig in range(0,growth_signals.shape[0]):
        ax.plot(time,growth_signals.iloc[sig,:][signal_columns],label='g')
    
    for sig in range(0,control_signals.shape[0]):
        ax.plot(time,control_signals.iloc[sig,:][signal_columns],label='control')



def plot_high_background_noise(bad_controls,plate_signals,path):

    time = np.linspace(0,48,193)
    for plate in bad_controls.index:
        signals = plate_signals.loc[(plate_signals['Plate']==bad_controls['Plate'][plate])&(plate_signals['Strain']==bad_controls['Strain'][plate])]
        control_signals = signals.loc[signals['Well']=='A01']
        fig,ax = plt.subplots(1,1)

        for well in wells[1:]:

            for i in range(0,control_signals.shape[0]):
                ax.plot(time,control_signals.iloc[i,8:],label='control',color='palevioletred')
            other_signals = signals.loc[signals['Well']==well]
            for i in range(0,other_signals.shape[0]):
                ax.plot(time,other_signals.iloc[i,8:],label=well,color='deepskyblue')
            ax.set_ylabel('time')
            ax.set_xlabel('signal')
            ax.set_title(well)
            ax.legend()
            
            if(not os.path.isdir(path)):
                os.mkdir(path)
             
            if(not os.path.isdir(path+'/'+bad_controls['Plate'][plate])):
                os.mkdir(path+'/'+bad_controls['Plate'][plate])
                        
            if(os.path.isdir(path+'/'+bad_controls['Plate'][plate]+'/'+bad_controls['Strain'][plate])):
                fig.savefig(path+'/'+bad_controls['Plate'][plate]+'/'+bad_controls['Strain'][plate]+'/'+well+'.png',dpi=500)
            else:
                os.mkdir(path+'/'+bad_controls['Plate'][plate]+'/'+bad_controls['Strain'][plate])
                fig.savefig(path+'/'+bad_controls['Plate'][plate]+'/'+bad_controls['Strain'][plate]+'/'+well+'.png',dpi=500)

            ax.clear()






