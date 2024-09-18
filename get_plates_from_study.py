#!/usr/bin/env python
# coding: utf-8

# In[19]:


import os
import pandas as pd

def get_plates_from_study():
    path = 'test'

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


# In[ ]:




