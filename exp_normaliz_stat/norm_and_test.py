# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 14:50:13 2019

@author: Dariia Vyshenska
"""

import pandas as pd

contr_g = 'NCD'

# NEEDS REWRITING: take file names from args and loop through them,
# also appending them to each other (skip append the first time loop runs)
df1 = pd.read_csv('mic1-pheno-w-cholest.csv', index_col = 0).T
map1 = pd.read_csv('map.csv', index_col = 0)
map1['expt'] = 1
df1_m = df1.join(map1, how='outer')
df1_m.set_index(['expt', 'Treatment'],append = True, inplace=True)
df1_m.index.names = ['sample', 'expt', 'group']
df1_m = df1_m.reorder_levels([1, 2, 0])


# second file. need to be transformed into a loop
df2 = pd.read_csv('mic3-pheno-w-cholest.csv', index_col = 0).T
map2 = pd.read_csv('map2.csv', index_col = 0)
map2['expt'] = 2
df2_m = df2.join(map2, how='outer')
df2_m.set_index(['expt', 'Treatment'],append = True, inplace=True)
df2_m.index.names = ['sample', 'expt', 'group']
df2_m = df2_m.reorder_levels([1, 2, 0])

# then merge this two experiments together
full_table = df1_m.append(df2_m)



# next - for each expt, only for treatment index do mean or median test 
# to get the parameter for each column
c_all = full_table.iloc[full_table.index.get_level_values('group')==contr_g]
m_normz = c_all.groupby(level='expt', axis=0).mean()
full_table.divide(m_normz, axis=1,level='expt') #check if normalization happens!


# apply normalization by the parameter to the respective experiment
# each experiment each column will have it's own parameter value
# but normalization by the parameter will be independent of the group