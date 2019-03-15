# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 14:50:13 2019

@author: Dariia Vyshenska
"""

import pandas as pd
from scipy.stats import ttest_ind
#import numpy as np
#%%
# FUNCTIONS
# gets either mean or median of the grouped table (specified by user)
def dispatch_if(operator, x):
    if operator == 'mean':
        return x.mean()
    elif operator == 'median':
        return x.median()
    else:
        return None
    
# converting table into a dict of dict
def table_to_dic(test_table):
    clean_list = lambda list: [x for x in list if str(x) != 'nan']
    test_table_dic = dict()
    for column in test_table:    
        sub_dic = dict()
        for g in groups:
            sub_dic[g] = clean_list(list(test_table[column][test_table.index==g]))
        test_table_dic[column] = sub_dic.copy()
        del sub_dic
    return test_table_dic

# getting p-values as separate table
def pval_calc(test_table_dic, groups):
    pval_table = pd.DataFrame(index=groups)
    for key in test_table_dic:
        p = list()
        for sub_key in groups:
            t_res = ttest_ind(test_table_dic[key][sub_key], \
                      test_table_dic[key][contr_g],\
                      equal_var = False)
            p.append(t_res[1])
        pval_table[key] = p
        del p
    return pval_table


#%%
# COMMAND LINE PARAMETERS
contr_g = 'NCD'
operator = 'median'

#%%
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

#%%

# next - for each expt, only for treatment index do mean or median test 
# to get the parameter for each column
c_all = full_table.loc[pd.IndexSlice[:, contr_g],:]
m_normz = dispatch_if(operator, c_all.groupby(level='expt', axis=0))
nfull_table = full_table.subtract(m_normz, axis=1,level='expt')

# export normalized table, full (can be customized with sorting by indexes)
nfull_table.T.to_csv (r'normalized_full_table.csv', index=True, header=True)

# next - do statistical comparison of pooled groups (each) to the control 
# pooled group
test_table = nfull_table.droplevel([0,2])
pvalues = pval_calc(table_to_dic(test_table), set(test_table.index))

# merging parameters table into second output table
table_medians = dispatch_if(operator,nfull_table.groupby(level='group',axis=0))


# next - export pvalues, fdrs, means and medians for each feature tested etc - 
# see the page with specs of second output file