# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 14:50:13 2019

@author: Dariia Vyshenska
"""
# FUNCTIONS
# takes files from arguments and collapses them into one table per experiment
def collapse_files(data, mapping):
    data_f = pd.read_csv(data, index_col = 0).T
    data_f.index.names = ['sample']
    map_f = pd.read_csv(mapping, header=None)
    map_f.columns = ['sample', 'group']
    map_f.set_index('sample', inplace=True)
    map_f['expt'] = 1
    dm = data_f.join(map_f, how='outer')
    dm.set_index(['expt', 'group'],append = True, inplace=True)
    dm = dm.reorder_levels([1, 2, 0])
    return(dm)

# gets either mean or median of the grouped table (specified by user)
def dispatch_if(operator, x):
    if operator == 'mean':
        return x.mean()
    elif operator == 'median':
        return x.median()
    else:
        return None

def dispatch_if_T(operatorT, x, y):
    if operatorT == 'parametric':
        return ttest_ind(x, y, equal_var = False)[1]
    elif operatorT == 'non-parametric':
        return mannwhitneyu(x,y, use_continuity=True, \
                            alternative='two-sided')[1]
    else:
        return None

# converting table into a dict of dict
def table_to_dic(test_table):
    clean_list = lambda list: [x for x in list if str(x) != 'nan']
    test_table_dic = dict()
    for column in test_table:    
        sub_dic = dict()
        for g in set(test_table.index):
            sub_dic[g]=clean_list(list(test_table[column][test_table.index==g]))
        test_table_dic[column] = sub_dic.copy()
        del sub_dic
    return test_table_dic

# getting p-values as separate table
def pval_calc(test_table_dic, groups):
    pval_table = pd.DataFrame(index=groups)
    for key in test_table_dic:
        p = list()
        for sub_key in groups:
            t_res= dispatch_if_T(operatorT, test_table_dic[key][sub_key], \
                                 test_table_dic[key][contr_g])
            p.append(t_res)
        pval_table[key] = p
        del p
    return pval_table

# adding parameter index and sorting the dataframe based on first index(group)
def r_name(table, parameter):
    new_table = (table.set_index([table.index, [parameter]* len(table)])
    .rename_axis(['group', 'parameter'])
    .sort_index(axis=0, level=0))
    return new_table

#%%

if __name__ == "__main__":
    
    try:
        import pandas as pd
        from scipy.stats import ttest_ind
        from scipy.stats import mannwhitneyu
        import argparse
        import sys

    except ModuleNotFoundError:
        print("This script requires packages: sys, pandas, argparse and \
        scipy.stats!\n")
    
    parser = argparse.ArgumentParser()

    parser.add_argument("--d_file", "-df", type=str, required=True, nargs='+',\
                        help="Files with data: csv file with header. "+ \
                        "Header - unique sample identifiers. " + \
                        "First column - name of the assay type.")
    parser.add_argument("--m_file", "-mf", type=str, required=True, nargs='+',\
                        help="Mapping files, same order as for data files. "+ \
                        "csv files, no header. " + \
                        "First column - unique sample identifiers. "+ \
                        "Second column - name of treatment group.")
    parser.add_argument("--control", "-C", type=str, required=True, nargs=1, \
                        help="Name of control treatment group.")
    parser.add_argument("--test_parameter", "-tP", type=str, required=True, \
                        nargs=1, help="Normalization parameter: " + \
                        "'mean' or 'median'.")
    
    parser.add_argument("--st_test", "-stT", type=str, required=True, \
                        nargs=1, help="T-test type: 'parametric' for Welch "+\
                        " t-test, or 'non-parametric' for Mann-Whitney test.")
    
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()


    if len(args.d_file) == len(args.m_file):
        contr_g = args.control[0]
        operator = args.test_parameter[0]
        operatorT = args.st_test[0]
        
        # combining mapping and data files into one table
        for count, d_file in enumerate(args.d_file):
            if count == 0:
                full_table = collapse_files(d_file, args.m_file[count])
            elif count > 0:
               t_table =  collapse_files(d_file, args.m_file[count])
               full_table = full_table.append(t_table)
        # next - for each expt, only for treatment index do mean or median test 
        # to get the parameter for each column
        c_all = full_table.loc[pd.IndexSlice[:, contr_g],:]
        m_normz = dispatch_if(operator, c_all.groupby(level='expt', axis=0))
        nfull_table = full_table.subtract(m_normz, axis=1,level='expt')
        
        # next - do statistical comparison of pooled groups (each) to the control 
        # pooled group
        test_table = nfull_table.droplevel([0,2])
        pvalues = pval_calc(table_to_dic(test_table), set(test_table.index))
        pvalues = r_name(table=pvalues, parameter='pvalue')
        
        # merging parameters table into second output table
        table_m = dispatch_if(operator,nfull_table.groupby(level='group', \
                                                           axis=0))
        table_m = r_name(table=table_m, parameter=operator)
        
        # export normalized table and table with calculated parameters & p-values
        nfull_table.T.to_csv (r'normalized_full_table.csv', index=True, \
                              header=True)
        table_m.append(pvalues).T.to_csv(r'parameters_table.csv',\
                                   index=True, header=True)

    else:
        print("Data and mapping files must be provided in pairs in " + \
              "respective order")
        quit()