# -*- coding: utf-8 -*-
"""
@author: Dariia Vyshenska

This script is created to combine multiple readouts from different experiments,
normalize the data experiment-wise 
(each_measurement minus control_group_mean_or_median), pool normalized data 
treatment-wise, and run t-test on all treatment combinations with the control 
treatment.
Output: table with normalized data, table with mean/median per treatment 
combined with t-test p-values and  fdr values. Both tables are in .csv format.

Build with: python 3.6.4, R 3.5.1
Require: sys, pandas(0.24.0), argparse(1.1), rpy2, R (require library stats 
version 3.5.1)
Require: r_statfun.R (to be located in the same folder)

usage: norm_and_test.py [-h] --d_file D_FILE [D_FILE ...] --m_file M_FILE
                        [M_FILE ...] --control CONTROL --test_statistic
                        TEST_STATISTIC --st_test ST_TEST

optional arguments:
  -h, --help            show this help message and exit
  --d_file D_FILE [D_FILE ...], -df D_FILE [D_FILE ...]
                        Files with data: csv file with header. Header - unique
                        sample identifiers. First column - name of the assay
                        type.
  --m_file M_FILE [M_FILE ...], -mf M_FILE [M_FILE ...]
                        Mapping files, same order as for data files. csv
                        files, no header. First column - unique sample
                        identifiers. Second column - name of treatment group.
  --control CONTROL, -C CONTROL
                        Name of control treatment group.
  --test_statistic TEST_STATISTIC, -tP TEST_STATISTIC
                        Normalization statistic: 'mean' or 'median'.
  --st_test ST_TEST, -stT ST_TEST
                        T-test type: 'parametric' for Welch two-sided t-test
                        (R t.test, default args), or 'non-parametric' for
                        Mann-Whitney two-sided test (R wilcox.test, default
                        args).

"""
# FUNCTIONS
# takes files from arguments and collapses them into one table per experiment
def collapse_files(data, mapping, count):
    try:
        data_f = pd.read_csv(data, index_col = 0).T
    except FileNotFoundError:
        print("ERROR: file ", data, " does not exist!\n")
        sys.exit()
    try:
        map_f = pd.read_csv(mapping, header=None)
    except FileNotFoundError:
        print("ERROR: file ", mapping, " does not exist!\n")
        sys.exit()
    data_f.index.names = ['sample']
    map_f.columns = ['sample', 'group']
    map_f.set_index('sample', inplace=True)
    map_f['expt'] = count+1
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

# performs either parametric or non-parametric two-tailed two sample t-test
def dispatch_if_T(operatorT, x, y):
    if operatorT == 'parametric':
        return list(r.r['wt'](FloatVector(x),FloatVector(y)))[0]
    elif operatorT == 'non-parametric':
        return list(r.r['wx'](FloatVector(x),FloatVector(y)))[0]
    else:
        return None

# converting table into a dict of dict-s
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

# getting p-values for all groups as a table
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

# adding statistic index and sorting the dataframe based on first index (group)
def r_name(table, statistic):
    new_table = (table.set_index([table.index, [statistic]* len(table)])
    .rename_axis(['group', 'statistic'])
    .sort_index(axis=0, level=0))
    return new_table

#%%

if __name__ == "__main__":
    
    try:
        import pandas as pd
        import argparse
        import sys
        from rpy2 import robjects as r
        from rpy2.robjects.vectors import FloatVector
        r.r.source('r_statfun.R')

    except ModuleNotFoundError:
        print("This script requires packages: sys, pandas, argparse and \
        scipy.stats.\n")

# parsing arguments
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
    parser.add_argument("--test_statistic", "-tP", type=str, required=True, \
                        nargs=1, help="Normalization statistic: " + \
                        "'mean' or 'median'.")
    
    parser.add_argument("--st_test", "-stT", type=str, required=True, \
                        nargs=1, help="T-test type: 'parametric' for Welch "+\
                        "two-sided t-test (R t.test, default args), " + \
                        "or 'non-parametric' for Mann-Whitney two-sided " + \
                        "test (R wilcox.test, default args).")
    
    # handling no arguments
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit()

    args = parser.parse_args()

    # checking if each data file has corresponding mapping file and vice versa
    if len(args.d_file) == len(args.m_file):
        
        # setting several operators to shorter names for shorter scripting
        contr_g = args.control[0]
        operator = args.test_statistic[0]
        operatorT = args.st_test[0]
        
        # combining mapping and data files into one table
        for count, d_file in enumerate(args.d_file):
            if count == 0:
                full_table = collapse_files(d_file, args.m_file[count], count)
            elif count > 0:
                t_table =  collapse_files(d_file, args.m_file[count], count)
                full_table = full_table.append(t_table)
        
        # next - for each expt, only for treatment group do mean or median test 
        # to get normalization statistic for each column
        c_all = full_table.loc[pd.IndexSlice[:, contr_g],:]
        m_normz = dispatch_if(operator, c_all.groupby(level='expt', axis=0))
        nfull_table = full_table.subtract(m_normz, axis=1,level='expt')
        
        # statistical comparison of pooled groups (each) to the  
        # control pooled group
        test_table = nfull_table.copy()
        test_table.index = test_table.index.droplevel([0,2])
        #test_table = nfull_table.droplevel([0,2])
        pvalues = pval_calc(table_to_dic(test_table), set(test_table.index))
        pvalues = r_name(table=pvalues, statistic='pvalue')
        
        # getting fdr values for each readout
        pvalues_noC = pvalues.drop(index=contr_g).sort_index()
        fdr_table = (pd.DataFrame(index=set(test_table.index))
        .drop(index=contr_g)
        .sort_index())
        fdr_table = (fdr_table
                     .set_index([fdr_table.index, ['fdr']*len(fdr_table)]))
        fdr_table.index.names = ['group', 'statistic']
        for i in pvalues_noC:
            fdr_table[i] = r.r['fdr'](FloatVector(pvalues_noC[i]))
            
        # getting statistics table for second output table (for all groups)
        table_m = dispatch_if(operator,nfull_table.groupby(level='group', \
                                                           axis=0))
        table_m = r_name(table=table_m, statistic=operator)
        
        # export normalized table and table with calculated statistics,
        # p-values & fdr values
        nfull_table.T.to_csv (r'normalized_full_table.csv', index=True, \
                              header=True)
        (table_m.append(pvalues).append(fdr_table).T
         .to_csv(r'statistics_table.csv', index=True, header=True))

    else:
        print('\n\n' + "Error: Data and mapping files must be provided " + \
              "in pairs in respective order." + '\n\n')
        parser.print_help(sys.stderr)
        sys.exit()