# -*- coding: utf-8 -*-
import scipy.stats as scst
import pandas as pd
import numpy as np
from itertools import chain
import multiprocessing as mp
from functools import partial
import time
import argparse
import sys

def dispatch_if(operator, x, y):
    """dispatch function for type of statistical test"""
    if operator == 'spearman':
        return scst.spearmanr(x, y)
    elif operator == 'pearson':
        return scst.pearsonr(x, y)
    else:
        return None

def cor_fun(operator, row1, row2):
    if sum(row1) == 0 or sum(row2) == 0:
        return "NA", "NA"
    else:
        coef, p = dispatch_if(operator, row1, row2)
        return coef, p

def calc_cor(operator, data_v2, data_v1):
    # getting table of correlation coefficients and their p-values
    df_corr = data_v1.apply(lambda y: \
        data_v2.apply(lambda x: cor_fun(operator, y, x), axis=1), axis=1)
    #df_corr.to_csv("corr_test.csv")
    
    # unpacking matrix into two lists
    val_as_list = list(chain.from_iterable(df_corr.values))
    cor, p = zip(*val_as_list)
    
# TRY TO IMPROVE THIS PART AND THEN MOVE ON (EVERYTHING ELSE IS GOOD TO GO):
    # creating comparison pair names (from rows and columns of df)
    r = np.array(df_corr.index)
    c = np.array(df_corr.columns)
    
    name1, name2, names = list(), list(), list()
    
    for x in r:
        for y in c:
            name1.append(x)
            name2.append(y)
            names.append(x + "<==>" + y)
#######################################################
            
    comb_output = list(zip(name1, name2, names, cor, p))
    return comb_output


# split one of the dataframes into chanks
def df_chunking(df, chunksize):
    """Splits df into chunks, drops data of original df inplace"""

    count = 0 # Counter for chunks
    while len(df):
        count += 1
        print('Preparing chunk {}'.format(count))
        # Return df chunk
        yield df.iloc[:chunksize].copy()
        # Delete data in place because it is no longer needed
        df.drop(df.index[:chunksize], inplace=True)

def run_corr_parallel(data_v1, data_v2, cores, operator):
    """Performs splitting initial table into subtables and doing 
    parallel correlation analysis"""

    # df_chunking is a generator. I can use it to give it to multipro or I can 
    # place chunks into a list like below. i may need to try both ways
    chunksize = data_v1.shape[0] // cores + 1
    test_list = list()
    for i in df_chunking(data_v1, chunksize):
        test_list.append(i)
        
    # figuring how to feed chanks and second dataframe to parallel python
    # (I have 4 cores)
    mp.freeze_support()

    
    ctx = mp.get_context('spawn')
    pool = ctx.Pool(cores)
    func = partial(calc_cor, operator, data_v2)
    results = pool.imap(func, test_list)
    pool.close()
    pool.join()
    
    joined_res = list()
    for res in results:
        joined_res += res
        
    return pd.DataFrame(joined_res)

def main():
    
    # parsing arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("--d_file", "-df", type=str, required=True, nargs=1,\
                        help="Files with data: csv file with a header. "+ \
                        "Header - unique sample identifiers. " + \
                        "First column - variable ID.")
    parser.add_argument("--v_file", "-vf", type=str, required=True, nargs='+', \
                        help="Files with variables of interest: csv files" + \
                        "without a header. Contain one column of variable" +\
                        " IDs from data file")
    parser.add_argument("--m_file", "-mf", type=str, required=True, nargs=1,\
                        help="Mapping file: csv file with a header.")
    parser.add_argument("--group", "-g", type=str, required=True, nargs=1,\
                        help="Name of the column in the mapping file" + \
                        "that contains group identifiers")
    parser.add_argument("--sample_id", "-s_id", type=str, required=True, 
                        nargs=1, help="Name of the column in the mapping " + \
                        "file that contains sample identifiers")
    parser.add_argument("--st_test", "-stT", type=str, required=True, \
                        nargs=1, help="Statistical test type: "+\
                        "'spearman' or 'pearson' for respectful " + \
                        "correlation type.  " + \
                        "Other types of tests can be implemented in future.")
    parser.add_argument("--cores", "-c", type=int, required=True, nargs=1,\
                        help="Number of cores to use, integer.")
    parser.add_argument("--out_name", "-oN", type=str, required=True, nargs=1,\
                        help="String to add to output file names.")
    
    parser.add_argument("--a_file", "-af", type=str, required=True, nargs=1,\
                        help="File with analysis info: csv file" + \
                        " without a header. See additional instructions on" + \
                        " how to construct one.")
    
    
    # handling no arguments
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit()

    args = parser.parse_args()

    arg_data = args.d_file[0]
    arg_mapfile = args.m_file[0]
    arg_v1 = args.v_file[0]
    arg_v2 = args.v_file[1]
    group_col = args.group[0]
    id_col = args.sample_id[0]
    operator = args.st_test[0]
    cores = args.cores[0]
    arg_analys = args.a_file[0]
    fout_n = args.out_name[0]

    # reading files
    data_table = pd.read_csv(arg_data, header=0, index_col=0)
    v1=np.array((pd.read_csv(arg_v1, header=None))[0])
    v2=np.array((pd.read_csv(arg_v2, header=None))[0])
    analysis_file = pd.read_csv(arg_analys, header=None)
    maptable = pd.read_csv(arg_mapfile, header=0)
    
    
    # getting names for samples filtering
    
    if "correlation" in np.array(analysis_file[3]):
        for cor_group in analysis_file[0][analysis_file[3] == "correlation"]:
            # creating two dataframes for correlations
            samples = maptable.loc[maptable[group_col] == cor_group, id_col]
            if len(v1) > len(v2):
                data_v1 = data_table.loc[v1,samples].copy()
                data_v2 = data_table.loc[v2, samples]
            else:
                data_v1 = data_table.loc[v2, samples].copy()
                data_v2 = data_table.loc[v1, samples]
                
            f_name = str(fout_n) +str(cor_group) + "_corrout.csv"
            
            if cores > 1:
                start = time.time()
                df_results = run_corr_parallel(data_v1, data_v2, cores, operator)
                end = time.time()
                print(end - start)
                print('\n')
                
            elif cores == 1: 
                start = time.time()
                df_results = pd.DataFrame(calc_cor(operator, data_v2, data_v1))
                end = time.time()
                print(end - start)
                print('\n')
            
            df_results.columns = ['name1', 'name2', 'pairName', operator, 'p-value']
            # removing self-loops
            no_loops = df_results.loc[df_results.name1 != df_results.name2].copy()
            print("DONE!")
            no_loops.to_csv(f_name, index=False)
                
    else:
        pass
        # here you can implement other analysis such as paired or unpaired
        # t-tests or non-parametric group comparisons
#      # getting sample names for two groups (in case I need to add t-tests or 
#      # something similar)
#        all_groups = (maptable[group_col]).dropna().unique()
#     group_id1 = maptable[id_col][maptable[group_col] == all_groups[0]]
#     if len(all_groups) > 1:
#         group_id2 = maptable[id_col][maptable[group_col] == all_groups[1]]
#     elif len(all_groups) ==1:
#            group_id2 = maptable[id_col][maptable[group_col] == all_groups[0]]
    


if __name__ == '__main__':
    main()






