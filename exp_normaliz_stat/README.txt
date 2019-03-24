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
                        [M_FILE ...] --control CONTROL --test_parameter
                        TEST_PARAMETER --st_test ST_TEST

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
  --test_parameter TEST_PARAMETER, -tP TEST_PARAMETER
                        Normalization parameter: 'mean' or 'median'.
  --st_test ST_TEST, -stT ST_TEST
                        T-test type: 'parametric' for Welch two-sided t-test
                        (R t.test, default args), or 'non-parametric' for
                        Mann-Whitney two-sided test (R wilcox.test, default
                        args).