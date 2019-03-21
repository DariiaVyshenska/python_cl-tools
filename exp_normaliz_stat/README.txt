This script is created to combine multiple readouts from different experiments, normalize the data experiment-wise
(each_measurement - control_group_mean/median), pool normalized data treatment-wise, and run t-test on all
treatment combinations with the control treatment.
Output: table with normalized data, table with mean/median per treatment combined with t-test p-value. Both tables 
are in .csv format.
Require: sys, pandas(0.24.0), argparse(1.1) and scipy.stats(1.0.0).

NOTE: srcipt will be complimented with fdr calculations in future.

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
                        (scipy.stats.ttest_ind), or 'non-parametric' for Mann-
                        Whitney two-sided test with no continuity correction
                        (scipy.stats.mannwhitneyu)
