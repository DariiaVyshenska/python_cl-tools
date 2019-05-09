########################################################################
#
# This script was written using Python 3.6.4
# It requires packages scikit-allel (1.2.0), requests (2.18.4), and sys.
#
# Goal: variant annotation
# Input:         .vcf file(s)
# Output:          two .csv files
#     1) .csv (same name as corresponding .vcf file) contains following 
#     information:
#         CHROM – chromosome
#         REF - Reference allele nucleotide(s)
#         ALT - Alternate allele nucleotide(s)
#         POS – position of alternate allele
#         TYPE - The type of allele, either snp, mnp, ins, del, or complex
#         AO - Alternate allele observation count
#         DP - Read Depth
#         DPRA - Alternate allele depth ratio
#         FREQ - Allele frequency of variant from Broad Institute ExAC 
#         Project API
# If a several alleles were detected at the position, first most 
# deleterious was chosen for information extraction (priority, with 1 as
# most and 5 as least deleterious: del - 1, ins - 2, complex - 3, 
# mnp - 4, snp - 5).
#
#     2) _extr.csv contains raw fields extracted from .vcf file, can be 
#    used for debugging.
#
# Usage:
# e.g. input for annotation:
#
#    python tempus_var_annotation.py input.vcf
#
#     input.vcf          :  The .vcf file(s) to be annotated, can be *.
#
# Additional options:
#     -h                 :  Help, displays this page
#
########################################################################



try:
    import sys
    import requests
    import allel

except ModuleNotFoundError:
    print("This script requires packages: sys, requests and \
    allel[full]!\n")

alt_prior_dic = {'del': 1, 'ins': 2, 'complex': 3, 'mnp':4, 'snp':5, \
'': 6}

def priority_index(in_data):
    """Requires alt_prior_dic dictionary with prioritized variant types.
    
    Function: input - a list of allele's types detected for a given 
    variant, output - index number (integer) of the type of allele 
    withing input list that has highest deleterious priority (priority
    is set in alt_prior_dic with the highest being 1, lowest being 6."""
    
    keys_in_data = dict((k, alt_prior_dic[k]) for k in in_data)
    priority_key = min(keys_in_data, key=keys_in_data.get)
    prior_i = in_data.index(priority_key)
    return(int(prior_i))


def allele_freq(var_info):
    """Function: for a given variant this function calls Broad Institute 
    ExAC Project API and returns information about the frequency of this
    variant in the database. input of function is a list of chromomosome
    number, position of variant, reference sequence, new alternative 
    sequence. Output - allele frequency. If allele frequency is not 
    available, function returns NA."""
    
    dbcall = "http://exac.hms.harvard.edu/rest/variant/variant/"
    url = dbcall + '-'.join(var_info)
    r = requests.get(url)
    response_dict = r.json()
    if "allele_freq" in response_dict:
        return response_dict['allele_freq']
    else:
        return "NA"


def conv_file(extr_file):
    """Function: reads _extr.csv file that contains extracted information
    from .vcf file, reduces number of observed allels per variation 
    to 1, gets allel frequency from Broad Institute ExAC Project that
    allel, outputs .csv file with chromocome, ref. sequence, alt. 
    sequence, position, type of variation, alternate allele observations,
     with partial observations recorded fractionally, total read depth 
     at the locus, alternate allele depth ratio, allele frequency."""
    
    outfile = extr_file.replace("_extr.csv", ".csv")
    fileO = open(outfile, "w")
    
    # for a _extr.csv file with extracted fields from .vcf, new .csv
    # is created (same name as .vcf file, different extension.
    with open(extr_file) as fileI:
        header = fileI.readline().rstrip().split(',')
        
        # extracting indexes for relevant infromation from the header
        chrom_i = header.index("CHROM")
        ref_i = header.index("REF")
        alt_i = [i for i, s in enumerate(header) if 'ALT' in s]
        pos_i = header.index("POS")
        type_i = [i for i, s in enumerate(header) if 'TYPE' in s]
        ao_i = [i for i, s in enumerate(header) if 'AO' in s]
        dp_i = header.index("DP")
        dpra_i = [i for i, s in enumerate(header) if 'DPRA' in s]
        
        # creating and writing to a file new header
        new_header = ['CHROM', 'REF', 'ALT', 'POS', 'TYPE', 'AO', 'DP',\
        'DPRA', 'FREQ', '\n']
        new_header = ','.join(str(e) for e in new_header)
        fileO.write(new_header)
        
        # reading the rest of the lines in the file
        lines = fileI.readlines()
        for line in lines:
            line_l = line.rstrip().split(',')
            
            # extracting infromation that is always present as a single
            # number
            chrom = line_l[chrom_i]
            ref = line_l[ref_i]
            pos = line_l[pos_i]
            dp = line_l[dp_i]
            
            # extract infromation about type of allel
            type_raw = [line_l[i] for i in type_i]
            # test if there is more then one type of allel detected
            # for the variant
            mult_test = len(list(filter(None, type_raw)))
            if mult_test == 1:
                # if only one type - extract the rest of information
                type_new = line_l[type_i[0]]
                alt_new = line_l[alt_i[0]]
                ao_new = line_l[ao_i[0]]
                dpra_new = line_l[dpra_i[0]]

            elif mult_test > 1:
                # if there is several - extract the rest of the 
                # information only for
                # the most deleterious mutation (or, if there are 
                # several with the highest priority - extracting the
                # first one
                type_priority = priority_index(type_raw)
                type_new = type_raw[type_priority]    
                alt_new = line_l[alt_i[0]+type_priority]
                ao_new = line_l[ao_i[0]+type_priority]
                dpra_new = line_l[dpra_i[0]+type_priority]

            else:
                print("Error: The variant type column is empty!")
                
            # getting the frequency for extracted variant type
            v = [chrom, pos, ref, alt_new]
            all_freq = allele_freq(v)
            
            # creating and writing to .csv file new line with 
            # extracted information
            new_line = [chrom, ref, alt_new, pos, type_new, ao_new, dp,\
            dpra_new, all_freq, '\n']
            new_line = ','.join(str(e) for e in new_line)
            fileO.write(new_line)

    fileI.close()


usage = ['Usage:',
  'e.g. input for annotation:',
  '',
  '   python tempus_var_annotation.py input.vcf',
  '',
  '    input.vcf          :  The .vcf file(s) to be annotated, can be *.',
    '',
    'Additional options:',
  '    -h                 :  Help, displays this page',
  ]

if __name__ == "__main__":

    inputL = sys.argv
    # extracting only .vcf files
    vcf_files = [f for f in inputL if ".vcf" in f]

    if "-h" in inputL or len(vcf_files)==0:
        print("\n".join(usage))
    else:
        # for each .vcf file provided to the script
        for vcf_file in vcf_files:
            try:
                # if input contains .vcf files, script extracts
                # chromosome number, reference and alternate sequence, 
                # position, type of variation, allele observations, read 
                # depth and  alternate allele depth ratio columns into
                # new _extr.csv file
                extr_file = vcf_file.replace('.vcf', '_extr.csv')
                ext_fields = ['CHROM', 'REF', 'ALT', 'POS','TYPE', 'AO',\
                'DP', 'DPRA']
                allel.vcf_to_csv(vcf_file, extr_file, fields=ext_fields)
                # after _extr.csv file is created, script reduces number
                # of alternate variations (to be equal 1) by set priority
                # and exports the results as new .cvs file
                conv_file(extr_file)
                print("File ", vcf_file, "is processed!")
            except FileNotFoundError:
                print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("Attention: file ", vcf_file, " does not exist!\n")
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
                
