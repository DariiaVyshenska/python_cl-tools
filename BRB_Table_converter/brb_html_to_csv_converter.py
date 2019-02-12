"""Script that transforms brb comparison .html table to .csv table.

From a command line provide a) absolute path/name of .html file(s) that need
to be converted to .csv files; b) --fcR parameter if fold change value
should be reversed; no --fcR argument if fold change should remain intact.
This code will accept only .html files, it will alter a header, remove
count column and substitute <1e-07 value with 0.0000001 value.
Example of the call: 
    brb_html_to_csv_converter.py * --fcR
"""

import re
import sys

# function that reads from .html file, transforms the text and writes
# transformed text into .csv file
# as arguments takes the file path/name and Fold change value - either
# forward or reverse
def transFile(infile):
    outfile = infile.replace(".html", ".csv")
    
    fileI = open(infile, "r")
    fileO = open(outfile, "w")
    lines = fileI.readlines()
 
    for line in lines:
        # getting names of two classes that were compared (to use later)
        if "<font size=2>Class 1:" in line:
            groups = line.split(";")    
            gr1 = re.sub('<font size=2>Class 1: <I>|</I></font>',\
            '', groups[0])
            gr2 = re.sub(' <font size=2>Class 2: <I>|</I></font>.\n',\
            '', groups[1])
        # getting and transforming the header
        if "<TR><TH>&nbsp;" in line:
            head = line.replace("&nbsp;  ", "Count ")
            head = head.replace(" <TH> ", ",")
            head = re.sub('<TR><TH>|<NOBR>|</NOBR>', '', head)
            head = re.sub('\s\s', ' ', head)
            splithead = head.split(",")
            # removing count column from the table
            splithead.pop(0)
            # adding description of which class geom. mean was devided
            # by which class geom. mean (simply - updating Fold-change
            # column name)
            FCindex = splithead.index("Fold-change")
            if fc_flag:
                splithead[FCindex] = splithead[FCindex] + "(" + gr2 + \
                "/" + gr1 + ")"

            elif not fc_flag:
                splithead[FCindex] = splithead[FCindex] + "(" + gr1 + \
                "/" + gr2 + ")"
                
            newhead = ",".join(splithead)
            # writing the header to output file
            fileO.write(newhead)
        # reading and writing the body of the table
        elif "<TR><TD>" in line:
            # cleaing lines
            lineW = line.replace("<BR><TD>", ",")
            lineW = re.sub('<\w*>|</.>|<A HREF=".*" TARGET=_blank>|\s', \
            "", lineW)
            splitlineW = lineW.split(",")
            # removing count column from body of the table
            splitlineW.pop(0)
            # for p-value and fdr column - substituting "<1e-07" string
            # for the number "0.0000001"
            if splitlineW[0] == "<1e-07":
                splitlineW[0] = "0.0000001"
            if splitlineW[1] == "<1e-07":
                splitlineW[1] = "0.0000001"
            # based on desired appearance of Fold-change - changing or
            # keeping intact the Fold-change column
            if fc_flag:
                FCvalue = splitlineW.pop(FCindex)
                newFCvalue = 1/float(FCvalue)
                splitlineW.insert(FCindex,str(newFCvalue))
                newlineW = ",".join(splitlineW)
            elif not fc_flag:
                newlineW = ",".join(splitlineW)
            newlineW+= "\n"
            # writing the line to a file
            fileO.write(newlineW)
            
            
    fileO.close()
    fileI.close()

if __name__ == "__main__":
    
    # reading arguments from the command line to: 1) get fc argument if
    # specified; identify all .html files.
    inputL = sys.argv
    
    fc_flag = False
    if ("--fcR" in inputL):
        fc_flag = True
    htmlfiles = [f for f in inputL if ".html" in f]
    
    # looping through all identified html files and applying transFile 
    # function to each of them.
    for htmlfile in htmlfiles:
        transFile(htmlfile)
    print("Finished!")

