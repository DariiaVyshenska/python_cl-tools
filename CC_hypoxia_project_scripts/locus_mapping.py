# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 13:43:25 2019

@author: vyshe
"""

import pandas as pd
import re
import matplotlib.pyplot as plt

#%%
def chab_freq_count(nw, nc):
    nc_probes = list()
    nc_cyt = list()
    nc_matchsum = list()
    
    for p in nc.locus:
        pattern = re.compile(r'^' + p + '[^0-9]')
        m_probe = list()
        m_cytoband = list()
        for count, band in enumerate(nw.cytoband):
            matchObj = re.match(pattern, band, flags=0)
            if matchObj:
    #            print("matching: ", p)
    #            print ("matchObj.group() : ", matchObj.group())
    #            print("full match cytoband is: ", band)                       
                m_probe.append(nw.Probe_Id[count])
                m_cytoband.append(band)         
    #    print(m_probe)
    #    print(m_cytoband)
    #    print("")
        
        if len(m_probe) > 0:
            nc_matchsum.append(len(m_probe))
            nc_probes.append(';'.join(m_probe))
            nc_cyt.append(';'.join(m_cytoband))
        if len(m_probe) == 0:
            nc_matchsum.append(0)
            nc_probes.append('NA')
            nc_cyt.append('NA')
        
        del m_probe, m_cytoband
    
    nc['in_nw_count'] = nc_matchsum
    nc['nw_probe'] = nc_probes
    nc['nw_cytoband'] = nc_cyt
    
    nc.to_csv("NC_NW_gainloss.csv", index=False)
    
    # vusualization
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.set_xlabel('locus')
    for tick in ax1.get_xticklabels():
        tick.set_rotation(90)
        tick.set_fontsize(5) 
    ax1.set_ylabel('frequency', color=color)
    ax1.plot(nc.locus, nc.gain, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.plot(nc.locus, nc.loss, color='tab:blue')
    ax1.legend(loc=2)
    ax2 = ax1.twinx()
    color = 'tab:grey'
    ax2.set_ylabel('gene count', color=color)
    ax2.plot(nc.locus, nc.in_nw_count, color=color, linewidth=0.5)
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout()
    plt.show()
    ax2.legend(loc=1)
    
    fig.savefig("plot.pdf", bbox_inches='tight')

def match_loc_cytoband(nw, nc):
    m_locus = [None] * len(nw.cytoband)
    
    for p in nc.locus:
        pattern = re.compile(r'^' + p + '[^0-9]')
        for count, band in enumerate(nw.cytoband):
            matchObj = re.match(pattern, band, flags=0)
            if matchObj:
    #            print("matching: ", p)
    #            print ("matchObj.group() : ", matchObj.group())
    #            print("full match cytoband is: ", band)
                if m_locus[count] is None:
                    m_locus[count] = p
                else:
                    print("You need to debug: count ", count)
    new_nw = nw.copy()
    new_nw['Matched locus'] = m_locus
    new_nw.to_csv("nw_ctyoband_locus.csv", index=False)
    
    print("Done!")
    return new_nw


#%%


# importing data
nw_table_fn = "nw_cytobands.csv"
nc_table_fn = "nc_cytobands.csv"
nw = pd.read_csv(nw_table_fn)
nc = pd.read_csv(nc_table_fn)

nw_table_fnPos = "nw_cytobands_Hpos.csv"
nwPos = pd.read_csv(nw_table_fnPos)

nw_table_fnNeg = "nw_cytobands_Hneg.csv"
nwNeg = pd.read_csv(nw_table_fnNeg)

# getting graphs of gene counts on gain-loss graph from NatCom
chab_freq_count(nw, nc)
chab_freq_count(nwPos, nc)
chab_freq_count(nwNeg, nc)

# matching locus from NatCom to cytobands of NW genes:
new_nw = match_loc_cytoband(nw, nc)






  
    















