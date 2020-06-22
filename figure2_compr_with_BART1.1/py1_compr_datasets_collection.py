import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"


def bar_plot(scores1,scores2,figname):
    plt.figure(figsize=(2.6,3))
    a = plt.bar(1,scores1[0],width=1,color='darkgrey')
    b = plt.bar(2,scores1[1],width=1,color='r')
    a = plt.bar(4,scores2[0],width=1,color='darkgrey')
    b = plt.bar(5,scores2[1],width=1,color='r')
#     a = plt.bar(7,scores3[0],width=1,color='k')
#     b = plt.bar(8,scores3[1],width=1,color='darkgrey')
    plt.xlim([0,6])
#     plt.ylim([0.5,1])
    plt.ylabel('Number of datasets')
    plt.legend([a,b],['BART1.1','BARTweb'],fontsize=12,\
    borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)#,bbox_to_anchor=[1.5,1]
    plt.axes().set_xticks([1.5,4.5])
    plt.axes().set_xticklabels(['Human','Mouse'],rotation=30, ha='center',fontsize=15,color='k')
#     plt.axes().set_yticks([0.5,0.75,1])
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()


def main():

    outdir='f1_compr_collection'
    os.makedirs(outdir,exist_ok=True)
    
    old_2016_hg=[3490,454]
    old_2016_mm=[3055,409]
    old_2017_hg=[7032,883]
    old_2017_mm=[5023,531]
    new_hg=[7968,918]
    new_mm=[5851,565]
    
    
    ii=0
    scores1 = [old_2016_hg[ii],new_hg[ii]]
    scores2 = [old_2016_mm[ii],new_mm[ii]]

#     scores1 = [old_2016_hg[ii],old_2016_mm[ii]]
#     scores2 = [old_2017_hg[ii],old_2017_mm[ii]]
#     scores3 = [new_hg[ii],new_mm[ii]]
    figname = outdir+os.sep+'compr_datasets.png'
    bar_plot(scores1,scores2,figname)




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
