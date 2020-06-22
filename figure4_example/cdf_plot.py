import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def cdf_plot(ID,target,background,figname):
    
    
    fig=plt.figure(figsize=(2.6,2.6))   
    dx = 0.01
    x = np.arange(0,1,dx)       
    by,ty = [],[]      
    for xi in x:
        by.append(sum(i< xi for i in background )/len(background))
        ty.append(sum(i< xi for i in target )/len(target))
    b = plt.plot(x,ty,'r-',label='{}'.format(ID)) 
    a = plt.plot(x,by,color='dimgrey',label='Background') 
    plt.legend(fontsize = 12,frameon=False,borderaxespad=0.,handletextpad=0.2,labelspacing=.1,handlelength=1,loc='upper left')
    #maxval = max(background)
    #minval = min(background)
    #plt.ylim([0,1])
    plt.xlim([0.0,1])
    plt.ylabel('Cumulative Fraction')
    plt.xlabel('AUC')
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()





def main_plot(infile,outdir):

    
    basename = os.path.basename(infile).split('.txt')[0]
    ID='POU5F1'
    lines = open(infile).readlines()
    target = []
    background = []
    for line in lines:
        score = float(line.strip().split('= ')[-1])
        background.append(score)
        if line.startswith(ID):
            target.append(score) 
    figname = outdir+os.sep+basename+'.pdf'
    cdf_plot(ID,target,background,figname)

    
def main(indir,outdir):
    
    os.makedirs(outdir,exist_ok=True)

    infiles = glob.glob('{}/*auc.txt'.format(indir))
    for infile in infiles:
        main_plot(infile,outdir)   




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
#     parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
#     parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>',default='./')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='figs')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.indir,args.outdir)
