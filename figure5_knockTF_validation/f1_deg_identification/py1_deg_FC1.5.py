import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
#sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
#sns.despine(offset=0, trim=True)




def main():

#     indir = 'ctcf_sites_output'
    outdir = 'deg_FC1.5_append'
    os.makedirs(outdir,exist_ok=True)

    factor_file="../data/human_factor.txt"
    factor_df = pd.read_csv(factor_file,index_col=0,sep='\t')
    factors = set(factor_df['Factor'])#;print(len(factors));exit()
    
    deg_file='../data/differential_expression_of_genes_in_all_datasets.txt'
    with open(deg_file) as inf:
        df = pd.read_csv(inf,sep='\t',index_col=0,low_memory=False)
    tf_total = set(df['TF']) 
    print('# total TF:',len(tf_total))  
    print('# w/ Cistrome:',len(tf_total.intersection(factors)))  
    print('# total gene list:',len(set(df.index)))  
    
    for geneset_label in set(df.index):
        geneset_df = df.loc[geneset_label]
        tf = geneset_df['TF'].iloc[0]
        if tf not in factors:
#         if 1:
            up_genes = geneset_df[geneset_df['up_down']==1]
            down_genes = geneset_df[geneset_df['up_down']==2]
            up_genes.to_csv(outdir+os.sep+'{}_{}_UpGenes.csv'.format(tf,geneset_label))
            down_genes.to_csv(outdir+os.sep+'{}_{}_DnGenes.csv'.format(tf,geneset_label))
            if up_genes.shape[0]>100:
                up_genes['Gene'].to_csv(outdir+os.sep+'{}_{}_UpGenes.txt'.format(tf,geneset_label),header=None,index=False)
            if down_genes.shape[0]>100:
                down_genes['Gene'].to_csv(outdir+os.sep+'{}_{}_DnGenes.txt'.format(tf,geneset_label),header=None,index=False)
        
        
    






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
