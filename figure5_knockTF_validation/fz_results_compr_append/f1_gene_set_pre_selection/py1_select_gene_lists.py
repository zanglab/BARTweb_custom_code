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

    
    deg_dir="/nv/vol190/zanglab/zw5j/since2019_projects/marge_bart_hic_lasso/f5_knockTF_validation/f1_deg_identification/"
    sub_dir_names = [['deg_FC1.5','deg_FC1.5_append'],['deg_FC1.5_p005','deg_FC1.5_p005_append']]
    
    for dir_name in sub_dir_names:
        outf = open(dir_name[0]+'_name_list.txt','w')
        total_names=set()
        for sub_dir_name in dir_name:
            deg_sub_dir = deg_dir+os.sep+sub_dir_name
            up_files = glob.glob(deg_sub_dir+"/*UpGenes.txt")#;print(deg_sub_dir)
            dn_files = glob.glob(deg_sub_dir+"/*DnGenes.txt");print(deg_sub_dir)
            up_basename = [os.path.basename(up_file).split('_UpGenes.txt')[0] for up_file in up_files]
            dn_basename = [os.path.basename(dn_file).split('_DnGenes.txt')[0] for dn_file in dn_files]
            basenames = set(up_basename+dn_basename)
            if len(total_names)==0:
                total_names = basenames
            else:
                total_names = total_names.union(basenames)

        outf.write('\n'.join(total_names)+'\n')
        outf.close()



           

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
