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


def return_index_rank(tf_name,df_index):

    for rank_i in np.arange(len(df_index)):
        if re.match(tf_name,df_index[rank_i].upper()):
            return rank_i+1
    return -1

def main():

    # change tool subname if needed
    results_dir="figure5_knockTF_validation"
    results_sub_dirs = ['f6_TFEA_ChIP/tfea_deg_FC1.5']
    method_name="TFEA"
    
    gene_list_dir="figure5_knockTF_validation/fz_results_compr_append/f1_gene_set_pre_selection/"
    gene_sub_names=['deg_FC1.5']
    
    for ii in np.arange(2):
        gene_list_file =  gene_list_dir+os.sep+gene_sub_names[ii]+'_name_list.txt'
        gene_lists = [i.strip() for i in open(gene_list_file).readlines()]#;print(len(gene_lists));exit()
        
        # collect results for each type of DEG
        out_df = pd.DataFrame(columns=['TF','Up_rank','Up_total','Up_pvalue','Dn_rank','Dn_total','Dn_pvalue','Final_rank','Final_pvalue'])
        for gene_list in gene_lists[:]:
            tf_name = gene_list.split('_DataSet')[0];print(gene_list)
            out_df.loc[gene_list]=[tf_name,-1,-1,-1,-1,-1,-1,-1,-1]
            up_result_file = results_dir+os.sep+results_sub_dirs[ii]+os.sep+'{}_UpGenes.txt'.format(gene_list)
            if not os.path.isfile(up_result_file):
                up_result_file = results_dir+os.sep+results_sub_dirs[ii]+"_append"+os.sep+'{}_UpGenes.txt'.format(gene_list)
            dn_result_file = results_dir+os.sep+results_sub_dirs[ii]+os.sep+'{}_DnGenes.txt'.format(gene_list)
            if not os.path.isfile(dn_result_file):
                dn_result_file = results_dir+os.sep+results_sub_dirs[ii]+"_append"+os.sep+'{}_DnGenes.txt'.format(gene_list)
            
            if os.path.isfile(up_result_file):
                up_df = pd.read_csv(up_result_file,sep=',',index_col=4)#;print(up_df);exit()
                if tf_name in up_df.index:
                    out_df.loc[gene_list,'Up_rank'] = list(up_df.index).index(tf_name)+1
                    out_df.loc[gene_list,'Up_pvalue'] = up_df.iloc[list(up_df.index).index(tf_name)]['adj.p.value']
#                 out_df.loc[gene_list,'Up_rank'] = return_index_rank(tf_name,up_df.index)
                out_df.loc[gene_list,'Up_total'] = len(up_df.index)
            if os.path.isfile(dn_result_file):
                dn_df = pd.read_csv(dn_result_file,sep=',',index_col=4)
                if tf_name in dn_df.index:
                    out_df.loc[gene_list,'Dn_rank'] = list(dn_df.index).index(tf_name)+1
                    out_df.loc[gene_list,'Dn_pvalue'] = dn_df.iloc[list(dn_df.index).index(tf_name)]['adj.p.value']
        
#                 out_df.loc[gene_list,'Dn_rank'] = return_index_rank(tf_name,dn_df.index)
                out_df.loc[gene_list,'Dn_total'] = len(dn_df.index)
            
        # select either up or down as final rank, which is smaller
        for index in out_df.index:
            if (out_df.loc[index,'Up_rank']>0 and out_df.loc[index,'Up_rank'] < out_df.loc[index,'Dn_rank']) or (out_df.loc[index,'Dn_rank']<0):
                out_df.loc[index,'Final_rank'] = out_df.loc[index,'Up_rank']
                out_df.loc[index,'Final_pvalue'] = out_df.loc[index,'Up_pvalue']
            elif (out_df.loc[index,'Dn_rank']>0 and out_df.loc[index,'Up_rank'] > out_df.loc[index,'Dn_rank']) or (out_df.loc[index,'Up_rank']<0):
                out_df.loc[index,'Final_rank'] = out_df.loc[index,'Dn_rank']
                out_df.loc[index,'Final_pvalue'] = out_df.loc[index,'Dn_pvalue']
       
        out_df = out_df.sort_values(by=['Final_pvalue'])
        out_df = out_df.replace(-1,np.nan)#.dropna(axis=0)
            
        out_df.to_csv('{}_{}_summary.csv'.format(method_name,gene_sub_names[ii]))
            
       



           

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
