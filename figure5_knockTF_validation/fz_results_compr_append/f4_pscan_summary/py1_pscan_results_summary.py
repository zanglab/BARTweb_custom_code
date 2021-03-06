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


tf_match_file="figure5_knockTF_validation/f4_deg_Pscan/results/terminal_pscan/tf_list.txt"
tf_match_df = pd.read_csv(tf_match_file,sep='\t',index_col=1,header=None)
tf_match_df.index = [i.split('.')[0] for i in tf_match_df.index]
tf_match_df.columns=['TF']
tf_match_df = tf_match_df.loc[~tf_match_df.index.duplicated()]

def return_index_rank(tf_name,df_index):

    for rank_i in np.arange(len(df_index)):
        if re.match(tf_name,df_index[rank_i].upper()):
            return rank_i+1
    return -1

def main():

    # change tool subname if needed
    results_dir="figure5_knockTF_validation/"
    results_sub_dirs = ['f4_deg_Pscan/results/refseq_to_fasta/deg_FC1.5_refseq_fa']
    method_name="Pscan"
    
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
            up_result_file = results_dir+os.sep+results_sub_dirs[ii]+os.sep+'{}_UpGenes.txt.fai.fa.res'.format(gene_list)
            if not os.path.isfile(up_result_file):
                up_result_file = results_dir+os.sep+results_sub_dirs[ii]+"_append"+os.sep+'{}_UpGenes.txt.fai.fa.res'.format(gene_list)
            dn_result_file = results_dir+os.sep+results_sub_dirs[ii]+os.sep+'{}_DnGenes.txt.fai.fa.res'.format(gene_list)
            if not os.path.isfile(dn_result_file):
                dn_result_file = results_dir+os.sep+results_sub_dirs[ii]+"_append"+os.sep+'{}_DnGenes.txt.fai.fa.res'.format(gene_list)
            
            if os.path.isfile(up_result_file):
                up_df = pd.read_csv(up_result_file,sep='\t',index_col=1)#;print(up_df)
                up_df = up_df.iloc[:719,:].join(tf_match_df).sort_values(by=['P_VALUE'])#;print(up_df)#;exit()
                rank = return_index_rank(tf_name,up_df['TF'])
                out_df.loc[gene_list,'Up_rank'] = rank
                out_df.loc[gene_list,'Up_total'] = len(up_df.index)
                if rank !=-1:
                    out_df.loc[gene_list,'Up_pvalue'] = up_df.iloc[rank-1]['P_VALUE']
            if os.path.isfile(dn_result_file):
                dn_df = pd.read_csv(dn_result_file,sep='\t',index_col=1)
                dn_df = dn_df.iloc[:719,:].join(tf_match_df).sort_values(by=['P_VALUE'])#;print(up_df);exit()
                rank = return_index_rank(tf_name,dn_df['TF'])
                out_df.loc[gene_list,'Dn_rank'] = rank
                out_df.loc[gene_list,'Dn_total'] = len(dn_df.index)
                if rank !=-1:
                    out_df.loc[gene_list,'Dn_pvalue'] = dn_df.iloc[rank-1]['P_VALUE']
        
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
