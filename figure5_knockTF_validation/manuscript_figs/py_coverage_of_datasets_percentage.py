import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
#sns.despine(offset=0, trim=True)




def main():

    outdir="figs_coverage_of_datasets"
    os.makedirs(outdir,exist_ok=True)

    method_match_names={'BART':'f2_bart_summary/BART',
                        'HOMER':'f3_homer_summary/HOMER',
                        'Pscan':'f4_pscan_summary/Pscan',
                        'TFEA.ChIP':'f6_tfea_summary/TFEA',
                        'ChEA3':'f7_chea3_summary/CHEA3',
                        'old_BART':'f8_old_bart_summary/old_BART',}
#                         'chip_atlas':'f9_chip_atlas/ChIP-Atlas',
#                         'lisa':'f10_lisa/Lisa'}

    method_match_colors={'BART':'r',
                        'HOMER':'goldenrod',
                        'Pscan':'g',
                        'TFEA.ChIP':'b',
                        'ChEA3':'purple',
                        'old_BART':'darkgrey',
                        'chip_atlas':'k',
                        'lisa':'blue'}
    
    lable_matchness=   {'BART':'BARTweb',
                        'HOMER':'HOMER',
                        'Pscan':'Pscan',
                        'TFEA.ChIP':'TFEA.ChIP',
                        'ChEA3':'ChEA3',
                        'old_BART':'BART1.1',
                        'chip_atlas':'ChIP-Atlas',
                        'lisa':'Lisa'}


    results_sub_names=['deg_FC1.5','deg_FC1.5_p005']
#     thres=[[0.1,0.01],[0.2,0.05],[0.2,0.01],[1,0.05],[1,0.01]] # %rank & p-value
    thres=[[0.1,0.01]] 
    
    
    for thre in thres:
        for ii in np.arange(1):
            # ==== plot compr tools
            plt.figure(figsize=(3,3))
            compr_df = pd.DataFrame()
#             with_results_index = set()
            for key in method_match_names.keys():
                result_file = '../fz_results_compr_append/{}_{}_summary.csv'.format(method_match_names[key],results_sub_names[ii])
                df = pd.read_csv(result_file,index_col=0,)#;print(df,result_file);exit()
                # ==== total number of datasets and TFs
                dataset_total = df.shape[0] # total number of datasets
                tf_total = df['Up_total'].dropna()[0]#;print(tf_total,key);exit()
                kept_NotNull_col = ['Final_rank','Final_pvalue']
                # ==== with true prediction
                df_true = df.loc[(df['Final_rank']/tf_total<thre[0])&(df['Final_pvalue']<thre[1])]
                
#                 kept_col=['Final_rank','Final_pvalue']
#                 df = df[kept_col].replace(-1,np.nan).dropna();with_results_index = with_results_index.union(df.index)#;print(key,len(with_results_index))
#                 df = df.dropna()
#                 with_results_index = with_results_index.union(df.index)#;print(key,len(with_results_index))
#                 none_results = total - df.shape[0]
                
                # == coverage of TF with prediction
                data_coverage = df[kept_NotNull_col].dropna().index
                data_coverage_true = df_true[kept_NotNull_col].dropna().index
                compr_df.loc[key,'data_coverage'] = len(data_coverage)
                compr_df.loc[key,'data_coverage_true'] = len(data_coverage_true)
                # plot to compr
                label_key=lable_matchness[key]
                plt.scatter(len(data_coverage),len(data_coverage_true)/len(data_coverage),label=label_key,color=method_match_colors[key])
                if key in ['Pscan','ChEA3']:
                    plt.text(len(data_coverage)-33,len(data_coverage_true)/len(data_coverage)-.03,label_key,fontsize=12,ha='left')
                elif key in ['BART']:
                    plt.text(len(data_coverage)-62,len(data_coverage_true)/len(data_coverage)+.01,label_key,fontsize=12,ha='left')
                elif key in ['TFEA.ChIP']:
                    plt.text(len(data_coverage)-55,len(data_coverage_true)/len(data_coverage)+.01,label_key,fontsize=12,ha='left')
                else:
                    plt.text(len(data_coverage)-5,len(data_coverage_true)/len(data_coverage)+.01,label_key,fontsize=12)
                
            plt.xlim([150,375])
            plt.axes().set_xticks([150,200,250,300,350])
            
            plt.ylim([0,.35])
#             plt.plot([0,1000],[0,280],lw=.7,c='k',ls='--')
            plt.xlabel('Dataset coverage',fontsize=17)
            plt.ylabel('Fraction of datasets \n with true prediction',fontsize=17)
            figname = outdir+os.sep+'compr_coverage_of_data_{}_rthre{}_pthre{}_percentage.pdf'.format(results_sub_names[ii],thre[0],thre[1])
            plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
            plt.close()
#             exit()
    
            compr_df['%']=compr_df['data_coverage_true']/compr_df['data_coverage']
            compr_df.to_csv(outdir+os.sep+'compr_coverage_of_data_{}_rthre{}_pthre{}.csv'.format(results_sub_names[ii],thre[0],thre[1]))
                            
        
        


           

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
