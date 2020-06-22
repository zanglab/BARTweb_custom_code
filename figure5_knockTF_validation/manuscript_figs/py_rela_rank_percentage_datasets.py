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


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])



def main():

    outdir="figs_rela_rank_percentage_datasets"
    os.makedirs(outdir,exist_ok=True)

    method_match_names={'BART':'f2_bart_summary/BART',
                        'HOMER':'f3_homer_summary/HOMER',
                        'Pscan':'f4_pscan_summary/Pscan',
                        'TFEA.ChIP':'f6_tfea_summary/TFEA',
                        'ChEA3':'f7_chea3_summary/CHEA3',
                        'old_BART':'f8_old_bart_summary/old_BART',
                        'chip_atlas':'f9_chip_atlas/ChIP-Atlas'}

    method_match_colors={'BART':'r',
                        'HOMER':'goldenrod',
                        'Pscan':'g',
                        'TFEA.ChIP':'b',
                        'ChEA3':'purple',
                        'old_BART':'gray',
                        'chip_atlas':'k'}
    
    lable_matchness=   {'BART':'BARTweb',
                        'HOMER':'HOMER',
                        'Pscan':'Pscan',
                        'TFEA.ChIP':'TFEA.ChIP',
                        'ChEA3':'ChEA3',
                        'old_BART':'BART1.1',
                        'chip_atlas':'ChIP-Atlas'}

    

    results_sub_names=['deg_FC1.5','deg_FC1.5_p005']
    thres=[[0.1,0.01]] 
    top_thres = [0.1,0.01,0.005]
    colors = ['lightgray','grey','k']
#     colors = ['lightskyblue','royalblue','navy']
    labels = ['Top 1%-10%','Top 0.5%-1%','Top 0.5%']
    for thre in thres:
        for ii in np.arange(1):
            # ==== plot compr tools
            plt.figure(figsize=(3,3))
            compr_df = pd.DataFrame()
            key_pos = 0
            keys=['BART','old_BART','TFEA.ChIP','ChEA3','Pscan','HOMER'][::-1]
            keys_labels=['BARTweb','BART1.1','TFEA.ChIP','ChEA3','Pscan','HOMER'][::-1]
#             keys=['BART','old_BART','TFEA.ChIP','ChEA3','Pscan','HOMER']
#             keys_labels=['BARTweb','BART1.1','TFEA.ChIP','ChEA3','Pscan','HOMER']
            for key in keys:
                result_file = '../fz_results_compr_append/{}_{}_summary.csv'.format(method_match_names[key],results_sub_names[ii])
                df = pd.read_csv(result_file,index_col=0,)#;print(df,result_file);exit()
                # ==== total number of datasets and TFs
                dataset_total = df.shape[0] # total number of datasets
                tf_total = df['Up_total'].dropna()[0]#;print(tf_total,key);exit()
                kept_NotNull_col = ['Final_rank','Final_pvalue']
                
                # ==== with true prediction
                df_true = df.loc[(df['Final_rank']/tf_total<thre[0])&(df['Final_pvalue']<thre[1])]
                # == coverage of TF with prediction
                data_coverage = df[kept_NotNull_col].dropna().index
                data_coverage_true = df_true[kept_NotNull_col].dropna().index
                compr_df.loc[key,'data_coverage'] = len(data_coverage)
                compr_df.loc[key,'data_coverage_true'] = len(data_coverage_true)
                for bar_pos in np.arange(len(top_thres)):
                    top_thre = top_thres[bar_pos]
                    top_val = len(df_true.loc[df_true['Final_rank']/tf_total<top_thre][kept_NotNull_col].dropna().index)
                    compr_df.loc[key,'true_top{}'.format(top_thre)] = top_val

                    if top_thre==0.005:
                        color = colors[bar_pos]
                        color = lighten_color(method_match_colors[key],1.33)
                        # color = method_match_colors[key]
                        alpha= 1
                        hatch=''
                    elif top_thre==0.01:
                        color = method_match_colors[key]
                        alpha= .66
                        hatch=''
                    else:
                        color = method_match_colors[key]
                        alpha= .25
                        hatch=''
                    print(thre,key,color,alpha)
                    
                    plt.bar(key_pos,top_val,color=color,lw=0,edgecolor='k',width=0.66,alpha=alpha,hatch=hatch)
                    if key=='old_BART':
                        plt.bar(key_pos,top_val,color=color,lw=0,edgecolor='k',width=0.66,label=labels[bar_pos],alpha=alpha,hatch=hatch)

                key_pos+=1
            

#             plt.ylim([0,110])
            plt.legend(fontsize=13,\
                borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper left",bbox_to_anchor=[0,1],frameon=False)
            plt.axes().set_xticks(np.arange(len(keys)))
            plt.axes().set_xticklabels(keys_labels,fontsize=15,rotation=45,ha='right')
            plt.ylabel('Number of datasets',fontsize=17)
            figname = outdir+os.sep+'compr_rela_rank_{}_rthre{}_pthre{}.pdf'.format(results_sub_names[ii],thre[0],thre[1])
            plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
            plt.close()

            compr_df.to_csv(outdir+os.sep+'compr_rela_rank_{}_rthre{}_pthre{}.csv'.format(results_sub_names[ii],thre[0],thre[1]))
                            
        
        


           

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
