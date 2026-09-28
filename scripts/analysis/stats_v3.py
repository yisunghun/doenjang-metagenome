import pandas as pd,numpy as np,json,re
from scipy.stats import kruskal,spearmanr,rankdata,t as T
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
m=pd.read_csv('metadata.tsv',sep='\t',index_col=0); S=list(m.index); O=['Me','0M','1M','3M','6M','12M']
order={'Me':0,'0M':1,'1M':2,'3M':3,'6M':4,'12M':5}
sm=pd.read_csv('sample_summary.tsv',sep='\t',index_col=0)
d=json.load(open('multiqc_general_stats.json'))['general_stats_table']['datasets'][0]['violin_value_by_sample_by_metric']
qc=pd.DataFrame({k.replace('fastp-',''):pd.Series(v) for k,v in d.items() if k.startswith('fastp')}).loc[S]
host=pd.DataFrame({h:{s:d['bowtie_2_hisat2-overall_alignment_rate'].get(f'{s}.{h}') for s in S} for h in ['glycine_max','homo_sapiens','hordeum_vulgare','oryza_sativa','triticum_aestivum']})
print('== QC'); print(qc.describe().round(2).T[['min','50%','max','mean']])
print('pairs(M) = reads/2:', (qc['filtering_result_passed_filter_reads']/2).describe().round(1).to_dict())
print('host soybean %',host.glycine_max.describe().round(2).to_dict()); print('human',host.homo_sapiens.describe().round(2).to_dict())
fp=pd.read_csv('fungi_2step_summary_parsed.tsv',sep='\t',index_col=0).loc[S]
mag=pd.read_csv('mag_relative_abundance.tsv',sep='\t',index_col=0)[S]
cls=pd.DataFrame({'bact':sm.loc[S,'bacteria_pct'],'fungi':fp.fungi_pct,'uncl1':sm.loc[S,'unclassified_pct'],'MAGmapped':100-mag.loc['unmapped'],'kropp':mag.loc['0M_5_p.bin.10']})
cls['kraken_total']=cls.bact+cls.fungi; cls['residual']=100-cls.bact-cls.fungi-cls.kropp
print('== classification'); print(cls.describe().round(1).T[['min','50%','max','mean']])
print(cls.join(m).groupby('stage')[['bact','fungi','uncl1','kropp','residual']].median().reindex(O).round(1))
# alpha
print('== alpha'); a=sm.loc[S,['shannon','observed_species_ge0.01pct']].join(m)
for c in ['shannon','observed_species_ge0.01pct']:
    print(c, a.groupby('stage')[c].median().reindex(O).round(2).to_dict(), a.groupby('region')[c].median().round(2).to_dict(),
      'KW stage p=%.3g region p=%.3g'%(kruskal(*[g[c] for _,g in a.groupby('stage')]).pvalue,kruskal(*[g[c] for _,g in a.groupby('region')]).pvalue))
# composition % of total reads
bg=pd.read_csv('bracken_genus_relab.tsv',sep='\t',index_col=0)[S]*sm.loc[S,'bacteria_pct']/100
fg=pd.read_csv('fungi_bracken_G_pct_of_total.tsv',sep='\t',index_col=0)[S]
tot=pd.concat([bg,fg.rename(index=lambda x:x+' (F)')]); tot.loc['Kroppenstedtia (MAG)']=mag.loc['0M_5_p.bin.10']
tot.to_csv('composition_pct_total_reads.tsv',sep='\t',float_format='%.4g')
key=['Bacillus','Enterococcus','Pediococcus','Tetragenococcus','Heyndrickxia','Weissella','Pseudomonas','Kroppenstedtia (MAG)','Aspergillus (F)','Penicillium (F)','Rhizomucor (F)','Pichia (F)']
T3=tot.loc[key].T.join(m).groupby('stage')[key].mean().reindex(O).T.round(1); print('== stage mean % total'); print(T3)
T3.to_csv('table_stage_mean.tsv',sep='\t')
R2=tot.loc[key].T.join(m).groupby('region')[key].mean().T.round(1); print(R2); R2.to_csv('table_region_mean.tsv',sep='\t')
# stage trend within region (rm-rank correlation with ordinal stage)
reg=m.loc[S,'region'].values; so=np.array([order[x] for x in m.loc[S,'stage']])
def rmr(v): r=pd.Series(rankdata(v),index=S); return (r-r.groupby(reg).transform('mean')).values
out=[]
for k in key:
    x=rmr(tot.loc[k,S].values); y=rmr(so); r=np.corrcoef(x,y)[0,1]; df=len(S)-6-1; tt=r*np.sqrt(df/(1-r*r)); out.append((k,round(r,2),2*T.sf(abs(tt),df)))
print('== within-region stage trend'); print(pd.DataFrame(out,columns=['taxon','r','p']).round(4))
pd.DataFrame(out,columns=['taxon','r','p']).to_csv('stage_trend.tsv',sep='\t',index=False)
