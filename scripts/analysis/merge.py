import pandas as pd, numpy as np, os, re
R=os.environ.get('RESULTS_DIR','/mnt/f/meta_hong/results')  # per-sample pipeline output folder
S=['Me_%d_p'%i for i in range(1,7)]+['0M_%d_p'%i for i in range(1,7)]+['1M_%d_p'%i for i in range(1,7)]+['3M_%d_p'%i for i in range(1,5)]+['6M_%d'%i for i in range(1,6)]+['12M_%d_p'%i for i in [1,2,3,5,6]]
tp={'Me':'Meju','0M':'0M','1M':'1M','3M':'3M','6M':'6M','12M':'12M'}
meta=pd.DataFrame([{'sample':s,'stage':s.split('_')[0],'region':int(s.split('_')[1])} for s in S]).set_index('sample')
# bracken species
sp={}; 
for s in S:
    b=pd.read_csv(f'{R}/{s}/05_taxonomy_reads_kraken/{s}.bracken_S.txt',sep='\t')
    sp[s]=b.set_index('name')['new_est_reads']
sp=pd.DataFrame(sp).fillna(0)
rel=sp/sp.sum()*100
gen=rel.groupby(rel.index.str.split().str[0]).sum()
# kraken report summary
rows=[]
for s in S:
    k=pd.read_csv(f'{R}/{s}/05_taxonomy_reads_kraken/{s}.kraken2.report',sep='\t',header=None,names=['pct','clade','direct','rank','taxid','name'])
    k['name']=k['name'].str.strip()
    g=lambda t: k.loc[k.taxid==t,'clade'].sum()
    tot=k.loc[k['rank'].isin(['U','R']),'clade'].sum()
    rows.append({'sample':s,'total_reads_kraken':tot,'unclassified_pct':100*g(0)/tot,'bacteria_pct':100*g(2)/tot,'fungi_pct':100*g(4751)/tot,'eukaryota_pct':100*g(2759)/tot,'viruses_pct':100*g(10239)/tot,'bracken_species_n':(sp[s]>0).sum()})
ks=pd.DataFrame(rows).set_index('sample')
# alpha diversity (species, bracken)
def shannon(x): p=x[x>0]/x.sum(); return -(p*np.log(p)).sum()
def simpson(x): p=x/x.sum(); return 1-(p**2).sum()
ks['shannon']=[shannon(rel[s]) for s in S]; ks['simpson']=[simpson(rel[s]) for s in S]
ks['observed_species_ge0.01pct']=[(rel[s]>=0.01).sum() for s in S]
# humann
pa={}
for s in S:
    h=pd.read_csv(f'{R}/{s}/11_function_community_humann/humann_output/{s}_pathabundance.tsv',sep='\t',index_col=0)
    pa[s]=h.iloc[:,0]
pa=pd.DataFrame(pa).fillna(0)
unstr=pa[~pa.index.str.contains(r'\|')]
ks['humann_unmapped_pct']=100*unstr.loc['UNMAPPED']/unstr.sum()
ks['humann_unintegrated_pct']=100*unstr.loc['UNINTEGRATED']/unstr.sum()
path=unstr.drop(['UNMAPPED','UNINTEGRATED'])
path_cpm=path/path.sum()*1e6
strat=pa[pa.index.str.contains(r'\|') & ~pa.index.str.startswith('UNINTEGRATED')]
ks['humann_pathways_n']=[(path[s]>0).sum() for s in S]
meta.join(ks).to_csv('sample_summary.tsv',sep='\t',float_format='%.4g')
rel.loc[rel.mean(1).sort_values(ascending=False).index].to_csv('bracken_species_relab.tsv',sep='\t',float_format='%.5g')
sp.to_csv('bracken_species_counts.tsv',sep='\t')
gen.loc[gen.mean(1).sort_values(ascending=False).index].to_csv('bracken_genus_relab.tsv',sep='\t',float_format='%.5g')
path_cpm.loc[path_cpm.mean(1).sort_values(ascending=False).index].to_csv('humann_pathabundance_cpm_unstratified.tsv',sep='\t',float_format='%.5g')
strat.to_csv('humann_pathabundance_stratified_raw.tsv',sep='\t',float_format='%.5g')
meta.to_csv('metadata.tsv',sep='\t')
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(meta.join(ks).round(2).to_string())
print(gen.mean(1).sort_values(ascending=False).head(15).round(2))
print(path.shape, strat.shape)
