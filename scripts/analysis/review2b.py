import pandas as pd,numpy as np
from scipy.stats import rankdata
exec(open(__file__.replace('review2b.py','correlate.py')).read().split("rows=[]")[0].replace("print(","(lambda *a,**k:None)("))
rng=np.random.default_rng(5); NP=9999
def bh(p):
    p=np.asarray(p);n=len(p);o=np.argsort(p);q=np.empty(n);prev=1
    for i in range(n-1,-1,-1): prev=min(prev,p[o[i]]*n/(i+1)); q[o[i]]=prev
    return q
regs=np.unique(reg); idx=[np.where(reg==b)[0] for b in regs]
def rc(M):  # rows=variables, cols=samples: rank within row, center within region
    R=np.apply_along_axis(rankdata,1,M).astype(float)
    for i in idx: R[:,i]-=R[:,i].mean(1,keepdims=True)
    return R/np.linalg.norm(R,axis=1,keepdims=True)
TX=rc(tax[S].values); MT=rc(met.loc[S].T.values)
obs=TX@MT.T   # 12x46
cnt=np.zeros_like(obs)
for _ in range(NP):
    perm=np.arange(len(S))
    for i in idx: perm[i]=i[rng.permutation(len(i))]
    cnt+=np.abs(rc(tax[S].values[:,perm])@MT.T)>=np.abs(obs)-1e-12
P=(cnt+1)/(NP+1)
R=pd.read_csv('taxa_metabolite_correlation.tsv',sep='\t')
pm=pd.DataFrame(P,index=tax.index,columns=met.columns); om=pd.DataFrame(obs,index=tax.index,columns=met.columns)
R['within_r_check']=[om.loc[a,b] for a,b in zip(R.taxon,R.metabolite)]
R['p_within_perm']=[pm.loc[a,b] for a,b in zip(R.taxon,R.metabolite)]
R['q_within_perm']=bh(R.p_within_perm.values)
print('max diff within r',np.abs(R.within_r_check-R.within_region_r).max())
print('within-region perm q<0.05:',(R.q_within_perm<0.05).sum(),' old t-based:',(R.q_within<0.05).sum(),' both simple q & perm q:',((R.q<0.05)&(R.q_within_perm<0.05)).sum())
R.drop(columns='within_r_check').to_csv('taxa_metabolite_correlation.tsv',sep='\t',index=False,float_format='%.4g')
T7=pd.read_csv('table5_sel.tsv',sep='\t').drop(columns=[c for c in ['p_within_perm','q_within_perm'] if c in pd.read_csv('table5_sel.tsv',sep='\t').columns])
T7=T7.merge(R[['taxon','metabolite','p_within_perm','q_within_perm']],on=['taxon','metabolite']); T7.to_csv('table5_sel.tsv',sep='\t',index=False)
print(T7[['taxon','metabolite','within_region_r','q_within','p_within_perm','q_within_perm']].round(4).to_string())
