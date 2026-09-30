# Alpha diversity (Shannon, richness) of the combined community (bacterial MAGs + fungal genera)
# and within/between-region Bray-Curtis dissimilarity; writes alpha_whole.tsv
import pandas as pd,numpy as np
from scipy.stats import rankdata
from scipy.spatial.distance import pdist,squareform
rng=np.random.default_rng(21)
m=pd.read_csv('metadata.tsv',sep='\t',index_col=0); S=list(m.index); reg=m.region.values; st=m.stage.values
def pw(x):
    y=x.copy()
    for b in np.unique(reg):
        i=np.where(reg==b)[0]; y[i]=x[rng.permutation(i)]
    return y
def design(fs,n):
    cols=[np.ones(n)]
    for f in fs:
        for l in np.unique(f)[1:]: cols.append((f==l).astype(float))
    return np.column_stack(cols)
n=len(S); Xr=design([reg],n); Hr=Xr@np.linalg.pinv(Xr)
def Fstage(y,stl):
    X=design([reg,stl],n); H=X@np.linalg.pinv(X)
    ssf=((np.eye(n)-H)@y)@((np.eye(n)-H)@y); ssr=((np.eye(n)-Hr)@y)@((np.eye(n)-Hr)@y)
    return ((ssr-ssf)/5)/(ssf/(n-np.linalg.matrix_rank(X)))
t=pd.read_csv('combined_bacteriaMAG_fungiGenus_relab.tsv',sep='\t',index_col=0)[S]
p=t/t.sum(); H=-(p*np.log(p.where(p>0))).sum(); obs=(p>1e-4).sum()
out=pd.DataFrame({'shannon_whole':H,'richness_whole':obs}).join(m)
O=['Me','0M','1M','3M','6M','12M']
for c in ['shannon_whole','richness_whole']:
    y=rankdata(out.loc[S,c].values); F=Fstage(y,st); NP=9999
    cnt=sum(Fstage(y,pw(st))>=F for _ in range(NP)); print(c,'stage (region-adjusted ranks, within-region perm) F=%.2f p=%.4f'%(F,(cnt+1)/(NP+1)))
    print('  median by stage:',out.groupby('stage')[c].median().reindex(O).round(2).to_dict())
    print('  median by region:',out.groupby('region')[c].median().round(2).to_dict())
D=squareform(pdist(t.T.values,'braycurtis')); iu=np.triu_indices(n,1); same=reg[iu[0]]==reg[iu[1]]
print('Bray-Curtis within region median %.2f (IQR %.2f-%.2f); between regions %.2f (%.2f-%.2f)'%(np.median(D[iu][same]),*np.percentile(D[iu][same],[25,75]),np.median(D[iu][~same]),*np.percentile(D[iu][~same],[25,75])))
out.to_csv('alpha_whole.tsv',sep='\t')
