import pandas as pd,numpy as np
from scipy.spatial.distance import pdist,squareform
from scipy.stats import kruskal, spearmanr
rel=pd.read_csv('bracken_species_relab.tsv',sep='\t',index_col=0); meta=pd.read_csv('metadata.tsv',sep='\t',index_col=0)
path=pd.read_csv('humann_pathabundance_cpm_unstratified.tsv',sep='\t',index_col=0)
rng=np.random.default_rng(1)
def ss_w(D2,g):
    s=0
    for k in np.unique(g):
        i=np.where(g==k)[0]; s+=D2[np.ix_(i,i)].sum()/2/len(i)
    return s
def perm_seq(D,factors,nperm=9999):
    # sequential SS (adonis type I) : factors list of arrays; returns R2,F,p per term
    D2=D**2; n=len(D); A=-0.5*D2; J=np.eye(n)-1/n; G=J@A@J
    def hat(X): return X@np.linalg.pinv(X)
    def design(fs):
        cols=[np.ones(n)]
        for f in fs:
            lv=np.unique(f)
            for l in lv[1:]: cols.append((f==l).astype(float))
        return np.column_stack(cols)
    def stats(fs):
        out=[];prev=hat(design([]));SSt=np.trace(G)
        Xf=design(fs);Hf=hat(Xf);SSres=np.trace((np.eye(n)-Hf)@G@(np.eye(n)-Hf));dfres=n-np.linalg.matrix_rank(Xf)
        for k in range(len(fs)):
            H=hat(design(fs[:k+1])); ss=np.trace(G@(H-prev)); df=np.linalg.matrix_rank(design(fs[:k+1]))-np.linalg.matrix_rank(design(fs[:k]))
            out.append((ss/SSt,(ss/df)/(SSres/dfres))); prev=H
        return out
    obs=stats(factors); cnt=np.zeros(len(factors))
    for _ in range(nperm):
        p=rng.permutation(n); fs=[f[p] for f in factors]
        st=stats(fs); cnt+=[st[k][1]>=obs[k][1] for k in range(len(factors))]
    return [(r2,F,(c+1)/(nperm+1)) for (r2,F),c in zip(obs,cnt)]
def run(tab,label):
    X=tab[meta.index].T.values
    D=squareform(pdist(X,'braycurtis'))
    r=meta.region.astype(str).values; s=meta.stage.values
    for name,fs in [('region then stage',[r,s]),('stage then region',[s,r])]:
        res=perm_seq(D,fs,1999)
        print(label,name,[(round(a,3),round(b,2),round(c,4)) for a,b,c in res])
run(rel,'Species BC')
run(path,'Pathway BC')
sm=pd.read_csv('sample_summary.tsv',sep='\t',index_col=0)
for m in ['shannon','observed_species_ge0.01pct','unclassified_pct']:
    print(m,'KW stage p=%.3g'%kruskal(*[g[m].values for _,g in sm.groupby('stage')]).pvalue,'KW region p=%.3g'%kruskal(*[g[m].values for _,g in sm.groupby('region')]).pvalue)
    print(sm.groupby('stage')[m].median().reindex(['Me','0M','1M','3M','6M','12M']).round(2).to_dict())
