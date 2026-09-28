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
sm=pd.read_csv('sample_summary.tsv',sep='\t',index_col=0)
for c in ['shannon','observed_species_ge0.01pct']:
    y=rankdata(sm.loc[S,c].values); F=Fstage(y,st); NP=9999
    cnt=sum(Fstage(y,pw(st))>=F for _ in range(NP)); print(c,'stage (region-adjusted ranks, within-region perm) F=%.2f p=%.4f'%(F,(cnt+1)/(NP+1)))
# PERMDISP stage restricted
def G_of(D): J=np.eye(len(D))-1/len(D); return J@(-0.5*D**2)@J
def dists(D,groups):
    G=G_of(D); w,v=np.linalg.eigh(G); pos=w>1e-10; neg=w<-1e-10
    Xp=v[:,pos]*np.sqrt(w[pos]); Xn=v[:,neg]*np.sqrt(-w[neg]); z=np.zeros(len(D))
    for g in np.unique(groups):
        i=groups==g; dp=((Xp[i]-Xp[i].mean(0))**2).sum(1); dn=((Xn[i]-Xn[i].mean(0))**2).sum(1) if neg.any() else 0
        z[i]=np.sqrt(np.clip(dp-dn,0,None))
    return z
def Fg(z,g):
    gm=z.mean(); L=np.unique(g); ssb=sum(len(z[g==l])*(z[g==l].mean()-gm)**2 for l in L); ssw=sum(((z[g==l]-z[g==l].mean())**2).sum() for l in L)
    return (ssb/(len(L)-1))/(ssw/(len(z)-len(L)))
tabs={'Bacteria species (Kraken2)':pd.read_csv('bracken_species_relab.tsv',sep='\t',index_col=0)[S]}
fs=pd.read_csv('fungi_bracken_S_relab_within_fungi.tsv',sep='\t',index_col=0)[S]; fs.loc['Aspergillus flavus/oryzae']=fs.loc[['Aspergillus oryzae','Aspergillus flavus']].sum(); tabs['Fungi species']=fs.drop(['Aspergillus oryzae','Aspergillus flavus'])
tabs['Fungi genus']=pd.read_csv('fungi_bracken_G_relab_within_fungi.tsv',sep='\t',index_col=0)[S]
tabs['MAG 63']=pd.read_csv('mag_relab_within_mapped.tsv',sep='\t',index_col=0)[S]
tabs['Bacteria MAG + fungi genus']=pd.read_csv('combined_bacteriaMAG_fungiGenus_relab.tsv',sep='\t',index_col=0)[S]
tabs['MetaCyc pathways']=pd.read_csv('humann_pathabundance_cpm_unstratified.tsv',sep='\t',index_col=0)[S]
out={}
for k,t in tabs.items():
    D=squareform(pdist(t.T.values,'braycurtis'))
    z=dists(D,st); F=Fg(z,st); NP=4999; cnt=0
    for _ in range(NP):
        sp=pw(st); cnt+=Fg(dists(D,sp),sp)>=F
    zr=dists(D,reg); out[k]=((cnt+1)/(NP+1), pd.Series(zr).groupby(reg).mean().round(2).to_dict())
    print(k,'PERMDISP stage restricted p=%.4f'%out[k][0],'region mean dist-to-centroid',out[k][1])
rv=pd.read_csv('review_permanova.tsv',sep='\t'); rv['PERMDISP_p_stage_restricted']=rv.data.map(lambda k:out[k][0])
rv['disp_region']=rv.data.map(lambda k:out[k][1]); rv.to_csv('review_permanova.tsv',sep='\t',index=False)
