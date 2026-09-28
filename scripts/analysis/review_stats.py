import pandas as pd,numpy as np
from scipy.spatial.distance import pdist,squareform
from scipy.stats import rankdata
rng=np.random.default_rng(7); NP=4999
m=pd.read_csv('metadata.tsv',sep='\t',index_col=0); S=list(m.index)
reg=m.region.astype(str).values; st=m.stage.values
sm=pd.read_csv('sample_summary.tsv',sep='\t',index_col=0)
def G_of(D): n=len(D); J=np.eye(n)-1/n; return J@(-0.5*D**2)@J
def design(fs,n):
    cols=[np.ones(n)]
    for f in fs:
        if f.dtype.kind in 'fc': cols.append(f.astype(float)); continue
        for l in np.unique(f)[1:]: cols.append((f==l).astype(float))
    return np.column_stack(cols)
hat=lambda X: X@np.linalg.pinv(X)
def ss_term(G,full,red,n):
    Hf=hat(design(full,n)); Hr=hat(design(red,n)); return np.trace(G@(Hf-Hr)), np.trace((np.eye(n)-Hf)@G@(np.eye(n)-Hf))
def perm_within(labels,blocks):
    out=labels.copy()
    for b in np.unique(blocks):
        i=np.where(blocks==b)[0]; out[i]=labels[rng.permutation(i)]
    return out
def marginal(D,cov=None):
    n=len(D); G=G_of(D); SSt=np.trace(G); base=[] if cov is None else [cov]
    res={}
    for name,term,other,block in [('region',reg,st,st),('stage',st,reg,reg)]:
        full=base+[reg,st]; red=base+[other]
        ss,sr=ss_term(G,full,red,n); dft=len(np.unique(term))-1; dfr=n-np.linalg.matrix_rank(design(full,n))
        F=(ss/dft)/(sr/dfr); cnt=0
        for _ in range(NP):
            tp=perm_within(term,block)   # 다른 요인 블록 내 제한 순열
            f2=[tp if x is term else x for x in full]
            s2,r2=ss_term(G,base+([tp,st] if name=='region' else [reg,tp]),red,n)
            cnt+=((s2/dft)/(r2/dfr))>=F
        res[name]=(ss/SSt,F,(cnt+1)/(NP+1))
    return res
def permdisp(D,groups):
    n=len(D); G=G_of(D); w,v=np.linalg.eigh(G); pos=w>1e-10; neg=w<-1e-10
    Xp=v[:,pos]*np.sqrt(w[pos]); Xn=v[:,neg]*np.sqrt(-w[neg])
    z=np.zeros(n)
    for g in np.unique(groups):
        i=groups==g; cp=Xp[i].mean(0); cn=Xn[i].mean(0) if neg.any() else 0
        dp=((Xp[i]-cp)**2).sum(1); dn=((Xn[i]-cn)**2).sum(1) if neg.any() else 0
        z[i]=np.sqrt(np.clip(dp-dn,0,None))
    def Fstat(z,g):
        gm=z.mean(); k=len(np.unique(g)); ssb=sum(len(z[g==l])*(z[g==l].mean()-gm)**2 for l in np.unique(g)); ssw=sum(((z[g==l]-z[g==l].mean())**2).sum() for l in np.unique(g))
        return (ssb/(k-1))/(ssw/(len(z)-k))
    F=Fstat(z,groups); cnt=sum(Fstat(rng.permutation(z),groups)>=F for _ in range(NP))
    return F,(cnt+1)/(NP+1), pd.Series(z).groupby(groups).mean().round(3).to_dict()
data={}
rel=pd.read_csv('bracken_species_relab.tsv',sep='\t',index_col=0)[S]; data['Bacteria species (Kraken2)']=rel
fs=pd.read_csv('fungi_bracken_S_relab_within_fungi.tsv',sep='\t',index_col=0)[S]
fs.loc['Aspergillus flavus/oryzae']=fs.loc[['Aspergillus oryzae','Aspergillus flavus']].sum(); fs=fs.drop(['Aspergillus oryzae','Aspergillus flavus']); data['Fungi species']=fs
data['Fungi genus']=pd.read_csv('fungi_bracken_G_relab_within_fungi.tsv',sep='\t',index_col=0)[S]
data['MAG 63']=pd.read_csv('mag_relab_within_mapped.tsv',sep='\t',index_col=0)[S]
data['Bacteria MAG + fungi genus']=pd.read_csv('combined_bacteriaMAG_fungiGenus_relab.tsv',sep='\t',index_col=0)[S]
data['MetaCyc pathways']=pd.read_csv('humann_pathabundance_cpm_unstratified.tsv',sep='\t',index_col=0)[S]
rows=[]
for k,tab in data.items():
    D=squareform(pdist(tab.T.values,'braycurtis'))
    r=marginal(D); pr=permdisp(D,reg); ps=permdisp(D,st)
    rows.append([k,*[round(x,4) for x in r['region']],*[round(x,4) for x in r['stage']],round(pr[0],2),round(pr[1],4),round(ps[0],2),round(ps[1],4)])
    print(k,rows[-1][1:],'disp region',pr[2])
T=pd.DataFrame(rows,columns=['data','R2_region_marg','F_region','p_region_restricted','R2_stage_marg','F_stage','p_stage_restricted','PERMDISP_F_region','PERMDISP_p_region','PERMDISP_F_stage','PERMDISP_p_stage'])
T.to_csv('review_permanova.tsv',sep='\t',index=False); print(T.to_string())
# HUMAnN sensitivity: unmapped% as covariate
D=squareform(pdist(data['MetaCyc pathways'].T.values,'braycurtis'))
un=sm.loc[S,'humann_unmapped_pct'].values
r=marginal(D,cov=un); print('Pathway with unmapped covariate:',{k:tuple(round(x,4) for x in v) for k,v in r.items()})
print('HUMAnN unmapped pct: min %.1f median %.1f max %.1f'%(un.min(),np.median(un),un.max()))
print(pd.Series(un,index=S).groupby(m.stage).median().round(1).to_dict(), pd.Series(un,index=S).groupby(m.region).median().round(1).to_dict())
keep=un<30; D2=squareform(pdist(data['MetaCyc pathways'].T.values[keep],'braycurtis'))
print('n kept(<30%)',keep.sum())
