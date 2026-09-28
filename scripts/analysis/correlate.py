import pandas as pd,numpy as np
from scipy.stats import spearmanr, rankdata
def bh(p):
    p=np.asarray(p); n=len(p); o=np.argsort(p); q=np.empty(n); prev=1
    for i in range(n-1,-1,-1):
        prev=min(prev,p[o[i]]*n/(i+1)); q[o[i]]=prev
    return q
d=pd.read_excel('HMT.xlsx',header=0); d['sample']=d['sample'].astype(str).str.strip()
conv=lambda s:(lambda p:f"{p[0].upper().replace('ME','Me')}_{p[1]}"+("" if p[0].upper()=='6M' else "_p"))(s.replace('-','_').split('_'))
d.index=d['sample'].map(conv); met=d.drop(columns=['숙성기간','sample']).apply(pd.to_numeric,errors='coerce')
m=pd.read_csv('metadata.tsv',sep='\t',index_col=0)
S=[s for s in m.index if s in met.index]; met=met.loc[S]
met=met.loc[:,(met>0).sum()>=int(0.5*len(S))]           # 절반 이상 검출된 대사체만
sm=pd.read_csv('sample_summary.tsv',sep='\t',index_col=0)
bg=pd.read_csv('bracken_genus_relab.tsv',sep='\t',index_col=0)[S]*sm.loc[S,'bacteria_pct']/100
fg=pd.read_csv('fungi_bracken_G_pct_of_total.tsv',sep='\t',index_col=0)[S]
mag=pd.read_csv('mag_relative_abundance.tsv',sep='\t',index_col=0)[S]
tax=pd.concat([bg.rename(index=lambda x:'B:'+x), fg.rename(index=lambda x:'F:'+x)])
tax.loc['B:Kroppenstedtia (MAG)']=mag.loc['0M_5_p.bin.10']
tax=tax[(tax.mean(1)>=1)|(tax.index=='B:Kroppenstedtia (MAG)')]
tax=tax.loc[tax.mean(1).sort_values(ascending=False).index]
print('taxa:',list(tax.index)); print('metabolites:',met.shape[1],' samples:',len(S))
reg=m.loc[S,'region'].values
def rm_rank(v):  # 지역 내 순위 중심화 (지역 효과 제거)
    r=pd.Series(rankdata(v),index=S); return (r-r.groupby(reg).transform('mean')).values
rows=[]
for t in tax.index:
    for c in met.columns:
        x=tax.loc[t,S].values; y=met[c].values
        rho,p=spearmanr(x,y)
        r2=np.corrcoef(rm_rank(x),rm_rank(y))[0,1]
        # rm 상관 p: t-분포, df = n - 1 - (지역수-1) - 1
        from scipy.stats import t as T
        df=len(S)-len(set(reg))-1; tt=r2*np.sqrt(df/(1-r2**2)) if abs(r2)<1 else np.inf
        p2=2*T.sf(abs(tt),df)
        rows.append((t,c,rho,p,r2,p2))
R=pd.DataFrame(rows,columns=['taxon','metabolite','spearman_rho','p','within_region_r','p_within']).dropna()
R['q']=bh(R.p.values); R['q_within']=bh(R.p_within.values)
R.to_csv('taxa_metabolite_correlation.tsv',sep='\t',index=False,float_format='%.4g')
print('pairs:',len(R),' Spearman q<0.05:',(R.q<0.05).sum(),' within-region q<0.05:',(R.q_within<0.05).sum(),' both:',((R.q<0.05)&(R.q_within<0.05)).sum())
key=['Tyr','Phe','Leu','Lys','Putrescine','Tyramine','Histamine','Cadaverine','Glu','GABA','Ile','Val']
key=[k for k in key if k in met.columns]; print('key metabolites present:',key)
pd.set_option('display.width',220)
K=R[R.metabolite.isin(key)&((R.q<0.05)|(R.q_within<0.05))].sort_values('q')
print(K.round(3).head(40).to_string(index=False))
top=R[(R.q<0.05)&(R.q_within<0.05)].sort_values('p_within')
print(top.round(3).head(25).to_string(index=False))
