import os; os.makedirs("fig", exist_ok=True)
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial.distance import pdist,squareform
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.edgecolor':'#8a8985','axes.linewidth':0.6,
  'axes.labelcolor':'#0b0b0b','xtick.color':'#52514e','ytick.color':'#52514e','axes.spines.top':False,'axes.spines.right':False,
  'savefig.dpi':300,'figure.facecolor':'white','legend.frameon':False})

import re as _re, matplotlib.figure as _mf, matplotlib.text as _mt
_GEN='Aspergillus|Bacillus|Brochothrix|Chromohalobacter|Enterococcus|Heyndrickxia|Kroppenstedtia|Lichtheimia|Mucor|Pediococcus|Penicillium|Pichia|Pseudomonas|Rhizomucor|Rhizopus|Saccharomyces|Staphylococcus|Tetragenococcus|Thermoascus|Wallemia|Weissella|Weizmannia|Leuconostoc|Lactobacillus|Caldifermentibacillus|Lactiplantibacillus'
_EPI='oryzae|flavus|licheniformis|sonorensis|thermosphacta|canadensis|faecalis|faecium|halophilus|confusa|coagulans|acidilactici|kudriavzevii|solitum|lactis|ichthyophaga|eburnea|subtilis|velezensis|pentosaceus'
_TRE=_re.compile(r'(?<![A-Za-z$])((?:'+_GEN+r')(?: (?:'+_EPI+r'))?|[A-Z]\. (?:'+_EPI+r'))(?![A-Za-z])')
def _it(m): return ' '.join(r'$\mathit{'+w+'}$' for w in m.group(1).split(' '))
_orig=_mf.Figure.savefig
def _save(self,*a,**k):
    for ax in self.axes:
        for get,setl in ((ax.get_yticklabels,ax.set_yticklabels),(ax.get_xticklabels,ax.set_xticklabels)):
            L=[t.get_text() for t in get()]
            if any(l and '$' not in l and _TRE.search(l) for l in L):
                setl([_TRE.sub(_it,l) if '$' not in l else l for l in L])
    for t in self.findobj(_mt.Text):
        s=t.get_text()
        if s and '$' not in s and _TRE.search(s): t.set_text(_TRE.sub(_it,s))
    return _orig(self,*a,**k)
_mf.Figure.savefig=_save
C=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#4a3aa7','#e34948']
GRAY='#b9b8b3'; LG='#e4e3df'; INK2='#52514e'
m=pd.read_csv('metadata.tsv',sep='\t',index_col=0); O=['Me','0M','1M','3M','6M','12M']; OL=['Meju','0M','1M','3M','6M','12M']
so={s:i for i,s in enumerate(O)}
S=sorted(m.index,key=lambda s:(m.region[s],so[m.stage[s]])); m=m.loc[S]
MK=['o','s','^','D','v','P']
RN={1:'Gangwon',2:'Gyeonggi',3:'Gyeongsang',4:'Jeolla',5:'Chungcheong',6:'Jeju'}
# ---- Fig1 PCoA
rel=pd.read_csv('bracken_species_relab.tsv',sep='\t',index_col=0)[S]
D=squareform(pdist(rel.T.values,'braycurtis')); n=len(D); J=np.eye(n)-1/n; B=-0.5*J@(D**2)@J
w,v=np.linalg.eigh(B); i=np.argsort(w)[::-1]; w,v=w[i],v[:,i]; X=v[:,:2]*np.sqrt(w[:2]); ev=w[:2]/w[w>0].sum()*100
fig,ax=plt.subplots(figsize=(4.6,3.8))
for k,r in enumerate(range(1,7)):
    idx=[j for j,s in enumerate(S) if m.region[s]==r]
    ax.plot(X[idx,0],X[idx,1],'-',color=C[k],lw=1,alpha=.6,zorder=1)
    ax.scatter(X[idx,0],X[idx,1],s=34,marker=MK[k],color=C[k],edgecolor='white',linewidth=.8,zorder=2,label=f'Region {r} ({RN[int(r)]})')
    j0=idx[0]; ax.annotate('Meju',(X[j0,0],X[j0,1]),xytext=(4,4),textcoords='offset points',fontsize=6,color=INK2)
ax.set_xlabel(f'PCoA1 ({ev[0]:.1f}%)'); ax.set_ylabel(f'PCoA2 ({ev[1]:.1f}%)')
ax.axhline(0,color=LG,lw=.6,zorder=0); ax.axvline(0,color=LG,lw=.6,zorder=0)
ax.legend(loc='center left',bbox_to_anchor=(1,.5),fontsize=7,handletextpad=.3)
ax.set_title('Bacterial community (species, Bray–Curtis)',fontsize=8.5,loc='left')
fig.tight_layout(); fig.savefig('fig/Fig1_PCoA.png'); plt.close()
print('PCoA ev',ev.round(1))
# ---- Fig2 stacked composition (% of total reads)
tot=pd.read_csv('composition_pct_total_reads.tsv',sep='\t',index_col=0)[S]
fung=tot[tot.index.str.endswith('(F)')].sum()
bact=tot[~tot.index.str.endswith('(F)')]
keys=['Bacillus','Enterococcus','Pediococcus','Tetragenococcus','Kroppenstedtia (MAG)']
other=bact.drop(keys).sum()
comp=pd.DataFrame({k:bact.loc[k] for k in keys}); comp['Other bacteria']=other; comp['Fungi']=fung
sc=comp.sum(1).clip(lower=100)/100; comp=comp.div(sc,axis=0); comp['Unclassified']=(100-comp.sum(1)).clip(lower=0)
cols=[C[0],C[1],C[2],C[3],C[6],GRAY,C[4],LG]
fig,ax=plt.subplots(figsize=(7,3.1)); x=[]; pos=0; ticks=[]
for r in range(1,7):
    ss=[s for s in S if m.region[s]==r]
    for s in ss: x.append(pos); ticks.append(OL[so[m.stage[s]]]); pos+=1
    pos+=0.8
x=np.array(x); bottom=np.zeros(len(S))
for c,col in zip(comp.columns,cols):
    ax.bar(x,comp[c].values,bottom=bottom,width=.82,color=col,edgecolor='white',linewidth=.5,label=c.replace(' (MAG)','*').replace('Unclassified','Unexplained (approx.)')); bottom+=comp[c].values
ax.set_xticks(x); ax.set_xticklabels(ticks,rotation=90,fontsize=6); ax.set_ylabel('% of total reads'); ax.set_ylim(0,100)
xs=0
for r in range(1,7):
    k=sum(m.region==r); ax.text(x[xs:xs+k].mean(),103,f'Region {r}\n({RN[int(r)]})',ha='center',fontsize=6.5,color='#0b0b0b'); xs+=k
ax.legend(ncol=8,loc='upper center',bbox_to_anchor=(.5,-.22),fontsize=6.5,handlelength=1,columnspacing=.8)
fig.tight_layout(); fig.savefig('fig/Fig2_composition.png'); plt.close()
comp.to_csv('fig2_data.tsv',sep='\t',float_format='%.2f')
# ---- Fig3 succession
kk=['Bacillus','Enterococcus','Tetragenococcus','Pediococcus','Kroppenstedtia (MAG)','Aspergillus (F)']
cc=[C[0],C[1],C[3],C[2],C[6],C[4]]
fig,ax=plt.subplots(figsize=(4.8,3.1))
for k,col in zip(kk,cc):
    g=tot.loc[k].groupby(m.stage); mu=g.mean().reindex(O); se=(g.std()/np.sqrt(g.count())).reindex(O)
    ax.plot(range(6),mu.values,'-o',color=col,lw=2,ms=4,mec='white',mew=.8)
    ax.fill_between(range(6),(mu-se).values,(mu+se).values,color=col,alpha=.12,lw=0)
    ax.text(5.15,mu.values[-1],k.replace(' (MAG)','*').replace(' (F)',' (fungi)'),color='#0b0b0b',fontsize=6.5,va='center')
ax.set_xticks(range(6)); ax.set_xticklabels(OL); ax.set_xlim(-.2,6.6); ax.set_ylabel('% of total reads (mean ± SE)')
ax.grid(axis='y',color=LG,lw=.5); ax.set_axisbelow(True)
ax.axvline(2.5,color=INK2,lw=.8,ls='--',zorder=1)
ax.text(2.5,1.01,'Soy sauce separation (day 60)',transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=6.5,color=INK2)
ax.text(1.0,.97,'Meju / tojang',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=6.5,color=INK2)
ax.text(4.0,.97,'Doenjang',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=6.5,color=INK2)
fig.tight_layout(); fig.savefig('fig/Fig4_succession.png'); plt.close()
# ---- Fig4 fungi
fp=pd.read_csv('fungi_2step_summary_parsed.tsv',sep='\t',index_col=0).loc[S]
fr=pd.read_csv('fungi_bracken_G_relab_within_fungi.tsv',sep='\t',index_col=0)[S]
fig,(a1,a2)=plt.subplots(1,2,figsize=(7,2.9),gridspec_kw={'width_ratios':[1,1.25]})
for k,r in enumerate(range(1,7)):
    ss=[s for s in S if m.region[s]==r]
    a1.scatter([so[m.stage[s]]+(k-2.5)*.07 for s in ss],fp.loc[ss,'fungi_pct'],s=18,marker=MK[k],color=C[k],edgecolor='white',linewidth=.6,label=f'Region {r} ({RN[int(r)]})',zorder=2)
med=fp.fungi_pct.groupby(m.stage).median().reindex(O)
a1.hlines(med.values,np.arange(6)-.3,np.arange(6)+.3,color='#0b0b0b',lw=1.2,zorder=3)
a1.set_xticks(range(6)); a1.set_xticklabels(OL); a1.set_ylabel('Fungal reads (% of total)'); a1.set_title('A  Fungal fraction',loc='left',fontsize=8.5)
a1.legend(fontsize=6,loc='upper right',handletextpad=.2); a1.grid(axis='y',color=LG,lw=.5); a1.set_axisbelow(True)
fk=['Aspergillus','Penicillium','Rhizomucor','Pichia','Mucor','Lichtheimia']
reg=fr.T.groupby(m.region).mean().T; rr=reg.loc[fk]; rr.loc['Other']=100-rr.sum()
bottom=np.zeros(6); fc=[C[0],C[1],C[2],C[3],C[4],C[6],GRAY]
for f,col in zip(rr.index,fc):
    a2.bar(range(6),rr.loc[f].values,bottom=bottom,color=col,width=.7,edgecolor='white',linewidth=.5,label=f); bottom+=rr.loc[f].values
a2.set_xticks(range(6)); a2.set_xticklabels([f'R{r} {RN[r]}' for r in range(1,7)],fontsize=7,rotation=35,ha='right',rotation_mode='anchor'); a2.set_ylabel('% of fungal reads'); a2.set_ylim(0,100)
a2.set_title('B  Fungal genera by region (mean of all stages)',loc='left',fontsize=8.5)
a2.legend(loc='center left',bbox_to_anchor=(1,.5),fontsize=6.5,handlelength=1)
fig.tight_layout(); fig.savefig('fig/Fig3_fungi.png'); plt.close()
# ---- Fig5 pathway heatmap (top 25 by region KW)
pw=pd.read_csv('humann_pathabundance_cpm_unstratified.tsv',sep='\t',index_col=0)[S]
kw=pd.read_csv('pathway_kw.tsv',sep='\t',index_col=0)
top=kw.index[:25]
Z=np.log10(pw.loc[top]+1).T.groupby(m.region).mean().T; Z=(Z.sub(Z.mean(1),axis=0)).div(Z.std(1)+1e-9,axis=0)
div=LinearSegmentedColormap.from_list('d',['#2a78d6','#f0efec','#e34948'])
fig,ax=plt.subplots(figsize=(5.6,5.2))
im=ax.imshow(Z.values,aspect='auto',cmap=div,vmin=-2,vmax=2)
ax.set_xticks(range(6)); ax.set_xticklabels([f'R{r} {RN[r]}' for r in range(1,7)],fontsize=7.5,rotation=40,ha='right',rotation_mode='anchor'); ax.set_yticks(range(len(top)))
import html as _h; ax.set_yticklabels([(lambda u:u if len(u)<62 else u[:59]+'…')(_h.unescape(t)) for t in top],fontsize=5.8)
for s in ax.spines.values(): s.set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.02); cb.set_label('Row z-score (log10 CPM)',fontsize=6.5); cb.ax.tick_params(labelsize=6)
fig.tight_layout(); fig.savefig('fig/FigS1_pathways.png'); plt.close()
# ---- Fig7 correlation heatmap
R=pd.read_csv('taxa_metabolite_correlation.tsv',sep='\t')
mets=['Tyr','Phe','Leu','Lys','Val','Ile','Glu','Thr','Ser','Ornithine','Hypoxanthine','Putrescine','Tyramine','GABA','Lactic acid','Citric acid']
mets=[x for x in mets if x in set(R.metabolite)]
tx=list(dict.fromkeys(R.taxon))
M=R.pivot(index='taxon',columns='metabolite',values='spearman_rho').loc[tx,mets]
Q1=R.pivot(index='taxon',columns='metabolite',values='q').loc[tx,mets]
Q2=R.pivot(index='taxon',columns='metabolite',values='q_within_perm').loc[tx,mets]
fig,ax=plt.subplots(figsize=(6.2,3.9))
im=ax.imshow(M.values,cmap=div,vmin=-1,vmax=1,aspect='auto')
for i in range(len(tx)):
    for j in range(len(mets)):
        if Q1.iat[i,j]<0.05 and Q2.iat[i,j]<0.05: ax.text(j,i,'●',ha='center',va='center',fontsize=6,color='#0b0b0b')
        elif Q1.iat[i,j]<0.05 or Q2.iat[i,j]<0.05: ax.text(j,i,'○',ha='center',va='center',fontsize=6,color='#0b0b0b')
ax.set_xticks(range(len(mets))); ax.set_xticklabels(mets,rotation=60,ha='right',fontsize=6.5)
ax.set_yticks(range(len(tx))); ax.set_yticklabels([t.replace('B:','').replace('F:','').replace(' (MAG)','*')+(' (fungi)' if t.startswith('F:') else '') for t in tx],fontsize=6.5)
for s in ax.spines.values(): s.set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.02); cb.set_label('Spearman ρ',fontsize=6.5); cb.ax.tick_params(labelsize=6)
fig.tight_layout(); fig.savefig('fig/Fig5_correlation.png'); plt.close()
print('done')
