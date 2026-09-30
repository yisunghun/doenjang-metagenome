import os; os.makedirs("fig", exist_ok=True)
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.edgecolor':'#8a8985','axes.linewidth':0.6,
  'axes.labelcolor':'#0b0b0b','xtick.color':'#52514e','ytick.color':'#52514e','axes.spines.top':False,'axes.spines.right':False,
  'savefig.dpi':300,'figure.facecolor':'white','legend.frameon':False})
C=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#4a3aa7']
d=pd.read_excel('HMT.xlsx'); d['st']=d['숙성기간'].ffill()
order=['Meju','담금직후','1개월','3개월','6 months']; lab=['Meju','0M','1M','3M','6M']
g=d.groupby('st'); m=g.mean(numeric_only=True).loc[order]; se=(g.std(numeric_only=True)/np.sqrt(g.size().values[:,None] if False else g.count())).loc[order]
fig,(a1,a2)=plt.subplots(1,2,figsize=(7.2,2.9))
x=np.arange(5)
for k,c in enumerate(['Tyr','Phe','Leu','Lys']):
    a1.errorbar(x,m[c],yerr=se[c],marker='o',ms=4,lw=1.6,capsize=2,color=C[k],label=c)
for k,c in enumerate(['Putrescine','Tyramine']):
    a2.errorbar(x,m[c],yerr=se[c],marker='s',ms=4,lw=1.6,capsize=2,color=C[k+4],label=c)
for a,t in ((a1,'A  Free amino acids'),(a2,'B  Biogenic amines')):
    a.set_xticks(x); a.set_xticklabels(lab); a.set_xlabel('Fermentation stage'); a.set_title(t,loc='left',fontsize=9)
    a.grid(axis='y',color='#e4e3df',lw=.6); a.set_axisbelow(True); a.legend(loc='upper left',fontsize=7)
a1.set_ylabel('Relative peak area (mean ± SE)')
fig.tight_layout(); fig.savefig('fig/FigS2_metabolites.png'); plt.close()
print(m[['Tyr','Lys','Putrescine','Tyramine']].round(0))
