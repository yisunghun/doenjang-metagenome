# 미생물-대사체 편순위상관 (지역·숙성단계 동시 통제) — correlate.py 실행 후
import pandas as pd, numpy as np
from scipy.stats import rankdata, t as T
exec(open(__file__.replace('partial_corr.py', 'correlate.py')).read().split("rows=[]")[0].replace("print(", "(lambda *a,**k:None)("))
st = m.loc[S, 'stage'].values
X = np.column_stack([np.ones(len(S))] + [(reg == r).astype(float) for r in sorted(set(reg))[1:]] + [(st == s).astype(float) for s in sorted(set(st))[1:]])
H = X @ np.linalg.pinv(X); dfres = len(S) - np.linalg.matrix_rank(X) - 1
res = lambda v: (lambda r: r - H @ r)(rankdata(v))
R = pd.read_csv('taxa_metabolite_correlation.tsv', sep='\t'); pr = []; pp = []
for _, row in R.iterrows():
    a = res(tax.loc[row.taxon, S].values); b = res(met[row.metabolite].values)
    r = np.corrcoef(a, b)[0, 1] if a.std() > 0 and b.std() > 0 else np.nan
    tt = r * np.sqrt(dfres / (1 - r ** 2)) if abs(r) < 1 else np.inf
    pr.append(r); pp.append(2 * T.sf(abs(tt), dfres))
R['partial_r_region_stage'] = pr; R['p_partial'] = pp
ok = R.p_partial.notna(); R.loc[ok, 'q_partial'] = bh(R.loc[ok, 'p_partial'].values)
R.to_csv('taxa_metabolite_correlation.tsv', sep='\t', index=False, float_format='%.4g')
sel = [('B:Tetragenococcus', x) for x in ['Tyr', 'Phe', 'Leu', 'Lys', 'Glu']] + [('B:Enterococcus', x) for x in ['Hypoxanthine', 'Lys', 'Putrescine', 'Tyramine']] + [('B:Bacillus', x) for x in ['Phe', 'Glu']] + [('F:Penicillium', 'GABA')]
R.set_index(['taxon', 'metabolite']).loc[sel][['spearman_rho', 'q', 'within_region_r', 'q_within', 'partial_r_region_stage', 'p_partial']].to_csv('table5_sel.tsv', sep='\t')
print('partial q<0.05:', (R.q_partial < 0.05).sum())
