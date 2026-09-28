# 숙성 단계 효과의 시계열 민감도 분석: 지역 내 시점 순서를 보존한 순환 이동(cyclic shift) 순열, 전체 21,600가지 열거
import pandas as pd, numpy as np, itertools
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata
m = pd.read_csv('metadata.tsv', sep='\t', index_col=0); S = list(m.index)
so = {'Me': 0, '0M': 1, '1M': 2, '3M': 3, '6M': 4, '12M': 5}
reg = m.region.values; st = m.stage.values; n = len(S)
blocks = [sorted(np.where(reg == r)[0], key=lambda i: so[st[i]]) for r in sorted(set(reg))]
def design(fs):
    cols = [np.ones(n)]
    for f in fs:
        for l in np.unique(f)[1:]: cols.append((f == l).astype(float))
    return np.column_stack(cols)
hat = lambda X: X @ np.linalg.pinv(X)
Hr = hat(design([reg]))
def all_shifts():
    for ks in itertools.product(*[range(len(b)) for b in blocks]):
        lab = st.copy()
        for b, k in zip(blocks, ks):
            seq = [st[i] for i in b]; rot = seq[k:] + seq[:k]
            for i, v in zip(b, rot): lab[i] = v
        yield lab
shifts = list(all_shifts()); print('permutations:', len(shifts))
def Fstage_G(G, lab):
    H = hat(design([reg, lab])); I = np.eye(n)
    ss = np.trace(G @ (H - Hr)); sr = np.trace((I - H) @ G @ (I - H))
    return (ss / 5) / (sr / (n - np.linalg.matrix_rank(design([reg, lab]))))
tabs = {'Bacteria species (Kraken2)': pd.read_csv('bracken_species_relab.tsv', sep='\t', index_col=0)[S]}
fs = pd.read_csv('fungi_bracken_S_relab_within_fungi.tsv', sep='\t', index_col=0)[S]
fs.loc['Aspergillus flavus/oryzae'] = fs.loc[['Aspergillus oryzae', 'Aspergillus flavus']].sum(); tabs['Fungi species'] = fs.drop(['Aspergillus oryzae', 'Aspergillus flavus'])
tabs['Fungi genus'] = pd.read_csv('fungi_bracken_G_relab_within_fungi.tsv', sep='\t', index_col=0)[S]
tabs['MAG 63'] = pd.read_csv('mag_relab_within_mapped.tsv', sep='\t', index_col=0)[S]
tabs['Bacteria MAG + fungi genus'] = pd.read_csv('combined_bacteriaMAG_fungiGenus_relab.tsv', sep='\t', index_col=0)[S]
tabs['MetaCyc pathways'] = pd.read_csv('humann_pathabundance_cpm_unstratified.tsv', sep='\t', index_col=0)[S]
rows = []
for k, t in tabs.items():
    D = squareform(pdist(t.T.values, 'braycurtis')); J = np.eye(n) - 1 / n; G = J @ (-0.5 * D ** 2) @ J
    F0 = Fstage_G(G, st); Fs = np.array([Fstage_G(G, lab) for lab in shifts])
    rows.append((k, F0, (Fs >= F0 - 1e-12).mean())); print(k, round(F0, 3), 'series p=%.4f' % rows[-1][2])
# 분류군 추세: 순위 변환 지역 내 상관, series p
tot = pd.read_csv('composition_pct_total_reads.tsv', sep='\t', index_col=0)[S]
def rc(v): r = pd.Series(rankdata(v), index=S); return (r - r.groupby(reg).transform('mean')).values
keys = ['Bacillus', 'Enterococcus', 'Pediococcus', 'Tetragenococcus', 'Heyndrickxia', 'Weissella', 'Pseudomonas', 'Kroppenstedtia (MAG)', 'Aspergillus (F)', 'Penicillium (F)', 'Rhizomucor (F)', 'Pichia (F)']
ordv = lambda lab: np.array([so[x] for x in lab]); Y = [rc(ordv(lab)) for lab in shifts]; y0 = rc(ordv(st))
tr = []
for k in keys:
    x = rc(tot.loc[k].values); r0 = np.corrcoef(x, y0)[0, 1]; rs = np.array([np.corrcoef(x, y)[0, 1] for y in Y])
    tr.append((k, r0, (np.abs(rs) >= abs(r0) - 1e-12).mean()))
T = pd.DataFrame(tr, columns=['taxon', 'rm_r', 'p_series']); print(T.round(4).to_string())
pd.DataFrame(rows, columns=['data', 'F_stage', 'p_series']).to_csv('series_perm_permanova.tsv', sep='\t', index=False)
T.to_csv('series_perm_trend.tsv', sep='\t', index=False)
