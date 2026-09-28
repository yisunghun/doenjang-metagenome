# 진균 2단계 Bracken 결과 병합 (입력: 각 샘플 05c_taxonomy_fungi_2step/*.fungi.bracken_{S,G}.txt, fungi_2step_summary.tsv)
import pandas as pd, re, os, glob
R = os.environ.get('RESULTS_DIR', '/mnt/f/meta_hong/results')
m = pd.read_csv('metadata.tsv', sep='\t', index_col=0); S = list(m.index)
rows = []
for s in S:
    for L in ['S', 'G']:
        f = f'{R}/{s}/05c_taxonomy_fungi_2step/{s}.fungi.bracken_{L}.txt'
        b = pd.read_csv(f, sep='\t'); b['s'] = s; b['lvl'] = L; rows.append(b)
L = pd.concat(rows)
fp = []
for l in open(f'{R}/v3_merged/fungi_2step_summary.tsv'):
    if l.strip():
        fp.append((l.split('\t')[0], float(re.search(r'fungi=\d+ \(([\d.]+)%', l).group(1)),
                   float(re.search(r'unclassified=\d+ \(([\d.]+)%', l).group(1))))
fp = pd.DataFrame(fp, columns=['sample', 'fungi_pct', 'step1_U_pct']).set_index('sample').loc[S]
fp['stage'] = m.stage; fp['region'] = m.region; fp['residual_U_pct'] = fp.step1_U_pct - fp.fungi_pct
fp[['stage', 'region', 'fungi_pct', 'step1_U_pct', 'residual_U_pct']].to_csv('fungi_2step_summary_parsed.tsv', sep='\t')
for lvl in ['G', 'S']:
    t = L[L.lvl == lvl].pivot_table(index='name', columns='s', values='fraction_total_reads', aggfunc='sum').fillna(0)[S] * 100
    t = t.loc[t.mean(axis=1).sort_values(ascending=False).index]
    t.to_csv(f'fungi_bracken_{lvl}_relab_within_fungi.tsv', sep='\t', float_format='%.4g')
    (t * fp.loc[S, 'fungi_pct'] / 100).to_csv(f'fungi_bracken_{lvl}_pct_of_total.tsv', sep='\t', float_format='%.4g')
