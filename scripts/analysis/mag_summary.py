# 대표 MAG 요약표 (CheckM 완성도/오염도는 MultiQC general stats에서, GTDB-Tk 분류는 reps_gtdb.tsv)
import pandas as pd, json
m = pd.read_csv('metadata.tsv', sep='\t', index_col=0); S = list(m.index)
g = pd.read_csv('reps_gtdb.tsv', sep='\t', header=None, names=['rep', 'cls'])
g['sp'] = g.cls.str.extract(r's__([^;]*)$')[0].fillna(''); g['gen'] = g.cls.str.extract(r'g__([^;]*)')[0]
cl = pd.read_csv('clusters.tsv', sep='\t', header=None, names=['rep', 'mem']); cl['rep'] = cl.rep.str.replace('.fa', '', regex=False)
d = json.load(open('multiqc_general_stats.json'))['general_stats_table']['datasets'][0]['violin_value_by_sample_by_metric']
g['comp'] = g.rep.map(d['checkm-Completeness']); g['cont'] = g.rep.map(d['checkm-Contamination']); g['n_members'] = g.rep.map(cl.groupby('rep').size())
mag = pd.read_csv('mag_relative_abundance.tsv', sep='\t', index_col=0)[S].drop('unmapped')
g['max_abund'] = g.rep.map(mag.max(axis=1)); g['prevalence'] = g.rep.map((mag > 0.1).sum(axis=1))
g.sort_values('max_abund', ascending=False).to_csv('mag_reps_summary.tsv', sep='\t', index=False, float_format='%.3g')
