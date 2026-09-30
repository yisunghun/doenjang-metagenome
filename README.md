# Shotgun metagenomics from meju through tojang to doenjang: a case study of six traditional producers

Data and code accompanying the manuscript
**"Microbial community dynamics and their associations with metabolite changes from meju through tojang to doenjang: a shotgun metagenomic case study of six traditional producers"**
(Yi SH†, Lee JE†, Chin YW, Lee SJ, Hong SP; submitted to *Food Science and Biotechnology*).

## Study design
- One certified traditional producer in each of six Korean provinces: regions 1–6 = Gangwon, Gyeonggi, Gyeongsang, Jeolla, Chungcheong, and Jeju. Each region represents a single producer.
- Process: meju → tojang (meju ripened with brine; 0 and 1 month) → soy sauce separation on day 60 → doenjang (3, 6, and 12 months).
- 32 shotgun metagenome libraries (missing: region 4 12M, region 5 3M, region 6 3M/6M); Illumina HiSeq 4000, 2 × 101 bp; samples and data obtained in 2016 (KFRI Main Research Program E0141500-03).
- Metabolome: HMT CE-TOF/MS and LC-TOF/MS, 77 compounds, meju–6 months (30 samples; 27 paired with metagenomes).

## Repository layout
```
data/
  metadata/sample_metadata.tsv          sample ID mapping (metagenome ↔ region/stage ↔ metabolome ID)  [Supplementary Table S8]
  metabolome/HMT_CE-LC-TOFMS_relative_quantification.csv   raw relative quantification (77 compounds × 30 samples)
  reference_db/fungal_refseq_674_assemblies.csv            RefSeq fungal assemblies used for the custom Kraken2 DB [Supplementary Data 1]
  processed/                             merged tables used for all statistics and figures (see below)
mags/                                    63 species-level representative MAGs (*.fa.gz) + GTDB-Tk classification
scripts/
  pipeline/                              per-sample processing (WSL2/Ubuntu, conda)
  analysis/                              merging, statistics, figures (Python 3)
figures/                                 manuscript Figs. 1–5 and Supplementary Figs. S1–S2
```

## Pipeline (per sample)
| Step | Tool (version) | Script |
|---|---|---|
| QC / trimming | fastp 0.24.1 | `03_run_shotgun_meta_hong.sh` |
| Host/plant read removal (soybean, human, barley, rice, wheat) | Bowtie2 2.5.4 | 〃 |
| Bacterial taxonomy | Kraken2 2.14 + Standard-16 index (2024-09), Bracken 3.1 (read length 100) | 〃 |
| Assembly / binning / QC / taxonomy | MEGAHIT 1.2.9, MetaBAT2 2.17, CheckM 1.2.3, GTDB-Tk 2.4.1 (GTDB R226) | 〃 |
| Functional profiling | HUMAnN 3.9 | 〃 |
| Fungal 2-step classification (unclassified reads → custom RefSeq fungi DB) | Kraken2 2.14, Bracken 3.1, dustmasker, seqkit | `05_fungal_kraken_2step_meta_hong.sh` |
| MAG catalog (ANI 95%) and read-mapping abundance | galah 0.5.2, CoverM 0.8.0 | `06_mag_abundance_meta_hong.sh` |

`05_…sh --build-db` rebuilds the fungal database from NCBI RefSeq (list in `data/reference_db`). Paths at the top of each script (`RESULTS_DIR`, `DB_ROOT`, `CONDA_ENV`) must be adapted.

## Analysis (reproducing tables and figures)
Run from `data/processed/` in this order (Python ≥3.10 with pandas, numpy, scipy, matplotlib, openpyxl):
```bash
cd data/processed
python ../../scripts/analysis/stats_v3.py      # QC/classification summaries, stage & region tables
python ../../scripts/analysis/permanova.py     # sequential PERMANOVA (initial)
python ../../scripts/analysis/review_stats.py  # marginal PERMANOVA with restricted permutations, PERMDISP, HUMAnN sensitivity
python ../../scripts/analysis/review2.py       # stage-only tests blocked by region (alpha diversity, PERMDISP)
python ../../scripts/analysis/correlate.py     # taxa–metabolite Spearman / within-region correlations
python ../../scripts/analysis/partial_corr.py  # partial rank correlation (region + stage)
python ../../scripts/analysis/review2b.py      # within-region block-permutation p-values
python ../../scripts/analysis/mag_summary.py   # representative MAG summary
python ../../scripts/analysis/series_perm.py   # series (cyclic-shift) permutation sensitivity analysis
python ../../scripts/analysis/figs.py          # Figs. 1–5 and Fig. S1 → ./fig
python ../../scripts/analysis/figs_metabolites.py  # Fig. S2 → ./fig
```
`merge.py` and `fungi_merge.py` rebuild the merged tables in `data/processed/` from the per-sample pipeline outputs (set `RESULTS_DIR`).
Permutation-based results use fixed random seeds; exact p-values may differ slightly with numpy versions.

### Figures
| File | Manuscript |
|---|---|
| `Fig1_PCoA.png` | Fig. 1 PCoA of bacterial communities |
| `Fig2_composition.png` | Fig. 2 Microbial composition of each sample |
| `Fig3_fungi.png` | Fig. 3 Fungal reads by stage and fungal genus composition |
| `Fig4_succession.png` | Fig. 4 Changes in major taxa (dashed line, soy sauce separation on day 60) |
| `Fig5_correlation.png` | Fig. 5 Taxon–metabolite correlations |
| `FigS1_pathways.png` | Supplementary Fig. S1 MetaCyc pathways by region |
| `FigS2_metabolites.png` | Supplementary Fig. S2 Free amino acids and biogenic amines |

### Key processed tables
| File | Content |
|---|---|
| `bracken_species_counts.tsv`, `bracken_species_relab.tsv`, `bracken_genus_relab.tsv` | Bacterial Bracken estimates (new_est_reads; relative abundance within Bracken total) |
| `fungi_bracken_{S,G}_relab_within_fungi.tsv`, `fungi_bracken_{S,G}_pct_of_total.tsv` | Fungal 2-step Bracken (within fungi; % of total input reads) |
| `composition_pct_total_reads.tsv` | Bacteria + fungi + Kroppenstedtia (MAG) as % of total input reads |
| `mag_relative_abundance.tsv`, `mag_reps_summary.tsv`, `reps_gtdb.tsv`, `clusters.tsv` | CoverM abundance of 63 representative MAGs; quality and GTDB-Tk; galah clusters |
| `humann_pathabundance_cpm_unstratified.tsv` | HUMAnN MetaCyc pathways (CPM, UNMAPPED/UNINTEGRATED removed) |
| `sample_summary.tsv` | Per-sample classification rates, diversity, HUMAnN UNMAPPED fraction |
| `review_permanova.tsv`, `stage_trend_perm.tsv`, `series_perm_permanova.tsv`, `series_perm_trend.tsv`, `taxa_metabolite_correlation.tsv`, `tableS1.tsv` | Statistical results reported in the manuscript |

## Raw sequence data
Raw FASTQ files are available from the corresponding author (Sang-Pil Hong, sphong@kfri.re.kr) upon reasonable request.

## Citation
Yi SH, Lee JE, Chin YW, Lee SJ, Hong SP. Microbial community dynamics and their associations with metabolite changes from meju through tojang to doenjang: a shotgun metagenomic case study of six traditional producers. Submitted to Food Science and Biotechnology.

Data and code archive: Zenodo, https://doi.org/10.5281/zenodo.23010290 (all versions; version 1.1: https://doi.org/10.5281/zenodo.23049995)

## License
Code in `scripts/` is released under the MIT License (see `LICENSE`). Data in `data/`, `mags/` and `figures/` are released under the Creative Commons Attribution 4.0 International License (CC BY 4.0; https://creativecommons.org/licenses/by/4.0/). Please cite the associated manuscript when using these data.
