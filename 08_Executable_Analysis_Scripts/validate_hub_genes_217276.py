"""
Cross-validation of PBM/MLV hub genes against GSE217276 (Rustenhoven et al.,
J Exp Med 2023) -- bulk RNA-seq of sorted meningeal lymphatic endothelial
cells (LECs), AAV-IFNgamma overexpression vs vehicle in young mice, n=3/arm.

This is NOT a light/PBM test. IFNgamma is the causal driver this paper
identifies for age-related meningeal lymphatic dysfunction (elevated in aged
dura via T cell accumulation; AAV-IFNg in young mice phenocopies the aged
CSF-drainage deficit). So this asks: are the hub genes altered by the actual
molecular driver of the age-related dysfunction PBM aims to reverse -- a
disease-relevance question, kept explicitly separate from the light-response
question already answered (null) in GSE245658.
"""
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 220)

df = pd.read_csv('counts.csv', index_col=0).set_index('ensembl_gene_id')
df = df[~df.index.duplicated(keep='first')]  # 455 rows are annotation-join dupes w/ identical counts
count_cols = ['IFN.S1','IFN.S2','IFN.S3','VEH.S1','VEH.S2','VEH.S3']
counts = df[count_cols].astype(int)
symbol = df['external_gene_name']

sample_meta = pd.DataFrame({'condition': ['IFN','IFN','IFN','VEH','VEH','VEH']},
                            index=count_cols)

gene_map = {
    'ACTB':['Actb'], 'CTNNB1':['Ctnnb1'], 'IL6':['Il6'], 'TGFB1':['Tgfb1'], 'GAPDH':['Gapdh'],
    'TNF':['Tnf'], 'TP53':['Trp53'], 'HIF1A':['Hif1a'], 'KDR':['Kdr'], 'PTEN':['Pten'],
    'HSP90AA1':['Hsp90aa1'], 'STAT3':['Stat3'], 'INS':['Ins1','Ins2'], 'SRC':['Src'], 'ALB':['Alb']
}
hub_rows = []
for human_sym, mouse_syms in gene_map.items():
    for ms in mouse_syms:
        for ens_id in symbol[symbol == ms].index:
            hub_rows.append({'human_symbol': human_sym, 'mouse_symbol': ms, 'ensembl_id': ens_id})
hub_df = pd.DataFrame(hub_rows).set_index('ensembl_id')
hub_df['total_counts'] = hub_df.index.map(lambda i: counts.loc[i].sum())
hub_df['detected'] = hub_df['total_counts'] > 10  # more than trace noise

keep = counts.sum(axis=1) >= 10
counts_f = counts[keep]
print(f"Genes retained after prefilter: {counts_f.shape[0]} / {counts.shape[0]}")

dds = DeseqDataSet(counts=counts_f.T, metadata=sample_meta, design="~condition",
                    refit_cooks=True, quiet=True)
dds.deseq2()
print("Design matrix coefficients:", dds.obsm['design_matrix'].columns.tolist())

stat = DeseqStats(dds, contrast=["condition","IFN","VEH"], quiet=True)
stat.summary()
res = stat.results_df.copy()
design_cols = dds.obsm['design_matrix'].columns.tolist()
if 'condition[T.IFN]' in design_cols:
    coeff, sign = 'condition[T.IFN]', 1.0
else:
    coeff, sign = 'condition[T.VEH]', -1.0  # model's coeff is VEH-vs-IFN; flip sign for IFN-vs-VEH
stat.lfc_shrink(coeff=coeff)
res['log2FoldChange'] = sign * stat.results_df['log2FoldChange']
res['lfcSE'] = stat.results_df['lfcSE']
res.columns = [f"IFN_vs_VEH_{c}" for c in res.columns]
res.to_csv('full_transcriptome_results_GSE217276.csv')

n_sig05 = (res['IFN_vs_VEH_padj'] < 0.05).sum()
print(f"\nGenome-wide: {n_sig05} / {len(res)} genes padj<0.05 for IFN vs VEH")

merged = hub_df.join(res, how='left')
for grp, cols in [('IFN', ['IFN.S1','IFN.S2','IFN.S3']), ('VEH', ['VEH.S1','VEH.S2','VEH.S3'])]:
    merged[f'mean_count_{grp}'] = counts.loc[merged.index, cols].mean(axis=1)

golden = pd.read_csv('/home/claude/geo/Final_Candidates_With_AlphaScore.csv')[['Gene','Golden_Score','Origin']]
merged = merged.reset_index().rename(columns={'index':'ensembl_id'})
merged = merged.merge(golden, left_on='human_symbol', right_on='Gene', how='left').drop(columns=['Gene'])
merged = merged.sort_values('Golden_Score', ascending=False)

cols_order = ['human_symbol','Golden_Score','Origin','mouse_symbol','ensembl_id','detected',
              'mean_count_VEH','mean_count_IFN',
              'IFN_vs_VEH_log2FoldChange','IFN_vs_VEH_pvalue','IFN_vs_VEH_padj']
merged = merged[cols_order]
merged.to_csv('hub_gene_validation_results_GSE217276.csv', index=False)
print("\n=== HUB GENE RESULTS: IFNg vs VEH in meningeal LECs (GSE217276) ===")
print(merged.to_string(index=False))
