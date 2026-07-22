import sys, pandas as pd, numpy as np
sd = sys.argv[1]
rc = pd.read_csv(f'{sd}/readcounts.tsv', sep='\t', index_col=0)
baf_in = pd.read_csv(f'{sd}/precomputed_baf.tsv', sep='\t')
regions = pd.read_csv(f'{sd}/filtered_regions.bed', sep='\t', header=None, names=['chrom','start','end'])
regions['bin_idx'] = range(len(regions))
lookup = {}
for _, row in regions.iterrows():
    lookup[(row['chrom'], int(row['start']))] = int(row['bin_idx'])
    lookup[(row['chrom'], int(row['start'])-1)] = int(row['bin_idx'])
baf_in['bin_idx'] = baf_in.apply(lambda r: lookup.get((r['chrom'], int(r['start'])), lookup.get((r['chrom'], int(r['start'])+1), None)), axis=1)
baf_in = baf_in.dropna(subset=['bin_idx'])
baf_in['bin_idx'] = baf_in['bin_idx'].astype(int)
baf_mat = baf_in.pivot(index='cell', columns='bin_idx', values='BAF')
baf_mat = baf_mat.reindex(index=rc.index.tolist(), columns=range(rc.shape[1]), fill_value=0.5)
baf_mat.to_csv(f'{sd}/BAF.tsv', sep='\t')
print(f'BAF.tsv written: {baf_mat.shape}')
