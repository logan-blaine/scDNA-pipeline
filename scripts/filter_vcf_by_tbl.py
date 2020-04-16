from sv_vcf_utils import PairedVcf, Pair, StrandedRecord
from pysam import VariantFile
import pandas as pd


vcf_path = snakemake.input[0]
table_path = snakemake.input[1]
output_vcf = snakemake.output[0]
min_count = 3
# win_size = snakemake.params['window_size']

# print("Reading "+vcf_path)
pc = PairedVcf(vcf_path)
pairs = pc.get_pairs()

chr_pos_str=['chr1', 'pos1', 'str1', 'chr2', 'pos2', 'str2']
id_cols = chr_pos_str+['sample']

# print("Reading "+table_path)
tbl = pd.read_csv(table_path)
tbl_deduped = tbl.groupby(id_cols).first().reset_index()
tbl_deduped['n']=tbl_deduped.groupby(chr_pos_str)['sample'].transform('count')
tbl_filtered=tbl_deduped.query(f'n==1 and count>={min_count} and hq_count>0')

vcf_head=VariantFile(vcf_path).header
vcf_head.add_line(f'##command=recount_svs.py {vcf_path} {table_path} > {output_vcf}')

# print("Writing "+output_vcf)
# vcf_out=VariantFile(output_vcf, mode='w', header=vcf_head)

with open(output_vcf, 'w') as f:
    f.write(str(vcf_head))
    for _, row in tbl_filtered.iterrows():
        hash = '\t'.join(map(str, row[chr_pos_str]))
        if hash in pairs:
            mate1, mate2 = pairs[hash].get_pair()
            f.write(str(mate1.rec))
            f.write(str(mate2.rec))
