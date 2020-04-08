from sv_vcf_utils import StrandedRecord, Pair, PairedVcf
from pysam import AlignmentFile
import pandas as pd
import re
import os
import multiprocessing

vcf_path = snakemake.input[0]
bam_paths = snakemake.input[1:]
output_txt = snakemake.output[0]
group = snakemake.wildcards['group']
threads = snakemake.threads

win_size = 2000
tmp_ext = '.tmp.csv'

# win_size = snakemake.params['window_size']

pc = PairedVcf(vcf_path)
pairs = pc.get_pairs()
output_dir = os.path.dirname(output_txt)
temp_dir = os.path.join(output_dir, group)
os.makedirs(temp_dir, exist_ok=True)


def recount_on_file(bam_path):
    bam = AlignmentFile(bam_path)
    sample, ext = os.path.splitext(os.path.basename(bam_path))
    assert(ext == '.bam')
    output_path = os.path.join(temp_dir, sample + tmp_ext)
    ret = ['\t'.join(['chr1', 'pos1', 'str1', 'chr2',
                      'pos2', 'str2', 'count', 'hq_count', 'sample'])]

    for pair in pairs.values():
        paired_rec1, paired_rec2 = pair.get_pair()
        rec1 = paired_rec1.rec
        rec2 = paired_rec2.rec
        interchrom = (rec1.chrom != rec1.chrom)
        str1 = str(paired_rec1.strand * 2 - 1)
        str2 = str(paired_rec2.strand * 2 - 1)
        loc1 = paired_rec1.get_upstream_region(win_size)
        loc2 = paired_rec2.get_upstream_region(win_size)

        # DEBUG ONLY
        # if rec1.chrom != 'chr5' or rec2.chrom != 'chr5':
        #     continue

        bam1 = list(bam.fetch(region=loc1, multiple_iterators=True))
        bam2 = list(bam.fetch(region=loc2, multiple_iterators=True))
        id1 = {rec.query_name for rec in bam1}
        id2 = {rec.query_name for rec in bam2}
        n_shared = len(id1.intersection(id2))

        hq_bam1 = {rec.query_name for rec in bam1 if rec.mapq >= 30
                   and not (rec.is_supplementary or rec.is_secondary)}
        hq_bam2 = {rec.query_name for rec in bam2 if rec.mapq >= 30
                   and not (rec.is_supplementary or rec.is_secondary)}
        n_shared_hq = len(hq_bam1.intersection(hq_bam2))

        if n_shared:
            tokens = [rec1.chrom, str(rec1.pos), str1,
                      rec2.chrom, str(rec2.pos), str2,
                      str(n_shared), str(n_shared_hq), sample]
            ret.append('\t'.join(tokens))

    bam.close()
    with open(output_path, 'w') as f:
        f.write('\n'.join(ret) + '\n')
    return output_path


with multiprocessing.Pool(threads) as p:
    out_files = p.map(recount_on_file, bam_paths)

tables = (pd.read_table(o, index_col=None) for o in out_files)

(
    pd.concat(tables)
    .sort_values(['chr1', 'pos1', 'str1', 'chr2', 'pos2', 'str2', 'sample'])
    .to_csv(output_txt, index_label='idx')
)

# For debugging only:
# os.chdir('/Users/pellmanlab/GenomeAnalysis/svbam')
# bam_path = 'MN_SI_181116_P1_22_1.bam'
# bam_paths = ['MN_SI_181116_P1_22_1.bam', 'MN_SI_181116_P1_22_2.bam']
# vcf_path = 'MN_SI_181116_P1_22.svaba.prefiltered.somatic.sv.vcf'
# output_txt = 'MN_SI_181116_P1_22_1.sv_count.csv'
# threads=2
# tbl['hq_count']=1
