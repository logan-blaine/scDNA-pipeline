from pysam import VariantFile, AlignmentFile, VariantRecord
import pandas as pd
import re
import os
import multiprocessing


class StrandedRecord:
    def __init__(self, rec, strand):
        self.strand = strand
        self.rec = rec

    def get_upstream_region(self, bases):
        start = self.rec.pos - (self.strand) * bases
        end = self.rec.pos + (1 - self.strand) * bases
        return f'{self.rec.chrom}:{start}-{end}'


class Pair:
    def __init__(self, id):
        self.mate1 = None
        self.mate2 = None
        self.id = id

    def is_paired(self):
        return self.mate1 and self.mate2

    def is_interchromosomal(self):
        return self.mate1.chrom != self.mate2.chrom

    def add_mate(self, rec):
        if not self.mate1:
            self.mate1 = rec
        elif not self.mate2:
            self.mate2 = rec
        else:
            raise(RuntimeError)

    def get_pair(self):
        assert(self.is_paired())
        str1 = int(self.mate1.alts[0][0] in 'actgACTG')
        str2 = int(self.mate2.alts[0][0] in 'actgACTG')
        rec1 = StrandedRecord(self.mate1, str1)
        rec2 = StrandedRecord(self.mate2, str2)
        return((rec1, rec2))


class PairedVcf:
    def __init__(self, vcf_path):
        self.pairs = {}
        vcf = VariantFile(vcf_path)
        for rec in vcf.fetch():
            if chr_regex.match(rec.chrom):
                self.add_entry(rec)

    def add_entry(self, rec):
        id = rec.id.split(':')[0]
        if id not in self.pairs:
            p = Pair(id)
            self.pairs[id] = p
        self.pairs[id].add_mate(rec)

    def get_pairs(self):
        return {p.get_pair() for p in pc.pairs.values() if p.is_paired()}


chr_regex = re.compile('chr[0-9XY]+')
win_size = 2000
tmp_ext = ".tmp.txt"

vcf_path = snakemake.input[0]
bam_paths = snakemake.input[1:]
output_path = snakemake.output[0]
# win_size = snakemake.params["window_size"]

pc = PairedVcf(vcf_path)
pairs = pc.get_pairs()


def recount_on_file(bam_path):
    bam = AlignmentFile(bam_path)
    dir_sample, ext = os.path.splitext(bam_path)
    sample = os.path.basename(dir_sample)
    assert(ext == ".bam")
    output_path = dir_sample + tmp_ext
    ret = ['\t'.join(['chr1', 'pos1', 'str1', 'chr2',
                      'pos2', 'str2', 'count', 'hq_count', 'sample'])]

    for paired_rec1, paired_rec2 in pairs:
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

        hq_bam1 = {rec.query_name for rec in bam1 if
                   rec.mapq >= 30 and not rec.is_supplementary}
        hq_bam2 = {rec.query_name for rec in bam2 if
                   rec.mapq >= 30 and not rec.is_supplementary}
        n_shared_hq = len(hq_bam1.intersection(hq_bam2))

        if n_shared:
            tokens = [rec1.chrom, str(rec1.pos), str1,
                      rec2.chrom, str(rec2.pos), str2,
                      str(n_shared), str(n_shared_hq), sample]
            ret.append('\t'.join(tokens))

    with open(output_path, 'w') as f:
        f.write('\n'.join(ret) + '\n')
    return output_path


with multiprocessing.Pool(snakemake.threads) as p:
    out_files = p.map(recount_on_file, bam_paths)

tables = (pd.read_table(o, index_col=None) for o in out_files)

(
    pd.concat(tables)
    .sort_values(["chr1", "pos1", "str1", "chr2", "pos2", "str2", "sample"])
    .to_csv(output_path, index_label="idx")
)

# For debugging only:
# os.chdir('/Users/pellmanlab/GenomeAnalysis/svbam')
# bam_path = 'MN_SI_181116_P1_22_1.bam'
# bam_paths = ['MN_SI_181116_P1_22_1.bam', 'MN_SI_181116_P1_22_2.bam']
# vcf_path = 'MN_SI_181116_P1_22.svaba.prefiltered.somatic.sv.vcf'
# output_path = 'MN_SI_181116_P1_22_1.sv_count.csv'
