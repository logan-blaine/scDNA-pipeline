from pysam import VariantFile, AlignmentFile, VariantRecord
import re


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

    def span(self):
        if not self.is_paired():
            return None
        elif self.is_interchromosomal():
            return float('inf')
        else:
            return abs(self.mate1.pos - self.mate2.pos)

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

    def get_pairs(self, min_span=1e5):
        return {p.get_pair() for id, p in pc.pairs.items()
                if p.is_paired() and p.span() > min_span}
        # if p.is_paired() and (p.is_interchromosomal() or p.span()>min_span)}


if __name__ == '__main__':
    chr_regex = re.compile('chr[0-9XY]+')
    # bam_path = 'MN_SI_181116_P1_22_1.bam'
    # vcf_path = 'MN_SI_181116_P1_22.svaba.prefiltered.somatic.sv.vcf'
    # output_path = 'MN_SI_181116_P1_22_1.sv_count.tsv'
    bam_path = snakemake.input[0]
    vcf_path = snakemake.input[1]
    output_path = snakemake.output[0]

    bam = AlignmentFile(bam_path)
    pc = PairedVcf(vcf_path)

    ret=['\t'.join(['chr1','pos1','str1','chr2','pos2','str2','TotalCount'])]
    for paired_rec1, paired_rec2 in pc.get_pairs():
        rec1 = paired_rec1.rec
        rec2 = paired_rec2.rec
        if rec1.chrom!='chr5' or rec2.chrom!='chr5':
            continue
        loc1 = paired_rec1.get_upstream_region(2000)
        loc2 = paired_rec2.get_upstream_region(2000)
        loc1_ids = {rec.query_name for rec in bam.fetch(
            region=loc1, multiple_iterators=True)}
        loc2_ids = {rec.query_name for rec in bam.fetch(
            region=loc2, multiple_iterators=True)}
        supp_reads = loc1_ids.intersection(loc2_ids)
        if not supp_reads:
            continue
        line = '\t'.join([rec1.chrom, str(rec1.pos), str(paired_rec1.strand*2-1),
                         rec2.chrom, str(rec2.pos), str(paired_rec2.strand*2-1),
                         str(len(supp_reads))])
        ret.append(line)
    with open(output_path, 'w') as f:
        f.write('\n'.join(ret)+'\n')

    # p = Pool(4)
    # res = [t for t in p.map(get_supporting_reads, range(10)) if t]
    # res
