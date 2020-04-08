from pysam import VariantFile
import re


class StrandedRecord:
    strand_chars = ['-1', '1']

    def __init__(self, rec, strand):
        self.strand = strand
        self.rec = rec

    def get_upstream_region(self, bases):
        start = self.rec.pos - (self.strand) * bases
        end = self.rec.pos + (1 - self.strand) * bases
        return f'{self.rec.chrom}:{start}-{end}'

    def __str__(self):
        strand = StrandedRecord.strand_chars[self.strand]
        return('\t'.join([self.rec.chrom, str(self.rec.pos), strand]))


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
        strand = int(rec.alts[0][0] in 'actgACTG')
        if not self.mate1:
            self.mate1 = StrandedRecord(rec, strand)
        elif not self.mate2:
            self.mate2 = StrandedRecord(rec, strand)
        else:
            raise(RuntimeError)

    def get_pair(self):
        assert(self.is_paired())
        return((self.mate1, self.mate2))

    def __str__(self):
        return('\t'.join([str(self.mate1), str(self.mate2)]))


class PairedVcf:
    def __init__(self, vcf_path):
        self.pairs = {}
        vcf = VariantFile(vcf_path)
        for rec in vcf.fetch():
            if re.match('^chr[0-9XY]+$', rec.chrom):
                self.add_entry(rec)

    def add_entry(self, rec):
        id = rec.id.split(':')[0]
        if id not in self.pairs:
            p = Pair(id)
            self.pairs[id] = p
        self.pairs[id].add_mate(rec)

    def get_pairs(self):
        return {str(p): p for p in self.pairs.values() if p.is_paired()}
