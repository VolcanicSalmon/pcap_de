import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from typing import Dict, Optional

class GFF3:
    def __init__(self, seqid: str, source: str, feat: str, start: int, end: int,
                 score: str = '.', strand: str = '.', phase: str = '.', attrs: str = '.'):
        self.seqid = seqid
        self.source = source
        self.feat = feat
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attrs = attrs
        self.parsed_attrs = None

    def parse_attrs(self) -> Dict[str, str]:
        if self.parsed_attrs is None:
            self.parsed_attrs = {}
            for attr in self.attrs.split(';'):
                if '=' in attr:
                    key, val = attr.strip().split('=', 1)
                    self.parsed_attrs[key] = val
        return self.parsed_attrs

    def get_attr(self, key: str, default: Optional[str] = None) -> Optional[str]:
        return self.parse_attrs().get(key, default)

    def extract_promoter(self, sequences: Dict[str, SeqRecord], promoter_size: int = 2000) -> Optional[Seq]:
        if self.seqid not in sequences:
            return None
        genome_seq = sequences[self.seqid].seq
        seq_len = len(genome_seq)
        if self.strand == '+':
            promoter_start = max(0, self.start - 1 - promoter_size)
            promoter_end = self.start - 1
        elif self.strand == '-':
            promoter_start = self.end
            promoter_end = min(seq_len, self.end + promoter_size)
        else:
            return None
        promoter_seq = genome_seq[promoter_start:promoter_end]
        if self.strand == '-':
            promoter_seq = promoter_seq.reverse_complement()
        return promoter_seq

if __name__ == "__main__":
    
    gff_path, fasta_path, output_path = sys.argv[1], sys.argv[2], sys.argv[3]
    promoter_size = 2000

    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    features = []
    with open(gff_path, "r") as gff_file:
        for line in gff_file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                feat = GFF3(*fields[:8], attrs=fields[8])
                features.append(feat)

    with open(output_path, "w") as out_fasta:
        for feat in features:
            promoter = feat.extract_promoter(genome, promoter_size)
            if promoter:
                feat_id = feat.get_attr('ID') or feat.get_attr('Parent') or f"{feat.feat}_{feat.start}_{feat.end}"
                record = SeqRecord(promoter, id=feat_id, description=f"Promoter for {feat.feat}")
                SeqIO.write(record, out_fasta, "fasta")
