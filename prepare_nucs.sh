#clean header from fungidb
sed -E '/^>/ s/>.*gene=/>/' *Proteins.fasta > renameprots.faa
sed '/^>/ s/|.*//' *prots.faa > nopipe.faa
for i in *names.txt; do seqtk subseq nopipe.faa $i > ${i//txt/subseq.faa}; done
EffectorP.py -i cl7downames.subseq.faa -o cl7downames.subseq.efp.txt -E cl7downames.subseq.E.fa
for i in cl*E.fa; do seqkit seq -n $i | awk '{print $1}' > ${i//.fa/_head.txt}; done
for i in *E_head.txt; do seqkit grep -n -r -f $i *Proteins.fasta > ${i//txt/prot.fa}; done
#activate entrez env
#first get the sequence identifiers for query
for i in *prot.fa
do
sed -n 's/^>//; s/[[:space:]].*$//p' $i > ${i//prot.fa/entrezq.txt}
done
#from entrezq protein ids, get bed files from the proteome sequence with 
IDS=$entrezq.txt
PROT=*Proteins.fasta
awk -v IDS="$IDS" '
BEGIN{
  # Load IDs (both with and without version suffix)
  while ((getline x < IDS) > 0) {
    sub(/\r$/,"",x)
    if (x != "") {
      ids[x]=1
      xb = x; sub(/\.[0-9]+$/, "", xb); ids[xb]=1
    }
  }
  close(IDS)
}
# Process only header lines
/^>/{
  hdr = $0

  # token just after ">" up to first space or pipe
  id = hdr; sub(/^>/,"",id)
  split(id, a, /[ \t|]/)
  idtok  = a[1]
  idbase = idtok; sub(/\.[0-9]+$/, "", idbase)

  # keep only if in your list
  if (!( (idtok in ids) || (idbase in ids) )) next

  # parse location=CHR:START-END(STRAND)
  if (match(hdr, /location=([^|]+)/, m)) {
    loc = m[1]
    split(loc, p, ":"); chr=p[1]; rest=p[2]

    str = "+"
    if (match(rest, /\(([+-])\)/, ms)) str = ms[1]
    gsub(/\([+-]\)/, "", rest)

    split(rest, q, "-"); s=q[1]+0; e=q[2]+0
    if (s > e) { t=s; s=e; e=t }   # normalize
    s0 = s - 1                     # BED: 0-based start

    # BED6: chrom  start0  end  name  score  strand
    print chr "\t" s0 "\t" e "\t" idtok "\t0\t" str
  }
}
' "$PROT" > coords.bed

The bed files are now usable for bedtools
samtools faidx *Genome.fasta
for i in *subseq.bed; do bedtools getfasta -fi *Genome.fasta -bed $i -s -name > ${i//bed/fa}; done
