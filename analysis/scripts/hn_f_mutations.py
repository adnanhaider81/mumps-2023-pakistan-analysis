#!/usr/bin/env python3
# Compare HN and F protein AAs with JL minor and major components
import argparse, os
from Bio import Entrez, SeqIO, pairwise2
from Bio.Seq import Seq

def fetch_gb(acc, email, api_key=None):
    Entrez.email = email
    if api_key: Entrez.api_key = api_key
    h = Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
    rec = SeqIO.read(h, 'genbank'); h.close()
    return rec

def find_cds_coords(rec, gene_name_substr):
    for f in rec.features:
        if f.type == 'CDS':
            gene = ' '.join(f.qualifiers.get('gene', [])).lower() + ' ' + ' '.join(f.qualifiers.get('product', [])).lower()
            if gene_name_substr in gene:
                return int(f.location.start), int(f.location.end), int(f.location.strand or 1)
    raise SystemExit(f'CDS not found for {gene_name_substr}')

ap = argparse.ArgumentParser(description='Report AA mutations in HN and F vs JL minor and major components')
ap.add_argument('--email', required=False)
ap.add_argument('--api_key', required=False)
ap.add_argument('--ref_acc', required=True, help='Genotype G reference with annotated HN and F, for example ON148331.1')
ap.add_argument('--sample', required=True, help='FASTA with one consensus')
ap.add_argument('--jl_minor', default='AAL83745.1', help='Protein accession of JL minor HN')
ap.add_argument('--jl_major', default='AAK83223.1', help='Protein accession of JL major HN')
ap.add_argument('--out_tsv', required=True)
a = ap.parse_args()

email = a.email or os.getenv('NCBI_EMAIL')
if not email: raise SystemExit('Set --email or env NCBI_EMAIL')
api_key = a.api_key or os.getenv('NCBI_API_KEY')

ref = fetch_gb(a.ref_acc, email, api_key)
cons = list(SeqIO.parse(a.sample, 'fasta'))[0].seq
aln = pairwise2.align.globalms(str(ref.seq), str(cons), 2, -1, -2, -1, one_alignment_only=True)[0]
ref_aln, cons_aln = aln.seqA, aln.seqB

def extract_region(ref_start, ref_end, strand=1):
    ref_pos = 0
    seg = []
    in_region = False
    for i, base in enumerate(ref_aln):
        if base != '-':
            ref_pos += 1
        if not in_region and ref_pos == ref_start + 1:
            in_region = True
        if in_region:
            seg.append(cons_aln[i])
            if ref_pos == ref_end:
                break
    seq = ''.join([b for b in seg if b != '-'])
    return Seq(seq) if strand == 1 else Seq(seq).reverse_complement()

hn_s, hn_e, hn_strand = find_cds_coords(ref, 'hemagglutinin-neuraminidase')
f_s, f_e, f_strand = find_cds_coords(ref, 'fusion')
hn_nt = extract_region(hn_s, hn_e, hn_strand)
f_nt = extract_region(f_s, f_e, f_strand)
hn_aa = hn_nt.translate()
f_aa = f_nt.translate()

with open(a.out_tsv, 'w') as out:
    out.write('gene	pos	ref	alt	against
')
    # Fetch JL proteins
    Entrez.email = email
    h1 = Entrez.efetch(db='protein', id=a.jl_minor, rettype='fasta', retmode='text')
    jl_minor_seq = ''.join(line.strip() for line in h1 if not line.startswith('>')); h1.close()
    h2 = Entrez.efetch(db='protein', id=a.jl_major, rettype='fasta', retmode='text')
    jl_major_seq = ''.join(line.strip() for line in h2 if not line.startswith('>')); h2.close()
    L = min(len(hn_aa), len(jl_minor_seq))
    for i in range(L):
        if hn_aa[i] != jl_minor_seq[i]:
            out.write(f'HN	{i+1}	{jl_minor_seq[i]}	{hn_aa[i]}	JL_minor
')
    L = min(len(hn_aa), len(jl_major_seq))
    for i in range(L):
        if hn_aa[i] != jl_major_seq[i]:
            out.write(f'HN	{i+1}	{jl_major_seq[i]}	{hn_aa[i]}	JL_major
')
    # F vs ref AA
    for f in ref.features:
        if f.type == 'CDS' and 'fusion' in ' '.join(f.qualifiers.get('product', [])).lower():
            ref_nt = f.extract(ref.seq); ref_f_aa = ref_nt.translate()
            L = min(len(f_aa), len(ref_f_aa))
            for i in range(L):
                if f_aa[i] != ref_f_aa[i]:
                    out.write(f'F	{i+1}	{ref_f_aa[i]}	{f_aa[i]}	Ref_{a.ref_acc}
')
            break
