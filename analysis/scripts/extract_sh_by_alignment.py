#!/usr/bin/env python3
# Extract 316 nt SH region from a sample consensus by aligning to a reference that has SH CDS annotated
import argparse, os
from pathlib import Path
from Bio import Entrez, SeqIO, pairwise2
from Bio.Seq import Seq

def fetch_ref_gb(acc, email, api_key=None):
    Entrez.email = email
    if api_key: Entrez.api_key = api_key
    h = Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
    rec = SeqIO.read(h, 'genbank'); h.close()
    return rec

def find_sh_coords(rec):
    for f in rec.features:
        if f.type == 'CDS':
            gene = ' '.join(f.qualifiers.get('gene', [])).upper()
            product = ' '.join(f.qualifiers.get('product', [])).lower()
            if 'SH' in gene or 'small hydrophobic' in product:
                return int(f.location.start), int(f.location.end), int(f.location.strand or 1)
    raise SystemExit('Could not find SH CDS in reference')

ap = argparse.ArgumentParser(description='Extract SH region by aligning sample consensus to a reference')
ap.add_argument('--email', required=False)
ap.add_argument('--api_key', required=False)
ap.add_argument('--ref_acc', required=True, help='GenBank accession with annotated SH CDS, for example AF280799.1')
ap.add_argument('--sample', required=True, help='FASTA with one consensus')
ap.add_argument('--out_fa', required=True)
a = ap.parse_args()

email = a.email or os.getenv('NCBI_EMAIL')
if not email: raise SystemExit('Set --email or env NCBI_EMAIL')
api_key = a.api_key or os.getenv('NCBI_API_KEY')

ref = fetch_ref_gb(a.ref_acc, email, api_key)
s, e, strand = find_sh_coords(ref)

cons = list(SeqIO.parse(a.sample, 'fasta'))[0].seq
# Align full genomes
aln = pairwise2.align.globalms(str(ref.seq), str(cons), 2, -1, -2, -1, one_alignment_only=True)[0]
ref_aln, cons_aln = aln.seqA, aln.seqB

# Map SH coordinates through alignment
ref_pos = 0
sh_start_idx = None
sh_end_idx = None
for i, base in enumerate(ref_aln):
    if base != '-':
        ref_pos += 1
    if ref_pos == s + 1 and sh_start_idx is None:
        sh_start_idx = i
    if ref_pos == e:
        sh_end_idx = i
        break
seg = []
for i in range(sh_start_idx, sh_end_idx + 1):
    if cons_aln[i] != '-':
        seg.append(cons_aln[i])
sh_seq = ''.join(seg)
with open(a.out_fa, 'w') as out:
    out.write(f">SH_from_{Path(a.sample).stem}\n{sh_seq}\n")
print('Wrote', a.out_fa)
