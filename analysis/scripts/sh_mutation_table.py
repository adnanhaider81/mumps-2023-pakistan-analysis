#!/usr/bin/env python3
# Compare SH AA vs Sheffield and JL
import argparse, os
from Bio import Entrez, SeqIO

ap = argparse.ArgumentParser(description='Report AA mutations in SH vs Sheffield and Jeryl Lynn')
ap.add_argument('--email', required=False)
ap.add_argument('--api_key', required=False)
ap.add_argument('--sheffield', default='ON148331.1')
ap.add_argument('--jl_sh', default='FN431985.1')
ap.add_argument('--sample_sh', required=True, help='FASTA with SH sequence only')
ap.add_argument('--out_tsv', required=True)
a = ap.parse_args()

email = a.email or os.getenv('NCBI_EMAIL')
if not email: raise SystemExit('Set --email or env NCBI_EMAIL')
api_key = a.api_key or os.getenv('NCBI_API_KEY')
Entrez.email = email
if api_key: Entrez.api_key = api_key

def fetch_nt(acc):
    from Bio import SeqIO, Entrez
    h = Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
    rec = SeqIO.read(h, 'genbank'); h.close()
    return rec

def get_sh_nt(rec):
    for f in rec.features:
        if f.type == 'CDS':
            gene = ' '.join(f.qualifiers.get('gene', [])).upper()
            product = ' '.join(f.qualifiers.get('product', [])).lower()
            if 'SH' in gene or 'small hydrophobic' in product:
                return f.extract(rec.seq)
    raise SystemExit('SH not found in reference')

shef = fetch_nt(a.sheffield)
jl = fetch_nt(a.jl_sh)
shef_sh = get_sh_nt(shef).translate()
jl_sh = get_sh_nt(jl).translate()
sample = list(SeqIO.parse(a.sample_sh, 'fasta'))[0].seq.translate()

with open(a.out_tsv, 'w') as out:
    out.write('against	pos	ref	alt
')
    L = min(len(sample), len(shef_sh))
    for i in range(L):
        if sample[i] != shef_sh[i]:
            out.write(f'Sheffield	{i+1}	{shef_sh[i]}	{sample[i]}
')
    L = min(len(sample), len(jl_sh))
    for i in range(L):
        if sample[i] != jl_sh[i]:
            out.write(f'JerylLynn	{i+1}	{jl_sh[i]}	{sample[i]}
')
