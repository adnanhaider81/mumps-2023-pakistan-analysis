#!/usr/bin/env python3
import argparse, os
from Bio import Entrez, SeqIO
ap = argparse.ArgumentParser(description='Fetch single GenBank FASTA by accession')
ap.add_argument('--email', required=False, help='Entrez email, or set env NCBI_EMAIL')
ap.add_argument('--api_key', required=False, help='NCBI API key, or set env NCBI_API_KEY')
ap.add_argument('--acc', required=True)
ap.add_argument('--out_fasta', required=True)
a = ap.parse_args()

email = a.email or os.getenv('NCBI_EMAIL')
if not email:
    raise SystemExit('Set --email or env NCBI_EMAIL')
api_key = a.api_key or os.getenv('NCBI_API_KEY')
Entrez.email = email
if api_key:
    Entrez.api_key = api_key

h = Entrez.efetch(db='nucleotide', id=a.acc, rettype='fasta', retmode='text')
open(a.out_fasta, 'w').write(h.read())
h.close()
print('Wrote', a.out_fasta)
