#!/usr/bin/env python3
import argparse, csv, sys, os
ap = argparse.ArgumentParser()
ap.add_argument('--blast_tsv', required=True, help='Output from blast_closest.py')
ap.add_argument('--preferred', required=False, help='Optional newline-separated accession whitelist to prefer')
ap.add_argument('--out_accession', required=True)
a = ap.parse_args()

preferred = set()
if a.preferred and os.path.exists(a.preferred):
    preferred = set(x.strip() for x in open(a.preferred) if x.strip())

best = None
best_score = -1.0
with open(a.blast_tsv) as f:
    for row in csv.reader(f, delimiter='\t'):
        if not row: 
            continue
        qseqid, sacc, pident, length, bitscore, stax, stitle = row[:7]
        score = float(bitscore)
        if preferred and sacc not in preferred:
            score -= 1.0
        if score > best_score:
            best_score = score
            best = sacc
if not best:
    sys.exit('No BLAST hits found')
open(a.out_accession, 'w').write(best + '\n')
print('Selected', best)
