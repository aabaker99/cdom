#!/usr/bin/env python
import sys, argparse
import os, os.path
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument('--protein2ipr', '-p', required=True, help="InterPro protein2ipr file")
parser.add_argument('--uniprot-human', '-u', required=True, help="Uniprot taxonomic divison, human database file")
parser.add_argument('--outdir', '-o', required=True, help="Directory to place results")
args = parser.parse_args()

parse_uniprot_ofp = os.path.join(args.outdir, 'uniprot_parsed.tsv')
if not os.path.exists(parse_uniprot_ofp):
  sp.check_call(['parse_uniprot.py', '-i', args.uniprot_human, '-o', args.outdir])

uniprot_human_set = set()
with open(parse_uniprot_ofp) as fh:
  for line in fh:
    line = line.rstrip()
    id_v, uniprot, length = line.split('\t')
    uniprot_human_set.add(uniprot)

ofp = os.path.join(args.outdir, 'protein2ipr_human_filtered.tsv')
with open(ofp, 'w') as ofh:
  with open(args.protein2ipr) as ifh:
    for line in ifh:
      line = line.rstrip()
      words = line.split('\t')
      uniprot = words[0]
      if uniprot in uniprot_human_set:
        ofh.write(line)
        ofh.write("\n")
