#!/usr/bin/env python
import os, os.path
import re
import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--infile', '-i', required=True)
parser.add_argument('--outdir', '-o', required=True)
args = parser.parse_args()

id_regexp = re.compile('^ID\s+(\S+)\s+\S+\s+(\d+)\sAA.*$')
ac_regexp = re.compile('^AC(.*)$')

rec = None
line_no = 0
with open(os.path.join(args.outdir, 'uniprot_parsed.tsv'), 'w') as ofh:
  ofh.write('{}\t{}\t{}\n'.format('id', 'uniprot', 'length'))
  with open(args.infile, 'r') as fh:
    for line in fh:
      line_no += 1
      line = line.rstrip()

      match_data = id_regexp.match(line)
      if match_data is not None:
        if rec is not None:
          # write the previous record
          if ('length' in rec and 'uniprots' in rec):
            for uniprot in rec['uniprots']:
              ofh.write('{}\t{}\t{}\n'.format(rec['id'], uniprot, rec['length']))
          else:
            sys.stderr.write('Cannot write record with id = {}\n'.format(rec['id']))

        # start new record
        rec = {}
        id_line = line
        rec['id'] = match_data.group(1)
        length = int(match_data.group(2))
        rec['length'] = length

      match_data = ac_regexp.match(line)
      if match_data is not None:
        ac_line = line
        uniprot_ids = list(filter(lambda x: len(x) != 0, map(lambda x: x.strip(), match_data.group(1).split(';'))))
        rec['uniprots'] = uniprot_ids
