#!/usr/bin/env python
import sys, argparse
import pandas as pd
import xml.etree.ElementTree as ET
import copy
import os, os.path
import itertools as it

# maftools-output:
# write.table(pNF.pfam$domainSummary, 'maftools_domain_summary.tsv', sep='\t', row.names=FALSE)

def transform_cdd_accession(merged_tbl):
  """
  CDD and InterPro treat PFAM and other identifiers differently. An example from CDD is "pfam00001"
  while an example from InterPro is "PF00051". The InterPro version is correct. Transform CDD
  identifiers to match.

  CHL is a NCBI Protein Clusters database identifier for curated chloroplasts (https://www.ncbi.nlm.nih.gov/genomes/static/proteinclusters_help.html)
  COG is Clusters of Orthologous Groups (COG) (https://www.ncbi.nlm.nih.gov/COG/)
  """
  rv = pd.DataFrame(data=merged_tbl, copy=True)
  for index, series in merged_tbl.iterrows():
    # skip NaN
    old_accession = series['cd-accession']
    if pd.isna(old_accession):
      continue

    new_accession = None
    if 'pfam' in old_accession:
      new_accession = old_accession.replace('pfam', 'PF')
    elif 'smart' in old_accession:
      new_accession = old_accession.replace('smart', 'SM')
    else:
      new_accession = copy.copy(old_accession)

    rv.loc[index, 'cd-accession'] = new_accession
  return rv

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--maftools-output", required=True, help='maftools protein domain summary table')
  parser.add_argument("--cdd-file", required=True, help='cddid_all.tbl from NCBI Conserved Domain Database')
  parser.add_argument('--interpro-file', required=True, help='interpro.xml from InterPro')
  parser.add_argument("--outdir", required=True, help="Output directory")
  args = parser.parse_args()

  maf_tbl = pd.read_csv(args.maftools_output, sep='\t')
  # rename columns as column names are not included as part of the CDD-distributed file
  cdd_tbl = pd.read_csv(args.cdd_file, sep='\t', names=['pssm-id', 'cd-accession', 'cd-name', 'cd-description', 'pssm-length'])

  # TODO this is a one-to-many mapping
  merged_tbl = maf_tbl.merge(cdd_tbl, left_on="DomainLabel", right_on="cd-name", how='left')
  merged_tbl = transform_cdd_accession(merged_tbl)
  merged_tbl.to_csv(os.path.join(args.outdir, 'cdd_join.csv'), index=False, sep='\t')

  # find cd-accession in interpro xml
  # construct mapping of member database id to interpro id
  mdb_to_ip = {}
  et = ET.parse(args.interpro_file)
  for interpro_elem in et.iterfind('interpro'):
    for xref_elem in interpro_elem.iterfind('member_list/db_xref'):
      tpl = mdb_to_ip.get(xref_elem.attrib['dbkey'])
      if tpl is not None:
        (ipd, ipd_count) = tpl
        sys.stderr.write("[warning] overwriting mapping for db_xref key = {} from {} to {}\n".format(kk, mdb_to_ip[ipd], interpro_elem.attrib['id']))
      mdb_to_ip[xref_elem.attrib['dbkey']] = (interpro_elem.attrib['id'], interpro_elem.attrib['protein_count'])

  # report missing identifier mapping
  missing = []
  for cd_accession in merged_tbl['cd-accession']:
    tpl = mdb_to_ip.get(cd_accession)
    if tpl is None:
      missing.append(str(cd_accession))
  if len(missing) != 0:
    sys.stderr.write("[warning] missing InterPro mapping for CDD accession:\n")
    sys.stderr.write('\n'.join(sorted(missing)))

  ll = []
  for k,v in mdb_to_ip.items():
    ll.append((k, v[0], v[1]))
  interpro_tbl = pd.DataFrame(ll, columns=['cd-accession', 'interpro-accession', 'interpro-protein-count'])
  merged_tbl = merged_tbl.merge(interpro_tbl, left_on="cd-accession", right_on="cd-accession", how="left")
  merged_tbl.to_csv(os.path.join(args.outdir, 'table.csv'), index=False)
  #merged_tbl.to_csv(os.path.join(args.outdir, 'table.csv'), columns=['DomainLabel', 'nMuts', 'nGenes', 'interpro-accession', 'interpro-protein-count'], index=False)
