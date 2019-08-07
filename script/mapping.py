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
  def transform(old):
    new = None
    if pd.isna(old):
      new = copy.copy(old)
    else:
      if 'pfam' in old:
        new = old.replace('pfam', 'PF')
      elif 'smart' in old:
        new = old.replace('smart', 'SM')
      else:
        new = copy.copy(old)
    return new
  for index, series in merged_tbl.iterrows():
    rv.loc[index, 'cd-accession'] = transform(series['cd-accession'])
    rv.loc[index, 'pfam'] = transform(series['pfam'])
  return rv

def main(args):
  protein_tbl = pd.read_csv(args.protein_summary, sep='\t')
  domain_tbl = pd.read_csv(args.domain_summary, sep='\t')
  # rename columns as column names are not included as part of the CDD-distributed file
  cdd_tbl = pd.read_csv(args.cdd_file, sep='\t', names=['pssm-id', 'cd-accession', 'cd-name', 'cd-description', 'pssm-length'])

  # TODO this is a one-to-many mapping
  merged_tbl = domain_tbl.merge(cdd_tbl, left_on="DomainLabel", right_on="cd-name", how='left')
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

  # save interpro table
  mdb_to_ip_items = sorted(mdb_to_ip.items(), key=lambda mdb_ip: mdb_ip[1][0])
  with open(os.path.join(args.outdir, 'mdb_to_ip.tsv'), 'w') as ofh:
    for mdb_ip in mdb_to_ip_items:
      ofh.write('\t'.join([mdb_ip[1][0], mdb_ip[0]]) + '\n')

  ll = []
  for k,v in mdb_to_ip.items():
    ll.append((k, v[0], v[1]))
  interpro_tbl = pd.DataFrame(ll, columns=['cd-accession', 'interpro-accession', 'interpro-protein-count'])
  merged_tbl = merged_tbl.merge(interpro_tbl, left_on="cd-accession", right_on="cd-accession", how="left")

  # try to fill in domains with the pfam column and the interpro table
  n_row, n_col = merged_tbl.shape
  for i in range(n_row):
    ip_tpl = mdb_to_ip.get(merged_tbl.loc[i, 'pfam'])
    if ip_tpl is not None:
      merged_tbl.loc[i, 'interpro-accession'] = ip_tpl[0]

  merged_tbl = merged_tbl.sort_values(['DomainLabel', 'interpro-accession'])
  merged_tbl = merged_tbl.drop_duplicates(subset='DomainLabel', keep='last')

  # report missing identifier mapping
  missing = []
  for cd_accession in merged_tbl['cd-accession']:
    tpl = mdb_to_ip.get(cd_accession)
    if tpl is None:
      missing.append(str(cd_accession))
  if len(missing) != 0:
    sys.stderr.write("[warning] missing InterPro mapping for CDD accession:\n")
    sys.stderr.write('\n'.join(sorted(missing)))

  # save table
  merged_tbl = merged_tbl.sort_values(['nMuts'], ascending=False)
  merged_tbl.to_csv(os.path.join(args.outdir, 'table.csv'), index=False, sep='\t', columns=['DomainLabel', 'interpro-accession', 'nMuts', 'interpro-protein-count'])
  #merged_tbl.to_csv(os.path.join(args.outdir, 'table.csv'), columns=['DomainLabel', 'nMuts', 'nGenes', 'interpro-accession', 'interpro-protein-count'], index=False)

  # map mutation to domain
  all_tbl = protein_tbl.merge(merged_tbl, left_on='DomainLabel', right_on='DomainLabel', how='left')

  return all_tbl

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  # TODO is there some reason pyreadr cant read our particular .RData file? maybe we can work around? probably just port this to R?
  parser.add_argument("--protein-summary", required=True, help='maftools protein domain summary table')
  parser.add_argument("--domain-summary", required=True, help='maftools protein domain summary table')
  parser.add_argument("--cdd-file", required=True, help='cddid_all.tbl from NCBI Conserved Domain Database')
  parser.add_argument('--interpro-file', required=True, help='interpro.xml from InterPro')
  parser.add_argument("--outdir", required=True, help="Output directory")
  args = parser.parse_args()
  all_tbl = main(args)
