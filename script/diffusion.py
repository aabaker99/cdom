#!/usr/bin/env python
import argparse, sys
import subprocess as su
import os, os.path
import networkx as nx
import numpy as np
import scipy.sparse as sp
import pandas as pd
from scipy.sparse.linalg import norm

SAMPLE_COLUMN = 'Tumor_Sample_Barcode'
PROTEIN_COLUMN = 'ENSP'
ALLELE_FREQ_COLUMN = 'ExAC_AF'
POLYPHEN_COLUMN = "PolyPhen"
POLYPHEN_DAMAGING = 'damaging'

def write_provenance(outdir):
  # TODO write version as a bash comment?
  fpath = os.path.join(outdir, 'provenance.sh')
  with open(fpath, 'w') as fh:
    fh.write(" ".join(sys.argv))

def download_stringdb(outdir):
  """
  Download STRING database file.

  Returns
  -------
  fpath : str
    filepath of STRINGdb file

  Raises
  ------
  TODO an error from sp if download fails
  """
  url = 'https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz'
  fname = os.path.basename(url)
  fpath = os.path.join(outdir, fname)
  fpath_unzip, ext = os.path.splitext(fpath)
  if not os.path.exists(fpath_unzip):
    if not os.path.exists(fpath):
      args = ['wget', '-O', fpath, url]
      su.check_call(args)

    args = ['gunzip', fpath]
    su.check_call(args)
  return fpath_unzip

def parse_string_fh(fh, threshold=0.0, score='combined_score', threshold_as_percentile=False):
  """
  Parameters
  ----------
  fh : file-like
    STRING db database file

  threshold : float
    edge weight confidence threshold expressed as a value in [0,1]
    include all edges with a confidence greater than or equal to the threshold

  score : str
    the name of a column in the string format to use as the edge confidence score
    default: combined_score
    Can be "neighborhood" "fusion" "cooccurence" "coexpression" "experimental" "database" 
    "textmining" or "combined_score"

  threshold_as_percentile : bool
    if True, interpret the threshold has a percentile value: if threshold is 0.9, include
    only the top 10% scoring interactions

  Returns
  -------
  G : nx.Graph
    protein-protein interaction network with Ensembl protein ids as node ids
  """
  # skip header
  line_no = 1
  header = fh.__next__()
  columns = header.strip().split()
  if score not in columns[2:]:
    # TODO exception class
    raise Exception("Invalid score column {}; must be one of {}".format(score, ", ".join(columns[2:])))
  score_index = columns.index(score)

  G = nx.Graph()
  weights = []
  if threshold_as_percentile:
    for line in fh:
      line_no += 1
      p1, p2, conf = parse_line_abc(line, score_index=score_index)
      G.add_edge(p1, p2, **{'weight': conf})
      weights.append(conf)

    weights = np.array(weights)
    weights_per = np.percentile(weights, threshold*100)
    to_remove = []
    for u, v, attrs in G.edges.data():
      weight = attrs['weight']
      if weight <= weights_per:
        to_remove.append((u,v))
    for u,v in to_remove:
      G.remove_edge(u,v)
    
  else:
    for line in fh:
      line_no += 1
      p1, p2, conf = parse_line_abc(line, score_index=score_index)
      if(conf >= threshold):
        G.add_edge(p1, p2, **{'weight': conf})
  return G

def parse_line_abc(line, score_index=-1):
  """
  Filter STRING data to just protein,protein,confidence and remove taxonomy code from protein identifiers
  """
  def trim_taxonomy_code(protein_id):
    return protein_id[5:]

  line = line.rstrip()
  words = line.split()
  p1 = words[0]
  p2 = words[1]
  # remove human taxonomy code prefix
  p1 = trim_taxonomy_code(p1)
  p2 = trim_taxonomy_code(p2)
  score = words[score_index]
  conf = int(score) / 1000.0
  return (p1, p2, conf)

def diffusion(M, adj, alpha=0.7, tol=10e-6):  # TODO equation, M, alpha
    """
    Network propagation iterative process

    Iterative algorithm for apply propagation using random walk on a network:
        Initialize::
            X1 = M

        Repeat::
            X2 = alpha * X1.A + (1-alpha) * M
            X1 = X2

        Until::
            norm(X2-X1) < tol

        Where::
            A : degree-normalized adjacency matrix

    Parameters
    ----------
    M : sparse matrix
        Data matrix to be diffused.

    adj : sparse matrix
        Adjacency matrice.

    alpha : float, default: 0.7
        Diffusion/propagation factor with 0 <= alpha <= 1.
        For alpha = 0 : no diffusion.
        For alpha = 1 :

    tol : float, default: 10e-6
        Convergence threshold.

    Returns
    -------
    X2 : sparse matrix
        Smoothed matrix.

    Notes
    -----
    Copied from the stratipy Python library
    """
    n = adj.shape[0]
    adj = adj+sp.eye(n)

    d = sp.dia_matrix((np.array(adj.sum(axis=0))**-1, [0]), shape=(n,  n))
    A = adj.dot(d)

    X1 = M
    X2 = alpha * X1.dot(A) + (1-alpha) * M
    i = 0
    while norm(X2-X1) > tol:
        X1 = X2
        X2 = alpha * X1.dot(A) + (1-alpha) * M
        i += 1
    return X2

def read_tsv_or_synapse_id(tsv_or_synapse_id, outdir, delimiter='\t'):
  df = None
  if tsv_or_synapse_id[:3] == 'syn':
    synapse_id = tsv_or_synapse_id
    ofp = os.path.join(outdir, '{}.tsv'.format(synapse_id))
    if os.path.exists(ofp):
      df = pd.read_csv(ofp, sep='\t')
    else:
      import synapseclient
      syn = synapseclient.Synapse()
      syn.login()
      results = syn.tableQuery('SELECT * FROM {}'.format(synapse_id))
      df = results.asDataFrame()
      df.to_csv(ofp, sep='\t') # TODO naming of tsv/csv with sep='\t' or ',', etc.
  else:
    tsv = tsv_or_synapse_id
    df = pd.read_csv(tsv, sep=delimiter)
  return df

def upload_to_synapse(outdir):
  pass

def embed_arr(all_col_names, some_col_names, arr):
  """
  Embed <arr> with len(some_col_names) number of columns in a larger array with len(all_col_names)
  number of columns. The larger array is constructed and returned.
  """
  m,n = arr.shape
  if len(some_col_names) != n:
    raise Exception('some_col_names != #columns of arr: {} != {}'.format(some_col_names, n))

  n2 = len(all_col_names)
  rv = np.zeros((m,n2))

  all_name_to_ind = {}
  for i,name in enumerate(all_col_names):
    all_name_to_ind[name] = i

  failed_cols = set()
  for i in range(m):
    for j, name in enumerate(some_col_names):
      j2 = all_name_to_ind.get(name)
      if j2 is None:
        failed_cols.add(name)
      else:
        rv[i,j2] = arr[i,j]

  return rv, failed_cols

# TODO 
def filter_variants_by_expression(df):
  pass

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="""
Diffuse genetic variant data over the STRINGdb protein-protein interaction network
""")
  parser.add_argument('--variant-data', '-d', help="Tabular data file or synapse id with maftools output", required=True)
  parser.add_argument('--expression-data', '-e', help="Tabular data or synapse id with RNA expression levels for variant filtering, e.g. syn20449214")
  parser.add_argument('--join-column', '-j', help='If --expression-data is provided, the column shared by both datasets to join on')
  parser.add_argument('--delimiter', default='\t', help='Field delimiter of <--data>; default = <TAB>')
  parser.add_argument('--outdir', '-o', help="Directory to place results", required=True)
  parser.add_argument('--synapse-id', '-s', help='synapse.org folder to upload results to')
  parser.add_argument('--threshold', '-t', help='Edge weight threshold as a percentile in decimal e.g. 0.9 for 90th percentile', type=float)
  args = parser.parse_args()
  
  write_provenance(args.outdir)

  sys.stderr.write("[status] Downloading and parsing STRINGdb... ")
  stringdb_fpath = download_stringdb(args.outdir)
  G = None
  with open(stringdb_fpath, 'r') as fh:
    if args.threshold is None:
      G = parse_string_fh(fh)
    else:
      G = parse_string_fh(fh, threshold=args.threshold, threshold_as_percentile=True)
  nodelist = sorted(G.nodes())
  sys.stderr.write("done reading {} nodes and {} edges\n".format(G.order(), G.size()))

  df = read_tsv_or_synapse_id(args.variant_data, args.outdir, args.delimiter)

  # variant filtering - {{
  # ignore variants that are not mapped to an Ensembl protein identifier
  n_variants = df.shape[0]
  df_filt = df[pd.notna(df[PROTEIN_COLUMN]) & df[PROTEIN_COLUMN].str.contains(PROTEIN_COLUMN)]
  n_variants_f1 = df_filt.shape[0]
  sys.stderr.write("[status] Filtered out {} variants without ENSP identifiers, {} remain\n".format(n_variants - n_variants_f1, n_variants_f1))

  # filter out common variants
  if ALLELE_FREQ_COLUMN not in df_filt.columns:
    sys.stderr.write("[error] Cannot filter out common variants because alelle frequency column {} is missing\n".format(ALLELE_FREQ_COLUMN))
    sys.exit(21)
  df_filt = df_filt[df_filt[ALLELE_FREQ_COLUMN] < 0.01]
  n_variants_f2 = df_filt.shape[0]
  sys.stderr.write("[status] Filtered out {} variants with allele frequency >= 0.01, {} remain\n".format(n_variants_f1 - n_variants_f2, n_variants_f2))

  # filter out benign variants according to PolyPhen (keep those for which we have no call from PolyPhen)
  df_filt = df_filt[pd.isna(df_filt[POLYPHEN_COLUMN]) | df_filt[POLYPHEN_COLUMN].str.contains(POLYPHEN_DAMAGING)]
  n_variants_f3 = df_filt.shape[0]
  sys.stderr.write("[status] Filtered out {} variants with benign PolyPhen, {} remain\n".format(n_variants_f2 - n_variants_f3, n_variants_f3))

  # TODO
  # filter out weakly expressed genes
  #if args.expression_data is not None:
    # TODO potentially different delimiter
    #df_exp = read_tsv_or_synapse_id(args.expression_data, args.outdir, args.delimiter)
    # map HGNC -> ENSP

  # }} - variant filtering

  # group by patient
  # count the number of mutations in each ENSP
  sample_set = sorted(set(df_filt[SAMPLE_COLUMN]))
  row_vecs = []
  all_failed_cols = set()
  for sample in sample_set:
    df_sample = df_filt[df_filt[SAMPLE_COLUMN] == sample][[SAMPLE_COLUMN, PROTEIN_COLUMN]] 
    df_sample_group_count = df_sample.groupby(PROTEIN_COLUMN).count()
    # TODO unclear how to rename the count column
    row_vec_sample = df_sample_group_count[SAMPLE_COLUMN]
    m = row_vec_sample.shape[0]
    row_vec_embed, failed_cols = embed_arr(nodelist, list(row_vec_sample.index), row_vec_sample.to_numpy().reshape(1,m))
    row_vecs.append(row_vec_embed)
    all_failed_cols = all_failed_cols.union(failed_cols)
  sys.stderr.write("Unable to diffuse the following proteins {} because they are absent from STRINGdb: ".format(len(all_failed_cols)) + ", ".join(sorted(all_failed_cols)) + "\n")

  mat = sp.csc_matrix(np.concatenate(row_vecs, axis=0))
  adj = nx.to_scipy_sparse_matrix(G, nodelist=nodelist, dtype=bool)

  # free up memory
  df = None
  df_filt = None
  row_vecs = None
  G = None

  mat_diffused = diffusion(mat, adj)

  bn = os.path.basename(args.variant_data)
  bn, ext = os.path.splitext(bn)
  ofp = os.path.join(args.outdir, '{}_diffused.tsv'.format(bn))
  df_out = pd.DataFrame(mat_diffused.todense(), index=sample_set, columns=nodelist)
  df_out.index.name = SAMPLE_COLUMN
  df_out.to_csv(ofp, sep='\t')

  # TODO
  #upload_to_synapse(args.outdir)
