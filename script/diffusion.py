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

def parse_string_fh(fh, threshold=0.0, score='combined_score'):
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

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--data', '-d', help="Tabular data file or synapse id with maftools output", required=True)
  parser.add_argument('--delimiter', default='\t', help='Field delimiter of <--data>; default = <TAB>')
  parser.add_argument('--outdir', '-o', help="Directory to place results", required=True)
  parser.add_argument('--synapse-id', '-s', help='synapse.org folder to upload results to')
  args = parser.parse_args()
  
  write_provenance(args.outdir)

  sys.stderr.write("[status] Downloading and parsing STRINGdb... ")
  stringdb_fpath = download_stringdb(args.outdir)
  G = None
  with open(stringdb_fpath, 'r') as fh:
    G = parse_string_fh(fh)
  nodelist = sorted(G.nodes())
  sys.stderr.write("done reading {} nodes and {} edges\n".format(G.order(), G.size()))

  df = read_tsv_or_synapse_id(args.data, args.outdir, args.delimiter)

  # ignore variants that are not mapped to an Ensembl protein identifier
  # TODO other row filters based on variant frequency
  df_filt = df[pd.notna(df['ENSP']) & df['ENSP'].str.contains('ENSP')]

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

  bn = os.path.basename(args.data)
  bn, ext = os.path.splitext(bn)
  ofp = os.path.join(args.outdir, '{}_diffused.tsv'.format(bn))
  df_out = pd.DataFrame(mat_diffused.todense(), index=sample_set, columns=nodelist)
  df_out.index.name = SAMPLE_COLUMN
  df_out.to_csv(ofp, sep='\t')

  # TODO
  #upload_to_synapse(args.outdir)
