#!/usr/bin/env python
import sys, argparse
import os, os.path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib.collections import PatchCollection
from matplotlib import gridspec
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

NF1_domains = [{
  'id': 'IPR001936',
  'start': 1187,
  'stop': 1557,
  'label': 'Ras GTPase-activating'
}, {
  'id': 'IPR001251',
  'start': 1580,
  'stop': 1738,
  'label': 'CAL-TRIO lipid binding'
}]

def get_n_samples(maf_df):
  return len(set(maf_df.loc[:,'Tumor_Sample_Barcode']))

def get_mutations(maf_df, hgnc_symbol='NF1'):
  """
  Parameters
  ----------
  maf_df : pd.DataFrame
    pandas dataframe parsed from a maftools tsv

  hgnc_symbol : str
    Gene symbol to create lollipop plot for
  """
  # protein coding mutations only
  muts = maf_df.loc[maf_df.loc[:,'Hugo_Symbol'] == hgnc_symbol]
  n_genetic_mutations = len(muts)
  muts = muts.loc[maf_df.loc[:, 'Protein_position'].notna()]
  n_protein_mutations = len(muts)
  
  length = None
  positions = []
  for mut_str in muts.loc[:,'Protein_position']:
    # e.g. 1590/2839
    # e.g. 2429-2431/2839
    # TODO update in mutated_domains.R for non-SNP mutations
    # maf_df$Protein_position_int = sapply(strsplit(maf_df$Protein_position, "/"), function(pair) { return(as.numeric(pair[1])) })
    words = mut_str.split('/')
    if length is None:
      length = int(words[1])
    pos = int(mut_str.split('/')[0].split('-')[0])
    positions.append(pos)

  return {
    'positions': positions,
    'length': length,
    'n_genetic_mutations': n_genetic_mutations,
    'n_protein_mutations': n_protein_mutations
  }

def plot(mutations_dict, domains=[], n_samples=None, outfile=None):
  """
  Parameters
  ----------
  mutations_dict : dict
    Return value of get_mutations
  """
  gs = gridspec.GridSpec(2, 1, height_ratios=[4,1])
  ax1 = plt.subplot(gs[0])
  ax2 = plt.subplot(gs[1])

  # create data
  y = np.zeros(mutations_dict['length'])
  for pos in mutations_dict['positions']:
    y[pos-1] += 1
  x = np.where(y != 0)[0]
  y = y[y != 0]

  # plot with no marker
  # values is a list with size equal to size of gene in bp
  ax1.stem(x,y)
   
  # plot domains
  TRACK_HEIGHT = 5
  TRACK_Y_ORIGIN = 0
  ax2.set_ylim(TRACK_Y_ORIGIN, TRACK_HEIGHT)
  rects = []
  for domain in domains:
    rect = Rectangle((domain['start'],TRACK_Y_ORIGIN), domain['stop']-domain['start']+1, TRACK_HEIGHT)
    rects.append(rect)

    # TODO annotations are messy
    #cx = (domain['start'] + domain['stop']) / 2
    #cy = (TRACK_Y_ORIGIN + TRACK_HEIGHT) / 2
    #ax2.annotate(domain['label'], (cx, cy), color='k', weight='bold', ha='center', va='center')
  pc = PatchCollection(rects) # facecolor=facecolor, alpha=alpha, edgecolor=edgecolor)
  ax2.add_collection(pc)

  # for shared axis
  ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
  ax2.set_xticks(ax1.get_xticks())
  ax2.set_xbound(ax1.get_xbound())

  # labels
  title_str = 'NF1 Protein Mutations'
  if n_samples is not None:
    title_str += ' (N = {})'.format(n_samples)
  ax1.set_title(title_str)
  ax1.set_ylabel('Number of mutations')
  ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
  ax2.set_xlabel('Position')

  if outfile is None:
    plt.savefig('out.png')
  else:
    plt.savefig(outfile)

def get_table(synapse_id):
  import synapseclient
  syn = synapseclient.Synapse()
  syn.login()
  results = syn.tableQuery("select * from %s" % synapse_id)
  results_df = results.asDataFrame()
  return results_df

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--infile', '-i', help='maftools .tsv file')
  parser.add_argument('--synapse-id', '-s', help="Synapse ID where a maftools table resides")
  parser.add_argument('--outdir', '-o', help="Directory to place results", required=True)
  args = parser.parse_args()

  if args.infile is None and args.synapse_id is None:
    sys.stderr.write("Either --infile or --synapse-id is required\n")
    sys.exit(2)
  if not args.infile is None and not args.synapse_id is None:
    sys.stderr.write("Exactly one of --infile or --synapse-id is required\n")
    sys.exit(3)

  maf_df = None
  if args.infile is not None:
    maf_df = pd.read_csv(args.infile, sep='\t')
  if args.synapse_id is not None:
    maf_df = get_table(args.synapse_id)
  if maf_df is None:
    sys.stderr.write("Error reading --infile or --synapse-id\n")
    sys.exit(4)

  # TODO get domain data for other genes to automate
  mutations = get_mutations(maf_df, hgnc_symbol="NF1")
  plot(mutations, NF1_domains, n_samples=get_n_samples(maf_df), outfile=os.path.join(args.outdir, 'lollipop.png'))
