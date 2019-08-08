#!/usr/bin/env python
import sys, argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib import gridspec
from matplotlib.patches import Rectangle

# NF1
length = 2839
domains = [{
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
    pos = int(str1.split('/')[0].split('-')[0])
    positions.append(pos)

  return {
    'positions': positions,
    'length': length,
    'n_genetic_mutations': n_genetic_mutations,
    'n_protein_mutations': n_protein_mutations
  }

def plot(mutations_dict):
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
  # array of SNP position
  values = np.random.uniform(size=40)

  # plot with no marker
  # values is a list with size equal to size of gene in bp
  ax1.stem(values)
   
  def plot_domain(domain_name, domain_start, domain_stop):
    TRACK_HEIGHT = 5
    TRACK_Y_ORIGIN = 0
    rect = Rectangle((domain_start,TRACK_Y_ORIGIN), domain_stop, TRACK_HEIGHT)
    cx = (domain_start + domain_stop) / 2
    cy = (TRACK_Y_ORIGIN + TRACK_HEIGHT) / 2

    pc = PatchCollection([rect]) # facecolor=facecolor, alpha=alpha, edgecolor=edgecolor)
    ax2.set_ylim(TRACK_Y_ORIGIN, TRACK_HEIGHT)
    ax2.add_collection(pc)
    ax2.annotate(domain_name, (cx, cy), color='w', weight='bold', ha='center', va='center')

  # for shared axis
  ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
  ax2.set_xticks(ax1.get_xticks())
  ax2.set_xbound(ax1.get_xbound())

  # labels
  ax1.set_title('NF1 Lollipop Plot')
  ax1.set_ylabel('Number of mutations')
  ax2.set_xlabel('Position')

  plt.savefig('out.png')

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

  mutations = get_mutations(maf_df, hgnc_symbol="NF1")
  plot(mutations)

