#!/bin/sh

USAGE="usage: download_1000genomes_vcf.sh <OUTDIR>"
if [ -z "$1" ]; then
  echo $USAGE >&2
  exit 1
fi

export OUTDIR="$1"
( 
  cd $OUTDIR
  wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/*.vcf.gz*'
  wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/README_GRCh38_liftover_20170504.txt'
)
