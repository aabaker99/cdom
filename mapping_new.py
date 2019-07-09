#!/usr/bin/env python
import sys, argparse

def main(args):

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--maf-object", "-p", help='maftools output @data data frame written as a tsv', required=True)
  parser.add_argument("--outdir", "-o", help="Output directory", required=True)
  args = parser.parse_args()
  main(args)
