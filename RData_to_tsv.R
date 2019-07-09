requireNamespace('argparse', quietly=TRUE)

main = function() {
  parser = argparse::ArgumentParser()
  parser$add_argument("--rdata")
  parser$add_argument("--outdir")
  args = parser$parse_args()
  # Loads jhu_maf_file_pNF and pNF.pfam
  load(args$rdata)

  # TODO error if not loaded
  print(ls())

  write.table(pNF.pfam$proteinSummary, row.names=FALSE, file=file.path(args$outdir, 'proteinSummary.tsv'), sep='\t')
  write.table(pNF.pfam$domainSummary, row.names=FALSE, file=file.path(args$outdir, 'domainSummary.tsv'), sep='\t')
}

results = main()
