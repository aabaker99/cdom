library(biomaRt)

# choice of mart limits results to human
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
rv = getBM(
  attributes = c("uniprotswissprot"),
  filters = c("with_uniprotswissprot"),
  values = c(TRUE),
  mart = ensembl
)
print(str(rv))
# TODO save; reports 19418 proteins
# however, current release of uniprot says 20431 human proteins
