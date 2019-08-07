library(biomaRt)
library(maftools)

args = commandArgs(trailingOnly=TRUE)
load(args[1])

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = c("ensembl_transcript_id")
attributes = c('ensembl_transcript_id', 'uniprotswissprot')
values = unique(jhu_maf_file_pNF@data$Transcript_ID)

rv = getBM(
  attributes = attributes,
  filters = filters,
  values = values,
  mart = ensembl
)

write.csv(rv, args[2])
