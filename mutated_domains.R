#!/usr/bin/env Rscript
library(argparse)
library(RSQLite)
library(biomaRt)

get_mutated_domains = function(maf_obj, db_file, outdir) {
  # maf_obj - output of maftools
  # db_file - sqlite3 database from protein2ipr data of InterPro
  #   CREATE TABLE protein2ipr (uniprotswissprot, interpro, interproname, sourceid, domainstart, domainstop);
  #   CREATE INDEX protein2ipr_uniprot on protein2ipr (uniprotswissprot);
  maf_df = maf_obj@data
  maf_df$Protein_position_int = sapply(strsplit(maf_df$Protein_position, "/"), function(pair) { return(as.numeric(pair[1])) })

  uniprot_ids = unique(maf_df$SWISSPROT)
  uniprot_ids_sql = paste0(sapply(uniprot_ids, function(x) { paste0('"', x, '"') }), collapse=",")

  con = dbConnect(RSQLite::SQLite(), db_file)
  res = dbSendQuery(con, sprintf('select * from protein2ipr where uniprotswissprot in (%s)', uniprot_ids_sql))
  res_fetch = dbFetch(res)

  merged_db = merge(x = maf_df, y = res_fetch, by.x = 'SWISSPROT', by.y = 'uniprotswissprot', all.x = TRUE)
  merge_fail = merged_db[is.na(merged_db$domainstart),]
  write(sprintf("Failed to merge the following Uniprot identifiers (no domain information): %s", paste(unique(merge_fail$SWISSPROT), collapse=", ")), stderr())
  merged_db = merged_db[!is.na(merged_db$domainstart),]
  merged_db = merged_db[merged_db$domainstart <= merged_db$Protein_position_int & merged_db$Protein_position_int <= merged_db$domainstop,]

  return(merged_db)
  #saveRDS(merged_db, 'merged_db.RDS')
}

get_uniprot_human = function(outdir) {
  # return a list of uniprot identifiers which exist in human
  rds_path = file.path(outdir, 'uniprot_human.RDS')
  uniprot_human = NA
  if (file.exists(rds_path)) {
    uniprot_human = readRDS(rds_path)
  } else {
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    uniprot_human = getBM(
      attributes = c("uniprotswissprot"),
      filters = c("with_uniprotswissprot"),
      values = c(TRUE),
      mart = ensembl
    )
    saveRDS(uniprot_human, rds_path)
  }
  return(uniprot_human)
}

get_multispecies_uniprot = function(interpro_ids, db_file) {
  # collect uniprot and interpro data in one dataframe for mutated domains
  # interpro_ids - mutated interpro domains
  # db_file - sqlite3 database from protein2ipr data of InterPro
  con = dbConnect(RSQLite::SQLite(), db_file)
  interpro_ids_sql = paste0(sapply(interpro_ids, function(x) { paste0('"', x, '"') }), collapse=",")
  res = dbSendQuery(con, sprintf('select * from protein2ipr where interpro in (%s)', interpro_ids_sql))
  res_fetch = dbFetch(res)
}

main = function(args) {
  # TODO currently using an RData with specific name
  #df = readRDS(args$maftools_object)
  load(args$maftools_object)
  rm(pNF.pfam)

  # get domains for which there is a mutation
  # (each row is a mutation which occurs in a domain the domains may be overlapping, and in
  # that case multiple rows may correspond to one InterPro identifier)
  merged_db = get_mutated_domains(jhu_maf_file_pNF, args$interpro_protein2ipr_db, args$outdir)
  mutated_interpro = unique(merged_db$interpro)
  # TODO
  saveRDS(merged_db, 'merged_db.RDS')

  # count the number of human proteins with that domain
  protein2ipr_mut_filt_df = get_multispecies_uniprot(mutated_interpro, args$interpro_protein2ipr_db)
  uniprot_human = get_uniprot_human(args$outdir)
  mutated_interpro_counts = sapply(mutated_interpro, function(interpro) {
    multispecies_df = protein2ipr_mut_filt_df[protein2ipr_mut_filt_df$interpro == interpro,]
    human_df = multispecies_df[multispecies_df$uniprotswissprot %in% uniprot_human,]
    return(dim(human_df)[1])
  })
  saveRDS(mutated_interpro, file.path(args$outdir, 'mutated_interpro.RDS'))
  saveRDS(mutated_interpro_counts, file.path(args$outdir, 'mutated_interpro_counts.RDS'))
  interpro_counts_df = data.frame(interpro = mutated_interpro, count = mutated_interpro_counts)
  saveRDS(interpro_counts_df, file.path(args$outdir, 'interpro_counts_df.RDS'))

  # TODO report "significantly enriched" domains
}

parser = ArgumentParser()
parser$add_argument("--interpro-protein2ipr-db", "-i", help="InterPro protein2ipr SQLite database")
parser$add_argument("--maftools-object", "-m", help="maftools object @data as a .Rds; TODO")
parser$add_argument("--outdir", "-o", help="Output directory")
args = parser$parse_args()
main(args)
