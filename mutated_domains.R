#!/usr/bin/env Rscript
library(argparse)
library(RSQLite)
library(biomaRt)

count_mutations_table = function(db_file) {
  # NOTE this is only valid if no records are deleted (however the use here does not include deletions)
  con = dbConnect(RSQLite::SQLite(), db_file)
  res = dbGetQuery(con, 'SELECT MAX(_ROWID_) as count FROM mutations LIMIT 1;')
  print(res)
  dbDisconnect(con)
  rv = res$count
  if (is.na(rv)) {
    # then there are no records yet
    rv = 0
  }
  return(rv)
}

write_mutations_table = function(maf_obj, db_file) {
  # db_file - sqlite3 database with
  #   create table if not exists mutations (uniprotswissprot PRIMARY KEY, position);
  maf_df = maf_obj@data
  maf_df$Protein_position_int = sapply(strsplit(maf_df$Protein_position, "/"), function(pair) { return(as.numeric(pair[1])) })

  con = dbConnect(RSQLite::SQLite(), db_file)
  my_df = data.frame(uniprotswissprot = maf_df$SWISSPROT, position = maf_df$Protein_position_int)
  dbWriteTable(con, 'mutations', my_df, append=TRUE)
  dbDisconnect(con)
}

get_mutated_domains = function(db_file) {
  # db_file - sqlite3 database from protein2ipr data of InterPro
  #   CREATE TABLE protein2ipr (uniprotswissprot, interpro, interproname, sourceid, domainstart, domainstop);
  #   CREATE INDEX protein2ipr_uniprot on protein2ipr (uniprotswissprot);
  con = dbConnect(RSQLite::SQLite(), db_file)
  res = dbGetQuery(con, 'SELECT * FROM protein2ipr,mutations WHERE protein2ipr.uniprotswissprot = mutations.uniprotswissprot AND CAST(mutations.position AS INTEGER) BETWEEN CAST(protein2ipr.domainstart AS INTEGER) AND CAST(protein2ipr.domainstop AS INTEGER)')
  dbDisconnect(con)
  return(res)
}

#get_uniprot_human = function(outdir) {
#  # return a list of uniprot identifiers which exist in human
#  rds_path = file.path(outdir, 'uniprot_human.RDS')
#  uniprot_human = NA
#  if (file.exists(rds_path)) {
#    uniprot_human = readRDS(rds_path)
#  } else {
#    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#    uniprot_human = getBM(
#      attributes = c("uniprotswissprot"),
#      filters = c("with_uniprotswissprot"),
#      values = c(TRUE),
#      mart = ensembl
#    )
#    saveRDS(uniprot_human, rds_path)
#  }
#  return(uniprot_human)
#}

get_human_uniprot_interpro_counts = function(interpro_ids, db_file) {
  # collect uniprot and interpro data in one dataframe for mutated domains
  # interpro_ids - mutated interpro domains
  # db_file - sqlite3 database from protein2ipr data of InterPro
  # return a data frame with interpro,count
  con = dbConnect(RSQLite::SQLite(), db_file)
  interpro_ids_sql = paste0(sapply(interpro_ids, function(x) { paste0('"', x, '"') }), collapse=",")
  # TODO update query after manual experimentation
  res = dbGetQuery(con, sprintf('select interpro,uniprotswissprot from protein2ipr,uniprot_human where interpro in (%s) and protein2ipr.uniprotswissprot = uniprot_human.uniprotswissprot', interpro_ids_sql))
  return(res)
}

get_human_uniprot_count = function(db_file) {
  con = dbConnect(RSQLite::SQLite(), db_file)
  res = dbGetQuery(con, sprintf('select count(*) from (select distinct id from uniprot_human);'))
  return(res)
}

get_pr_dgp_sql = function() {
  return("
SELECT uniprot_human.id as id, protein2ipr.uniprotswissprot as uniprot, interpro, interproname, 
  min(domainstart) as domainstart_min, max(domainstop) as domainstop_max, uniprot_human.length as length, 
  (domainstop_max - domainstart_min + 1) / length as pr_dgp
FROM protein2ipr,uniprot_human
WHERE protein2ipr.uniprotswissprot = uniprot_human.uniprot
GROUP BY interpro;")
}
get_pr_dgp = function(db_file) {
  con = dbConnect(RSQLite::SQLite(), db_file)
  res = dbGetQuery(con, get_pr_dgp_sql())
  return(res)
}

get_pr_p_sql = function() {
  return("
SELECT id, uniprot, length / t1.total_length as pr_p
FROM uniprot_human, (
  SELECT sum(length) as total_length FROM uniprot_human
) as t1;")
}
get_pr_p = function(db_file) {
  con = dbConnect(RSQLite::SQLite(), db_file)
  res = dbGetQuery(con, get_pr_p_sql())
  return(res)
}

get_pr_d = function(db_file) {
  con = dbConnect(RSQLite::SQLite(), db_file)
  res = dbGetQuery(con, sprintf("
WITH pr_dgp_tbl AS ( %s ),
pr_p_tbl AS ( %s )
SELECT pr_dgp_tbl.interpro, pr_dgb_tbl.interproname, sum(pr_dgp_tbl.pr_dgp * pr_p_tbl.pr_p) as pr_d
FROM pr_dgp_tbl, pr_p_tbl
WHERE pr_dgp_tbl.id = pr_p_tbl.id
GROUP BY pr_dgb_tbl.interpro
", get_pr_dgp_sql(), get_pr_p_sql()))
  return(res)
}

main = function(args) {
  # TODO currently using an RData with specific name
  #df = readRDS(args$maftools_object)
  maf_obj = NA
  db_file = args$interpro_protein2ipr_db
  n_mutations_sql = count_mutations_table(db_file)
  print(n_mutations_sql)
  if (n_mutations_sql == 0) {
    load(args$maftools_object)
    rm(pNF.pfam)
    maf_obj = jhu_maf_file_pNF
    write_mutations_table(maf_obj, db_file)
  }

  # get domains for which there is a mutation
  # (each row is a mutation which occurs in a domain the domains may be overlapping, and in
  # that case multiple rows may correspond to one InterPro identifier)
  merged_db = get_mutated_domains(db_file)
  saveRDS(merged_db, 'merged_db.RDS')

  #n_mutations_in_domains = 0

  # count the number of human proteins with that domain
  mutated_interpro = unique(merged_db$interpro)
  interpro_count_df = get_human_uniprot_interpro_counts(mutated_interpro, db_file)

  # TODO report "significantly enriched" domains
  count_sum = sum(interpro_count_df$count)
  interpro_count_df$freq = interpro_count_df$count / count_sum

  
}

# TODO rename interpro-protein2ipr-db (more tables now)
parser = ArgumentParser()
parser$add_argument("--interpro-protein2ipr-db", "-i", help="InterPro protein2ipr SQLite database")
parser$add_argument("--maftools-object", "-m", help="maftools object @data as a .Rds; TODO")
parser$add_argument("--uniprot-human", "-u", help="File containing a newline-delimited list of Uniprot SwissProt identifiers which occur in humans")
parser$add_argument("--outdir", "-o", help="Output directory")
args = parser$parse_args()
main(args)
