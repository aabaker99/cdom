.bail on
create table protein2ipr (uniprotswissprot, interpro, interproname, sourceid, domainstart INTEGER, domainstop INTEGER);
create index protein2ipr_uniprot on protein2ipr (uniprotswissprot);
create index protein2ipr_interpro on protein2ipr (interpro);
.separator "\t"
.import ../data/protein2ipr.dat protein2ipr
