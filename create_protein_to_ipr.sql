create table if not exists protein2ipr (uniprotswissprot, interpro, interproname, sourceid, domainstart INTEGER, domainstop INTEGER);
create index if not exists protein2ipr_uniprot on protein2ipr (uniprotswissprot);
create index if not exists protein2ipr_interpro on protein2ipr (interpro);
.separator "\t"
.import protein2ipr.dat protein2ipr
