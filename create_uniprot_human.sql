.bail on
create table uniprot_human(id, uniprot, length);
.separator "\t"
.import uniprot_parsed.tsv uniprot_human
create index uniprot_human_idx on uniprot_human(uniprot);
delete from uniprot_human where id = 'id';
