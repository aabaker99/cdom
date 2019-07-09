create table if not exists uniprot_human (uniprotswissprot PRIMARY KEY);
.separator "\t"
.import uniprot_sprot_human.txt uniprot_human
