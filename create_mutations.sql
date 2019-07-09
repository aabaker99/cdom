create table if not exists mutations (uniprotswissprot, position);
create index if not exists mutations_key on mutations (uniprotswissprot, position);
