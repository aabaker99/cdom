#!/bin/sh

sqlite3 ../data/cdom.db < create_protein_to_ipr.sql
sqlite3 ../data/cdom.db < create_uniprot_human.sql
sqlite3 ../data/cdom.db < create_mutations.sql
