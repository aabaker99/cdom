#!/bin/sh
cat uniprot_sprot_human.dat | grep -E '^AC' | sed 's/^AC\W\W*//' | awk -F';' '{ for (i=1;i<NF;i++) { print $i }}' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' | sort | uniq > uniprot_sprot_human.txt
