# cdom
Identify commonly mutated protein domains in a given [Mutation Annotation Format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) file.

The main script is mutated\_domains.R.

## Data
- [uniprot\_sprot\_human](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz)
- [protein2ipr](ftp://ftp.ebi.ac.uk/pub/databases/interpro/74.0/protein2ipr.dat.gz)

## Documentation
- [uniprot](https://web.expasy.org/docs/userman.html)

## Dependencies
Scripts to install (some) of the dependencies are located in the *install* directory.

- Python 3.x
  - pyreadr
  - pandas
  - argparse
- R
