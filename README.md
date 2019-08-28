# cdom
Identify commonly mutated protein domains in a given [Mutation Annotation Format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) file.

## scripts/mutated\_domains.R
Analyzes mutation data in the MAF format.
This script is applied to analyze mutations in patients with NF and to analyze mutations in the general population according to the 1000 genomes project.
In either case, it reports the most frequently mutated protein domains.

## scripts/diffusion.py
Diffuse genetic mutations over a protein interaction network as done in [Network-Based Stratification](https://doi.org/10.1038/nmeth.2651).
Genetic variants are filtered to include only those which have a rare allele frequency and a damaging PolyPhen prediction.

## scripts/filter_protein_to_ipr.py
Make the InterPro protein2ipr data more managble by excluding all non-human proteins.

## scripts/parse_uniprot.py
Called by filter_protein_to_ipr.py

## scripts/nf_lollipop.py
Plot the frequency of amino acid mutations for a protein of interest.
Currently only analyzes the NF1 protein.

## Data
- uniprot\_sprot\_human at ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz
- protein2ipr at ftp://ftp.ebi.ac.uk/pub/databases/interpro/74.0/protein2ipr.dat.gz

## Documentation
- [uniprot](https://web.expasy.org/docs/userman.html)

## Dependencies
Scripts to install (some) of the dependencies are located in the *install* directory.

- Python 3.x
  - pyreadr
  - pandas
  - argparse
- R

## TODO
- Finalize application of mutated_domains.R to 1000 genomes data
- Program to synthesize frequently mutated domains in disease and control populations (which domains are mutated in disease state significantly more than the control?)
- Check if subsetting InterPro's protein2ipr data to include only human proteins is a sufficient reduction in the size of data to be able to perform in-memory joins of InterPro and mutation data in R (an obviate the need for SQLite).
- Generalize nf_lollipop.py to genes other than NF1
