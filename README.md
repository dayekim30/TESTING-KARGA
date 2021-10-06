# TESTING-KARGA


### Description
K-mer-based Antibiotic Resistance Gene Variant Analyzer (KARGVA). KARGVA is a Multi-platform Toolkit for Identification of Antibiotic Resistance from Sequencing Data Conferred by Point Mutations in Bacterial Genes.

### How to install
Copy the git-link on Visual Studio or other IDEs and run this program with the release version.

### How to run
1. Before running the program, locate the fasta and fastq files in the same directory where this program class files exist.

2. This program takes 3 inputs for running in the order: 1. the number of K for K-mer  2. The name of FASTA file   3. The name of FASTQ file.

3. Put 1 if the detailed information about the depth of hit is required, otherwise, 0 will should be the input for the following question "Do you need mapped genes information in detail? Yes -> 1, No ->0" 

4. The result data will be stored in the same directory where this program class files exist.

### Output
KARGVA_mappedGenes.csv

KARGVA_mappedReads.csv
