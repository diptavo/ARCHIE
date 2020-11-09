# ARCHIE (Aggregative trans-association analysis to identify trait-specific target genes)
Manuscript: Novel Aggregative trans-eQTL Association Analysis of Known Genetic Variants Detect Trait-specific Target Gene-sets. Diptavo Dutta, Yuan He, Ashis Saha, Marios Arvanitis, Alexis Battle, Nilanjan Chatterjee.

Preprint can be found on [medrxiv](https://www.medrxiv.org/content/10.1101/2020.09.29.20204388v2)

# Installation and pre-requisites
The codes in this repository are writted in R. To run them, please download and install the following R packages: `optparse` and `data.table`.

# Running instructions

We have provided an example run of ARCHIE for demonstration purposes. The individual application and customization would critically depend on the particular application. We encourage the user to closely follow the outlined algorithm in the [preprint](https://www.medrxiv.org/content/10.1101/2020.09.29.20204388v2) to adapt it to their cases.

The repo can be cloned as:
```
git clone https://github.com/diptavo/ARCHIE/
```
We have provided an example dataset to demonstrate the functionality of ARCHIE. The simulated dataset contains 100 independent (no LD) SNPs that are significantly (GWAS) associated with a phenotype and 1000 independent genes. The simulation setting was such that SNPs 1 to 5 are associated with genes 1 to 10; SNPs 11 to 15 are associated with genes 11 to 20. Three example input files for $`\Sigma_{GE}`$, $`\Sigma_{GG}`$ and $`\Sigma_{EE}`$ are [provided](https://github.com/diptavo/ARCHIE/files/).
By applying ARCHIE, we aim to map them to their distal targets.

```shell

Rscript ~/archie/codes/run_archie.R --sigmage ~/archie/files/sigmage.txt --sigmagg ~/archie/files/sigmagg.txt --sigmaee ~/archie/files/sigmaee.txt

```
