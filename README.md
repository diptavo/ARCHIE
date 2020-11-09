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
We have provided an example dataset to demonstrate the functionality of ARCHIE. The simulated dataset contains 100 independent (no LD) SNPs that are significantly (GWAS) associated with a phenotype and 1000 independent genes. The simulation setting was such that SNPs 1 to 5 are associated with genes 1 to 10; SNPs 11 to 15 are associated with genes 11 to 20. Three example input files for Sigma<sub>GE</sub>, Sigma<sub>GG</sub> and Sigma<sub>EE</sub> are [provided](https://github.com/diptavo/ARCHIE/tree/master/files).

Run ARCHIE using the following example command:

```shell

Rscript ~/archie/codes/run_archie.R --sigmage ~/archie/files/sigmage.txt --sigmagg ~/archie/files/sigmagg.txt --sigmaee ~/archie/files/sigmaee.txt

```
This will produce `out.rds`.

```R
a1 = readRDS("out.rds");
a1$q
0.07748791 0.02815841

head(a1$us,20)
      [,1] [,2]
 [1,]    0    1
 [2,]    0    1
 [3,]    0    1
 [4,]    0    1
 [5,]    0    1
 [6,]    0    0
 [7,]    0    0
 [8,]    0    0
 [9,]    0    0
[10,]    0    0
[11,]    1    0
[12,]    1    0
[13,]    1    0
[14,]    1    0
[15,]    1    0
[16,]    0    0
[17,]    0    0
[18,]    0    0
[19,]    0    0
[20,]    0    0

head(a1$vs,20)
      [,1] [,2]
 [1,]    0    1
 [2,]    0    1
 [3,]    0    1
 [4,]    0    1
 [5,]    0    1
 [6,]    0    1
 [7,]    0    1
 [8,]    0    1
 [9,]    0    1
[10,]    0    1
[11,]    1    0
[12,]    1    0
[13,]    1    0
[14,]    1    0
[15,]    1    0
[16,]    1    0
[17,]    1    0
[18,]    1    0
[19,]    1    0
[20,]    1    0



```
