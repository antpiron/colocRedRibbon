# Description

__colocRedRibbon__ is a [RedRibbon](https://github.com/antpiron/RedRibbon) based colocalization of large genome-wide association studies (GWAS) and expression quantitative trait locus (eQTL) analyses. colocRedRibbon was developped to pinpoint variants associated both with increased risk for a disease and gene expression variation. The method uses novel pre-filtering steps, shortlisting variants before applying colocalization analysis. 


# Installation

This procedure has been tested on debian/ubuntu but should work on any linux distribution.

## Directly from Github

In R, just run

```R
devtools::install_github("antpiron/colocRedRibbon")
```

or for `dev` branch,

```R
devtools::install_github("antpiron/colocRedRibbon", ref = "dev")
```

# Documentation

A [R vignette with fully reproducible examples is available here](https://antpiron.github.io/colocRedRibbon.html). You can also open it from R with

```R
vignette("colocRedRibbon")
```

## Short summary

```R
library(colocRedRibbon)
data("th", package = "colocRedRibbon")

## th.dt is a data.frame containing
## 
##        rsid pval.GWAS n.GWAS eaf.GWAS   or.GWAS ea.GWAS nea.GWAS pval.eQTL       pos n.eQTL zscore.eQTL ea.eQTL nea.eQTL eaf.eQTL
##   rs1003483     0.350 231420    0.510 1.0063199       T        G   0.68710   2167543    404       0.403       T        G    0.510
##   rs1003484     0.640 231420    0.260 0.9965061       A        G   0.92180   2167618    404      -0.098       A        G    0.260
##   rs1003889     0.990 187126    0.011 0.9993002       T        G   0.58720   1970108    317       0.543       T        G    0.011
## ...

## Create a new column indicating the direction of GWAS and assign the logarithms of the odds ratio of GWAS
th.dt$dir.GWAS <- log(th.dt$or.GWAS)

## Construct a colocRedRibbon S3 object
rrc.dec <- RedRibbonColoc(th.dt, risk="a", effect=`<=`,
                          columns=c(id="rsid", position="pos", a.type="cc", a="pval.GWAS", b="pval.eQTL",
                                    a.n="n.GWAS", a.eaf="eaf.GWAS", a.dir="dir.GWAS", 
                                    b.n="n.eQTL", b.eaf="eaf.eQTL", b.dir="zscore.eQTL"))
```

The plots are generated with

```R
ggRedRibbonColoc(rrc.dec, shortid = "TH")
```


# Citation

Please cite our [medXriv preprint](https://doi.org/10.1101/2024.10.19.24315808) and optionnaly, the [RedRibbon Life Science Alliance publication](https://doi.org/10.26508/lsa.202302203),

```text
Anthony Piron, Florian Szymczak, Lise Folon, Daniel J. M. Crouch, Theodora Papadopoulou, Maria Inês Alvelos, Maikel L. Colli, Xiaoyan Yi, Marcin Pekalski, Type 2 Diabetes Global Genomics Initiative, Matthieu Defrance, John A. Todd, Décio L. Eizirik, Josep M. Mercader, Miriam Cnop.
Identification of novel type 1 and type 2 diabetes genes by co-localization of human islet eQTL and GWAS variants with colocRedRibbon.
medRxiv 2024.10.19.24315808; doi: https://doi.org/10.1101/2024.10.19.24315808 

Anthony Piron, Florian Szymczak, Theodora Papadopoulou, Maria Inês Alvelos, Matthieu Defrance, Tom Lenaerts, Décio L Eizirik, Miriam Cnop.
RedRibbon: A new rank-rank hypergeometric overlap for gene and transcript expression signatures.
Life Science Alliance. 2023 Dec 8;7(2):e202302203. doi: 10.26508/lsa.202302203. PMID: 38081640; PMCID: PMC10709657.
```
