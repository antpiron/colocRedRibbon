# Description

colocRedRibbon: a [RedRibbon](https://github.com/antpiron/RedRibbon) based colocalization of GWAS and eQTL.

# Installation

This procedure have been tested on debian/ubuntu but should work on any linux distribution.

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

## Create colocRedRibbon S3 object
rrc.dec <- RedRibbonColoc(th.dt, risk="a", effect=`<=`,
                          columns=c(id="rsid", position="pos", a="pval.GWAS", b="pval.eQTL",
                                    a.n="n.GWAS", a.eaf="eaf.GWAS", a.or="or.GWAS", 
                                    b.n="n.eQTL", b.eaf="eaf.eQTL", b.beta="zscore.eQTL"))
		
## Run C. Wallace coloc()
rrc.dec <- coloc(rrc.dec)

ggRedRibbonColoc(rrc.dec, shortid = "TH")

```


# Citation

Please cite our [Life Science Alliance publication](https://doi.org/10.26508/lsa.202302203),

```text
Anthony Piron, Florian Szymczak, Theodora Papadopoulou, Maria Inês Alvelos, Matthieu Defrance, Tom Lenaerts, Décio L Eizirik, Miriam Cnop.
RedRibbon: A new rank-rank hypergeometric overlap for gene and transcript expression signatures.
Life Science Alliance. 2023 Dec 8;7(2):e202302203. doi: 10.26508/lsa.202302203. PMID: 38081640; PMCID: PMC10709657.
```
