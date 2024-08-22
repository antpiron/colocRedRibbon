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
Here, we use the real dataset "th", which is composed of a GWAS for type 2 diabetes (deriving from https://www.diagram-consortium.org/) and an eQTL colocalization study for the gene TH, encoding tyrosine hydroxylase (originating from http://tiger.bsc.es/). 

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
S3 object is built on a data.frame, here "th.dt", requiring the columns _id_, _a_, _b_ and _position_, representing the name of the SNP, the p-value of the first analysis, 
the p-value of the second analysis, and the position on the chromosome, respectively. Additionnal columns are used for computing the colocalisation.

The specified parameters are: <br/> 
*  ```risk``` which can be either _NULL_, _'a'_ or _'b'_. It represents the odd-ratios for the analysis we assigned to _a_ or _b_ column, (e.g. or.GWAS or or.eQTL). __Is it correct? in the code you wrote: risk the GWAS dataset with an odd-ratio (e.g. a.or or b.or) either NULL, 'a' or 'b'__
*  ```effect``` is an operator like `>=` or `<=` indicating the effect direction for __maybe compared to?__ the non-risk dataset __which is the non-risk dataset?__
*  ```columns``` is a named vector including the column names in the provided data.frame.
*  When ```shortlist``` is _TRUE_ it first runs _RedRibbon_ and then _coloc_ (default = _TRUE_).

To compute the colocalization __self-explanatory?__ only the colocRedRibbon object is required.

```R
## Run C. Wallace coloc(). What is Wallace coloc?
rrc.dec <- coloc(rrc.dec)
```
Optional parameters can be specified: <br/> 
*  ```n.reduce``` is a function to reduce the number of sample columns for _a_ and _b_ into a number (default = _max_). This is needed as some eQTL/GWAS toolchains output an effective number of samples by SNP. _a.n_ or _b.n_ parameters are used if not _NULL_. __not clear to me__
*  ```a.n``` is the number of __what?__ for _a_ (default = _NULL_)
*   ```b.n``` is the number of __what?__ for _b_ (default = _NULL_)
*   ```a.type``` can either be _quant_ or _cc_ to denote quantitative or case-control mode for _a_ (default = _quant_)
*   ```b.type``` can either be _quant_ or _cc_ to denote quantitative or case-control mode for _b_ (default = quant)
*   ```region.mode``` can be set to _"IQR"_, to fed coloc with the interquantile range region from RedRibbon overlap. Or it can be set to _"below"_ __then what?__ (default = _NULL_)

To extract the best SNP __what does best mean here?__ :

```R
coloc.res["bestSnp"]
## Is this correct command?
## Anything else that we can extract?
```
## Plotting the colocalisation 
The colocalisation can be plotted with the helper function _ggRedRibbonColoc_.

```R
ggRedRibbonColoc(rrc.dec, shortid = "TH")
```

```rrc.dec``` is the colocRedRibbon object and ```shortid``` indicates the name of the gene of interest. <br/>
Additional parameters that can be used include: <br/>
*  ```plot.order``` is a vector specifying the plot order (default = 1:4, where 1 = RedRibbon plot, 2 =  manhantan plot for _a_, 3 = plot for _a_ and _b_, 4 = manhantan plot for _b_)
*  ```show.title``` if TRUE shows the title of the plot (default = _TRUE_)
*  ```title``` the title of the plot
*  ```labels``` the axis labels (default = _NULL_)
*  ```tss``` is a numerical value used to indicate the position of transcription start site (default = _NULL_)
*  ```highlight``` provides a list of SNPs to highlight
*  ```.log10``` is used to output the log10 pval (default = FALSE)

# Citation

Please cite our [Life Science Alliance publication](https://doi.org/10.26508/lsa.202302203),

```text
Anthony Piron, Florian Szymczak, Theodora Papadopoulou, Maria Inês Alvelos, Matthieu Defrance, Tom Lenaerts, Décio L Eizirik, Miriam Cnop.
RedRibbon: A new rank-rank hypergeometric overlap for gene and transcript expression signatures.
Life Science Alliance. 2023 Dec 8;7(2):e202302203. doi: 10.26508/lsa.202302203. PMID: 38081640; PMCID: PMC10709657.
```
