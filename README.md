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

vignette("RedRibbon")



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
These parameters are specified as follows: <br/> 
*  ```risk```: This parameter can be either _'a'_ or _'b'_. It represents the odds-ratios for the analysis assigned to the _a_ or _b_ column, such as or.GWAS or or.eQTL. __Is it correct?__
*  ```effect```: Specifies an operator like `>=` or `<=` to indicate the effect direction for __maybe compared to?__ the non-risk dataset __which is the non-risk dataset?__
*  ```columns```: A named vector that includes the column names present in the provided data.table.
*  ```shortlist```: When this parameter is set to _TRUE_ the function first runs _RedRibbon_ and then _coloc_ (default = _TRUE_).

```R
## Run C. Wallace coloc(). 
rrc.dec <- coloc(rrc.dec)
```
Several optional parameters can be specified to customize the behavior of the _colocRedRibbon_ object: <br/> 
*  ```n.reduce```: This parameter accepts a function that reduces the number of sample columns for _a_ and _b_ to a single value (default = _max_). This is useful because some eQTL/GWAS toolchains report an effective number of samples per SNP. If _a.n_ or _b.n_ are provided and not _NULL_, they will be used instead.  __not clear to me__
*  ```a.n```: Represents the number of samples for _a_. The default value is NULL, meaning it will use the effective number of samples if not specified. __correct?__
*   ```b.n```: Represents the number of samples for _b_. Similar to _a.n_, the default value is NULL.
*   ```a.type```: This parameter can be set to either _quant_ or _cc_ to denote whether the analysis for _a_ is in quantitative mode or case-control mode, respectively. The default value is _quant_.
*   ```b.type```: This parameter functions similarly to _a.type_, allowing you to specify _quant_ or _cc_ for the analysis of _b_. The default value is also _quant_.
*   ```region.mode```: This parameter can be set to _"IQR"_ to feed coloc with the interquartile range region derived from the overlap in _RedRibbon_. Alternatively, it can be set to _"below"_ (default is _NULL_), which specifies a different region mode. __then what?__ 

```R
#Extract the best SNP
coloc.res["bestSnp"]
## Is this command correct?
## Anything else that we can extract?
```
## Plotting the colocalization 

```R
ggRedRibbonColoc(rrc.dec, shortid = "TH")
```

_ggRedRibbonColoc_ requires to specify the colocRedRibbon object containing the colocalization data, here _rrc.dec_, and to indicate the name of the gene of interest by assigning it to the parameter ```shortid```. <br/>

Additional optional parameters to customize the visualization of colocalization results: <br/>

*  ```plot.order```: A vector that determines the order of the plots (default is 1:4, where 1 represents the RedRibbon plot, 2 is the Manhattan plot for _a_, 3 is the plot for _a_ and _b_, and 4 is the Manhattan plot for _b_) 
*  ```show.title```: If set to _TRUE_, displays the title of the plot (default = _TRUE_)
*  ```title```: Specifies the title of the plot
*  ```labels```: Provides axis labels (default = _NULL_)
*  ```tss```: A numerical value used to indicate the position of the transcription start site (default = _NULL_)
*  ```highlight```: Accepts a list of SNPs to be highlighted in the plot.
*  ```.log10```: If set to _TRUE_, outputs the log10 p-values (default = _FALSE_)

# colocRedRibbon workflow
colocRedRibbon is a method designed to identify common causal candidates by examining the colocalization of GWAS and eQTL. This approach aims to pinpoint variants that are linked to both disease risk and variation in gene expression. <br/>

The method employs a two-step approach for shortlisting variants: 
-  __Risk allele effect step:__  In this step, variants are categorized into two distinct groups based on their direction of effect on gene expression, i.e., down- or upregulating. The upregulating variant set comprises the variants whose risk alleles increase gene expression, and the downregulating set risk alleles that decrease gene expression. Each of these variant sets is analyzed independently in subsequent steps.
-  __RedRibbon Overlap Step:__ In this step, the RedRibbon rank-rank hypergeometric overlap method is applied on both GWAS and eQTL variants, which are ranked according to their P-values. This analysis examines the potential overlap between the two ranked lists. If a significant overlap is detected by RedRibbon, these shortlisted SNPs are further analyzed by the coloc package. If no significant overlap is found, the coloc method is still applied to the two effect sets without the preliminary overlap shortlisting. <br/>
<p align="center">
<img src="https://github.com/user-attachments/assets/6c1e4cb5-d331-428e-9002-fec2e2e6976f" width="2500" height="250">
</p>

# FAQ

## When to use the IQR mode?
The area between the 25th percentile (first quartile) and the 75th percentile (third quartile) of the chromosomal positions of the core set is referred to as the interquartile range (IQR).
Using the IQR mode allows for a more focused analysis by including only those variants that are relevant to the central distribution of the core set. This mode is particularly beneficial when you want to minimize the influence of outliers and concentrate on the most relevant variants in your analysis. <br/>

<p align="center">
<img src="https://github.com/user-attachments/assets/dcb95536-9a13-4de4-85b4-31d2a88a2c1e" width="700" height="200">
</p>

# Citation

Please cite our [Life Science Alliance publication](https://doi.org/10.26508/lsa.202302203),

```text
Anthony Piron, Florian Szymczak, Theodora Papadopoulou, Maria Inês Alvelos, Matthieu Defrance, Tom Lenaerts, Décio L Eizirik, Miriam Cnop.
RedRibbon: A new rank-rank hypergeometric overlap for gene and transcript expression signatures.
Life Science Alliance. 2023 Dec 8;7(2):e202302203. doi: 10.26508/lsa.202302203. PMID: 38081640; PMCID: PMC10709657.
```
