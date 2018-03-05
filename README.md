---
title: "README"
author: "Jingw Yang"
date: "2018/01/30"
output: md_document
---

# AnceTran

*AnceTran* is an *R* package that performs analyses of transcriptome evolution
based on *RNA-seq* expression data or *ChIP-seq* TF binding data, including optimized input formatting, normalization
and pair-wise distance evaluation, transcriptome character tree construction and ancestral
transcriptome state inference.

*AnceTran* package is under active developing, current developing version 2.0 is available at <https://github.com/jingwyang/AnceTran>.

A convenient way to install package from github is through *devtools* package:

```{r, eval=FALSE}
install.packages('devtools')
devtools::install_github("jingwyang/AnceTran")
```

Users can also download *AnceTran* package and install locally through:

```{r, eval=FALSE}
install.packages("filePath/AnceTran.2.0.tar.gz", repos = NUll, type = "source")
```
