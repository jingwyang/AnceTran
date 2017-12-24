---
title: "README"
author: "Jingw Yang"
date: "2017/12/24"
output: md_document
---

# AnceTran

*AnceTran* is an *R* package that performs analyses of TF-binding evolution
from *ChIP-seq* data, including optimized input formatting, normalization
and pair-wise distance evaluation, TF-binding character tree construction and ancestral
TF binding state inference

*AnceTran* package is under active developing, current developing version 1.0 is available at <https://github.com/jingwyang/AnceTran>.

A convenient way to install package from github is through *devtools* package:

```{r, eval=FALSE}
install.packages('devtools')
devtools::install_github("jingwyang/AnceTran")
```

Users can also download *AnceTran* package and install locally through:

```{r, eval=FALSE}
install.packages("filePath/AnceTran.1.0.tar.gz", repos = NUll, type = "source")
```
