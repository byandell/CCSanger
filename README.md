---
title: "Process Sanger data for CC Founders"
author: "Brian S Yandell"
date: "8 December 2016"
output: html_document
---

These files are getting features for use with `DOQTL` data. These are derived from
code developed by Karl Broman.

* `R/get_mgi_features.R`: get MGI features (genes, exons, ...) for Sanger file
* `inst/CreatSQL/vcf_snp_2db.R`: (re)create `cc_foundersnps.sqlite` (update of Karl Broman's `R/0_vcf2db.R`)
* `inst/CreatSQL/vcf_indel_2db.R`: create `cc_founderindels.sqlite` (new)
* `inst/CreatSQL/svs.Rmd`: create `svs8_*.rds` files (new)

The `vcf` and `svs` files have hard-wired `dirpath` that needs to be locally edited.

### Required packages

For MGI features:

```
library(dplyr)
library(readr)
source("check_interval.R")
```

For VCF:

```
library(VariantAnnotation)
library(RSQLite)
```

For SVS:

```
library(dplyr)
```

### SQLite Files

The MGI features are in their own SQLite. This is much smaller. Not clear that it is needed this way.

I modified Karl's SNP SQLite code to include consequence and Ensembl IDs.
The InDel SQLite code is a minor change to get indels. Note that there is a column labelled `allele`, which is the name used in creation of the VCF, but this is pretty close to the column `alleles`, which we use in DOQTL stuff.

### SVS Files

The structural variants were small, but unclear what users might want. I mainly use the `svs8_len.rds` file.
These are in a different format from VCF. Look carefully at this. Not totally happy, but it works.
