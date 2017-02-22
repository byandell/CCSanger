---
title: "Derived Data for R/doqtl2"
author: "Brian S Yandell"
date: "February 22, 2017"
output: html_document
---

These files are for converting the `DOQTL` data (or rather, that curated by Karl Broman with `R/qtl2` packages) for use with the `R/qtl2` suite. There are some conversions and some new files.

* [vcf_snp_2db.R](inst/derived_data/vcf_snp_2db.R): (re)create `cc_foundersnps.sqlite` (update of `R/0_vcf2db.R`)
* [vcf_indel_2db.R](inst/derived_data/vcf_indel_2db.R): create `cc_founderindels.sqlite` (new)
* [svs.Rmd](inst/derived_data/vcf_svs_2csv.R): create `svs8_*.rds` files (new)

The derived data now go into a folder called `DerivedData`. This is in a
hard-wired folder that is, for now, machine specific. On my laptop it is:

```
dirpath <- file.path("~/Documents/Research/attie_alan/DO", "data")
datapath <- file.path(dirpath, "DerivedData")
```

The files in this `derived_data` folder have this hard-wired. The only other place where it is hard-wired is for the shiny server, which is in
`inst/shiny/setup.R`. These will have to be changed for local installation.
It would be nice to have this in only one place. Next iteration.

### RDS Files

RDS files are apparently an improvement on RData files. They hold one R object and are the same size, but are quicker to use. See for instance [R bloggers](http://www.r-bloggers.com/a-better-way-of-saving-and-loading-objects-in-r/).

### SQLite Files

I modified the SNP SQLite code to include consequence and Ensembl IDs.
The InDel SQLite code is a minor change to get indels. Note that there is a column labelled `allele`, which is the name used in creation of the VCF, but this is pretty close to the column `alleles`, which we use in DOQTL stuff.

### SVS Files

These are in a different format. Look carefully at this. Not totally happy, but it works.
