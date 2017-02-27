---
title: "Derived Data for R/doqtl2"
author: "Brian S Yandell"
date: "June 12, 2016"
output: html_document
---

These files are for converting the `DOQTL` data (or rather, that curated by Karl Broman with `R/qtl2` packages) for use with the `R/doqtl2` package. There are some conversions and some new files.

* [convertRDS.Rmd](inst/derived_data/convertRDS.Rmd): convert DOQTL RData and CSV to RDS

The derived data go into a folder called `DerivedData`. This is in a
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
