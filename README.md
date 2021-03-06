----------
# fgga-lnc
----------
FGGA-lnc is a graph-based Machine Learning approach for the automatic and consistent GO annotation from lncRNA sequences.

---------
# Download data used in paper
---------

Data is hosted on Datasets file 

--------------
# Requirements
-------------
The following software and libraries must be installed on your machine:

[R](https://cran.r-project.org/) Rscript: tested with version 3.6.3

[mfold](http://www.unafold.org/) package: tested with version 2.3

[jsonlite](https://cran.r-project.org/web/packages/jsonlite/) library: tested with version 1.6.1

[seqinr](https://cran.r-project.org/web/packages/seqinr/) library: tested with version 3.6-1

------------
# How to use
------------

```
lncRNA_SS2("path","name-seq")
```

Where name-seq is the name of the sequence such as URS000195C635_9615 and path is the directory where the predicted secondary structure is stored such as /home/user/temp

----------
# Authors
----------
Flavio E. Spetale

Javier Murillo
