# 1. Introduction
dbCNV is a tool to predict the pathogenicity for five-tier classification and binary classification of CNVs based on the deleterious significance of features. The quantitative evaluation of features was based on their pathogenicity levels in CNVs of different classifications.  According to the deleterious significance, we formulated quantitative methods for features, which fall into two categories: the first is variable type, including maximum, minimum and mean; the second is attribute type, which is measured by numerical sum. The reference genome version used by dbCNV is GRCh37/hg19.
# 2. Requirements
The local version dbCNV requires Linux (version: 3.10.0-514.el7.x86_64) environment, perl language (version:5.28.3) and R language (version: 3.6.0). Four perl packages (List::Util qw/sum min max/, List::MoreUtils qw/uniq/, Statistics::Sequences::Runs 0.10, Statistics::Distributions) and five R packages (xgboost, caret, ggplot2, lattice, Matrix). This pipeline runs dbCNV.pl in the same directory as Database.
# 3. Installation

```
git clone https://github.com/lllllv-1/dbCNV
```

# 4. Usage and example
# Usage:

```
perl dbCNV.pl -i prefix.txt -n 2 -r /usr/bin/Rscript
perl dbCNV.pl -i prefix.txt -n 5 -r /usr/bin/Rscript
```

  -i File    CNV region file (Chr Start End Type)
  
  -n Number  the number of disease type (2 or 5)
  
  -r File    the file of Rscript path
  
  -h Help    help information

 
# Example:

```
perl dbCNV.pl -i example.txt -n 2 -r /usr/bin/Rscript
```

The results can be seen in the gain_2_predication_result.txt and loss_2_predication_result.txt.
 
# 5. Input & output
Input file format (the header is required):

Chr Start End Type

chr10	260000	700000	gain

Column 1: The chromosome

Column 2: Start

Column 3: End

Column 4: CNV type (gain or loss)

The output file has 5 columns including chromosome, Start, End, CNV type and pathogenicity of CNV.

Chr     Start   End     Type    Disease

chr10	260000	700000	gain  benign
