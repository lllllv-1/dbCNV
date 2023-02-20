# 1. Introduction
dbCNV is a tool to predict the pathogenicity for five-tier classification and binary classification of CNVs based on the deleterious significance of features. The quantitative evaluation of features was based on their pathogenicity levels in CNVs of different classifications.  According to the deleterious significance, we formulated quantitative methods for features, which fall into two categories: the first is variable type, including maximum, minimum and mean; the second is attribute type, which is measured by numerical sum. The reference genome version used by dbCNV is GRCh37/hg19.
# 2. Requirements
The local version dbCNV requires Linux (version: 3.10.0-514.el7.x86_64) environment, perl language (version:5.28.3) and R language (version: 3.6.0). This pipeline runs dbCNV.pl in the same directory as Database.
# 3. Installation

```
git clone https://github.com/kbvstmd/XCNV.git
```

# 4. Usage and example
# Usage:

```
perl dbCNV.pl -i prefix.txt -n 2
perl dbCNV.pl -i prefix.txt -n 5
```

 -i FILE    CNV region file (Chr Start End Tye)
 
 -n number  the number of pathogenic classification  (2 or 5)
 
 -h HELP    help information
 
# Example:

```
perl dbCNV.pl -i example.txt -n 2
```

The results can be seen in the gain_2_predication_result.txt and loss_2_predication_result.txt.
 
# 5. Input & output
Input file format (The columns are separated by TAB key and the header is required):

chr13 48611882  49054207	gain

Column 1: The chromosome

Column 2: Start

Column 3: End

Column 4: CNV type (gain or loss)

The output file has 5 columns including chromosome, Start, End, CNV type and pathogenicity of CNV.
