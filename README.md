# Automatic calculation of fixed, shared and private snp via a shiny app

##Input data
Your data must have a specific format, a tab separated table with different columns: 
scaffold, position, reference_nucleotide, alternative_nucleotide, allele_count and number_of_sample

To produce this kind of file from a vcf.gz you can use the following bash command: 
```bash
bcftools view -S population1_individual_list vcf_file.vcf.gz |bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' > population1_summary.txt 
```
The word before _summary.txt will be used as a name for each population so name your files accordingly

##How to run the shiny app
To run the shiny, you can install the **shiny** package in R, and
use the function `runGitHub()`. For example:

```R
if (!require('shiny')) install.packages("shiny")
shiny::runGitHub("fixed_shared_private_count", "BastienBennetot",ref="main")
```

Or you can clone or download this repository, and use `runApp()` function.
