# Automatic calculation of fixed, shared and private snp via a shiny app

## Input data
Your data must have a specific format, a tab separated table with different columns: 
scaffold, position, reference_nucleotide, alternative_nucleotide, allele_count and number_of_sample

To produce this kind of file from a vcf.gz you can use the following bash command: 
```bash
bcftools view -S population1_individual_list vcf_file.vcf.gz |bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' > population1_summary.txt 
```
bcftools view allows you to subset your vcf by a list of sample name and bcftools query reformat the vcf into a table.

The word before _summary.txt will be used as a name for each population so name your files accordingly. You will need one of these file for each population you have.

---

## How to run the shiny app
To run the shiny, you can install the **shiny** package in R, and
use the function `runGitHub()`. For example:

```R
if (!require('shiny')) install.packages("shiny")
shiny::runGitHub("fixed_shared_private_count", "BastienBennetot",ref="main")
```

Or you can clone or download this repository, and use `runApp()` function.

## How SNP are computed ?
There is two information when we compare alleles from population A and B 
The SNP is said to be "fixed" if all samples from a population have the same allele 
The SNP is said "snp" if at least one sample have a different allele within a population 
When we compare alleles present in both population, they can have same allele content (same ref or alt or same combination of multiple alleles) then we declare them "same" 
If allele content is different between population, then they are declared "diff"

## We can obtain these combinations
fixed fixed	same	=> IDENTIC 
snp snp	same	=>SHARED 
fixed snp	same	=>NOT POSSIBLE 
snp fixed	same	=>NOT POSSIBLE 
fixed fixed	diff	=>FIXED 
snp snp	diff	=>SHARED 
fixed snp	diff	=>PRIVATE 
snp fixed	diff	=>PRIVATE
