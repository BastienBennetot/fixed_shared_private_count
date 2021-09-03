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
Within population:  
The SNP is said to be **"fixed"** if all samples from a population have the same allele  
The SNP is said **"snp"** if at least one sample have a different allele within a population 
Pairwise population:  
When we compare alleles present in both populations, they can have same allele content (same ref or alt or same combination of multiple alleles) then we declare them **"same"**  
If allele content is different between populations, then they are declared **"diff"**  

## We can obtain these combinations 
fixed fixed same  &#10132;  IDENTIC (means this site is FIXED in population A and FIXED in population B with the SAME nucleotide)  
snp snp same  &#10132;  SHARED  
fixed snp same  &#10132;  NOT POSSIBLE  
snp fixed same  &#10132;  NOT POSSIBLE  
fixed fixed diff  &#10132;  FIXED  
snp snp diff  &#10132;  SHARED  
fixed snp diff  &#10132;  PRIVATE  
snp fixed diff  &#10132;  PRIVATE

## Download the computed data frame
On the third tab you can download the big dataframe containing all information about which snp is shared fixed or private in any pairwise population comparison  
It has different columns :  
scaffold: the name of the scaffold of the site considered  
position: the position of the site considered  
**Since there is a pairwise comparison between population x and population y some columns are duplicated (column.x belongs to population x)**
ref.x: the ref nucleotide 
alt.x: the alt nucleotide (can be *)
sites.x: number of samples in each	alt genotype separated by a comma  
individual_known.x: number of sample that genotype is known  
source.x: name of the "x" population  
number_of_ref.x: number of samples that have the same genotype as the reference  
sites_fixed.x: a binary code of genotype present in the population with the first number equal to the reference and others are ALT genotype in order (of alt.x column). 1=present, 0=missing  
fixed.x: Does all samples have the same genotype, if yes it writes "fixed"  

diff: Does allele content (genotypes) are different between population x and population y  
status: Status of this site for this pairwise comparison  
pairwise: a pairwise comparison name in the format of "x vs y"
