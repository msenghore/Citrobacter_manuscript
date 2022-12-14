Data and code for “Beyond consensus sequence: a quantitative scheme for inferring transmission using deep sequencing in a bacterial transmission model”


These scripts are compatible with the following system requirements:
R version 4.2.2, RStudio 2022.07.1+554, Perl v5.30.3

Description of scripts used to analyse the data

| Script       | Description|
| ------------- |:-------------:| 
|Variant_calling_steps.sh| A bash file showing the steps taken to map sequence reads, call variants and generate a multi sample vcf file|
| GVCF_snp_extract_to_longtbl.pl| Perl script that takes multi-sample gvcf and extracts allelic ratio and allele calls for each sample|
| Crod_prepare_distance_tables.R| An R scripts that takes the table of allelic frequencies, a table of pairwise SNP distances and calculates changes in allelic frequencies for each pairwise comparison, the likelihood of transmission and the number of shared within host variants|


Description of data files generated in the analysis

| Data file     | Description     | 
| ------------- |:-------------:| 
| All_snps_clean.vcf.gz| A multifile gvcf file generated using bwa, samtools, gatk, bedtools and vcfutils as described in methods| 
| All_snps_summary_long_clean.csv| A table showing the allele. Frequency at each variable site for each sample in the multi-sample vcf| 
| Citrobacter_distances_table.csv| A table showing the pairwise SNP distance between each pair of isolates (equivalent of a distance matrix) | 

