This repository contains a Snakemake (Mölder et al 2021) pipeline for performing HLA imputation on a plink2 genotype dataset (.pgen + .pvar + .psam) using the softwares HIBAG (Zheng et al 2014) and CookHLA (Cook et al 2021). To run the pipeline, you must clone the CookHLA repository (https://github.com/WansonChoi/CookHLA.git) to a separate folder and obtain the InfiniumOmniExpress-24-European-HLA4-hg19.RData HIBAG model file. Then update the following parameters in /snakefiles/parameters/config.yaml

- input_trunk: the path to the input plink2 file set without the file endings
- release: an appropriate name for the release (e.g. the date you start running the pipeline)
- output_folder: the path to the folder where you want the output files to be stored
- hibag_model: the path to the InfiniumOmniExpress-24-European-HLA4-hg19.RData HIBAG model file
- CookHLA_path: the path to the folder where you cloned the CookHLA repository

You can now run the file /snakefiles/Snakefile using Snakemake (see https://snakemake.readthedocs.io/). The final ouput of the pipeline consists of the files alleles_merged and report.md. The former gives the the imputed alleles (hibag1, hibag2, cook1, cook2) for each gene and each sample, as well as relevant statistics output from HIBAG and CookHLA (prefixed with "hibag_" and "cook_"). hibag_mendel and cook_mendel (for HIBAG and CookHLA, respectively) are TRUE if the sample is from a child in a trio where the imputed alleles are a possible combination of the parents' alleles, FALSE if the sample is from a child in a trio where the imputed alleles are *not* a possible combination of the parents' alleles (i.e. if there is a Mendelian error) and NA if the sample is not from a child in a trio. The column 'inconsistencies' gives the number of inconsistencies between HIBAG and CookHLA for the given sample and gene (0 if the two softwares give the same result, 1 if they disagree on one allele, 2 if they disagree on both alleles.) The file report.md includes plots and statistics indicating the quality of the imputation from the two softwares, as follows:

- The number and rate of Mendelian errors for each software and gene. These are also found in a machine-readable format in the file /docs/{release}/mendelian_error_rates.

- Number of samples with 0, 1 and 2 inconsistent alleles between HIBAG and CookHLA, per gene.

- Single inconsistency counts per gene and allele. These tables list the number of times HIBAG and CookHLA disagree on each particular pair of alleles, where the rows are the HIBAG alleles and the columns the CookHLA alleles. For example the row 02:01, column 01:01 in the HLA-A table gives the number of times HIBAG suggests allele 02:01 while CookHLA instead suggest allele 01:01, for HLA-A. These tables are also reproduced in a machine-readable format in the folder /docs/{release}/inconsistencies_aggregated.

- Allele frequencies and probability density plots. The allele frequency plots compare the frequencies of the imputed alleles to the frequencies in the Norwegian population, published in Lande et al (2018). The probability density plots show the densities of the output 'prob' from HIBAG and 'confidence' from CookHLA (ideally these should show a sharp peak close to 1).

The pipeline also outputs lists of samples with single and double inconsistencies in the folder {output_folder}/{release}/inconsistencies_samples/

The [docs for the current release](https://github.com/people-of-norway/HLA-imputation/tree/main/docs/2025.09.16) show the report for the [beta release](https://github.com/fhi-beta/mobaGenetics-qc) of the genotype data from the [Norwegian Mother, Father and Child Cohort Study (MoBa)](https://www.fhi.no/en/studies/moba/).

Øystein Kapperud and Kinnie Le Roy have contributed to the work in this repository. It is currently maintained by Øystein Kapperud.

References:

Cook, S., Choi, W., Lim, H., Luo, Y., Kim, K., Jia, X., ... & Han, B. (2021). Accurate imputation of human leukocyte antigens with CookHLA. Nature communications, 12(1), 1264.

Lande, A., Andersen, I., Egeland, T., Lie, B. A., & Viken, M. K. (2018). HLA-A,-C,-B,-DRB1,-DQB1 and-DPB1 allele and haplotype frequencies in 4514 healthy Norwegians. Human immunology, 79(7), 527-529.

Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., (2021). Sustainable data analysis with Snakemake. F1000Res 10, 33.

Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir B (2014). “HIBAG – HLA Genotype Imputation with Attribute Bagging.” The Pharmacogenomics Journal. https://www.nature.com/articles/tpj201318.