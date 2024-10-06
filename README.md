# Scripts for processing qSIP samples.

## Description of files
Filename | Description
------------ | -------------
main.sh | main entry point to run the entire pipeline from scratch.
biom2LongMatrix.R | Reformat a biom file as a long matrix.
concatenate-files.sh | Join files from each sample type into a single table.
dada2-pipeline.sh | Generic 16S rRNA amplicon pipeline over several studies. Includes OP3-specific classification.
format-inputs.R | Format metadata files for each isotope input.
p-vs-np.R | Execute comparison of in-group to out-of-group differences in assimilation.
plot-soil-qsip.R | Generate plots.
qiime2.sh | Execute qiime2 pipeline on a single sample-folder.
run-ranked-qSIP.R | Execute comparison of in-taxon to out-of-group differences in assimilation for OP3.
tax_glom_fast.R | Much faster version of phyloseq's tax_glom function.
