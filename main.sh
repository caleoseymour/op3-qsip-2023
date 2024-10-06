## Cale Seymour
## 2021
## Process Bri Finley's qSIP data.

## Run Qiime2 on all of the datasets.
export WD='whole_DNA'
bash qiime2.sh
export WD='SEQ_plate_G'
bash qiime2.sh
export WD='SEQ_plates_A-B'
bash qiime2.sh
export WD='SEQ_plates_C-F'
bash qiime2.sh
export WD='2014-07-02_Bri_16S'
bash qiime2.sh
export WD='2014-07-10_Bri_16S_run2_rerun'
bash qiime2.sh
export WD='2014-08-07_Bri_16S_plate3'
bash qiime2.sh
export WD='2014-08-14_Bri_16S_plate4'
bash qiime2.sh
export WD='2014-08-27_Bri_Mau_16S'
bash qiime2.sh

## Concatenate results and process in R.
sh concatenate-files.sh
Rscript format-inputs.R

## Run qSIP.
Rscript p-vs-np.R
