source include/q2tils.sh

if [ ! -f $WD/table-"$1"f-"$2"r-dada2.qza ];
then
    logger "Denoise..."
     qiime dada2 denoise-paired \
        --i-demultiplexed-seqs $WD/paired-end-demux.qza \
        --p-trunc-len-f $1 \
        --p-trunc-len-r $2 \
        --o-table $WD/table-"$1"f-"$2"r-dada2.qza \
        --o-representative-sequences $WD/rep-seqs-"$1"f-"$2"r-dada2.qza \
        --o-denoising-stats $WD/dada2-"$1"f-"$2"r-stats.qza \
        --verbose
    stop_on_error $?
    logger "Denoise complete."
fi


if [ ! -f $WD/rep-seqs-"$1"f-"$2"r-dada2-taxonomy.qza ];
then
    logger "Classifying via SILVA..."
    classifierDir='/mnt/d/Ubuntu/Archive/Qiime2/emp/deblur/qiime-emp/'
    qiime feature-classifier classify-sklearn \
        --i-classifier $classifierDir/naiive-bayes-classifier-q2_2020.8-515F-806R-silva-99nr-138.qza \
        --i-reads $WD/rep-seqs-"$1"f-"$2"r-dada2.qza \
        --o-classification $WD/rep-seqs-"$1"f-"$2"r-dada2-taxonomy.qza \
        --p-reads-per-batch 1000
    stop_on_error $?
fi

#if [ ! -f $WD/rep-seqs-"$1"f-"$2"r-dada2-op3-taxonomy-2.qza ];
#then
    logger "Classifying via OP3 classifier..."
    classifierDir='/mnt/d/Ubuntu/science/OP3/try-5/marker-trees/16S/'
    qiime feature-classifier classify-sklearn \
        --i-classifier $classifierDir/omnitrophota-classifier.515F-806R.qza \
        --i-reads $WD/rep-seqs-"$1"f-"$2"r-dada2.qza \
        --o-classification $WD/rep-seqs-"$1"f-"$2"r-dada2-op3-taxonomy-2.qza
    stop_on_error $?
#fi

q2export $WD/table-"$1"f-"$2"r-dada2.qza feature-table.biom
biom convert --to-tsv -i $WD/table-"$1"f-"$2"r-dada2.biom -o $WD/table-"$1"f-"$2"r-dada2.tsv
Rscript biom2LongMatrix.R $WD/table-"$1"f-"$2"r-dada2.tsv $WD/table-"$1"f-"$2"r-dada2-long.tsv

q2export $WD/rep-seqs-"$1"f-"$2"r-dada2-op3-taxonomy-2.qza taxonomy.tsv
q2export $WD/rep-seqs-"$1"f-"$2"r-dada2-taxonomy.qza taxonomy.tsv

q2export $WD/rep-seqs-"$1"f-"$2"r-dada2.qza dna-sequences.fasta
