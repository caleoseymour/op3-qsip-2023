source include/q2tils.sh

logger "Working within $WD"

mkdir $WD/reads
mkdir $WD/temp

## Find the read files.
if [ ! -f $WD/reads/forward.fastq.gz ];
then
    FRAW=$(find $WD -iregex '.+_S0_L001_R1_001\.fastq\.gz')
    logger "Copying reads from $FRAW to $WD/reads/forward.fastq.gz"
    cp $FRAW $WD/reads/forward.fastq.gz
    stop_on_error $?
fi

if [ ! -f $WD/reads/reverse.fastq.gz ];
then
    RRAW=$(find $WD -iregex '.+_S0_L001_R2_001\.fastq\.gz')
    logger "Copying reads from $RRAW to $WD/reads/reverse.fastq.gz"
    cp $RRAW $WD/reads/reverse.fastq.gz
    stop_on_error $?

fi

## Find the barcode files, deal with them.
if [ ! -f $WD/reads/barcodes.fastq.gz ];
then
    I1RAW=$(find $WD -iregex '.+_S0_L001_I1_001\.fastq\.gz')
    I2RAW=$(find $WD -iregex '.+_S0_L001_I2_001\.fastq\.gz')
    if [ "$I2RAW" != "" ];
    then
        
        logger "Concatenating dual-index barcodes from $I1RAW and $I2RAW to $WD/reads/barcodes.fastq.gz"
        gunzip -ckq $I1RAW > $WD/temp/I1.fastq
        gunzip -ckq $I2RAW > $WD/temp/I2.fastq
        
        ## Stole this command from akutils, concatenates the barcode reads.
        ## Never done this before, but it looks like that's what qiime1 does to handle dual-index barcodes as well.
        paste -d '' <(echo; sed -n '1,${n;p;}' $WD/temp/I1.fastq | sed G) $WD/temp/I2.fastq | sed '/^$/d' > $WD/temp/barcodes.fastq
        stop_on_error $?
        rm $WD/temp/I1.fastq
        rm $WD/temp/I2.fastq
        
        gzip -cq $WD/temp/barcodes.fastq > $WD/reads/barcodes.fastq.gz
        rm $WD/temp/barcodes.fastq
    else
        logger "Copying single-index barcodes from $I1RAW to $WD/reads/barcodes.fastq.gz"
        cp $I1RAW $WD/reads/barcodes.fastq.gz
    fi
fi 

## Find the mapping file.
MAPPING=$(find $WD -iregex '.+map.*.txt')
logger "Mapping file chosen: $MAPPING"

## Import into qiime2.
if [ ! -f $WD/paired-end-seqs.qza ];
then
    logger "Importing sequences..."
    qiime tools import \
      --type 'EMPPairedEndSequences' \
      --input-path $WD/reads \
      --output-path $WD/paired-end-seqs.qza
    stop_on_error $?
fi

if [ ! -f $WD/paired-end-demux.qza ];
then
    logger "Demultiplexing sequences..."
    qiime demux emp-paired \
      --i-seqs $WD/paired-end-seqs.qza \
      --m-barcodes-file $MAPPING \
      --m-barcodes-column BarcodeSequence \
      --p-no-golay-error-correction \
      --o-per-sample-sequences $WD/paired-end-demux.qza \
      --o-error-correction-details $WD/demux-ec-details.qza
    stop_on_error $?
fi

## Run deblur piplines, truncating at 150 and 250 nt
#bash deblur-pipeline.sh 150
#bash deblur-pipeline.sh 250

## Run dada2 pipeline, setting read trimming at 150f 150r.
#logger "Executing 150f 150r dada2 pipeline"
#bash dada2-pipeline.sh 150 150

if [ "$WD" == "SEQ_plates_A-B" ];
then
    logger "Executing 150f 144r dada2 pipeline"
    bash dada2-pipeline.sh 150 144
else
    logger "Executing 150f 150r dada2 pipeline"
    bash dada2-pipeline.sh 150 150
fi


