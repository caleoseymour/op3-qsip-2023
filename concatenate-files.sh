cat SEQ_plate_G/rep-seqs-150f-150r-dada2.fasta \
    SEQ_plates_A-B/rep-seqs-150f-144r-dada2.fasta \
    SEQ_plates_C-F/rep-seqs-150f-150r-dada2.fasta \
    whole_DNA/rep-seqs-150f-150r-dada2.fasta \
    > concatenated-rep-seqs-13C.fasta
    
cat SEQ_plate_G/table-150f-150r-dada2-long.tsv \
    SEQ_plates_A-B/table-150f-144r-dada2-long.tsv \
    SEQ_plates_C-F/table-150f-150r-dada2-long.tsv \
    whole_DNA/table-150f-150r-dada2-long.tsv \
    > concatenated-table-13C-long.tsv

cat SEQ_plate_G/rep-seqs-150f-150r-dada2-op3-taxonomy-2.tsv > concatenated-table-13C-op3-taxonomy-2.tsv
    tail -n+1 SEQ_plates_A-B/rep-seqs-150f-144r-dada2-op3-taxonomy-2.tsv >> concatenated-table-13C-op3-taxonomy-2.tsv
    tail -n+1 SEQ_plates_C-F/rep-seqs-150f-150r-dada2-op3-taxonomy-2.tsv >> concatenated-table-13C-op3-taxonomy-2.tsv
    tail -n+1 whole_DNA/rep-seqs-150f-150r-dada2-op3-taxonomy-2.tsv >> concatenated-table-13C-op3-taxonomy-2.tsv
    
# cat SEQ_plate_G/rep-seqs-150f-150r-dada2-taxonomy.tsv > concatenated-table-13C-taxonomy.tsv
    # tail -n+1 SEQ_plates_A-B/rep-seqs-150f-144r-dada2-taxonomy.tsv >> concatenated-table-13C-taxonomy.tsv
    # tail -n+1 SEQ_plates_C-F/rep-seqs-150f-150r-dada2-taxonomy.tsv >> concatenated-table-13C-taxonomy.tsv
    # tail -n+1 whole_DNA/rep-seqs-150f-150r-dada2-taxonomy.tsv >> concatenated-table-13C-taxonomy.tsv
    
    
cat 2014-07-02_Bri_16S/rep-seqs-150f-150r-dada2.fasta \
    2014-07-10_Bri_16S_run2_rerun/rep-seqs-150f-150r-dada2.fasta \
    2014-08-07_Bri_16S_plate3/rep-seqs-150f-150r-dada2.fasta \
    2014-08-14_Bri_16S_plate4/rep-seqs-150f-150r-dada2.fasta \
    2014-08-27_Bri_Mau_16S/rep-seqs-150f-150r-dada2.fasta \
    > concatenated-rep-seqs-18O.fasta
    
cat 2014-07-02_Bri_16S/table-150f-150r-dada2-long.tsv \
    2014-07-10_Bri_16S_run2_rerun/table-150f-150r-dada2-long.tsv \
    2014-08-07_Bri_16S_plate3/table-150f-150r-dada2-long.tsv \
    2014-08-14_Bri_16S_plate4/table-150f-150r-dada2-long.tsv \
    2014-08-27_Bri_Mau_16S/table-150f-150r-dada2-long.tsv \
    > concatenated-table-18O-long.tsv


tail -n+1 2014-07-02_Bri_16S/rep-seqs-150f-150r-dada2-op3-taxonomy-2.tsv > concatenated-table-18O-op3-taxonomy-2.tsv
    tail -n+1 2014-07-10_Bri_16S_run2_rerun/rep-seqs-150f-150r-dada2-op3-taxonomy-2.tsv >> concatenated-table-18O-op3-taxonomy-2.tsv
    tail -n+1 2014-08-07_Bri_16S_plate3/rep-seqs-150f-150r-dada2-op3-taxonomy-2.tsv >> concatenated-table-18O-op3-taxonomy-2.tsv
    tail -n+1 2014-08-14_Bri_16S_plate4/rep-seqs-150f-150r-dada2-op3-taxonomy-2.tsv >> concatenated-table-18O-op3-taxonomy-2.tsv
    tail -n+1 2014-08-27_Bri_Mau_16S/rep-seqs-150f-150r-dada2-op3-taxonomy-2.tsv >> concatenated-table-18O-op3-taxonomy-2.tsv
    
# tail -n+1 2014-07-02_Bri_16S/rep-seqs-150f-150r-dada2-taxonomy.tsv >> concatenated-table-18O-taxonomy.tsv
    # tail -n+1 2014-07-10_Bri_16S_run2_rerun/rep-seqs-150f-150r-dada2-taxonomy.tsv >> concatenated-table-18O-taxonomy.tsv
    # tail -n+1 2014-08-07_Bri_16S_plate3/rep-seqs-150f-150r-dada2-taxonomy.tsv >> concatenated-table-18O-taxonomy.tsv
    # tail -n+1 2014-08-14_Bri_16S_plate4/rep-seqs-150f-150r-dada2-taxonomy.tsv >> concatenated-table-18O-taxonomy.tsv
    # tail -n+1 2014-08-27_Bri_Mau_16S/rep-seqs-150f-150r-dada2-taxonomy.tsv >> concatenated-table-18O-taxonomy.tsv
    