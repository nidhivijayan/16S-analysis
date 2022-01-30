# 16S-analysis
## QIIME2 pipeline


qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path /home/FCAM/nvijayan/(your input folder here) \
 --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path your_demuxed_filename.qza

qiime demux summarize --i-data your_demuxed_filename.qza --o-visualization your_demuxed_filename.qzv

qiime dada2 denoise-paired --i-demultiplexed-seqs your_demuxed_filename.qza  \
--o-table table_demux_dada.qza --o-representative-sequences rep-set.qza \
--o-denoising-stats denoising_stat.qza --p-trunc-len-f 250 --p-trunc-len-r 250

TO FILTER OUT THE CONTROL SEQS
qiime feature-table filter-samples \
  --i-table table_demux_dada.qza \
 --m-metadata-file metadata_squid_Jan2022.txt \
  --p-where '"#SampleID" IN ("Control1", "Control2","Control3")' \
  --o-filtered-table table-blanks.qza

qiime feature-table summarize \
  --i-table table-blanks.qza \
  --o-visualization table-blanks.qzv

# I donwloaded the "Frequency per feature detail in csv, which is the "feature-frequency-detail.csv" file.
echo 'Feature ID\tFrequency' | cat - feature-frequency-detail.csv | tr "," "\\t" > features-to-filter.tsv

qiime feature-table filter-samples \
  --i-table table_demux_dada.qza \
  --m-metadata-file metadata_squid_Jan2022.txt \
  --p-where '"#SampleID" IN ("Control1", "Control2","Control3")' \
  --p-exclude-ids \
  --o-filtered-table table-sans-blanks.qza

qiime feature-table filter-features \
  --i-table table-sans-blanks.qza \
  --m-metadata-file features-to-filter.tsv \
  --p-exclude-ids \
  --o-filtered-table filtered-table.qza
  
  qiime feature-table summarize \
  --i-table filtered-table.qza \
  --o-visualization filtered-table.qzv
  
# Convert SILVA database to QIIME classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy 99_ref-taxonomy.qza \
--o-classifier gg_99_classifier_2.qza

# Add taxonomy
qiime feature-classifier classify-sklearn --i-classifier gg_99_classifier.qza --i-reads rep-set.qza --o-classification taxonomy.qza 

qiime taxa barplot --i-table filtered-table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata_squid_Jan2022.txt  --o-visualization taxa-bar-plots.qzv
