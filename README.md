# 16S-analysis
## QIIME2 pipeline

qiime/2-2020.6

```ruby
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path /home/FCAM/nvijayan/(your input folder here) \
 --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path your_demuxed_filename.qza
```

```ruby
qiime demux summarize --i-data your_demuxed_filename.qza --o-visualization your_demuxed_filename.qzv
```

```ruby
qiime dada2 denoise-paired --i-demultiplexed-seqs your_demuxed_filename.qza  \
--o-table table_demux_dada.qza --o-representative-sequences rep-set.qza \
--o-denoising-stats denoising_stat.qza --p-trunc-len-f 250 --p-trunc-len-r 250
```

### TO FILTER OUT THE CONTROL SEQS
```ruby
qiime feature-table filter-samples \
  --i-table table_demux_dada.qza \
 --m-metadata-file metadata_squid_Jan2022.txt \
  --p-where '"#SampleID" IN ("Control1", "Control2","Control3")' \
  --o-filtered-table table-blanks.qza

qiime feature-table summarize \
  --i-table table-blanks.qza \
  --o-visualization table-blanks.qzv
```
### I donwloaded the "Frequency per feature detail in csv, which is the "feature-frequency-detail.csv" file.
echo 'Feature ID\tFrequency' | cat - feature-frequency-detail.csv | tr "," "\\t" > features-to-filter.tsv

```ruby
qiime feature-table filter-samples \
  --i-table table_demux_dada.qza \
  --m-metadata-file metadata_squid_Jan2022.txt \
  --p-where '"#SampleID" IN ("Control1", "Control2","Control3")' \
  --p-exclude-ids \
  --o-filtered-table table-sans-blanks.qza
```

```ruby
qiime feature-table filter-features \
  --i-table table-sans-blanks.qza \
  --m-metadata-file features-to-filter.tsv \
  --p-exclude-ids \
  --o-filtered-table filtered-table.qza
  
  qiime feature-table summarize \
  --i-table filtered-table.qza \
  --o-visualization filtered-table.qzv
  ```
  
### Convert SILVA database to QIIME classifier
```ruby
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy 99_ref-taxonomy.qza \
--o-classifier gg_99_classifier_2.qza
```

### Add taxonomy
```ruby
qiime feature-classifier classify-sklearn --i-classifier gg_99_classifier.qza --i-reads rep-set.qza --o-classification taxonomy.qza 
```

```ruby
qiime taxa barplot --i-table filtered-table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata_squid_Jan2022.txt  --o-visualization taxa-bar-plots.qzv
```
### To convert table to biom file
```ruby
qiime tools export  --input-path filtered-table.qza --output-path exported-feature-table
```
### To convert biom to tsv
```ruby
biom convert -i feature-table.biom -o table.txt --to-tsv
```

## The remaining analyses with ASV and taxonomy table can be merged and processed in R
