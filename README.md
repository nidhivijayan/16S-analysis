# 16S-analysis
## QIIME2 pipeline

conda activate qiime/2-2020.6 (Or module load if using local computer)

### Convert fastq to q2 artifact
```ruby
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path /home/FCAM/nvijayan/(your input folder here) \
 --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path your_demuxed_filename.qza
```

```ruby
qiime demux summarize --i-data your_demuxed_filename.qza --o-visualization your_demuxed_filename.qzv
```

### Denoise reads with DADA2
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
### To convert table to biom and taxonomy file
```ruby
qiime tools export  --input-path filtered-table.qza --output-path exported-feature-table

qiime tools export --input-path /Users/nidhivijayan/Documents/QIIME/all_squids_data/dada_qiime/taxonomy_demx_paired.qza --output-path taxonomy_phyloseq

```
### To convert biom to tsv
```ruby
biom convert -i feature-table.biom -o table.txt --to-tsv
```

## The remaining analyses with ASV and taxonomy table can be merged and processed in R

# In R
# Load packages
library(qiime2R)
library(tidyverse)
library("phyloseq", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("ggpubr")
library("vegan")
library("RColorbrewer")
library("scales")
library("reshape2")

### Convert qiime objects to R
physeq_all<-qza_to_phyloseq(
  features="table-all.qza",
  tree="rooted_tree_all.qza",
  "taxonomy_all.qza",
  metadata = "metadata_squid_All.txt"
)

### To normalize
physeq_norm = transform_sample_counts(physeq_all, function(x) 1538 * x/sum(x))

### To rarefy
ps2.rarefied = rarefy_even_depth(physeq_all, rngseed=1, sample.size=min(sample_sums(physeq_all)), replace=F)

### Beta diversity
#### For weighted Unifrac
wunifrac_norm_all = phyloseq::distance(physeq_norm, method="unifrac", weighted=T) 

#### For unweighted Unifrac
unifrac_norm_all = phyloseq::distance(physeq_norm, method="unifrac", weighted=F)

### Bray-Curtis distance
bdist = phyloseq::distance(physeq_norm,"bray")

### To make Ordination plot
ordination_all_norm = ordinate(physeq_norm, method="PCoA", distance=wunifrac_norm_all) 
method= PCoA or NMDS
distance= beta diversity distance

p_pcoa=plot_ordination(physeq_norm, ordination_pcoa, color="Family")+geom_point(size=5)+geom_point(aes(fill=Family), colour="black", size=5,shape=21)+stat_ellipse(linetype = 2)+stat_ellipse(geom = "polygon",aes(fill = Family), alpha = 0.25)+scale_fill_manual(values =dark2)+scale_color_manual(values =dark2)+theme(legend.title = element_text(size=20),legend.text = element_text(size=18),axis.text = element_text(size=14),legend.key.size = unit(2,"line"))

* For the color palette "dark2"
dark2 <- brewer.pal(8, "Dark2") 
dark2 <- colorRampPalette(coul)(14)

color3<- c("black","black","black","black","black")

### For Alpha diversity
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "ns", "ns"))

comparisons <- list(c("a_Sepiidae","b_Idiosepiidae"),c("a_Sepiidae","c_Sepiolidae"), c("a_Sepiidae", "d_Loliginidae"),c ("b_Idiosepiidae","c_Sepiolidae"), c("d_Loliginidae", "c_Sepiolidae"))

p1_alpha=plot_richness(physeq_norm, x="Family", measures=c("Observed","Shannon","InvSimpson"), color = "Family")+
   geom_boxplot(aes(fill=Family))+scale_fill_manual(values = dark2)+
   theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+scale_color_manual(values=color3)+
   stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", symnum.args = symnum.args)

rich = estimate_richness(physeq_norm,measures = "Observed")
*measures can be different alpha diversity measures

#### To perform wilcoxon test
pairwise.wilcox.test(rich$Oberseved, sample_data(physeq_norm)$Family)

#### To perform Dunn's test
alphaChao = estimate_richness(physeq_norm, measures="Observed")
alpha.stats <- cbind(alphaChao, sample_data(physeq_norm))
dunn.test(alpha.stats$Observed, alpha.stats$Description, method="bonferroni")







