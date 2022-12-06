# 16S-analysis
## QIIME2 pipeline

conda activate qiime/2-2020.6 (Or module load if using local computer)

### Convert fastq to q2 artifact
```ruby
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path /home/FCAM/nvijayan/(your input folder here) \
 --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux.qza
```

```ruby
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
```

### Denoise reads with DADA2
```ruby
qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza  \
--o-table table_demux_dada.qza --o-representative-sequences rep-set.qza \
--o-denoising-stats denoising_stat.qza --p-trunc-len-f 250 --p-trunc-len-r 250
```

  
### Convert SILVA database to QIIME classifier
```ruby
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy 99_ref-taxonomy.qza \
--o-classifier gg_99_classifier.qza
```

### Add taxonomy
```ruby
qiime feature-classifier classify-sklearn --i-classifier gg_99_classifier.qza --i-reads rep-set.qza --o-classification taxonomy.qza 
```

```ruby
qiime taxa barplot --i-table filtered-table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata_squid_Jan2022.txt  --o-visualization taxa-bar-plots.qzv
```
```ruby
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-set.qza \
  --output-dir phylogeny-align-to-tree-mafft-fasttree
 ```
 
 ### Diversity metrics
 ```ruby
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny phylogeny-align-to-tree-mafft-fasttree/rooted_tree.qza \
  --i-table Eup-table-no-mitochondria-no-chloroplast.qza \
  --p-sampling-depth 6728 \
  --m-metadata-file metadata_squid.txt \
  --output-dir diversity-core-metrics-phylogenetic
  
qiime diversity beta-group-significance \
  --i-distance-matrix diversity-core-metrics-phylogenetic/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata_squid.txt \
  --m-metadata-column Description \
  --p-pairwise \
  --o-visualization weighted-unifrac-group-significance.qzv
  ```

### To convert table to biom and taxonomy file
```ruby
qiime tools export  --input-path filtered-table.qza --output-path exported-feature-table

qiime tools export --input-path taxonomy_demx_paired.qza --output-path taxonomy_phyloseq

```
### To convert biom to tsv
```ruby
biom convert -i feature-table.biom -o table.txt --to-tsv
```

## The remaining analyses with ASV and taxonomy table can be merged and processed in R

# In R
# Load packages
```
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
```

### Convert qiime objects to R
```
physeq_all<-qza_to_phyloseq(
  features="table-all.qza",
  tree="rooted_tree_all.qza",
  "taxonomy_all.qza",
  metadata = "metadata_squid_All.txt"
)
```

### To normalize
```
physeq_norm = transform_sample_counts(physeq_all, function(x) 1538 * x/sum(x))
```
### To rarefy
```
ps2.rarefied = rarefy_even_depth(physeq_all, rngseed=1, sample.size=min(sample_sums(physeq_all)), replace=F)
```
### Beta diversity
#### For weighted Unifrac
```
wunifrac_norm_all = phyloseq::distance(physeq_norm, method="unifrac", weighted=T) 
```
#### For unweighted Unifrac
```
unifrac_norm_all = phyloseq::distance(physeq_norm, method="unifrac", weighted=F)
```
### Bray-Curtis distance
```
bdist = phyloseq::distance(physeq_norm,"bray")
```
### To make Ordination plot
```
ordination_all_norm = ordinate(physeq_norm, method="PCoA", distance=wunifrac_norm_all) 
```
method= PCoA or NMDS
distance= beta diversity distance
```
p_pcoa=plot_ordination(physeq_norm, ordination_pcoa, color="Family")+geom_point(size=5)+geom_point(aes(fill=Family), colour="black", size=5,shape=21)+stat_ellipse(linetype = 2)+stat_ellipse(geom = "polygon",aes(fill = Family), alpha = 0.25)+scale_fill_manual(values =dark2)+scale_color_manual(values =dark2)+theme(legend.title = element_text(size=20),legend.text = element_text(size=18),axis.text = element_text(size=14),legend.key.size = unit(2,"line"))
```
### For the color palette "dark2"
```
dark2 <- brewer.pal(8, "Dark2") 
dark2 <- colorRampPalette(coul)(14)

paired<-brewer.pal(12, "Paired") 
paired <- colorRampPalette(paired)(15)

color3<- c("black","black","black","black","black")
```
### For Alpha diversity
```
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "ns", "ns"))

comparisons <- list(c("a_Sepiidae","b_Idiosepiidae"),c("a_Sepiidae","c_Sepiolidae"), c("a_Sepiidae", "d_Loliginidae"),c ("b_Idiosepiidae","c_Sepiolidae"), c("d_Loliginidae", "c_Sepiolidae"))

p1_alpha=plot_richness(physeq_norm, x="Family", measures=c("Observed","Shannon","InvSimpson"), color = "Family")+
   geom_boxplot(aes(fill=Family))+scale_fill_manual(values = dark2)+
   theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+scale_color_manual(values=color3)+
   stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", symnum.args = symnum.args)

rich = estimate_richness(physeq_norm,measures = "Observed")
```
*measures can be different alpha diversity measures

#### To perform wilcoxon test
```
pairwise.wilcox.test(rich$Oberseved, sample_data(physeq_norm)$Family)
```
#### To perform Dunn's test
```
alphaChao = estimate_richness(physeq_norm, measures="Observed")
alpha.stats <- cbind(alphaChao, sample_data(physeq_norm))
dunn.test(alpha.stats$Observed, alpha.stats$Description, method="bonferroni")
```

### PERMANOVA
```
adonis(wunifrac_norm_all ~ sample_data(physeq_norm)$Description,permutations = 10000)
dispr <- vegan::betadisper(wunifrac_norm_all, phyloseq::sample_data(physeq_norm)$Description)
```
## To run mantel test
### First upload host distance matrix
```
distance_mya_m=read.delim("distance_mya_mean.txt",row.names = 1)
distance_mya<-as.matrix(distance_mya_m)
```
### Create UPGMA microbiome distance (Trevelline et al. 2020)
```
physeq.for.upgma.table = as(otu_table(physeq.rarefied), "matrix")
physeq.for.upgma.table <- t(physeq.for.upgma.table)
physeq.for.upgma.table <- as.data.frame(physeq.for.upgma.table)

physeq.for.upgma.table$SampleID <- row.names(physeq.for.upgma.table)
physeq.for.upgma.table$Species <- metadata2$Species[match(metadata2$Description, physeq.for.upgma.table$SampleID)]
```
### FORMAT DATA
```
physeq.for.upgma.table.molten <- melt(physeq.for.upgma.table, id.var = c("SampleID", "Species"), var = 'OTUID')
```
### COMPUTE AVERAGES
```
physeq.for.upgma.table.mean <- aggregate(value~Species * OTUID, FUN = 'mean', data = physeq.for.upgma.table.molten)
physeq.for.upgma.table.mean$value <- ceiling(physeq.for.upgma.table.mean$value) #Round values up (mean-ceiling)
physeq.for.upgma.table.mean.casted <- dcast(physeq.for.upgma.table.mean, OTUID ~ Species)
write.csv(physeq.for.upgma.table.mean.casted, "physeq.for.upgma.table.mean.casted.csv")

otu.table2 = read.csv("physeq.for.upgma.table.mean.casted.csv", sep = ',', row.names = 1) # must delete first column of numbers before import
head(otu.table)
otu.table2 = as.matrix(otu.table2)

metadata.by.species = read.csv("metadata_by_species.csv", sep = ",", row.names = 1)
head(metadata.by.species)

OTU2 = otu_table(otu.table2, taxa_are_rows = TRUE) #Import physeq.for.upgma.mean.casted as Phyloseq object
META2a = sample_data(metadata.by.species) #Import metadata as Phyloseq object
phy.tree2 <- read_tree("pruned_tree.nwk") #Import rooted tree from QIIME. Must be trimmed of all exluded OTUs.

physeq.for.upgma.mean = phyloseq(OTU2, TAX, META2a, phy.tree2) #Merge all files into a single Phyloseq object.
```
### Bray-Curtis
```
physeq.for.upgma.mean.rarefied.bray.dist <- distance(physeq.for.upgma.mean, method = "bray") #Generate distance matrix
physeq.for.upgma.mean.rarefied.bray.dist.UPGMA <- hclust(physeq.for.upgma.mean.rarefied.bray.dist, method = "average")
plot(physeq.for.upgma.mean.rarefied.bray.dist.UPGMA)
physeq.for.upgma.mean.rarefied.bray.dist.UPGMA <- as.phylo(physeq.for.upgma.mean.rarefied.bray.dist.UPGMA)
write.tree(physeq.for.upgma.mean.rarefied.bray.dist.UPGMA, "physeq.rarefied.bray.dist.UPGMA.newick")
```
### Unweighted
```
physeq.for.upgma.mean.rarefied.unweighted.dist <- distance(physeq.for.upgma.mean, method = "uunifrac") #Generate distance matrix
physeq.for.upgma.mean.rarefied.unweighted.dist.UPGMA <- hclust(physeq.for.upgma.mean.rarefied.unweighted.dist, method = "average")
plot(physeq.for.upgma.mean.rarefied.unweighted.dist.UPGMA)
physeq.for.upgma.mean.rarefied.unweighted.dist.UPGMA <- as.phylo(physeq.for.upgma.mean.rarefied.unweighted.dist.UPGMA)
write.tree(physeq.for.upgma.mean.rarefied.unweighted.dist.UPGMA, "physeq.rarefied.unweighted.dist.UPGMA.newick")
```
### Weighted
```
physeq.for.upgma.mean.rarefied.weighted.dist <- distance(physeq.for.upgma.mean, method = "wunifrac") #Generate distance matrix
physeq.for.upgma.mean.rarefied.weighted.dist.UPGMA <- hclust(physeq.for.upgma.mean.rarefied.weighted.dist, method = "average")
plot(physeq.for.upgma.mean.rarefied.weighted.dist.UPGMA)
physeq.for.upgma.mean.rarefied.weighted.dist.UPGMA <- as.phylo(physeq.for.upgma.mean.rarefied.weighted.dist.UPGMA)
write.tree(physeq.for.upgma.mean.rarefied.weighted.dist.UPGMA, "physeq.rarefied.weighted.dist.UPGMA_avg.newick")
```
### Run Mantel with vegan package
```
mantel(physeq.for.upgma.mean.rarefied.weighted.dist,distance_mya,method = "spearman",permutations=999)
mantel(physeq.for.upgma.mean.rarefied.unweighted.dist,distance_mya,method = "pearson",permutations=999)
mantel(physeq.for.upgma.mean.rarefied.bray.dist,distance_mya,method = "pearson",permutations=999)
```

### To make distance decay plot
```
d_mya=as.dist(distance_mya)
plot(d_mya,physeq.for.upgma.mean.rarefied.weighted.dist,ylim=c(0,1),xlim = c(0,max(d_mya)))
decay.exp<-decay.model(physeq.for.upgma.mean.rarefied.weighted.dist,d_mya,y.type="dissim",model.type="exp",perm=100)
plot.decay(decay.exp, col="magenta", remove.dots=TRUE, add=TRUE)
```
## To make Barplots
### Make each species its own physeq object
```
speciesList<-tapply(sample_names(physeq_norm), get_variable(physeq_norm, "Description"), c)

speciesList_ID<-tapply(sample_names(physeq_norm), get_variable(physeq_norm, "Sample"), c)

speciesPhyseq <- lapply(speciesList_ID, prune_samples, physeq_norm)

speciesOTUtable <- lapply(speciesPhyseq,otu_table)
speciesAvg <- lapply(speciesOTUtable,rowMeans)
pooledOTUtable = t(do.call(rbind,speciesAvg))
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable)
View(pooledOTUtable)

TT = tax_table(ps2.rarefied)
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]

tdf = data.frame(TT, OTU = taxa_names(ps2.rarefied))
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")

pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:71]))
write.csv(pOTU,"pOTU.csv") *Inspect the table*
pOTU_avg=read.csv("pOTU.csv",header=TRUE) 
pOTU.avg.class =pOTU_avg[,c(2:13,16)]
melt.class.avg = melt(pOTU.avg.class,id.vars="Class") #Library(reshape2)
colnames(melt.class.avg)[2]="species"
agg.class.avg=aggregate(.~Class+species,melt.class.avg,sum)

bar=ggplot(agg.class.avg,aes(x=species,y=value,fill=Class)) +
       geom_bar(stat="identity",position="fill",color="black") +
        scale_y_continuous(labels = percent_format())+
         xlab("Species") +
        ylab("Relative Abudance")+
             theme(axis.text.x = element_text(angle=30,vjust=1, hjust = 1, size=12),strip.text.x = element_text(size=15),legend.text = element_text(size = 13),legend.key.size =  unit(2,"line"),axis.text.y=element_text(size=12),axis.title.y = element_text(size=15))+ scale_fill_manual(values = paired)+
     facet_grid(.~family+species2,drop = TRUE,scales="free",space="free")
```

### To prune tree to keep some otus only 
```
t1=read.tree("tree.nwk")
otuid2=read.csv("otuid2.csv",header=FALSE) #made table with just otus to keep
o2<-as.character(otuid2$V1) #Makes the column to character list
pruned.tree<-drop.tip(t1,t1$tip.label[-match(o2, t1$tip.label)])
write.tree(pruned.tree,"pruned_tree.nwk")
```

