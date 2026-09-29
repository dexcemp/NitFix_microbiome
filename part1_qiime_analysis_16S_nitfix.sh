
# activate conda virtual environment
conda activate qiime2-amplicon-2024.2

# Make a manifest file of this format
######################
#sample-id	forward-absolute-filepath	reverse-absolute-filepath
#A40	$PWD/test_seqs/A40_P1.fastq.gz	$PWD/test_seqs/A40_P2.fastq.gz
######################


# import sequences
# Qiime must have gzipped input files in a folder
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.txt --output-path import.qza --input-format PairedEndFastqManifestPhred33V2

# Trimming is to remove primers, set to 16S V3/V4 primers 515F/806R
qiime dada2 denoise-paired --i-demultiplexed-seqs import.qza --p-trim-left-f 19 --p-trim-left-r 20 --p-trunc-len-f 250 --p-trunc-len-r 200 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza
# If there is a generic error asking you to look at log files, and the log file says something like this: "Mismatched forward and reverse sequence files: 45896, 30375"
# Generally this means two very similar sample/sequence names that differ by _2 or something. I have been choosing the best of the "duplicates" and removing the other one from the manifest file.

# Visualize sequencing stats
qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

# Taxonomic analysis
qiime feature-classifier classify-sklearn --i-classifier ./greengenes_13_8_otus/classifier_gg_16S.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

# Remove chloroplast and mitochondrial 16S
qiime taxa filter-table --p-exclude c__Chloroplast,f__mitochondria --i-table table.qza --i-taxonomy taxonomy.qza --o-filtered-table table.hostremoved.qza

qiime taxa filter-seqs --p-exclude c__Chloroplast,f__mitochondria --i-sequences rep-seqs.qza --i-taxonomy taxonomy.qza --o-filtered-sequences rep-seqs.hostremoved.qza



# Taxa barplots, no host
qiime taxa barplot --i-table table.hostremoved.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.qzv


## Build GreenGenes database for 16S (this was already done on BotBot) -- these steps take a LONG time -- about 2 hours
## Format greengenes database (this was already done on BotBot)
#qiime tools import --type 'FeatureData[Sequence]' --input-path ./greengenes_13_8_otus/rep_set/97_otus.fasta --output-path 97_otus.qza
## Format greengenes database (this was already done on BotBot)
#qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path ./greengenes_13_8_otus/taxonomy/97_otu_taxonomy.txt --output-path 97_ref-taxonomy.qza
## Extract representative reads for OTUS from greengenes; set to 16S V3/V4 primers 515F/806R
## Already done on BotBot but could change if settings or primers change -- consider redoing but it takes a long time
#qiime feature-classifier extract-reads --i-sequences 97_otus.qza   --p-f-primer GTGCCAGCMGCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT --p-trunc-len 250 --p-min-length 200 --p-max-length 400 --o-reads 97_ref-seqs.qza
## Bayesian approach to taxonomic classifications on greengenes (explained in QIIME documentation -- simply taking the closest sequence has issues
#qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads 97_ref-seqs.qza --i-reference-taxonomy 97_ref-taxonomy.qza --o-classifier classifier_gg_16S.qza
#mv 97_*qza classifier_gg_16S.qza greengenes_13_8_otus/
