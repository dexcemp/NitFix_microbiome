
# activate conda virtual environment
conda activate qiime2-amplicon-2024.2

# Make a manifest file of this format
######################
#sample-id	forward-absolute-filepath	reverse-absolute-filepath
#SA40	$PWD/test_seqs/A40_P1.fastq.gz	$PWD/test_seqs/A40_P2.fastq.gz
######################


# import sequences
# Qiime must have gzipped input files in a folder
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.txt --output-path import.qza --input-format PairedEndFastqManifestPhred33V2

# Trimming is to remove primers, set to ITS primers ITS1FI2/ITS2
qiime dada2 denoise-paired --i-demultiplexed-seqs import.qza --p-trim-left-f 18 --p-trim-left-r 20 --p-trunc-len-f 250 --p-trunc-len-r 200 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza

# Visualize sequencing stats
qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

# Taxonomic analysis
qiime feature-classifier classify-sklearn --i-classifier ./UNITE/unite_nitfix_99-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

# Remove chloroplast and mitochondrial 16S
qiime taxa filter-table --p-exclude k__Viridiplantae --i-table table.qza --i-taxonomy taxonomy.qza --o-filtered-table table.hostremoved.qza

qiime taxa filter-seqs --p-exclude k__Viridiplantae --i-sequences rep-seqs.qza --i-taxonomy taxonomy.qza --o-filtered-sequences rep-seqs.hostremoved.qza



# Taxa barplots, no host
qiime taxa barplot --i-table table.hostremoved.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.qzv


## Build UNITE database for fungal ITS (this was already done on BotBot) -- these steps take a LONG time -- about 2 hours
#conda activate qiime2-amplicon-2024.2
#cd UNITE
#awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' developer/sh_refs_qiime_ver8_99_02.02.2019_dev.fasta  | tr -d ' ' > developer/sh_refs_qiime_ver8_99_02.02.2019_dev_uppercase.fasta
#cat developer/sh_refs_qiime_ver8_99_02.02.2019_dev_uppercase.fasta ITS_phlawd_dec2019.fasta > unite_nitfix_99.fasta
#qiime tools import --type FeatureData[Sequence] --input-path unite_nitfix_99.fasta --output-path unite_nitfix_99.qza
#cat developer/sh_taxonomy_qiime_ver8_99_02.02.2019_dev.txt ITS_phlawd_dec2019.csv > unite_nitfix_99.txt
#qiime tools import --type FeatureData[Taxonomy] --input-path unite_nitfix_99.txt --output-path unite_nitfix_99-tax.qza --input-format HeaderlessTSVTaxonomyFormat
## Skip the reference trimming, unlike 16S, per Qiime documentation
#qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads unite_nitfix_99.qza --i-reference-taxonomy unite_nitfix_99-tax.qza --o-classifier unite_nitfix_99-classifier.qza



