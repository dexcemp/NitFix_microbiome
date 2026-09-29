#activate conda environment
conda activate qiime2-amplicon-2024.2

#Run standard feature tables
qiime feature-table summarize --i-table table.hostremoved.qza --o-visualization table.qzv --m-sample-metadata-file sample-metadata.tsv
qiime feature-table tabulate-seqs --i-data rep-seqs.hostremoved.qza --o-visualization rep-seqs.qzv
qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.hostremoved.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza


#Filter just nodule samples
qiime feature-table filter-samples --i-table table.hostremoved.qza --m-metadata-file sample-metadata.tsv --p-where "[sample_type]='nodule'" --o-filtered-table nodule-table.qza

#Filter just soil samples
qiime feature-table filter-samples --i-table table.hostremoved.qza --m-metadata-file sample-metadata.tsv --p-where "[sample_type]='soil'" --o-filtered-table soil-table.qza

#Filter just rhizosphere samples
qiime feature-table filter-samples --i-table table.hostremoved.qza --m-metadata-file sample-metadata.tsv --p-where "[sample_type]='rhizosphere'" --o-filtered-table rhizosphere-table.qza

#Filter just root samples
qiime feature-table filter-samples --i-table table.hostremoved.qza --m-metadata-file sample-metadata.tsv --p-where "[sample_type]='root'" --o-filtered-table root-table.qza

#Create and view summary of feature table; helps inform what sampling depth to consider
qiime feature-table summarize --i-table nodule-table.qza --m-sample-metadata-file sample-metadata.tsv --o-visualization nodule-table.qzv
qiime feature-table summarize --i-table soil-table.qza --m-sample-metadata-file sample-metadata.tsv --o-visualization soil-table.qzv
qiime feature-table summarize --i-table rhizosphere-table.qza --m-sample-metadata-file sample-metadata.tsv --o-visualization rhizosphere-table.qzv
qiime feature-table summarize --i-table root-table.qza --m-sample-metadata-file sample-metadata.tsv --o-visualization root-table.qzv

# Rarefaction plots
# p-sampling-depth: change per experiment; look across your samples for feature counts -- choose a number as high as possible but lower thant the lowest sample
qiime diversity alpha-rarefaction --i-table table.hostremoved.qza --i-phylogeny rooted-tree.qza --p-max-depth 100 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction_100.qzv
qiime diversity alpha-rarefaction --i-table table.hostremoved.qza --i-phylogeny rooted-tree.qza --p-max-depth 1000 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction_1000.qzv
qiime diversity alpha-rarefaction --i-table table.hostremoved.qza --i-phylogeny rooted-tree.qza --p-max-depth 10000 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction_10000.qzv


# Diversity statistics and significance
# p-sampling-depth: change per experiment; look across your samples for feature counts by looking at table.qzv -- choose a number as high as possible but lower than the lowest sample
# Will not work for too few samples
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.hostremoved.qza --p-sampling-depth 100 --output-dir n_core-metrics-results --m-metadata-file sample-metadata.tsv 
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table nodule-table.qza --p-sampling-depth 100 --output-dir nodule_core-metrics-results --m-metadata-file sample-metadata.tsv 

qiime diversity alpha-group-significance --i-alpha-diversity n_core-metrics-results/faith_pd_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization n_core-metrics-results/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity nodule_core-metrics-results/faith_pd_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization nodule_core-metrics-results/faith-pd-group-significance.qzv

#Calculate Pielou evenness manually and plug in the result to alpha-group-significance
qiime diversity-lib pielou-evenness --i-table nodule-table.qza --p-drop-undefined-samples --o-vector n_core-metrics-results/noduleevenness.qza
qiime diversity alpha-group-significance --i-alpha-diversity n_core-metrics-results/noduleevenness.qza --m-metadata-file sample-metadata.tsv --o-visualization n_core-metrics-results/evenness-group-significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity n_core-metrics-results/noduleevenness.qza --m-metadata-file sample-metadata.tsv --o-visualization nodule_core-metrics-results/evenness-group-significance.qzv

#Shannon diversity group significance
qiime diversity alpha-group-significance --i-alpha-diversity n_core-metrics-results/shannon_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization n_core-metrics-results/shannon-group-significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity nodule_core-metrics-results/shannon_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization nodule_core-metrics-results/shannon-group-significance.qzv


#Export tsv
qiime tools export --input-path nodule_core-metrics-results/faith_pd_vector.qza --output-path export
mv export/alpha-diversity.tsv export/nodule_faith_pd_vector.tsv
# Add sample column header
sed -i 's/^\t/sample-id\t/g' export/*.tsv
qiime tools export --input-path nodule_core-metrics-results/shannon_vector.qza --output-path export
mv export/alpha-diversity.tsv export/nodule_shannon_vector.tsv
# Add sample column header
sed -i 's/^\t/sample-id\t/g' export/*.tsv


# Due to the longer time needed for beta diversity significance, a specific metadata column is specified to reduce computational time
# Can only run if there are both multiple categories and at least sometimes multiple individuals within each category (all can't be unique)
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column sample_type --o-visualization core-metrics-results/unweighted-unifrac-sample_type-significance.qzv --p-pairwise

qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column tribe --o-visualization core-metrics-results/unweighted-unifrac-tribe-significance.qzv --p-pairwise
qiime diversity beta-group-significance --i-distance-matrix nodule_core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column tribe --o-visualization nodule_core-metrics-results/unweighted-unifrac-tribe-significance.qzv --p-pairwise

qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column native --o-visualization core-metrics-results/unweighted-unifrac-native-significance.qzv --p-pairwise
qiime diversity beta-group-significance --i-distance-matrix nodule_core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column native --o-visualization nodule_core-metrics-results/unweighted-unifrac-native-significance.qzv --p-pairwise

qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column symbiont --o-visualization core-metrics-results/unweighted-unifrac-symbiont-significance.qzv --p-pairwise
qiime diversity beta-group-significance --i-distance-matrix nodule_core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column symbiont --o-visualization nodule_core-metrics-results/unweighted-unifrac-symbiont-significance.qzv --p-pairwise


## Ordination plots as an example -- requires numeric values to actually run -- regular PCoA is available in core metrics above
#qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes sample_type --o-visualization core-metrics-results/unweighted-unifrac-emperor-sample_type.qzv
#qiime emperor plot --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes sample_type --o-visualization core-metrics-results/bray-curtis-emperor-sample_type.qzv


# ANCOM -- differential abundance analysis
# Helps find which genera, families, etc. differ across metadata categories
# Imputation method for zero-abundance samples
# ANCOM with taxonomic level collapsing -- here level 6 genus
# Can run on original table without collapsing as well
# Some of these analyses cannot run with only one sample per category
qiime taxa collapse --i-table table.hostremoved.qza --i-taxonomy taxonomy.qza --p-level 6 --o-collapsed-table table.hostremoved-l6.qza
qiime composition add-pseudocount --i-table table.hostremoved.qza --o-composition-table comp-table.hostremoved-l6.qza

qiime taxa collapse --i-table nodule-table.qza --i-taxonomy taxonomy.qza --p-level 6 --o-collapsed-table nodule-table-l6.qza
qiime composition add-pseudocount --i-table nodule-table-l6.qza --o-composition-table comp-nodule-table-l6.qza

qiime composition ancom --i-table comp-table.hostremoved-l6.qza --m-metadata-file sample-metadata.tsv --m-metadata-column sample_type --o-visualization l6-ancom-sample_type.qzv
qiime composition ancom --i-table comp-nodule-table-l6.qza --m-metadata-file sample-metadata.tsv --m-metadata-column native --o-visualization l6-ancom-native.qzv
qiime composition ancom --i-table comp-nodule-table-l6.qza --m-metadata-file sample-metadata.tsv --m-metadata-column symbiont --o-visualization l6-ancom-symbiont.qzv
qiime composition ancom --i-table comp-nodule-table-l6.qza --m-metadata-file sample-metadata.tsv --m-metadata-column tribe --o-visualization l6-ancom-tribe.qzv
qiime composition ancom --i-table comp-nodule-table-l6.qza --m-metadata-file sample-metadata.tsv --m-metadata-column habitat --o-visualization l6-ancom-habitat.qzv
qiime composition ancom --i-table comp-nodule-table-l6.qza --m-metadata-file sample-metadata.tsv --m-metadata-column genus --o-visualization l6-ancom-genus.qzv


# Analyse how much was thrown away as host DNA using new taxonomic barplots
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.withhost.qzv
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.withhost.qzv
qiime taxa barplot --i-table nodule-table.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.withhost.noduleonly.qzv

#####
# Taxa barplots, no host, nodules only
qiime taxa barplot --i-table nodule-table.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.noduleonly.qzv


######
# Collapse taxon barplots by metadata
#Sample type
qiime feature-table group --i-table table.hostremoved.qza --m-metadata-file sample-metadata.tsv --m-metadata-column sample_type --p-axis sample --p-mode mean-ceiling --o-grouped-table table.grouped_sampletype # note lack of extension
qiime taxa barplot --i-table table.grouped_sampletype.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.grouped_sampletype.tsv --o-visualization taxa-bar-plots.grouped_sampletype.qzv

#state
qiime feature-table group --i-table table.hostremoved.qza --m-metadata-file sample-metadata.tsv --m-metadata-column state --p-axis sample --p-mode mean-ceiling --o-grouped-table table.grouped_state # note lack of extension
qiime taxa barplot --i-table table.grouped_state.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.grouped_state.tsv --o-visualization taxa-bar-plots.grouped_state.qzv

qiime feature-table group --i-table nodule-table.qza --m-metadata-file sample-metadata.tsv --m-metadata-column state --p-axis sample --p-mode mean-ceiling --o-grouped-table nodule-table.grouped_state # note lack of extension
qiime taxa barplot --i-table nodule-table.grouped_state.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.grouped_state.tsv --o-visualization taxa-bar-plots.grouped_state.noduleonly.qzv

#native
qiime feature-table group --i-table table.hostremoved.qza --m-metadata-file sample-metadata.tsv --m-metadata-column native --p-axis sample --p-mode mean-ceiling --o-grouped-table table.grouped_native # note lack of extension
qiime taxa barplot --i-table table.grouped_native.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.grouped_native.tsv --o-visualization taxa-bar-plots.grouped_native.qzv

qiime feature-table group --i-table nodule-table.qza --m-metadata-file sample-metadata.tsv --m-metadata-column native --p-axis sample --p-mode mean-ceiling --o-grouped-table nodule-table.grouped_native # note lack of extension
qiime taxa barplot --i-table nodule-table.grouped_native.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.grouped_native.tsv --o-visualization taxa-bar-plots.grouped_native.noduleonly.qzv

#habitat
qiime feature-table group --i-table table.hostremoved.qza --m-metadata-file sample-metadata.tsv --m-metadata-column habitat --p-axis sample --p-mode mean-ceiling --o-grouped-table table.grouped_habitat # note lack of extension
qiime taxa barplot --i-table table.grouped_habitat.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.grouped_habitat.tsv --o-visualization taxa-bar-plots.grouped_habitat.qzv

qiime feature-table group --i-table nodule-table.qza --m-metadata-file sample-metadata.tsv --m-metadata-column habitat --p-axis sample --p-mode mean-ceiling --o-grouped-table nodule-table.grouped_habitat # note lack of extension
qiime taxa barplot --i-table nodule-table.grouped_habitat.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.grouped_habitat.tsv --o-visualization taxa-bar-plots.grouped_habitat.noduleonly.qzv

#host tribe
qiime feature-table group --i-table nodule-table.qza --m-metadata-file sample-metadata.tsv --m-metadata-column tribe --p-axis sample --p-mode mean-ceiling --o-grouped-table nodule-table.grouped_tribe # note lack of extension
qiime taxa barplot --i-table nodule-table.grouped_tribe.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.grouped_tribe.tsv --o-visualization taxa-bar-plots.grouped_tribe.noduleonly.qzv


# must make a mock metadata table with the categories as header sample-id
