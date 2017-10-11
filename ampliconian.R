# @author:      Carl-Eric Wegner
# @affiliation: Küsel Lab - Aquatic Geomicrobiology
#               Friedrich Schiller University of Jena
#               carl-eric.wegner@uni-jena.de
#
# Exemplary amplicon data analysis using phyloseq and respective dependencies.
# This pipeline is built on ampliconian.sh
#
# What do we need?
#               ** an OTU table (legacy QIIME .biom format, that is .biom v. 1.0.0)
#               ** a metadata table
#                  --> Linker and primer columns are not necessary, however
#                      phyloseq expects 6 columns

# Load necessary modules.
library("phyloseq")
library("ape")
library("DESeq2")
library("ggplot2")

# Import OTU table and sample metadata (aka a QIIME1 mapping file).
setwd("/home/calle/Storage/Data/becca_peat_cores/")
OTUtab <- import_biom("otutab_norm_silva.biom")
SAMPLEmeta <- import_qiime_sample_data("peat_core_mapping.txt")

### 1. Data import
# Quick look at the loaded data.
OTUtab
SAMPLEmeta

# Merge both into a phyloseq object.
physeq <- merge_phyloseq(OTUtab, SAMPLEmeta)

# Qick look at our phyloseq object which should contain three objects now:
#             ** an OTU table
#             ** a Taxonomy Table
#             ** a Sample Data Table
physeq

# A comment about normalization:
# ==============================
# The authors of phyloseq started a debate about proper sample normalization, with the main outcome being that
# sampling to an equal sequencing depth (= normalization based on the sample with the lowest sequence count) 
# is wrong. I refer here to the respective publication:
# McMurdie and Holmes (2014)
# --> http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531
#
# I absolutely agree with that in the context of differential abundance analysis. With respect to "general"
# or community analysis I can not see any harm in rarefying. From experience, overall community profiles are not
# tremendously affected by rarefaction.
#
# ampliconian.sh normalizes datasets to the one with the lowest sequence count.

### 2. Alpha (intra-sample) diversity analysis
# Observed refers to the number of identified OTUs, while Chao1 extrapolates detected OTU
# numbers to estimate the "true" number of OTUs. Chao1 put's special emphasis on singletons
# and doubletons. Shannon gives us an idea about sample evenness.
plot_richness(physeq, color = "X.SampleID", measures = c("Observed", "Chao1", "Shannon"))

### 3. Data filtering
# The below function filters out OTUs represented by less than three counts and OTUs that are present
# in less than 10% of all samples. This values need to individually adjusted. This filtering should
# not be done before alpha diversity analysis as it skewes up metrics such as Chao1
physeq_filt <- filter_taxa(physeq, function(x) sum(x > 3) > (0.1*length(x)), TRUE)
physeq_filt_norm <- transform_sample_counts(physeq_filt, function(OTU) OTU/sum(OTU))
any(taxa_sums(physeq_filt_norm) == 0)

### 4. Basic community profiling
# Summarize data on phylum level, plot the 20 most abundant phyla.
phyla_sum = tax_glom(physeq_filt_norm, "Rank2")
top_phyla <- prune_taxa(names(sort(taxa_sums(phyla_sum), TRUE)[1:20]), phyla_sum)
plot_bar(top_phyla, x = "X.SampleID", fill = "Rank2")

### 5. Beta diversity analysis
# Principal coordinates analysis, using the Jensen Shannon Divergence (JSD). As the JSD is not
# that common have a look for details here:
# --> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123278/
# --> http://stm.sciencemag.org/content/4/132/132ra52.full
# One thing I like is that it is independent from normalization. 
physeq_filt_norm.ord <- ordinate(physeq_filt_norm, "PCoA", "jsd")
p = plot_ordination(physeq_filt_norm, physeq_filt_norm.ord, color = "X.SampleID")
p = p + geom_point(size = 6, alpha = 0.7)
p