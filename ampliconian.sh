#!/bin/bash
# @author:      Carl-Eric Wegner
# @affiliation: Küsel Lab - Aquatic Geomicrobiology
#               Friedrich Schiller University of Jena
#               carl-eric.wegner@uni-jena.de
#
# The outlined pipeline is loosely based on exemplary workflows
# outlined by Robert Edgar outlined on his page.
# --> https://www.drive5.com/usearch/manual/pipe_examples.html
# The primary objective is speed and highest possible accuracy.
# Output data is formatted for downstream analysis using QIIME1/2
# and R (e.g. phyloseq).
#
# Given that USEARCH unfortunately follows a "freemium" model one goal
# will be to find open-source alternatives and modify the workflow accordingly.

### 1. Merged paired-end reads
# - sample names are listed in "samples.txt" to faciliate iterating
# - sequences are renamed according to the respective sample
# Note to myself: incorporate vsearch, code syntax should be the same
# except for an optional -reverse flag (maybe).
for sample in `awk '{print $1}' samples.txt`
do
    usearch10 -fastq_mergepairs *"$sample"_R1_clipped.fastq \
	    -fastqout "$sample"_merged.fastq \
	    -relabel "$sample".
done

sleep 5

# Write the number of sequences for all files into an array, and determine
# the sample with the lowest sequence count. This information will be later
# used for normalization.
for merged in *_merged.fastq
do
    arr[$i]=`awk '{s++}END{print s/4}' "$merged"`
    i=`expr $i + 1`
done
# show the array
# echo ${arr[*]}
min=${arr[0]}
for v in "${arr[@]}"; do
    if (( "$v" < "$min" )); then min="$v"; fi; 
done

# Concatenate all paired-end assembled files.
cat *_merged.fastq > all_merged.fastq

### 2. Remove primers, the output has to be still in .fastq format
# as QC is was not yet done
# - 341F (17 nt long), 785R (19 nt long)
# -information about most primers available from probebase
# --> http://probebase.csb.univie.ac.at
# Note to myself: could be done by a simple python script 
# with python speed could be an issue.
usearch10 -fastx_truncate all_merged.fastq \
    -stripleft 17 \
    -stripright 19 \
    -fastqout all_stripped.fastq

### 3. Quality control
# I adapted and I believe that the expected error metrics are superior
# to anything else (sliding windows, every Q-Score) commonly used.
# - details can be found here:
# --> https://www.drive5.com/usearch/manual/expected_errors.html
# --> https://www.drive5.com/usearch/manual/avgq.html
# The output is saved as .fasta to save storage space and speed up
# subsequent steps.
# Note to myself: incorporate vsearch, code syntax should be the same.
usearch10 -fastq_filter all_stripped.fastq \
    -fastq_maxee 1.0 \
    -fastaout all_filtered.fasta \
    -relabel Filtered

### 4. Identify unique sequences and determine their abundances
# - abundances are written to the output due to the sizeout flag
# Note to myself: incorporate vsearch, should be possible by dereplication,
# adding the vsearch sizeout flag should lead to compatible output.
usearch10 -fastx_uniques all_filtered.fasta -sizeout -relabel Unique -fastaout unique_seqs.fasta

### 5. OTU clustering and chimera removal
# The standard 97% sequence identity threshold is applied.
# The script will be updated to consider in addition the recently proposed
# concept of ASV.
# to be added:
# **deblur
# **dada2 
# functionality
usearch10 -cluster_otus unique_seqs.fasta -otus otus.fasta -relabel OTU
usearch10 -unoise3 unique_seqs.fasta -zotus zotus.fasta

### 6. Make OTU table and normalize to lowest number of reads per sample
usearch10 -otutab all_merged.fastq -otus otus.fasta -otutabout otutab.txt -mapout otu_map.txt
usearch10 -otutab all_merged.fastq -zotus zotus.fasta -otutabout zotutab.txt -mapout zotu_map.txt
usearch10 -otutab_norm otutab.txt -sample_size "$min" -output otutab_norm.txt

### 7. Taxonomic assignment and summary with SINTAX
# and replcace empty entries by "d:__unknown__"
# Taxonomic assignments are done against SILVA LTP and the RDP classifier training set,
# these are both rather small reference datasets, please have a look here:
# --> https://www.drive5.com/usearch/manual/faq_tax_db.html
# to understand my rationale.
# SINTAX taxonomy databases are available here:
# --> https://www.drive5.com/usearch/manual/sintax_downloads.html
# The 0.8 cutoff is chosen to guarantee an accuracy roughly equivalent to the RDP classifier with
# an 80% bootstrap cutoff.
# The overestimation of errors is much lower for SINTAX in comparison to the RDP classifier.
# --> https://www.drive5.com/usearch/manual/tax_err.html
usearch10 -sintax otus.fasta -db /path/to/SINTAX_dbs/ltp.fa \
    -strand both \
    -tabbedout sintax_ltp.txt \
    -sintax_cutoff 0.8
usearch10 -sintax otus.fasta -db /path/to/SINTAX_dbs/rdp.fa \
    -strand both \
    -tabbedout sintax_rdp.txt \
    -sintax_cutoff 0.8
# The sintax output table needs to be fixed, I follow the recommendation:
# https://www.drive5.com/usearch/manual/bugs.html
awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "d:__unknown__" }; 1' sintax_ltp.txt > sintax_ltp_fixed.txt
awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "d:__unknown__" }; 1' sintax_rdp.txt > sintax_rdp_fixed.txt
# Generate summaries on phylum, family and genus level.
usearch10 -sintax_summary sintax_rdp_fixed.txt \
    -otutabin otutab_norm.txt \
    -rank p \
    -output rdp_phylum_summary.txt
usearch10 -sintax_summary sintax_rdp_fixed.txt \
    -otutabin otutab_norm.txt \
    -rank f \
    -output rdp_family_summary.txt
usearch10 -sintax_summary sintax_rdp_fixed.txt \
    -otutabin otutab_norm.txt \
    -rank g \
    -output rdp_genus_summary.txt

### 8. Generate QIIME1/2, phyloseq compatible output
# Convert tab separated OTU table into a .biom table (JSON, .biom 1.0.0).
# Add taxonomic assignments.
biom convert -i otutab_norm.txt -o otutab_norm.biom --table-type="OTU table" --to-json
biom add-metadata -i otutab_norm.biom -o otutab_norm_rdp.biom \
    --observation-metadata-fp sintax_rdp_fixed.txt \
    --observation-header ID,taxonomy \
    --sc-separated taxonomy
biom add-metadata -i otutab_norm.biom -o otutab_norm_ltp.biom \
    --observation-metadata-fp sintax_ltp_fixed.txt \
    --observation-header ID,taxonomy \
    --sc-separated taxonomy
biom convert -i otutab_norm_rdp.biom -o otutab_norm_rdp.txt --to-tsv --header-key taxonomy
biom convert -i otutab_norm_ltp.biom -o otutab_norm_ltp.txt --to-tsv --header-key taxonomy
sed "s/[(][^)]*[)]//g" otutab_norm_rdp.txt > otutab_norm_rdp_fixed.txt
sed "s/[(][^)]*[)]//g" otutab_norm_ltp.txt > otutab_norm_ltp_fixed.txt

### 9. Taxonomic assignment with QIIME1.9 and SILVA
qiime19
assign_taxonomy.py -r /path/toSILVA123_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -t /path/to/SILVA123_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_all_levels.txt -i otus.fasta -o uclust_silva123
biom add-metadata -i otutab_norm.biom -o otutab_norm_silva.biom --observation-metadata-fp uclust_silva123/otus_tax_assignments.txt --observation-header ID,taxonomy --sc-separated taxonomy
summarize_taxa_through_plots.py -i otutab_norm_silva.biom -o tax_summary_silva