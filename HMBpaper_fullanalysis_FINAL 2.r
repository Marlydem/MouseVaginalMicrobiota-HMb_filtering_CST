HMBpaper_fullanalysis_FINAL

%%%%%%%%%%%%%%%%%%%%%%%%%% Qiime2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Study files UCSD Dataset
gzip -d [filename.fastq.gz] #shows and opens everything
less [filename.fastq.gz] #opens in sectioned off window but gives you a new command line to go on
less [head(-8)] #not tried yet but something like this is in my notes from Sara

#started with directory for Dec 1-3 but decided to uploaded all at once
mkdir Dec1-3_emp-paired-end-sequences
  #download raw files from qiita
  #convert names(look to end of file name LIBR_38_41_44_45_S1_L001_R1_001.fastq)
  #R1= forward.fastq.gz
  #R2= revers.fastq.gz
  #I1= barcodes.fastq.gz
qiime tools import \
--type EMPPairedEndSequences \
--input-path Dec1-3_emp-paired-end-sequences \
--output-path Dec1-3_emp-paired-end-sequences.qza

qiime tools import \
--type EMPPairedEndSequences \
--input-path RerunMoMISeed1-2_emp-paired-end-sequences \
--output-path RerunMoMISeed1-2_emp-paired-end-sequences.qza

qiime tools import \
--type EMPPairedEndSequences \
--input-path UpdatedMoMISeed1-3rerun_second_set_emp-paired-end-sequences \
--output-path UpdatedMoMISeed1-3rerun_second_set_emp-paired-end-sequences.qza

qiime tools peek Dec1-3_emp-paired-end-sequences.qza
#Read out
"""
  UUID:        159ec267-dc8b-45ec-94dd-fd6d53433e45
  Type:        EMPPairedEndSequences
  Data format: EMPPairedEndDirFmt
"""

#DEMULTIPLEX -------------------------------------------------------------
#download "prep info" from qiita to get sample-metadata with barcodes
#rename column barcode-sequence (may be able to change below to "barcode" but the above comment worked for me)

%%CORRECTED CODE for GOLAY 12n, rev compl barcodes & rev compl mapping barcodes%%
  
%%UpdatedMoMISeed1-3rerun_second_set_emp-paired-end-sequences.qza%%
qiime demux emp-paired \
--i-seqs UpdatedMoMISeed1-3rerun_second_set_emp-paired-end-sequences.qza \
--m-barcodes-file sample-metadata1.txt \
--m-barcodes-column barcode-sequence \
--p-golay-error-correction TRUE \
--p-rev-comp-barcodes TRUE \
--p-rev-comp-mapping-barcodes TRUE \
--o-per-sample-sequences demux-full_1.qza \
--o-error-correction-details demux-details_1.qza
#metadata read as txt!! yay!

#just to view
qiime demux summarize \
--i-data demux-full_1.qza \
--o-visualization demux-full_1.qzv

qiime tools view demux-full_1.qzv

#READ-JOINING TO MATCH END-PAIRS
qiime vsearch join-pairs \
--i-demultiplexed-seqs demux-full_1.qza \
--o-joined-sequences qiime2-read-joining-tutorial/demux-full_1-joined.qza

qiime demux summarize \
--i-data qiime2-read-joining-tutorial/demux-full_1-joined.qza \
--o-visualization qiime2-read-joining-tutorial/demux-full_1-joined.qzv

#QUALITY FILTERING (referenced Qiita processing info on split libraries)
qiime quality-filter q-score-joined \
--i-demux qiime2-read-joining-tutorial/demux-full_1-joined.qza \
--p-min-quality 3 \
--p-quality-window 3 \
--p-min-length-fraction 0.75 \
--p-max-ambiguous 0 \
--o-filtered-sequences qiime2-read-joining-tutorial/demux-full_1-joined-filtered.qza \
--o-filter-stats qiime2-read-joining-tutorial/demux-full_1-joined-filter-stats.qza

qiime deblur denoise-16S \
--i-demultiplexed-seqs qiime2-read-joining-tutorial/demux-full_1-joined-filtered.qza \
--p-trim-length 150 \
--p-sample-stats \
--p-mean-error 0.005 \
--p-indel-prob 0.01 \
--p-indel-max 3 \
--p-min-reads 0 \
--p-min-size 2 \
--p-jobs-to-start 5 \
--o-representative-sequences rep-seqs_1.qza \
--o-table table_1.qza \
--o-stats deblur-stats_1.qza

#Visualize deblur
qiime metadata tabulate \
--m-input-file qiime2-read-joining-tutorial/demux-full_1-joined-filter-stats.qza \
--o-visualization demux-filter-stats_1.qzv
qiime deblur visualize-stats \
--i-deblur-stats deblur-stats_1.qza \
--o-visualization deblur-stats_1.qzv

------------

%%RerunMoMISeed1-2_emp-paired-end-sequences.qza%%
  qiime demux emp-paired \
--i-seqs RerunMoMISeed1-2_emp-paired-end-sequences.qza \
--m-barcodes-file sample-metadata2.txt \
--m-barcodes-column barcode-sequence \
--p-golay-error-correction TRUE \
--p-rev-comp-barcodes TRUE \
--p-rev-comp-mapping-barcodes TRUE \
--o-per-sample-sequences demux-full_2.qza \
--o-error-correction-details demux-details_2.qza
#metadata read as txt!! yay!

#just to view
qiime demux summarize \
--i-data demux-full_2.qza \
--o-visualization demux-full_2.qzv

qiime tools view demux-full_2.qzv

#READ-JOINING TO MATCH END-PAIRS
qiime vsearch join-pairs \
--i-demultiplexed-seqs demux-full_2.qza \
--o-joined-sequences qiime2-read-joining-tutorial/demux-full_2-joined.qza

qiime demux summarize \
--i-data qiime2-read-joining-tutorial/demux-full_2-joined.qza \
--o-visualization qiime2-read-joining-tutorial/demux-full_2-joined.qzv

#QUALITY FILTERING (referenced Qiita processing info on split libraries)
qiime quality-filter q-score-joined \
--i-demux qiime2-read-joining-tutorial/demux-full_2-joined.qza \
--p-min-quality 3 \
--p-quality-window 3 \
--p-min-length-fraction 0.75 \
--p-max-ambiguous 0 \
--o-filtered-sequences qiime2-read-joining-tutorial/demux-full_2-joined-filtered.qza \
--o-filter-stats qiime2-read-joining-tutorial/demux-full_2-joined-filter-stats.qza

qiime deblur denoise-16S \
--i-demultiplexed-seqs qiime2-read-joining-tutorial/demux-full_2-joined-filtered.qza \
--p-trim-length 150 \
--p-sample-stats \
--p-mean-error 0.005 \
--p-indel-prob 0.01 \
--p-indel-max 3 \
--p-min-reads 0 \
--p-min-size 2 \
--p-jobs-to-start 5 \
--o-representative-sequences rep-seqs_2.qza \
--o-table table_2.qza \
--o-stats deblur-stats_2.qza

#Visualize deblur
qiime metadata tabulate \
--m-input-file qiime2-read-joining-tutorial/demux-full_2-joined-filter-stats.qza \
--o-visualization demux-filter-stats_2.qzv
qiime deblur visualize-stats \
--i-deblur-stats deblur-stats_2.qza \
--o-visualization deblur-stats_2.qzv

------------

%%Dec1-3_emp-paired-end-sequences.qza%%
qiime demux emp-paired \
--i-seqs Dec1-3_emp-paired-end-sequences.qza \
--m-barcodes-file sample-metadata3.txt \
--m-barcodes-column barcode-sequence \
--p-golay-error-correction FALSE \
--o-per-sample-sequences demux-full_3.qza \
--o-error-correction-details demux_3-details.qza
  # Saved SampleData[PairedEndSequencesWithQuality] to: demux-full_3.qza
  # Saved ErrorCorrectionDetails to: demux_3-details.qza

qiime demux summarize \
--i-data demux-full_3.qza \
--o-visualization demux-full_3.qzv
  #Saved Visualization to: demux-full_3.qzv
qiime tools view demux-full_3.qzv


#READ-JOINING TO MATCH END-PAIRS
#preparation before deblur (dada2 would have included this step)
qiime vsearch join-pairs \
--i-demultiplexed-seqs demux-full_3.qza \
--o-joined-sequences qiime2-read-joining-tutorial/demux-full_3-joined.qza
  #Saved SampleData[JoinedSequencesWithQuality] to: qiime2-read-joining-tutorial/demux-full_3-joined.qza
qiime demux summarize \
--i-data qiime2-read-joining-tutorial/demux-full_3-joined.qza \
--o-visualization qiime2-read-joining-tutorial/demux-full_3-joined.qzv
  #Saved Visualization to: qiime2-read-joining-tutorial/demux-full_3-joined.qzv

#QUALITY FILTERING (referenced Qiita processing info on split libraries)
qiime quality-filter q-score-joined \
--i-demux qiime2-read-joining-tutorial/demux-full_3-joined.qza \
--p-min-quality 3 \
--p-quality-window 3 \
--p-min-length-fraction 0.75 \
--p-max-ambiguous 0 \
--o-filtered-sequences qiime2-read-joining-tutorial/demux-full_3-joined-filtered.qza \
--o-filter-stats qiime2-read-joining-tutorial/demux-full_3-joined-filter-stats.qza
  # Saved SampleData[JoinedSequencesWithQuality] to: qiime2-read-joining-tutorial/demux-full_3-joined-filtered.qza
  # Saved QualityFilterStats to: qiime2-read-joining-tutorial/demux-full_3-joined-filter-stats.qza

qiime deblur denoise-16S \
--i-demultiplexed-seqs qiime2-read-joining-tutorial/demux-full_3-joined-filtered.qza \
--p-trim-length 150 \
--p-sample-stats \
--p-mean-error 0.005 \
--p-indel-prob 0.01 \
--p-indel-max 3 \
--p-min-reads 0 \
--p-min-size 2 \
--p-jobs-to-start 5 \
--o-representative-sequences rep-seqs_3.qza \
--o-table table_3.qza \
--o-stats deblur-stats_3.qza
  # Saved FeatureTable[Frequency] to: table_3.qza
  # Saved FeatureData[Sequence] to: rep-seqs_3.qza
  # Saved DeblurStats to: deblur-stats_3.qza


#Visualize deblur
qiime metadata tabulate \
--m-input-file qiime2-read-joining-tutorial/demux-full_3-joined-filter-stats.qza \
--o-visualization demux-filter-stats_3.qzv
  #Saved Visualization to: demux-filter-stats_3.qzv
qiime deblur visualize-stats \
--i-deblur-stats deblur-stats_3.qza \
--o-visualization deblur-stats_3.qzv
  #Saved Visualization to: deblur-stats_3.qzv


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### HMB Dataset
# All of the below commands assume you are in the following directory
source activate qiime2-2020.8

cd ~/Desktop/Patras_Lab/Humanized_Models/HMB_BrittonMice/RawSequences/Patras_x
#prep the manifest with sample_id/absolute-filepath/direction

cd /Users/owner/Desktop/Patras_Lab/BCM_16Sseq_and_metadata/RawSequencesZip

#assuming all the folders have already been opened/unzipped from prior individual analyses
# convert to .gz
for f in `ls RawSequences_HMB_estrousGBS2022`; do
bzcat RawSequences_HMB_estrousGBS2022/$f/*.1.fq.bz2 | gzip -c > RawSequences_HMB_estrousGBS2022/$f/$f.1.fq.gz;
bzcat RawSequences_HMB_estrousGBS2022/$f/*.2.fq.bz2 | gzip -c > RawSequences_HMB_estrousGBS2022/$f/$f.2.fq.gz;
done

#HUVAMI_1 samples are going to be pulled from their processed filepath in  Humanized_Models/HMB_BrittonMice
#skip converting from .bz2 to gz since already done - jump to manifest - included in manifest where the extension will be changed in the following block of code

# Remove the bz2 compressed files, to reduce duplication
# this command will give you a list of all the files that will be found with the "find" command. These are the files that will then be pipped into the remove command
find RawSequences_HMB_estrousGBS2022 -name "*.bz2" | less
# Use Ctrl+Z or Ctrl+C to exit list back to command line
#Actual remove command (only do if you are sure the list being pipped has the correct files). "xargs" is a tool that can read a list of e.g. filenames and do the same command on them at once - kind of like a loop. https://ss64.com/osx/xargs.html
find RawSequences_HMB_estrousGBS2022 -name "*.bz2" | xargs rm


# update your manifest file (“sed” is the commandline search and replace)
cd /Users/owner/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/seqanalysis_final
cat HMBpaper_finalmanifest_V1.txt | sed "s/bz2/gz/g" > HMBpaper_finalmanifest_V2.txt

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path HMBpaper_finalmanifest_V2.txt \
--output-path paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2
#Imported murine-manifest-V2.txt as PairedEndFastqManifestPhred33V2 to paired-end-demux.qza

qiime tools peek paired-end-demux.qza
# UUID:        c11565f3-3ae6-46b9-b31d-8ed650ceeb4f
# Type:        SampleData[PairedEndSequencesWithQuality]
# Data format: SingleLanePerSamplePairedEndFastqDirFmt

#just to view
mkdir Output
mkdir Visualization

qiime demux summarize \
--i-data paired-end-demux.qza \
--o-visualization Visualization/demux.qzv
#Saved Visualization to: Visualization/demux.qzv

qiime tools view Visualization/demux.qzv
mkdir qiime2-read-joining-tutorial

#READ-JOINING TO MATCH END-PAIRS
qiime vsearch join-pairs \
--i-demultiplexed-seqs paired-end-demux.qza \
--o-joined-sequences qiime2-read-joining-tutorial/demux-joined.qza
#Saved SampleData[JoinedSequencesWithQuality] to: qiime2-read-joining-tutorial/demux-joined.qza

qiime demux summarize \
--i-data qiime2-read-joining-tutorial/demux-joined.qza \
--o-visualization qiime2-read-joining-tutorial/demux-joined.qzv

#QUALITY FILTERING (referenced Qiita processing info on split libraries)
qiime quality-filter
from qiime2.plugins import quality_filter
#changed from quality-filterq-score-joined

qiime quality-filter q-score \
--i-demux qiime2-read-joining-tutorial/demux-joined.qza \
--p-min-quality 3 \
--p-quality-window 3 \
--p-min-length-fraction 0.75 \
--p-max-ambiguous 0 \
--o-filtered-sequences qiime2-read-joining-tutorial/demux-joined-filtered.qza \
--o-filter-stats qiime2-read-joining-tutorial/demux-joined-filter-stats.qza
# Saved SampleData[JoinedSequencesWithQuality] to: qiime2-read-joining-tutorial/demux-joined-filtered.qza
# Saved QualityFilterStats to: qiime2-read-joining-tutorial/demux-joined-filter-stats.qza

qiime deblur denoise-16S \
--i-demultiplexed-seqs qiime2-read-joining-tutorial/demux-joined-filtered.qza \
--p-trim-length 150 \
--p-sample-stats \
--p-mean-error 0.005 \
--p-indel-prob 0.01 \
--p-indel-max 3 \
--p-min-reads 0 \
--p-min-size 2 \
--p-jobs-to-start 5 \
--o-representative-sequences rep-seqs_HMB_final.qza \
--o-table table_HMB_final.qza \
--o-stats deblur-stats_HMB_final.qza
# Saved FeatureTable[Frequency] to: table_GBSPBS.qza
# Saved FeatureData[Sequence] to: rep-seqs_GBSPBS.qza
# Saved DeblurStats to: deblur-stats.qza

#Visualize demux joined filter stats and deblur
qiime metadata tabulate \
--m-input-file qiime2-read-joining-tutorial/demux-joined-filter-stats.qza \
--o-visualization Visualization/demux-joined-filter-stats.qzv

qiime deblur visualize-stats \
--i-deblur-stats deblur-stats_HMB_final.qza \
--o-visualization Visualization/deblur-stats.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqs_HMB_final.qza \
--o-visualization Visualization/rep-seqs_HMB_final.qzv


##################### MERGE SEQUENCES ##run 06/22/22
cd Desktop/Patras_Lab #changed so I can reach into older files not in my paper directory -aka from HMO trials
mkdir Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported
mkdir Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output
mkdir Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization

#merging all three studies
qiime feature-table merge \
--i-tables Study_files16S/table_1.qza \
--i-tables Study_files16S/table_2.qza \
--i-tables Study_files16S/table_3.qza \
--i-tables Proposals/Manuscripts/Mine/2022_HMB_characterization/seqanalysis_final/table_HMB_final.qza \
--o-merged-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_feattable

#merging rep-seqs
qiime feature-table merge-seqs \
--i-data Study_files16S/rep-seqs_1.qza \
--i-data Study_files16S/rep-seqs_2.qza \
--i-data Study_files16S/rep-seqs_3.qza \
--i-data Proposals/Manuscripts/Mine/2022_HMB_characterization/seqanalysis_final/rep-seqs_HMB_final.qza \
--o-merged-data Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_seqs.qza

qiime tools export \
--input-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_seqs.qza \
--output-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/merged_rep-seqs

####Re-train classifier
mkdir training-feature-classifiers
cd training-feature-classifiers

#used 99_otu_taxonomy.txt and 99_otus.fasta (instead of aligned seq.fasta.qz) from downloaded file (from online)
#open in Excel, add headings "Feature ID" and "Taxon", then move to 16_Files_Study3/Qiime folder

qiime tools import  \
--type FeatureData[Sequence] \
--input-path 99_otus.fasta \
--output-path 99_otus.qza
#Imported 99_otus.fasta as DNASequencesDirectoryFormat to 99_otus.qza

qiime tools import \
--type FeatureData[Taxonomy] \
--input-path 99_otu_taxonomy.txt \
--output-path ref-taxonomy.qza
#Imported 99_otu_taxonomy.txt as TSVTaxonomyDirectoryFormat to ref-taxonomy.qza

qiime feature-classifier extract-reads \
--i-sequences 99_otus.qza -\
-p-f-primer GTGCCAGCMGCCGCGGTAA \
--p-r-primer GGACTACHVGGGTWTCTAAT \
--o-reads ref-seqs.qza
#Saved FeatureData[Sequence] to: ref-seqs.qza

#TRAIN THE CLASSIFIER!
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy ref-taxonomy.qza \
--o-classifier gg-13-8-99-classifier.qza
#Saved TaxonomicClassifier to: gg-13-8-99-classifier.qza

cd Patras_Lab #cd ../

#Generate taxonomy table
qiime feature-classifier classify-sklearn \
--i-classifier training-feature-classifiers/gg-13-8-99-classifier.qza \
--i-reads Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_seqs.qza \
--o-classification Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/mergedstudy_taxonomy.qza

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_seqs.qza \
--o-alignment Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/training-feature-classifiersaligned-mergedstudy_seqs.qza \
--o-masked-alignment Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/masked-aligned-mergedstudy_seqs.qza \
--o-tree Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/unrooted-tree_mergedstudy.qza \
--o-rooted-tree Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/rooted-tree_mergedstudy.qza
  # Saved FeatureData[AlignedSequence] to: aligned-rep-seqs_3.qza
  # Saved FeatureData[AlignedSequence] to: masked-aligned-rep-seqs_3.qza
  # Saved Phylogeny[Unrooted] to: unrooted-tree_3.qza
  # Saved Phylogeny[Rooted] to: rooted-tree_3.qza
      #rooted tree is same as insertion tree? Will use as such for now

qiime tools export \
--input-path rooted-tree_mergedstudy.qza \
--output-path exported/rooted-tree_mergedstudy
  #Exported rooted-tree_3.qza as NewickDirectoryFormat to directory exported/rooted-tree


##EXPORTING##
qiime tools export \
--input-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_feattable.qza \
--output-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/merged_feattable

biom convert -i Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/merged_feattable/feature-table.biom -o Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/merged_feattable/OTU_table.txt --to-tsv
#added taxonomy in txt and am converting back
biom convert -i Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/merged_feattable/OTU_table.txt -o Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/merged_feattable/OTU_mergedfeattable_json.biom --table-type="OTU table" --to-json

qiime tools export \
--input-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_seqs.qza \
--output-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/merged_rep-seqs
      
      qiime feature-table tabulate-seqs \
      --i-data mergedstudy_seqs.qza \
      --o-visualization Visualization/mergedstudy_seqs.qzv

qiime tools export \
--input-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/mergedstudy_taxonomy.qza \
--output-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/mergedstudy_taxonomy

qiime metadata tabulate \
--m-input-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/mergedstudy_taxonomy.qza \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/mergedstudy_taxonomy

##         ##

cd Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete

qiime feature-table filter-samples \
--i-table mergedstudy_feattable.qza \
--m-metadata-file samples-to-keep_mergedstudy.txt \
--o-filtered-table subset_mergedstudy_feattable.qza

#back to working on the feature table
qiime feature-table summarize \
--i-table subset_mergedstudy_feattable.qza \
--o-visualization Visualization/mergedstudy_Final_feattable.qzv \
--m-sample-metadata-file HMBpaper_metadata_FINAL.txt

#samples to keep for figures will come later, post filter

########## FILTER by min samples ##########
mkdir exported/feattable_min7samples

#filter table for feat IDs that appear in 7+ samples. Then will proceed to add taxonomic names
qiime feature-table filter-features \
--i-table subset_mergedstudy_feattable.qza \
--p-min-samples 7 \
--o-filtered-table subset_mergedstudy_feattable_min7samples.qza

#OVERALL peek into table to see what taxa are in blanks and need to be removed
qiime tools export \
--input-path mergedstudy_feattable.qza \
--output-path exported/feattable
biom convert -i exported/feattable/feature-table.biom -o exported/feattable/feat_table.txt --to-tsv

#peek into new SUBSET table to see what OTUs are present
qiime tools export \
--input-path subset_mergedstudy_feattable.qza \
--output-path exported/SUBSET_feattable
biom convert -i exported/SUBSET_feattable/feature-table.biom -o exported/SUBSET_feattable/feat_table.txt --to-tsv

#export reduced table based on features in 3+ samples
qiime tools export \
--input-path subset_mergedstudy_feattable_min7samples.qza \
--output-path exported/feattable_min7samples

biom convert -i exported/feattable_min7samples/feature-table.biom -o exported/feattable_min7samples/feat_table.txt --to-tsv --header-key taxonomy
#added taxonomy in txt and am converting back
#biom convert -i exported/feattable_min7samples/feat_table_taxa.txt -o exported/feattable_min7samples/feat_json.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy

biom convert -i feat_taxatable.txt -o feat_taxatable_json.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy


#------------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%% R studio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
######### FILTER OUT CONTAMINANTS USING DECONTAM #########
  
## ----loadPS----------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")

setwd("~")
#upload the matrix into phyloseq format
#Go in and add prefix "ID_" to all feature ID names in the taxonomy file and in the OTU tables and such. This way they can be matched later for taxa removal or merging
#feature table read with feature IDs and not taxonomy names
#make sure "taxonomy" column (now last row) is deleted
tFeatureTable <- read.csv("Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/feattable_min7samples/feat_table_matrix.csv", header=TRUE, sep = ",", row.names="sample_name")
metadata_BC <- read.csv("Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL.csv", header=TRUE, sep = ",", row.names="sample_name")
taxonomy <- read.csv("Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/mergedstudy_taxonomy/taxonomy_ID.csv", header=TRUE, sep = ",", row.names="Feature_ID")
tTaxonomy <- t(taxonomy)

#create directory "Decontam_process"
setwd("Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process")
dir.create("decontam_items_final")
setwd("~/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final")
OTUmatrixfilt = phyloseq(otu_table(tFeatureTable,taxa_are_rows = FALSE), sample_data(metadata_BC))


## ----see-meta-table--------------------------------------------------------
head(sample_data(OTUmatrixfilt))

dffilt <- as.data.frame(sample_data(OTUmatrixfilt)) # Put sample_data into a ggplot-friendly data.frame
dffilt$LibrarySize <- sample_sums(OTUmatrixfilt)
dffilt <- dffilt[order(dffilt$LibrarySize),]
dffilt$Index <- seq(nrow(dffilt))
ggplot(data=dffilt, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()


sample_data(OTUmatrixfilt)$is.neg <- sample_data(OTUmatrixfilt)$Sample_or_Control == "Control Sample"
# Make phyloseq object of presence-absence in negative controls and true samples
OTUmatrixfilt.pa <- transform_sample_counts(OTUmatrixfilt, function(abund) 1*(abund>0))
OTUmatrixfilt.pa.neg <- prune_samples(sample_data(OTUmatrixfilt.pa)$Sample_or_Control == "Control Sample", OTUmatrixfilt.pa)
OTUmatrixfilt.pa.pos <- prune_samples(sample_data(OTUmatrixfilt.pa)$Sample_or_Control == "True Sample", OTUmatrixfilt.pa)


## ----prevalence (0.1?)------------------------------------------------------------
contamdffilt.prev <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg")
table(contamdffilt.prev$contaminant)
head(which(contamdffilt.prev$contaminant))
# Make data.frame of prevalence in positive and negative samples
dffilt.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                              contaminant=contamdffilt.prev$contaminant)
ggplot(data=dffilt.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#withbatch
contamdffilt.prev01b <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.1, batch= "batch", batch.combine = "minimum")
table(contamdffilt.prev01b$contaminant)
dffilt01b.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                contaminant=contamdffilt.prev01b$contaminant)
ggplot(data=dffilt01b.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.1b)") + ylab("Prevalence (True Samples)")

## ----see-prev-.02-----------------------------------------------------------
contamdffilt.prev02 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.2, batch= "batch", batch.combine = "minimum")
table(contamdffilt.prev02$contaminant)
dffilt02.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                contaminant=contamdffilt.prev02$contaminant)
ggplot(data=dffilt02.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.2)") + ylab("Prevalence (True Samples)")

## ----see-prev-.03-----------------------------------------------------------
contamdffilt.prev03 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.3, batch= "batch", batch.combine = "minimum")
table(contamdffilt.prev03$contaminant)
dffilt03.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                contaminant=contamdffilt.prev03$contaminant)
ggplot(data=dffilt03.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.3)") + ylab("Prevalence (True Samples)")

## ----see-prev-.04-----------------------------------------------------------
contamdffilt.prev04 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.4, batch= "batch", batch.combine = "minimum")
table(contamdffilt.prev04$contaminant)
dffilt04.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                contaminant=contamdffilt.prev04$contaminant)
ggplot(data=dffilt04.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.4)") + ylab("Prevalence (True Samples)")

## ----see-prev-.45-----------------------------------------------------------
contamdffilt.prev04_5 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.45, batch= "batch", batch.combine = "minimum")
table(contamdffilt.prev04_5$contaminant)
dffilt04_5.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                  contaminant=contamdffilt.prev04_5$contaminant)
ggplot(data=dffilt04_5.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.45)") + ylab("Prevalence (True Samples)")

## ----prevalence-.05---------------------------------------------------------
contamdffilt.prev05 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.5, batch= "batch", batch.combine = "minimum")
table(contamdffilt.prev05$contaminant)
dffilt05.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                  contaminant=contamdffilt.prev05$contaminant)
ggplot(data=dffilt05.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.5)") + ylab("Prevalence (True Samples)")


write.csv(contamdffilt.prev01b,"contamdffilt.prev01b.csv", row.names = TRUE)
write.csv(contamdffilt.prev,"contamdffilt.prev01.csv", row.names = TRUE)
write.csv(contamdffilt.prev02,"contamdffilt.prev02.csv", row.names = TRUE)
write.csv(contamdffilt.prev03,"contamdffilt.prev03.csv", row.names = TRUE)
write.csv(contamdffilt.prev04,"contamdffilt.prev04.csv", row.names = TRUE)
write.csv(contamdffilt.prev04_5,"contamdffilt.prev04_5.csv", row.names = TRUE)
write.csv(contamdffilt.prev05,"contamdffilt.prev05.csv", row.names = TRUE)

write.csv(dffilt01b.pa,"dffilt01b.pa.csv", row.names = TRUE)
write.csv(dffilt.pa,"dffilt01.pa.csv", row.names = TRUE)
write.csv(dffilt02.pa,"dffilt02.pa.csv", row.names = TRUE)
write.csv(dffilt03.pa,"dffilt03.pa.csv", row.names = TRUE)
write.csv(dffilt04.pa,"dffilt04.pa.csv", row.names = TRUE)
write.csv(dffilt04_5.pa,"dffilt04_5.pa.csv", row.names = TRUE)
write.csv(dffilt05.pa,"dffilt05.pa.csv", row.names = TRUE)



#combine above tables to match everything by Feature ID
taxa_dffilt05 <- merge(dffilt05.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt05)[1] <- "Feature_ID"
taxa_contam05 <- merge(contamdffilt.prev05, taxonomy, by= 'row.names')
colnames(taxa_contam05)[1] <- "Feature_ID"
taxa_contam_dffilt05 <- merge(taxa_dffilt05, taxa_contam05, by= 'Feature_ID')
write.csv(taxa_contam_dffilt05,"taxa_contam_dffilt05.csv", row.names = FALSE)

taxa_dffilt04_5 <- merge(dffilt04_5.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt04_5)[1] <- "Feature_ID"
taxa_contam04_5 <- merge(contamdffilt.prev04_5, taxonomy, by= 'row.names')
colnames(taxa_contam04_5)[1] <- "Feature_ID"
taxa_contam_dffilt04_5 <- merge(taxa_dffilt04_5, taxa_contam04_5, by= 'Feature_ID')
write.csv(taxa_contam_dffilt04_5,"taxa_contam_dffilt04_5.csv", row.names = FALSE)

taxa_dffilt04 <- merge(dffilt04.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt04)[1] <- "Feature_ID"
taxa_contam04 <- merge(contamdffilt.prev04, taxonomy, by= 'row.names')
colnames(taxa_contam04)[1] <- "Feature_ID"
taxa_contam_dffilt04 <- merge(taxa_dffilt04, taxa_contam04, by= 'Feature_ID')
write.csv(taxa_contam_dffilt04,"taxa_contam_dffilt04.csv", row.names = FALSE)

taxa_dffilt03 <- merge(dffilt03.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt03)[1] <- "Feature_ID"
taxa_contam03 <- merge(contamdffilt.prev03, taxonomy, by= 'row.names')
colnames(taxa_contam03)[1] <- "Feature_ID"
taxa_contam_dffilt03 <- merge(taxa_dffilt03, taxa_contam03, by= 'Feature_ID')
write.csv(taxa_contam_dffilt03,"taxa_contam_dffilt03.csv", row.names = FALSE)

taxa_dffilt02 <- merge(dffilt02.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt02)[1] <- "Feature_ID"
taxa_contam02 <- merge(contamdffilt.prev02, taxonomy, by= 'row.names')
colnames(taxa_contam02)[1] <- "Feature_ID"
taxa_contam_dffilt02 <- merge(taxa_dffilt02, taxa_contam02, by= 'Feature_ID')
write.csv(taxa_contam_dffilt02,"taxa_contam_dffilt02.csv", row.names = FALSE)

taxa_dffilt <- merge(dffilt.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt)[1] <- "Feature_ID"
taxa_contam <- merge(contamdffilt.prev, taxonomy, by= 'row.names')
colnames(taxa_contam)[1] <- "Feature_ID"
taxa_contam_dffilt <- merge(taxa_dffilt, taxa_contam, by= 'Feature_ID')
write.csv(taxa_contam_dffilt,"taxa_contam_dffilt.csv", row.names = FALSE)

taxa_dffilt01b <- merge(dffilt01b.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt01b)[1] <- "Feature_ID"
taxa_contam01b <- merge(contamdffilt.prev01b, taxonomy, by= 'row.names')
colnames(taxa_contam01b)[1] <- "Feature_ID"
taxa_contam_dffilt01b <- merge(taxa_dffilt01b, taxa_contam01b, by= 'Feature_ID')
write.csv(taxa_contam_dffilt01b,"taxa_contam_dffilt01b.csv", row.names = FALSE)


#Selected threshold- remove these feature IDs from the matrix and then write a transposed version so that samples are rows
OTUmatrixfilt
OTUmatrixfilt.noncontam02 <- prune_taxa(!contamdffilt.prev02$contaminant, OTUmatrixfilt)
OTUmatrixfilt.noncontam02

OTUmatrixfiltDecontam02 <-t(as(otu_table(OTUmatrixfilt.noncontam02),"matrix"))
write.csv(OTUmatrixfiltDecontam02,"OTUmatrixfiltDecontam02.csv", row.names=TRUE)
#open file and add column name "#OTU ID" and use replace function to remove "ID_"




#------------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%% Qiime2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
########## FILTER OUT SPECIFIC TAXA ###########
cd /Users/owner/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete

#Convert to qza file
biom convert -i Decontam_process/decontam_items_final/OTUmatrixfiltDecontam02_usethis.txt -o Decontam_process/decontam_items_final/OTUmatrixfiltDecontam02_usethis.biom --table-type="OTU table" --to-json

mkdir Decontam_process/decontam_items_final/Output
    
qiime tools import \
--input-path Decontam_process/decontam_items_final/OTUmatrixfiltDecontam02_usethis.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path Decontam_process/decontam_items_final/Output/OTU_filtdecontam.qza

#re-import new table to qiime2 for filtering of known contaminants
qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/OTU_filtdecontam.qza \
--i-taxonomy Output/mergedstudy_taxonomy.qza \
--p-mode contains \
--p-exclude "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__Pseudomonas; s__veronii" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza \
--i-taxonomy Output/mergedstudy_taxonomy.qza \
--p-mode contains \
--p-exclude "k__Bacteria; p__Cyanobacteria; c__Chloroplast; o__Streptophyta;" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza \
--i-taxonomy Output/mergedstudy_taxonomy.qza \
--p-mode contains \
--p-exclude "k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Geobacillus; s__" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza \
--i-taxonomy Output/mergedstudy_taxonomy.qza \
--p-mode contains \
--p-exclude "Thermus" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza \
--i-taxonomy Output/mergedstudy_taxonomy.qza \
--p-mode contains \
--p-exclude "f__Phyllobacteriaceae" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza \
--i-taxonomy Output/mergedstudy_taxonomy.qza \
--p-mode contains \
--p-exclude "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Bradyrhizobiaceae; g__Bradyrhizobium" \
--o-filtered-table Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt.qza

#--p-exclude "k__Bacteria; p__Cyanobacteria; c__Chloroplast; o__Streptophyta; f__; g__; s__" \
#--p-exclude "k__Bacteria; p__Chloroflexi; c__Chloroflexi; o__Chloroflexales; f__Chloroflexaceae; g__Chloronema; s__" \

qiime tools export \
--input-path Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt.qza \
--output-path Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt
biom convert -i Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt/feature-table.biom -o Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt/feature-table.biom.txt --to-tsv
biom convert -i Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt/feature-taxatable.biom.txt -o Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt/feature-taxatable_json.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy


mkdir subset_GBS
mkdir subset_estrousstaging
mkdir subset_longitudinal
mkdir subset_vivarium


mkdir Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results
qiime diversity core-metrics-phylogenetic \
--i-phylogeny Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/rooted-tree_mergedstudy.qza \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt.qza \
--p-sampling-depth '1' \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_code.txt \
--o-unweighted-unifrac-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/unweighted_unifrac_distancematrix.qza \
--o-weighted-unifrac-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_distancematrix.qza \
--o-unweighted-unifrac-pcoa-results Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/unweighted_unifrac_pcoa.qza \
--o-weighted-unifrac-pcoa-results Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_pcoa.qza \
--o-unweighted-unifrac-emperor Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/unweighted_unifrac_emperor.qza \
--o-weighted-unifrac-emperor Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_emperor.qza \
--output-dir Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results_unspecified

qiime tools export \
--input-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_distancematrix.qza \
--output-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_distancematrix.txt

qiime diversity beta-group-significance \
--i-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_distancematrix.qza \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_code.txt \
--m-metadata-column site_line \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/HMBfinal_WeightedUnifrac_Beta_diversity_permanova_siteline
#representation, birthplace, mouse_line, CST_categorical, estrous_stage, vivarium

    # qiime diversity beta-group-significance \
    # --i-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_distancematrix.qza \
    # --m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_GBSCSTCFU_TissueSubset.txt \
    # --m-metadata-column CST_categorical \
    # --p-method permanova \
    # --p-pairwise True \
    # --p-permutations 999 \
    # --o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/HMBfinal_WeightedUnifrac_Beta_diversity_permanova_GBSCFUCST_TissueSubset_categorical
    # 
    # qiime diversity beta-group-significance \
    # --i-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_distancematrix.qza \
    # --m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
    # --m-metadata-column estrous_stage \
    # --p-method permanova \
    # --p-pairwise True \
    # --p-permutations 999 \
    # --o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/HMBfinal_WeightedUnifrac_Beta_diversity_permanova_estrousstaging
    # 
    # qiime diversity beta-group-significance \
    # --i-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_distancematrix.qza \
    # --m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_tissue_sites/metadata_FPVS.txt \
    # --m-metadata-column representation \
    # --p-method permanova \
    # --p-pairwise True \
    # --p-permutations 999 \
    # --o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/HMBfinal_WeightedUnifrac_Beta_diversity_permanova_representation
    # 
    # qiime diversity beta-group-significance \
    # --i-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_distancematrix.qza \
    # --m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_VivariumHMB/metadata_VivariumHMB.txt \
    # --m-metadata-column site_line \
    # --p-method permanova \
    # --p-pairwise True \
    # --p-permutations 999 \
    # --o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/HMBfinal_WeightedUnifrac_Beta_diversity_permanova_vivariumHMB
    # 
    # qiime diversity beta-group-significance \
    # --i-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_distancematrix.qza \
    # --m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_VivariumHMB/metadata_VivariumHMB.txt \
    # --m-metadata-column vivarium \
    # --p-method permanova \
    # --p-pairwise True \
    # --p-permutations 999 \
    # --o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/HMBfinal_WeightedUnifrac_Beta_diversity_permanova_vivariumHMB
    # 
    # qiime diversity beta-group-significance \
    # --i-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/core-metrics-results/weighted_unifrac_distancematrix.qza \
    # --m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_VivariumHMB/metadata_VivariumHMB.txt \
    # --m-metadata-column site_line \
    # --p-method permanova \
    # --p-pairwise True \
    # --p-permutations 999 \
    # --o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/HMBfinal_WeightedUnifrac_Beta_diversity_permanova_vivariumHMB
        #doesn't allow because of missing IDs in the metadata
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/Beta_diversity_permanova.qzv


#####create tables per study/figure type
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt.qza \
--m-metadata-file samples-to-keep_ESTROUSGBSd0_HMB.txt \
--o-filtered-table subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0.qza
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt.qza \
--m-metadata-file samples-to-keep_ESTROUSsubset.txt \
--o-filtered-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt.qza \
--m-metadata-file samples-to-keep_GBSsubset.txt \
--o-filtered-table subset_GBS/OTU_filtdecontamfilt_GBS.qza
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt.qza \
--m-metadata-file samples-to-keep_LONGITUDINALsubset.txt \
--o-filtered-table subset_longitudinal/OTU_filtdecontamfilt_longitudinal.qza
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt.qza \
--m-metadata-file samples-to-keep_VIVARIUMsubset.txt \
--o-filtered-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza #adjusted to keep all parity but only day 3, and only study 3 P11-25 and T11-25 since no GBS in run

%%bash -e
qiime feature-table relative-frequency \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0.qza \
--o-relative-frequency-table subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0_relfreq
#Saved FeatureTable[RelativeFrequency] to: subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq.qza
%%bash -e
qiime feature-table relative-frequency \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--o-relative-frequency-table subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq
#Saved FeatureTable[RelativeFrequency] to: subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq.qza
%%bash -e
qiime feature-table relative-frequency \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--o-relative-frequency-table subset_GBS/OTU_filtdecontamfilt_GBS_relfreq
#Saved FeatureTable[RelativeFrequency] to: subset_GBS/OTU_filtdecontamfilt_GBS_relfreq.qza
%%bash -e
qiime feature-table relative-frequency \
--i-table subset_longitudinal/OTU_filtdecontamfilt_longitudinal.qza \
--o-relative-frequency-table subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq
#Saved FeatureTable[RelativeFrequency] to: subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq.qza

qiime feature-table relative-frequency \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--o-relative-frequency-table subset_vivarium/OTU_filtdecontamfilt_vivarium_relfreq
#Saved FeatureTable[RelativeFrequency] to: subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq.qza


#export feature tables (reads and rel_frequency)
qiime tools export \
--input-path subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0.qza \
--output-path subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0
biom convert -i subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0/feature-table.biom -o subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0/filtfeat_table_estrousstaging_d0.txt --to-tsv

qiime tools export \
--input-path subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--output-path subset_estrousstaging/OTU_filtdecontamfilt_estrous
biom convert -i subset_estrousstaging/OTU_filtdecontamfilt_estrous/feature-table.biom -o subset_estrousstaging/OTU_filtdecontamfilt_estrous/filtfeat_table_estrousstaging.txt --to-tsv

qiime tools export \
--input-path subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--output-path subset_GBS/OTU_filtdecontamfilt_GBS
biom convert -i subset_GBS/OTU_filtdecontamfilt_GBS/feature-table.biom -o subset_GBS/OTU_filtdecontamfilt_GBS/filtfeat_table_GBS.txt --to-tsv

qiime tools export \
--input-path subset_longitudinal/OTU_filtdecontamfilt_longitudinal.qza \
--output-path subset_longitudinal/OTU_filtdecontamfilt_longitudinal
biom convert -i subset_longitudinal/OTU_filtdecontamfilt_longitudinal/feature-table.biom -o subset_longitudinal/OTU_filtdecontamfilt_longitudinal/filtfeat_table_longitudinal.txt --to-tsv

qiime tools export \
--input-path subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--output-path subset_vivarium/OTU_filtdecontamfilt_vivarium
biom convert -i subset_vivarium/OTU_filtdecontamfilt_vivarium/feature-table.biom -o subset_vivarium/OTU_filtdecontamfilt_vivarium/filtfeat_table_vivarium.txt --to-tsv



qiime tools export \
--input-path subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0_relfreq.qza \
--output-path subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0_relfreq
biom convert -i subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0_relfreq/feature-table.biom -o subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0_relfreq/filtfeat_table_estrousstaging_d0_relfreq.txt --to-tsv
##
qiime tools export \
--input-path subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq.qza \
--output-path subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq
biom convert -i subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq/feature-table.biom -o subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq/filtfeat_table_estrousstaging_relfreq.txt --to-tsv
#Saved FeatureTable[RelativeFrequency] to: subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq.qza
#Exported subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq.qza as BIOMV210DirFmt to directory subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq

qiime tools export \
--input-path subset_GBS/OTU_filtdecontamfilt_GBS_relfreq.qza \
--output-path subset_GBS/OTU_filtdecontamfilt_GBS_relfreq
biom convert -i subset_GBS/OTU_filtdecontamfilt_GBS_relfreq/feature-table.biom -o subset_GBS/OTU_filtdecontamfilt_GBS_relfreq/filtfeat_table_GBS_relfreq.txt --to-tsv
##
qiime tools export \
--input-path subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq.qza \
--output-path subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq
biom convert -i subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq/feature-table.biom -o subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq/filtfeat_table_longitudinal_relfreq.txt --to-tsv

qiime tools export \
--input-path subset_vivarium/OTU_filtdecontamfilt_vivarium_relfreq.qza \
--output-path subset_vivarium/OTU_filtdecontamfilt_vivarium_relfreq
biom convert -i subset_vivarium/OTU_filtdecontamfilt_vivarium_relfreq/feature-table.biom -o subset_vivarium/OTU_filtdecontamfilt_vivarium_relfreq/filtfeat_table_vivarium_relfreq.txt --to-tsv


#to convert back into biom
biom convert -i exported/filt_HMBestrous_final_relfreq/feature-table.biom -o exported/filt_HMBestrous_final_relfreq/filtfeat_table.txt --to-tsv --header-key taxonomy
biom convert -i subset_vivarium/OTU_filtdecontamfilt_vivarium/filtfeat_taxatable_vivarium.txt -o subset_vivarium/OTU_filtdecontamfilt_vivarium/filtfeat_taxatable_vivarium_json.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy



######## STOPPED CODE HERE ###############
######## CONTINUE DIVERSITY ANALYSES BELOW OR SKIP TO CST USING EXPORTED FILES ###############
######## BEGIN ESTROUS HEATMAP GENERATION in R ############


%%%%%%%%other diversity analyses%%%%%%%%
  #need to change file names accordingly
  
mkdir subset_GBS/Output
mkdir subset_GBS/Visualization

  
qiime feature-table summarize \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--m-sample-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--o-visualization subset_GBS/Visualization/Summary_HMBfinal_GBS_featuretable_filt.qzv
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/Summary_featuretableMM.qzv


%%bash -e
qiime taxa barplot \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--i-taxonomy Output/mergedstudy_taxonomy.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_barplot
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/_barplot.qzv

%%%%beta%%%%
  
  %%bash -e
qiime diversity beta \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix
#Saved DistanceMatrix to: subset_GBS/Output/_BCDmatrix.qza

%%bash -e
qiime diversity pcoa \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--p-number-of-dimensions 3 \
--o-pcoa subset_GBS/Output/HMBfinal_GBS__PCoA_BCD
#Saved PCoAResults to: subset_GBS/Output/_PCoA_BCD.qza

%%bash -e
qiime emperor plot \
--i-pcoa subset_GBS/Output/HMBfinal_GBS__PCoA_BCD.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--p-ignore-missing-samples False \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS__PCoA_BCD
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/_PCoA_BCD.qzv

%%bash -e
qiime feature-table relative-frequency \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--o-relative-frequency-table subset_GBS/Output/HMBfinal_GBS__relative
#Saved FeatureTable[RelativeFrequency] to: subset_GBS/Output/_relative.qza

%%bash -e
qiime diversity pcoa-biplot \
--i-pcoa subset_GBS/Output/HMBfinal_GBS__PCoA_BCD.qza \
--i-features subset_GBS/Output/HMBfinal_GBS__relative.qza \
--o-biplot subset_GBS/Output/HMBfinal_GBS_pcoa_BCD_biplot 
#Saved PCoAResults % Properties('biplot') to: subset_GBS/Output/pcoa_BCD_biplot.qza

%%bash -e
qiime emperor biplot \
--i-biplot subset_GBS/Output/HMBfinal_GBS_pcoa_BCD_biplot.qza \
--m-sample-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--m-feature-metadata-file Output/mergedstudy_taxonomy.qza \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_pcoaBCD_biplot
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/pcoaBCD_biplot.qzv

%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--m-metadata-column CST \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_Beta_diversity_permanova
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/Beta_diversity_permanova.qzv

%%%%alpha%%%%
  
  %%bash -e
qiime diversity alpha \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--p-metric observed_features \
--o-alpha-diversity subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_otus
#Saved SampleData[AlphaDiversity] to: subset_GBS/Output/HMBfinal_GBS_alpha_otus.qza
#no longer --p-metric observed_otus because updated qiime2

%%bash -e
qiime diversity alpha-correlation \
--i-alpha-diversity subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_otus.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_HMBfinal_GBS_alpha_otus
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/HMBfinal_GBS_alpha_otus.qzv

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_otus.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_HMBfinal_GBS_alpha_otus_usethis
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/HMBfinal_GBS_alpha_otus_usethis.qzv
qiime tools export \
--input-path subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_otus.qza \
--output-path subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_otus



%%bash -e
qiime diversity alpha \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--p-metric shannon \
--o-alpha-diversity subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_shannon
#Saved SampleData[AlphaDiversity] to: subset_GBS/Output/HMBfinal_GBS_alpha_shannon.qza

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_shannon.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_HMBfinal_GBS_alpha_shannon
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/HMBfinal_GBS_alpha_shannon.qzv


# qiime feature-table filter-samples \
# --i-table subset_GBS/Output/filt_idtaxa_feature-table_3.qza \
# --m-metadata-file samples-to-keep_3KO.txt \
# --o-filtered-table subset_GBS/Output/filt_idtaxa-table_3.qza


  %%%%longitudinal%%%%
  
%%bash -e
qiime diversity filter-distance-matrix \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--p-where "mouse_line='Conventional'" \
--o-filtered-distance-matrix subset_GBS/Output/HMBfinal_GBS_filtered-distance-matrix.qza
#Saved DistanceMatrix to: filtered-distance-matrix.qza

%%bash -e
qiime longitudinal first-distances \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS_filtered-distance-matrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--p-state-column timepoint \
--p-individual-id-column mouse_namevf \
--p-replicate-handling error \
--o-first-distances subset_GBS/Output/HMBfinal_GBS__BC_FirstDistances 
#Saved SampleData[FirstDifferences] to: _BC_FirstDistances.qza

qiime tools export \
--input-path subset_GBS/Output/HMBfinal_GBS__BC_FirstDistances.qza \
--output-path exported//HMBfinal_GBS__BC_FirstDistances


%%%%%%%%%%
  
  %%bash -e
qiime longitudinal pairwise-distances \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS_filtered-distance-matrix.qza \
--m-metadata-file murine-metadata_time_HVM.txt \
--p-group-column group_study \
--p-state-column timepoint \
--p-state-1 14 \
--p-state-2 49 \
--p-individual-id-column mouse_namevf \
--p-no-parametric True \
--p-palette Set1 \
--p-replicate-handling error \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_Pairwise_Distances_BC_week2-7
#Saved subset_GBS/Visualization to: Pairwise_Distances_BC_day1-2.qzv

%%bash -e
qiime longitudinal pairwise-distances \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS_filtered-distance-matrix.qza \
--m-metadata-file murine-metadata_time_HVM.txt \
--p-group-column group_study \
--p-state-column timepoint \
--p-state-1 14 \
--p-state-2 77 \
--p-individual-id-column mouse_namevf \
--p-no-parametric True \
--p-palette Set1 \
--p-replicate-handling error \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_Pairwise_Distances_BC_week2-11
#Saved subset_GBS/Visualization to: Pairwise_Distances_BC_day1.qzv

%%bash -e
qiime longitudinal pairwise-distances \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS_filtered-distance-matrix.qza \
--m-metadata-file murine-metadata_time_HVM.txt \
--p-group-column group_study \
--p-state-column timepoint \
--p-state-1 49 \
--p-state-2 77 \
--p-individual-id-column mouse_namevf \
--p-no-parametric True \
--p-palette Set1 \
--p-replicate-handling error \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_Pairwise_Distances_BC_week7-11
#Saved subset_GBS/Visualization to: Pairwise_Distances_BC_day1.qzv

qiime tools export \
--input-path subset_GBS/Output/HMBfinal_GBS__relative_frequency.qza \
--output-path exported/HMBfinal_GBS__relative_frequency

(For abundance plots)

%%bash -e
qiime feature-table relative-frequency \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--o-relative-frequency-table subset_GBS/Output/HMBfinal_GBS__relative_frequency
#Saved FeatureTable[RelativeFrequency] to: _1261_relative_frequency.qza


%%bash -e
qiime tools export \
--input-path subset_GBS/Output/HMBfinal_GBS__relative_frequency.qza \
--output-path exported/HMBfinal_GBS__relative_frequency
#Exported subset_GBS/Output/_1_relative_frequency.qza as BIOMV210DirFmt to directory exported/Study1/_1_relative_frequency

#reach into the exported/HMBfinal_GBS__relative_frequency directory
#then input the taxonomy to the end using advanced filter on excel
%%bash -e
biom convert \
-i feature-table.biom \
-o otu_table.tsv --to-tsv --header-key taxonomy


# qiime taxa barplot \
# --i-table murinecomp_feattable.qza \
# --i-taxonomy subset_GBS/Output/murinecomp_taxonomy.qza \
# --m-metadata-file murinecomp_meta.txt \
# --o-visualization subset_GBS/Visualization/murinecomp_taxonomy_barplot
# this one is not rarefied, rarefied alsready made


#move back to Study 16S
%%bash -e
qiime tools export \
--input-path subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_otus.qza \
--output-path exported/HMBfinal_GBS_HMBfinal_GBS_alpha_otus
#Exported HMBfinal_GBS_alpha_otus.qza as AlphaDiversityDirectoryFormat to directory HMBfinal_GBS_alpha_otus

#move into the exported/HMBfinal_GBS_HMBfinal_GBS_alpha_otus directory
%%bash -e
biom convert \
-i alpha-diversity.tsv \
-o otu_table.tsv --to-tsv --header-key taxonomy

# #move back to Study 16S
# %%bash -e
# qiime tools export \
# --input-path subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_otus.qza \
# --output-path exported/HMBfinal_GBS_HMBfinal_GBS_alpha_otus/alpha-diversity
# #Exported subset_GBS/Output/HMBfinal_GBS_alpha_otus_1.qza as AlphaDiversityDirectoryFormat to directory exported/Study1/alpha-diversity

#move back to Study_16S
%%bash -e
qiime tools export \
--input-path subset_GBS/Output/murinecomp_taxonomy.qza \
--output-path exported/murinecomp_taxonomy/taxonomy
#Exported subset_GBS/Output/taxonomy_1.qza as TSVTaxonomyDirectoryFormat to directory export/Study1/taxonomy

qiime tools export \
--input-path subset_GBS/Output/HMBfinal_GBS_pcoa_BCD_biplot.qza \
--output-path exported/HMBfinal_GBS_pcoa_BCD_biplot.qza

qiime tools export \
--input-path subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--output-path exported/HMBfinal_GBS__BCDmatrix.qza





cd /Users/owner/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete 

qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/OTU_filtdecontamfilt.qza \
--m-metadata-file HMBpaper_metadata_FINAL.txt \
--p-where "[study_topic]='HMB_GBSchallenge'" \
--o-filtered-table HMB_GBSuterine-table.qza

qiime feature-table filter-samples \
--i-table HMB_GBSuterine-table.qza \
--m-metadata-file HMBpaper_metadata_FINAL.txt \
--p-where "[body_site]='vaginal tract'" \
--o-filtered-table HMB_GBSuterinee-VTtable.qza

#might need to run... add HMB before.qza
qiime feature-table filter-samples \
--i-table HMB_GBSuterine-table.qza \
--m-metadata-file HMBpaper_metadata_FINAL.txt \
--p-where "[mouse_line]='HMB'" \
--o-filtered-table HMB_GBSuterine-VTtableHMB.qza

qiime taxa collapse \
--i-table HMB_GBSuterine-VTtableHMB.qza \
--i-taxonomy Output/mergedstudy_taxonomy.qza \
--p-level 6 \
--o-collapsed-table HMB_GBSuterine-collapsedl6tableHMB.qza

qiime composition add-pseudocount \
--i-table HMB_GBSuterine-collapsedl6tableHMB.qza \
--o-composition-table comp-HMB_GBSuterine-collapsedl6tableHMB.qza

qiime composition ancom \
--i-table comp-HMB_GBSuterine-collapsedl6tableHMB.qza \
--m-metadata-file HMBpaper_metadata_FINAL.txt \
--m-metadata-column Uterine_clearance \
--o-visualization l6-ancom-Uterine_ColonizeClear.qzv




#Phylogenetic analyses

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza






#------------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%% R studio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


######################CST-HEATMAP!!!
"Community State Types Baseline"
#perform in R
library(factoextra) # clustering visualization
library(phyloseq) #OTU table handling
library(ggplot2)
#package stats is automatically loaded in R, used for kmeans, wss, and heatmapping

setwd("~/Desktop")

%%%%%%%%%%%%%%% on estrous dataset %%%%%%%%%%%%%%%
#########################################################3
#use relative frequency file, consolidated by taxa and transposed to have sample_name column and taxa headings  
HMBfinalEstrous <- read.csv("Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/OTU_filtdecontamfilt_estrous_d0_relfreq/HMB_estrousstaging_taxahead_allincluded.csv", header=TRUE, sep = ",", row.names="sample_name")
meta_HMBfinalEstrousstages2 <- read.csv("Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/HMBpaper_metadata_estrous.csv", header=TRUE, sep = ",", row.names="sample_name")

setwd("~/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging")
setwd("~/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/EnvisionedFigs/prelimCST") #for quick look without estrous data


dHMBfinalEs <- dist(HMBfinalEstrous, method= "euclidean")
fitHMBfinalEs <- hclust(dHMBfinalEs, method="ward.D")
quartz()
plot(fitHMBfinalEs)

# Determine number of clusters
wssHMBfinalEs <- (nrow(HMBfinalEstrous-1)*sum(apply(HMBfinalEstrous,2,var)))
for (i in 2:15) wssHMBfinalEs[i] <- sum(kmeans(HMBfinalEstrous,
                                              centers=i)$withinss)

plot(1:15, wssHMBfinalEs, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#Return to main code
group3HMBfinalEs <- cutree(fitHMBfinalEs, k=3)
group4HMBfinalEs <- cutree(fitHMBfinalEs, k=4)
group5HMBfinalEs <- cutree(fitHMBfinalEs, k=5)
group6HMBfinalEs <- cutree(fitHMBfinalEs, k=6)
group7HMBfinalEs <- cutree(fitHMBfinalEs, k=7)
group8HMBfinalEs <- cutree(fitHMBfinalEs, k=8)
table(group3HMBfinalEs) # Number of members in each cluster
table(group4HMBfinalEs)
table(group5HMBfinalEs)
table(group6HMBfinalEs)
table(group7HMBfinalEs)
table(group8HMBfinalEs)

plot(fitHMBfinalEs)
rect.hclust(fitHMBfinalEs, k=3, border="yellow")
rect.hclust(fitHMBfinalEs, k=4, border="red")
rect.hclust(fitHMBfinalEs, k=5, border="green")
rect.hclust(fitHMBfinalEs, k=6, border="Blue")
rect.hclust(fitHMBfinalEs, k=7, border="purple")
rect.hclust(fitHMBfinalEs, k=8, border="orange")

print(group5HMBfinalEs)   
write.csv(group5HMBfinalEs,"newCST5_HMBfinalEs_mapping.csv")

#merge the two tables- HMBfinalEstrous and CST5_HMBinitEs_mapping.csv ,- read in as such
#Manually add first column heading "sample_name"
#printed and wrote CST 5 and CST 6 just to see who would be categorized differently
meta_HMBfinalEstrousstudy <- read.csv("newCST5_HMBfinalEs_mapping.csv", header=TRUE, sep = ",", row.names="sample_name")

#below is a metadata file I made one column at a time, but I already loaded the one above in the same order, so I don't need to recreate a metadata sheet.... I don't think
#CST_HMBinitEsstudy_metadata <- read.csv("CST_metadata_HMBinitEsstudy.csv", header=TRUE, sep = ",", row.names="sample_name")
#This is to select the top taxa and then I rearrange/ consolidate to "Other", but I have already done this
#   topN = 15
#   HMBfinalEstrousT <- t(HMBfinalEstrous)
#   HMBfinalEstrous_taxa <- otu_table(HMBfinalEstrousT, taxa_are_rows = TRUE)
#   most_abundant_taxaHMBinitEsstudy = sort(taxa_sums(HMBfinalEstrous_taxa), TRUE)[1:topN]
#   print(most_abundant_taxaHMBinitEsstudy)
#   write.csv(most_abundant_taxaHMBinitEsstudy,"top15HMBinitEsst_taxa.csv")
#   
#   most_abundant_taxa20HMBinitEsst = sort(taxa_sums(HMBfinalEstrous_taxa), TRUE)[1:20]
#   print(most_abundant_taxa20HMBinitEsst)
#   write.csv(most_abundant_taxa20HMBinitEsst,"top20GHMBinitEsst_taxa.csv")
#   
#   OTU_HMBinitEsre_ordername <- read.csv("exported/RelFreq_1900_HMBinitEsstudy/RelFreq_top24taxahead_HMBinitEsstudy.csv", header=TRUE, sep = ",", row.names="sample_name")

#Use this to see what indices I need to remove or combine into "other"

%%Rename the CST values accordingly%%
  # CST I: Staphylococcus
  # CST II: Staph-Enterococcus
  # CST III: Enterococcus
  #CST IV: Lactobacillus
  #CST V: High alpha diversity, heterogenous bacteria
  #CST VI: HMBinitEs
  
  %%HEATMAP%%
  
install.packages("RColorBrewer")
library("RColorBrewer")
colMainHMBfinalEs <- colorRampPalette(brewer.pal(9, "BuPu"))(25)

OTU_clustHMBfinalEs_matrix <- as.matrix(HMBfinalEstrous)

#for CST
my_group_HMBfinalEs <- as.numeric(as.factor(substr(meta_HMBfinalEstrousstudy$CST,1,1)))
#mycolHMBfinalEs <- c("tomato","yellow","orange","lightseagreen","royalblue")
mycolHMBfinalEs <- c("tomato","lightseagreen","yellow","orange","royalblue")
colSideHMBfinalEs <- mycolHMBfinalEs[my_group_HMBfinalEs]
heatmap(OTU_clustHMBfinalEs_matrix, cexCol = 0.8, scale= "none", RowSideColors=colSideHMBfinalEs, col=colMainHMBfinalEs, Colv= NA, distfun = function(x) dist(x, method="euclidean"), hclustfun = function(x) hclust(x, method="ward.D"))
legend(x="right",legend)
# 
# my_group_6HMBfinalEs <- as.numeric(as.factor(substr(meta_6HMBfinalEstrousstudy$CST,1,1)))
# mycol6HMBfinalEs <- c("tomato","yellow","orange","lightseagreen","royalblue", "pink")
# colSide6HMBfinalEs <- mycol6HMBfinalEs[my_group_6HMBfinalEs]
# heatmap(OTU_clustHMBfinalEs_matrix, cexCol = 0.8, scale= "none", RowSideColors=colSide6HMBfinalEs, col=colMainHMBfinalEs, Colv= NA, distfun = function(x) dist(x, method="euclidean"), hclustfun = function(x) hclust(x, method="ward.D"))

#need files to be only subject ID and parameter, so 2 columns total. But check if (1,1) means something about array column length and variable column

#for estrous stage
my_group_HMBfinalEstrx <- as.numeric(as.factor(substr(meta_HMBfinalEstrousstages2$estrous_stage,1,1)))
esstcol <- c("#FFFFFF","#7c4c5c","#969190","#e1ddce","#c4b444") #pink, gray, cream, yellow
esstcol_black <- c("#000000","#7c4c5c","#969190","#e1ddce","#c4b444") #pink, gray, cream, yellow
esstcol_other <- c("#FFFFFF","#E08C77","#CB7863","#4C839C","#D6CEA9","#C6A61F") #pink, gray, cream, yellow
esstcol_blues <- c("#FFFFFF", "#6f8da6", "#E8E1A6", "#254f82", "#c1d8d1") # (#FFFFFF", "#c1d8d1", "#6f8da6", "#254f82", "#d4caa9" = white, dark blue, light blue, seafoam, yellow
#("#2a4d28","#abc289","#1b260d","#42805c")
#("#b9aa22","#add279","#1c321a","#19865c")
#yellow, light mint, green, teal
#evergreen, faint mint, green #233111, teal
#
#esstcol2 <- c("#c4b444","#969190","#e1ddce","#7c4c5c") #yellow, gray, cream, pink

colSideHMBfinalEstrx <- esstcol[my_group_HMBfinalEstrx]
colSideHMBfinalEstrxBlack <- esstcol_black[my_group_HMBfinalEstrx]
colSideHMBfinalEstrxother <- esstcol_other[my_group_HMBfinalEstrx]
colSideHMBfinalEstrxblues <- esstcol_blues[my_group_HMBfinalEstrx]

heatmap(OTU_clustHMBfinalEs_matrix, cexCol = 0.8, scale= "none", RowSideColors=colSideHMBfinalEstrxblues, col=colMainHMBfinalEs,Colv= NA, distfun = function(x) dist(x, method="euclidean"), hclustfun = function(x) hclust(x, method="ward.D"))

#color palette from HMO treatments
#"#ecebe2","#9fc3cc","#466566","#09496c"

my_group_HMBfinalEstrx
my_group_HMBfinalEs
colSideHMBfinalEs
colSideHMBfinalEstrx
write.csv(my_group_HMBfinalEstrx,"new_my_group_HMBfinalEstrx_USEME.csv")
write.csv(my_group_HMBfinalEs,"new_my_group_HMBfinalEs_USEME.csv")
write.csv(colSideHMBfinalEs,"new_colSideHMBfinalEs_USEME.csv")
write.csv(colSideHMBfinalEstrxblues,"new_colSideHMBfinalEstrxblues_USEME.csv")



#saving images - example
# tiff(file="Heatmap_HMBfinalEs_CST_bupu.tif")
# heatmap(OTU_clustHMBfinalEs_matrix, cexCol = 0.8, scale= "none", RowSideColors=colSideHMBfinalEs, col=colMainHMBfinalEs, Colv= NA, distfun = function(x) dist(x, method="euclidean"), hclustfun = function(x) hclust(x, method="ward.D"))
# dev.off()


#image legend scale for body of heatmap
image(1:nrow(OTU_clustHMBfinalEs_matrix), 1, as.matrix(1:nrow(OTU_clustHMBfinalEs_matrix)), 
      col=colorRampPalette(brewer.pal(9, "BuPu"))(25),
      xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
