HMBpaper_fullanalysis_FINAL 3
#Combined Study16S (Jackson and UCSD from UCSD) and BCM (Jackson, CCM, and HMb)

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

#deblur for each UCSD run - took previously generated demux-full_1-joined-filter-stats.qza from old folder into UCSDseqs
qiime deblur denoise-16S \
--i-demultiplexed-seqs demux-full_1-joined-filtered.qza \
--p-left-trim-len 12 \
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
--m-input-file demux-full_1-joined-filter-stats.qza \
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
--i-demultiplexed-seqs demux-full_2-joined-filtered.qza \
--p-left-trim-len 12 \
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
--m-input-file demux-full_2-joined-filter-stats.qza \
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
--o-visualization demux-full_3.qzv#Saved Visualization to: demux-full_3.qzv
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
--i-demultiplexed-seqs demux-full_3-joined-filtered.qza \
--p-left-trim-len 12 \
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
--m-input-file demux-full_3-joined-filter-stats.qza \
--o-visualization demux-filter-stats_3.qzv
  #Saved Visualization to: demux-filter-stats_3.qzv
qiime deblur visualize-stats \
--i-deblur-stats deblur-stats_3.qza \
--o-visualization deblur-stats_3.qzv
  #Saved Visualization to: deblur-stats_3.qzv


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### HMB Dataset
# All of the below commands assume you are in the following directory
source activate qiime2-2022.8

cd ~/Desktop/Patras_Lab/Humanized_Models/HMB_BrittonMice/RawSequences/Patras_x
#prep the manifest with sample_id/absolute-filepath/direction

cd /Users/marlydmejia/Desktop/Patras_Lab/BCM_16Sseq_and_metadata/RawSequencesZip
cd /Users/marlydmejia/Desktop/Patras\ Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/OneDrive_1_7-14-2023 


#assuming all the folders have already been opened/unzipped from prior individual analyses
# convert to .gz
#repeat for all folders that need to be prepped
for f in `ls RawSequences_exp1`; do
bzcat RawSequences_exp1/$f/*.1.fq.bz2 | gzip -c > RawSequences_exp1/$f/$f.1.fq.gz;
bzcat RawSequences_exp1/$f/*.2.fq.bz2 | gzip -c > RawSequences_exp1/$f/$f.2.fq.gz;
done
for f in `ls RawSequences_exp2`; do
bzcat RawSequences_exp2/$f/*.1.fq.bz2 | gzip -c > RawSequences_exp2/$f/$f.1.fq.gz;
bzcat RawSequences_exp2/$f/*.2.fq.bz2 | gzip -c > RawSequences_exp2/$f/$f.2.fq.gz;
done
for f in `ls RawSequences_exp3`; do
bzcat RawSequences_exp3/$f/*.1.fq.bz2 | gzip -c > RawSequences_exp3/$f/$f.1.fq.gz;
bzcat RawSequences_exp3/$f/*.2.fq.bz2 | gzip -c > RawSequences_exp3/$f/$f.2.fq.gz;
done

#HUVAMI_1 samples are going to be pulled from their processed filepath in  Humanized_Models/HMB_BrittonMice
#skip converting from .bz2 to gz since already done - jump to manifest - included in manifest where the extension will be changed in the following block of code

# Remove the bz2 compressed files, to reduce duplication
# this command will give you a list of all the files that will be found with the "find" command. These are the files that will then be pipped into the remove command
find RawSequences_exp2 -name "*.bz2" | less
# Use Ctrl+Z or Ctrl+C to exit list back to command line
#Actual remove command (only do if you are sure the list being pipped has the correct files). "xargs" is a tool that can read a list of e.g. filenames and do the same command on them at once - kind of like a loop. https://ss64.com/osx/xargs.html
find RawSequences_exp2 -name "*.bz2" | xargs rm


# update your manifest file (“sed” is the commandline search and replace)
cd /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete
cat HMBpaper_finalmanifest_V1_updated.txt | sed "s/bz2/gz/g" > HMBpaper_finalmanifest_V2.txt

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path HMBpaper_finalmanifest_V2.txt \
--output-path paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2
#Imported HMBpaper_finalmanifest_V2.txt as PairedEndFastqManifestPhred33V2 to paired-end-demux.qza

qiime tools peek paired-end-demux.qza
# UUID:        3900adc7-2765-4b40-8956-57371588b1e1
# Type:        SampleData[PairedEndSequencesWithQuality]
# Data format: SingleLanePerSamplePairedEndFastqDirFmt

#second run since first showed 0 reads for ccxxxxxx and seqs were re-downloaded
# UUID:        fb2a4604-0ee8-41ac-912b-a4a58dcbcfd2
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
#Saved SampleData[JoinedSequencesWithQuality] to: qiime2-read-joining-tutorial/demux-joined-filtered.qza
#Saved QualityFilterStats to: qiime2-read-joining-tutorial/demux-joined-filter-stats.qza

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
--o-representative-sequences rep-seqs_HMB.qza \
--o-table table_HMB.qza \
--o-stats deblur-stats_HMB.qza
# Saved FeatureTable[Frequency] to: table_HMB.qza
# Saved FeatureData[Sequence] to: rep-seqs_HMB.qza
# Saved DeblurStats to: deblur-stats_HMB.qza

#Visualize demux joined filter stats and deblur
qiime metadata tabulate \
--m-input-file qiime2-read-joining-tutorial/demux-joined-filter-stats.qza \
--o-visualization Visualization/demux-joined-filter-stats.qzv
# Saved Visualization to: Visualization/demux-joined-filter-stats.qzv

qiime deblur visualize-stats \
--i-deblur-stats deblur-stats_HMB.qza \
--o-visualization Visualization/deblur-stats.qzv
# Saved Visualization to: Visualization/deblur-stats.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqs_HMB.qza \
--o-visualization Visualization/rep-seqs_HMB.qzv
# Saved Visualization to: Visualization/rep-seqs_HMB_final.qzv

##################### MERGE SEQUENCES ##run 06/22/22
cd /Users/marlydmejia/Desktop/Patras_Lab #changed so I can reach into older files not in my paper directory -aka from HMO trials
mkdir Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported
mkdir Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output
mkdir Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization

#merging all three studies
qiime feature-table merge \
--i-tables Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/UCSDseqs/table_1.qza \
--i-tables Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/UCSDseqs/table_2.qza \
--i-tables Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/UCSDseqs/table_3.qza \
--i-tables Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/table_HMB.qza \
--o-merged-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_feattable

#merging rep-seqs
qiime feature-table merge-seqs \
--i-data Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/UCSDseqs/rep-seqs_1.qza \
--i-data Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/UCSDseqs/rep-seqs_2.qza \
--i-data Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/UCSDseqs/rep-seqs_3.qza \
--i-data Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/rep-seqs_HMB.qza \
--o-merged-data Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_seqs.qza
#Saved FeatureData[Sequence] to: Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_seqs.qza

qiime tools export \
--input-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_seqs.qza \
--output-path Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/merged_rep-seqs
#Exported Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy_seqs.qza as DNASequencesDirectoryFormat to directory Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/merged_rep-seqs

#### integrating taxonomy and phylogeny

#using new Greengenes2 reference. Install plugin. re-cache qiime2
#cite McDonald et al bioRxiv 2022
pip install q2-greengenes2
conda activate qiime2-2022.8
wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza
cd Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete

#non-V4 because I have paired end reads vs forward read fragments alone that match the V4 pipeline
qiime greengenes2 non-v4-16s \
--i-table mergedstudy_feattable.qza \
--i-sequences mergedstudy_seqs.qza \
--i-backbone 2022.10.backbone.full-length.fna.qza \
--o-mapped-table mergedstudy_feattable_gg2.qza \
--o-representatives mergedstudy_seqs_gg2.qza
#Saved FeatureTable[Frequency] to: mergedstudy_feattable_gg2.qza
#Saved FeatureData[Sequence] to: mergedstudy_seqs_gg2.qza

qiime greengenes2 taxonomy-from-table \
--i-reference-taxonomy 2022.10.taxonomy.md5.nwk.qza \
--i-table mergedstudy_feattable_gg2.qza \
--o-classification mergedstudy.gg2.tabletaxonomymd5.qza
#Saved FeatureData[Taxonomy] to: mergedstudy.gg2.taxonomy.qza

wget http://ftp.microbio.me/greengenes_release/current/2022.10.phylogeny.asv.nwk.qza




##EXPORTING##
qiime tools export \
--input-path mergedstudy.gg2.tabletaxonomy.qza \
--output-path exported/mergedstudy_tabletaxonomy

qiime tools export \
--input-path mergedstudy_feattable.qza \
--output-path exported/merged_feattable

qiime tools export \
--input-path mergedstudy_feattable_gg2.qza \
--output-path exported/merged_feattable_gg2

#OVERALL peek into table to see what taxa are in blanks and need to be removed
biom convert -i exported/merged_feattable/feature-table.biom -o exported/merged_feattable/OTU_table.txt --to-tsv
biom convert -i exported/merged_feattable_gg2/feature-table.biom -o exported/merged_feattable_gg2/OTU_table_gg2.txt --to-tsv --header-key taxonomy

biom convert -i exported/merged_feattable_gg2/OTU_table_gg2_taxatable.txt -o exported/merged_feattable_gg2/OTU_mergedfeattable_gg2TAXA_json.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy

qiime tools export \
--input-path mergedstudy_seqs.qza \
--output-path exported/merged_rep-seqs

qiime feature-table tabulate-seqs \
--i-data mergedstudy_seqs.qza \
--o-visualization Visualization/mergedstudy_seqs.qzv

qiime tools export \
--input-path mergedstudy_seqs_gg2.qza \
--output-path exported/merged_rep-seqs_gg2

qiime feature-table tabulate-seqs \
--i-data mergedstudy_seqs_gg2.qza \
--o-visualization Visualization/mergedstudy_seqs_gg2.qzv

qiime metadata tabulate \
--m-input-file mergedstudy.gg2.tabletaxonomy.qza \
--o-visualization Visualization/mergedstudy_taxonomy

##         ##
#before proceeding with usual code, manually moved non gg2 filtered tables, seqs, and taxonomy to folder pre_gg2-filter
#added seqs to metadata files and samples-to-keep where needed
qiime feature-table filter-samples \
--i-table mergedstudy_feattable_gg2.qza \
--m-metadata-file samples-to-keep_mergedstudy.txt \
--o-filtered-table subset_mergedstudy_feattable.qza

#back to making working feature table
qiime feature-table summarize \
--i-table subset_mergedstudy_feattable.qza \
--o-visualization Visualization/mergedstudy_Final_feattable.qzv \
--m-sample-metadata-file HMBpaper_metadata_FINAL_addseqs.txt

#peek into new SUBSET table to see what OTUs are present
qiime tools export \
--input-path subset_mergedstudy_feattable.qza \
--output-path exported/SUBSET_feattable
biom convert -i exported/SUBSET_feattable/feature-table.biom -o exported/SUBSET_feattable/feat_table.txt --to-tsv

#samples to keep for figures will come later, post filter

########## FILTER by min samples ##########
mkdir exported/feattable_min5samples

#filter table for feat IDs that appear in 5+ samples (based on orig gg2 feattable). Then will proceed to add taxonomic names
qiime feature-table filter-features \
--i-table subset_mergedstudy_feattable.qza \
--p-min-samples 5 \
--o-filtered-table subset_mergedstudy_feattable_min5samples.qza

#export reduced table based on features in 5+ samples
qiime tools export \
--input-path subset_mergedstudy_feattable_min5samples.qza \
--output-path exported/feattable_min5samples

biom convert -i exported/feattable_min5samples/feature-table.biom -o exported/feattable_min5samples/feat_table.txt --to-tsv --header-key taxonomy
#added taxonomy in txt and am converting back
#biom convert -i exported/feattable_min5samples/feat_table_taxa.txt -o exported/feattable_min5samples/feat_json.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy

biom convert -i exported/feattable_min5samples/feat_taxatable.txt -o exported/feattable_min5samples/feat_taxatable_json.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy


#------------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%% R studio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
######### FILTER OUT CONTAMINANTS USING DECONTAM #########
  
## ----loadPS----------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")
"1.44.0"

library(ggplot2); packageVersion("ggplot2")
"3.4.2"

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")
library(decontam); packageVersion("decontam")
"1.20.0"

setwd("~")
#upload the matrix into phyloseq format
#transpose feat table and change first col title to "sample_name"
#Go in and add prefix "ID_" to all feature ID names in the taxonomy file and in the OTU tables and such. This way they can be matched later for taxa removal or merging
#feature table read with feature IDs and not taxonomy names
#phyloseq converts "-" to periods, so need to replace characters in original Taxonomy_IDs file prior to import
#make sure "taxonomy" column (now last row) is deleted
#add col in metadata "Sample_or_Control" and indicate Control or True Sample
tFeatureTable <- read.csv("Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/feattable_min5samples/feat_table_matrix.csv", header=TRUE, sep = ",", row.names="sample_name")
metadata_BC <- read.csv("Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.csv", header=TRUE, sep = ",", row.names="sample_name")
taxonomy <- read.csv("Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/mergedstudy_tabletaxonomy/taxonomy_ID.csv", header=TRUE, sep = ",", row.names="Feature_ID")
tTaxonomy <- t(taxonomy)

#create directory "Decontam_process"
dir.create("Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process")
setwd("Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process")
dir.create("decontam_items_final")
setwd("~/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final")
OTUmatrixfilt = phyloseq(otu_table(tFeatureTable,taxa_are_rows = FALSE), sample_data(metadata_BC))
OTUmatrixdecontam = phyloseq(otu_table(tFeatureTable,taxa_are_rows = FALSE))


## ----see-meta-table--------------------------------------------------------
head(sample_data(OTUmatrixfilt))
quartz()

dffilt <- as.data.frame(sample_data(OTUmatrixfilt)) # Put sample_data into a ggplot-friendly data.frame
dffilt$LibrarySize <- sample_sums(OTUmatrixfilt)
dffilt <- dffilt[order(dffilt$LibrarySize),]
dffilt$Index <- seq(nrow(dffilt))
ggplot(data=dffilt, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()


sample_data(OTUmatrixfilt)$is.neg <- sample_data(OTUmatrixfilt)$Sample_or_Control == "Control"
# Make phyloseq object of presence-absence in negative controls and true samples
OTUmatrixfilt.pa <- transform_sample_counts(OTUmatrixfilt, function(abund) 1*(abund>0))
OTUmatrixfilt.pa.neg <- prune_samples(sample_data(OTUmatrixfilt.pa)$Sample_or_Control == "Control", OTUmatrixfilt.pa)
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
contamdffilt.prev01b <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.1, batch= "submission_date", batch.combine = "minimum")
table(contamdffilt.prev01b$contaminant)
dffilt01b.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                contaminant=contamdffilt.prev01b$contaminant)
ggplot(data=dffilt01b.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.1b)") + ylab("Prevalence (True Samples)")

## ----see-prev-.2-----------------------------------------------------------
contamdffilt.prev02 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.2, batch= "submission_date", batch.combine = "minimum")
table(contamdffilt.prev02$contaminant)
dffilt02.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                contaminant=contamdffilt.prev02$contaminant)
ggplot(data=dffilt02.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.2)") + ylab("Prevalence (True Samples)")

## ----see-prev-.3-----------------------------------------------------------
contamdffilt.prev03 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.3, batch= "submission_date", batch.combine = "minimum")
table(contamdffilt.prev03$contaminant)
dffilt03.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                contaminant=contamdffilt.prev03$contaminant)
ggplot(data=dffilt03.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.3)") + ylab("Prevalence (True Samples)")

## ----see-prev-.4-----------------------------------------------------------
contamdffilt.prev04 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.4, batch= "submission_date", batch.combine = "minimum")
table(contamdffilt.prev04$contaminant)
dffilt04.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                contaminant=contamdffilt.prev04$contaminant)
ggplot(data=dffilt04.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.4)") + ylab("Prevalence (True Samples)")

## ----see-prev-.45-----------------------------------------------------------
contamdffilt.prev04_5 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.45, batch= "submission_date", batch.combine = "minimum")
table(contamdffilt.prev04_5$contaminant)
dffilt04_5.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                  contaminant=contamdffilt.prev04_5$contaminant)
ggplot(data=dffilt04_5.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.45)") + ylab("Prevalence (True Samples)")

## ----prevalence-.5---------------------------------------------------------
contamdffilt.prev05 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.5, batch= "submission_date", batch.combine = "minimum")
table(contamdffilt.prev05$contaminant)
dffilt05.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                                  contaminant=contamdffilt.prev05$contaminant)
ggplot(data=dffilt05.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.5)") + ylab("Prevalence (True Samples)")


## ----prevalence-.05---------------------------------------------------------
contamdffilt.prev005 <- isContaminant(OTUmatrixfilt, method="prevalence", neg="is.neg", threshold=0.05, batch= "submission_date", batch.combine = "minimum")
table(contamdffilt.prev005$contaminant)
dffilt005.pa <- data.frame(pa.pos=taxa_sums(OTUmatrixfilt.pa.pos), pa.neg=taxa_sums(OTUmatrixfilt.pa.neg),
                          contaminant=contamdffilt.prev005$contaminant)
ggplot(data=dffilt005.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls threshold 0.05)") + ylab("Prevalence (True Samples)")


write.csv(contamdffilt.prev01b,"contamdffilt.prev01b.csv", row.names = TRUE)
write.csv(contamdffilt.prev,"contamdffilt.prev01.csv", row.names = TRUE)
write.csv(contamdffilt.prev02,"contamdffilt.prev02.csv", row.names = TRUE)
write.csv(contamdffilt.prev03,"contamdffilt.prev03.csv", row.names = TRUE)
write.csv(contamdffilt.prev04,"contamdffilt.prev04.csv", row.names = TRUE)
write.csv(contamdffilt.prev04_5,"contamdffilt.prev04_5.csv", row.names = TRUE)
write.csv(contamdffilt.prev05,"contamdffilt.prev05.csv", row.names = TRUE)
write.csv(contamdffilt.prev005,"contamdffilt.prev005.csv", row.names = TRUE)

write.csv(dffilt01b.pa,"dffilt01b.pa.csv", row.names = TRUE)
write.csv(dffilt.pa,"dffilt01.pa.csv", row.names = TRUE)
write.csv(dffilt02.pa,"dffilt02.pa.csv", row.names = TRUE)
write.csv(dffilt03.pa,"dffilt03.pa.csv", row.names = TRUE)
write.csv(dffilt04.pa,"dffilt04.pa.csv", row.names = TRUE)
write.csv(dffilt04_5.pa,"dffilt04_5.pa.csv", row.names = TRUE)
write.csv(dffilt05.pa,"dffilt05.pa.csv", row.names = TRUE)
write.csv(dffilt005.pa,"dffilt005.pa.csv", row.names = TRUE)



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

taxa_dffilt005 <- merge(dffilt005.pa, taxonomy, by= 'row.names')
colnames(taxa_dffilt005)[1] <- "Feature_ID"
taxa_contam005 <- merge(contamdffilt.prev005, taxonomy, by= 'row.names')
colnames(taxa_contam005)[1] <- "Feature_ID"
taxa_contam_dffilt005 <- merge(taxa_dffilt005, taxa_contam005, by= 'Feature_ID')
write.csv(taxa_contam_dffilt005,"taxa_contam_dffilt005.csv", row.names = FALSE)

# remove_decontam <- read.csv("taxa_contam_dffilt_remove.csv", header=FALSE, sep = ",")
# names(remove_decontam) <- as.matrix(remove_decontam[1, ])
# remove_decontam <- remove_decontam[-1, ]
# remove_decontamTAXA <- remove_decontam[,-1]
# rownames(remove_decontamTAXA) <- remove_decontam[,1]
# 
# #Selected threshold- remove these feature IDs from the matrix and then write a transposed version so that samples are rows
# OTUmatrixfilt
# OTUmatrixfilt.decontamSel <- prune_taxa(!remove_decontamTAXA$contaminant, OTUmatrixfilt)
# OTUmatrixfilt.decontamSel

  #manually removing columns of decontam-flagged taxa minus the 5 (code removed wrong taxa. uploading other csv returned incorrect variables - all characters)
  #also, taxonomy uploaded where '-' became '.' so not understood is uploaded back to qiime2 
#------------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%% Qiime2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
########## FILTER OUT SPECIFIC TAXA ###########
cd /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete

#Convert to qza file
biom convert -i Decontam_process/decontam_items_final/feat_table_decontam.txt -o Decontam_process/decontam_items_final/feat_table_decontam.biom --table-type="OTU table" --to-json

mkdir Decontam_process/decontam_items_final/Output
    
qiime tools import \
--input-path Decontam_process/decontam_items_final/feat_table_decontam.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path Decontam_process/decontam_items_final/Output/OTU_filtdecontam.qza

#re-import new table to qiime2 for filtering of known contaminants
qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Geodermatophilaceae; g__Blastococcus; s__Blastococcus aggregatus" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Geodermatophilaceae; g__Blastococcus; s__Blastococcus litoris" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt2_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt2_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Geodermatophilaceae; g__Geodermatophilus_A; s__Geodermatophilus_A marinus" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt3_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt3_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Aquificota; c__Aquificae; o__Aquificales; f__Aquificaceae; g__Thermothrix; s__Thermothrix azorensis" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt4_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt4_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Deinococcota; c__Deinococci; o__Deinococcales; f__Deinococcaceae; g__Deinococcus_B; s__Deinococcus_B geothermalis" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt5_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt5_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Deinococcota; c__Deinococci; o__Deinococcales; f__Thermaceae_405955; g__Thermus_A; s__Thermus_A arciformis" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt6_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt6_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Firmicutes_A; c__Thermoanaerobacteria; o__Caldicellulosiruptorales; f__Caldicellulosiruptoraceae; g__Caldicellulosiruptor; s__Caldicellulosiruptor acetigenus" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt7_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt7_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Firmicutes_A; c__Thermoanaerobacteria; o__Thermoanaerobacterales_43504; f__Thermoanaerobacteraceae_43502; g__Thermoanaerobacterium; s__Thermoanaerobacterium butyriciformans" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt8_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt8_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Thermotogota; c__Thermotogae; o__Thermotogales; f__Fervidobacteriaceae; g__Fervidobacterium_A; s__Fervidobacterium_A pennivorans_14473" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt9_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt9_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "p__Cyanobacteria" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt10_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt10_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales_650611; f__Pseudomonadaceae; g__Pseudomonas_E_647464; s__Pseudomonas_E_647464 salomonii" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt11_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt11_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Deinococcota; c__Deinococci; o__Deinococcales; f__Thermaceae_405955; g__Meiothermus_B_405754; s__Meiothermus_B_405754 silvanus" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt12_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt12_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "o__Rhizobiales" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt13_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt13_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales_650611; f__Pseudomonadaceae; g__Pseudomonas_F; s__Pseudomonas_F furukawaii" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filt14_OTU_filtdecontam.qza

qiime taxa filter-table \
--i-table Decontam_process/decontam_items_final/Output/filt14_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-mode contains \
--p-exclude "d__; p__; c__; o__; f__; g__; s__" \
--o-filtered-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza

qiime feature-table relative-frequency \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--o-relative-frequency-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam_relfreq

qiime tools export \
--input-path Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--output-path Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam
biom convert -i Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam/feature-table.biom -o Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam/feature-table.txt --to-tsv
biom convert -i Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam/feature-taxatable.txt -o Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam/feature-taxatable_json.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy

qiime tools export \
--input-path Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam_relfreq.qza \
--output-path Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam_relfreq
biom convert -i Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam_relfreq/feature-table.biom -o Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam_relfreq/feature-table.txt --to-tsv

mkdir subset_GBS
mkdir subset_estrousstaging
mkdir subset_longitudinal
mkdir subset_vivarium
mkdir subset_FPVS
mkdir subset_vivariumHMB


wget http://ftp.microbio.me/greengenes_release/current/2022.10.phylogeny.asv.nwk.qza
wget http://ftp.microbio.me/greengenes_release/current/2022.10.phylogeny.md5.nwk.qza
#Saved Phylogeny[Rooted] and didn't have to generate it myself

qiime tools export \
--input-path 2022.10.phylogeny.md5.nwk.qza \
--output-path exported/rooted-tree_mergedstudy
#Exported 2022.10.phylogeny.asv.md5.qza as NewickDirectoryFormat to directory exported/rooted-tree_mergestudy

mkdir core-metrics-results
cd ../../../../..

qiime metadata tabulate \
--m-input-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/rooted-tree.qzv
# There was an issue with viewing the artifact Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza as QIIME 2 Metadata:
#   
#   Artifact <artifact: Phylogeny[Rooted] uuid: 1d6fd745-9191-448c-9066-6b754e53a272> cannot be viewed as QIIME 2 Metadata.
#alpha rarefaction is in the "16S Analysis -> old" folder since included the one extra in UCSD mice and the one HMBpt2d6m11

qiime diversity alpha-rarefaction \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--i-phylogeny Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-max-depth 1500 \
--p-metrics 'observed_features' \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-min-depth 5 \
--p-steps 100 \
--p-iterations 100 \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_otus_allseq.qzv \
 "11:19pm-2am"

#since overwrote otu.qza, saved all csv combinations. running again with less steps *STILL NEED TO DO
qiime diversity alpha-rarefaction \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--i-phylogeny Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-max-depth 1500 \
--p-metrics 'observed_features' \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-min-depth 5 \
--p-steps 20 \
--p-iterations 50 \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_otus_broad_allseq.qzv \

qiime diversity alpha-rarefaction \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--i-phylogeny Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-max-depth 100000 \
--p-metrics 'observed_features' \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-min-depth 5 \
--p-steps 20 \
--p-iterations 50 \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_otus_max_allseq.qzv \

qiime diversity alpha-rarefaction \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--i-phylogeny Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-max-depth 200 \
--p-metrics 'observed_features' \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-min-depth 5 \
--p-steps 20 \
--p-iterations 100 \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_OTU_fine_allseq.qzv \

qiime diversity alpha-rarefaction \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--i-phylogeny Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-max-depth 1500 \
--p-metrics 'shannon' \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-min-depth 5 \
--p-steps 20 \
--p-iterations 50 \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_shannon_broad_allseq.qzv \
--output-dir Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Jackknife_alpha

qiime diversity alpha-rarefaction \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--i-phylogeny Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-max-depth 100000 \
--p-metrics 'shannon' \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-min-depth 5 \
--p-steps 20 \
--p-iterations 50 \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_shannon__maxallseq.qzv \

qiime diversity alpha-rarefaction \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--i-phylogeny Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-max-depth 200 \
--p-metrics 'shannon' \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-min-depth 5 \
--p-steps 20 \
--p-iterations 100 \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_shannon_fine_allseq.qzv \

#Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.xlsx > Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.tsv
#export UNIFRAC_USE_GPU=N

cd Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete

#####create tables per study/figure type
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--m-metadata-file samples-to-keep_LONGITUDINAL_BASELINEsubset.txt \
--o-filtered-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza
  #Longitudinal and Estrous_GBSd0 are now the same content, so covered in subset above
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--m-metadata-file samples-to-keep_ESTROUSsubset.txt \
--o-filtered-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--m-metadata-file samples-to-keep_GBSsubset.txt \
--o-filtered-table subset_GBS/OTU_filtdecontamfilt_GBS.qza
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--m-metadata-file samples-to-keep_VIVARIUMsubset.txt \
--o-filtered-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza
#adjusted to keep all parity but only day 3 *added 1 and 5*, and only study 3 P11-25 and T11-25 since no GBS in run. All conventional BCM
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--m-metadata-file samples-to-keep_FPVS.txt \
--o-filtered-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza
#qiime feature-table filter-samples \
#--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
#--m-metadata-file samples-to-keep_VivariumHMB.txt \
#--o-filtered-table subset_VivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza
qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--m-metadata-file samples-to-keep_VIVARIUM_mixHMB.txt \
--o-filtered-table subset_VivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza

%%bash -e
qiime feature-table relative-frequency \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--o-relative-frequency-table subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq
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

qiime feature-table relative-frequency \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--o-relative-frequency-table subset_vivarium/OTU_filtdecontamfilt_vivarium_relfreq
#Saved FeatureTable[RelativeFrequency] to: subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq.qza

qiime feature-table relative-frequency \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--o-relative-frequency-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB_relfreq

qiime feature-table relative-frequency \
--i-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--o-relative-frequency-table subset_FPVS/OTU_filtdecontamfilt_FPVS_relfreq

#export feature tables (reads and rel_frequency)
qiime tools export \
--input-path subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--output-path subset_longitudinal/OTU_filtdecontamfilt_long_baseline
biom convert -i subset_longitudinal/OTU_filtdecontamfilt_long_baseline/feature-table.biom -o subset_longitudinal/OTU_filtdecontamfilt_long_baseline/filtfeat_table_long_baseline.txt --to-tsv

qiime tools export \
--input-path subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--output-path subset_estrousstaging/OTU_filtdecontamfilt_estrous
biom convert -i subset_estrousstaging/OTU_filtdecontamfilt_estrous/feature-table.biom -o subset_estrousstaging/OTU_filtdecontamfilt_estrous/filtfeat_table_estrousstaging.txt --to-tsv

qiime tools export \
--input-path subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--output-path subset_GBS/OTU_filtdecontamfilt_GBS
biom convert -i subset_GBS/OTU_filtdecontamfilt_GBS/feature-table.biom -o subset_GBS/OTU_filtdecontamfilt_GBS/filtfeat_table_GBS.txt --to-tsv

qiime tools export \
--input-path subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--output-path subset_vivarium/OTU_filtdecontamfilt_vivarium
biom convert -i subset_vivarium/OTU_filtdecontamfilt_vivarium/feature-table.biom -o subset_vivarium/OTU_filtdecontamfilt_vivarium/filtfeat_table_vivarium.txt --to-tsv

qiime tools export \
--input-path subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--output-path subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB
biom convert -i subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB/feature-table.biom -o subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB/filtfeat_table_vivariumHMB.txt --to-tsv

qiime tools export \
--input-path subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--output-path subset_FPVS/OTU_filtdecontamfilt_FPVS
biom convert -i subset_FPVS/OTU_filtdecontamfilt_FPVS/feature-table.biom -o subset_FPVS/OTU_filtdecontamfilt_FPVS/filtfeat_table_FPVS.txt --to-tsv



qiime tools export \
--input-path subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq.qza \
--output-path subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq
biom convert -i subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq/feature-table.biom -o subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq/filtfeat_table_longitudinal_relfreq.txt --to-tsv

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
--input-path subset_vivarium/OTU_filtdecontamfilt_vivarium_relfreq.qza \
--output-path subset_vivarium/OTU_filtdecontamfilt_vivarium_relfreq
biom convert -i subset_vivarium/OTU_filtdecontamfilt_vivarium_relfreq/feature-table.biom -o subset_vivarium/OTU_filtdecontamfilt_vivarium_relfreq/filtfeat_table_vivarium_relfreq.txt --to-tsv

qiime tools export \
--input-path subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB_relfreq.qza \
--output-path subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB_relfreq
biom convert -i subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB_relfreq/feature-table.biom -o subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB_relfreq/filtfeat_table_vivariumHMB_relfreq.txt --to-tsv

qiime tools export \
--input-path subset_FPVS/OTU_filtdecontamfilt_FPVS_relfreq.qza \
--output-path subset_FPVS/OTU_filtdecontamfilt_FPVS_relfreq
biom convert -i subset_FPVS/OTU_filtdecontamfilt_FPVS_relfreq/feature-table.biom -o subset_FPVS/OTU_filtdecontamfilt_FPVS_relfreq/filtfeat_table_FPVS_relfreq.txt --to-tsv


#tease apart Jackson sites
qiime diversity beta-group-significance \
--i-distance-matrix subset_vivariumHMB/Output/HMBfinal_vivariumHMB_BCDmatrix.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column house_siteline \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_vivariumHMB_byJAX_Beta_diversity_permanova
qiime diversity beta-group-significance \
--i-distance-matrix subset_vivariumHMB/Output/HMBfinal_vivariumHMB_BCDmatrix.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column house_siteline \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_vivariumHMB_byJAX_Beta_diversity_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_vivariumHMB.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column house_siteline \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_wnu_pair-wise/HMBfinal_vivariumHMB_byJAX_wnu_permanova
qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_vivariumHMB.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column house_siteline \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_wnu_pair-wise/HMBfinal_vivariumHMB_byJAX_wnu_permdisp


qiime diversity beta-group-significance \
--i-distance-matrix subset_vivariumHMB/Output/HMBfinal_vivariumHMB_BCDmatrix.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column house_siteline_bcmJAX \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_vivariumHMB_bybcmJAX_Beta_diversity_permanova
qiime diversity beta-group-significance \
--i-distance-matrix subset_vivariumHMB/Output/HMBfinal_vivariumHMB_BCDmatrix.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column house_siteline_bcmJAX \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_vivariumHMB_bybcmJAX_Beta_diversity_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_vivariumHMB.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column house_siteline_bcmJAX \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_wnu_pair-wise/HMBfinal_vivariumHMB_bybcmJAX_wnu_permanova
qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_vivariumHMB.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column house_siteline_bcmJAX \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_wnu_pair-wise/HMBfinal_vivariumHMB_bybcmJAX_wnu_permdisp
#house_siteline_Jax
#house_siteline_UCSDJax


#confirm UCSD subset sufficiency
qiime taxa barplot \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--o-visualization Visualization/HMBfinal_fulltable_barplot

qiime diversity alpha \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--p-metric observed_features \
--o-alpha-diversity Output/HMBfinal_fulltable_alpha_otus
#Saved SampleData[AlphaDiversity] to: subset_GBS/Output/HMBfinal_GBS_alpha_otus.qza
#no longer --p-metric observed_otus because updated qiime2

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity Output/HMBfinal_fulltable_alpha_otus.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--o-visualization Visualization/HMBfinal_full_alpha_otus_usethis

%%bash -e
qiime diversity alpha \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--p-metric shannon \
--o-alpha-diversity Output/HMBfinal_full_alpha_otus_usethis

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity Output/HMBfinal_full_alpha_otus_usethis.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--o-visualization Visualization/HMBfinal_full_alpha_shannon

qiime diversity beta \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix Output/HMBfinal_full__BCDmatrix

qiime diversity beta-phylogenetic \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-variance-adjusted FALSE \
--p-bypass-tips FALSE \
--o-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_full

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_full.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column house_siteline_split_UCSD_Jax \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_wnu_pair-wise/HMBfinal_full_splitUCSDJAX_wnu_permanova
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/Beta_diversity_permanova.qzv
#house_siteline_UCSDJax

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_full.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column house_siteline_split_UCSD_Jax \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_wnu_pair-wise/HMBfinal_full_splitUCSDJAX_wnu_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix Output/HMBfinal_full__BCDmatrix.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column house_siteline_split_UCSD_Jax \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_full_splitUCSDJAX_bc_permanova
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix Output/HMBfinal_full__BCDmatrix.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column house_siteline_split_UCSD_Jax \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_full_splitUCSDJAX_bc_permdisp

qiime diversity beta-rarefaction \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme RdGy_r \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_w_norm_unif_rarefaction_full.qzv \

qiime diversity beta-rarefaction \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--p-metric braycurtis \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme RdGy_r \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/JK_bc_rarefaction_full.qzv \

qiime feature-table filter-samples \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-where "[site_line]='Jackson Labs_Conventional'" \
--o-filtered-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/ancom_JAX_full.qza

qiime taxa collapse \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/ancom_JAX_full.qza \
--i-taxonomy /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy.gg2.tabletaxonomy.qza \
--p-level 7 \
--o-collapsed-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/ancom_JAX_full_collapsedl7table.qza

qiime composition add-pseudocount \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/ancom_JAX_full_collapsedl7table.qza \
--o-composition-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/comp-ancom_JAX_full_collapsedl7table.qza

qiime composition ancom \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/comp-ancom_JAX_full_collapsedl7table.qza \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column house_siteline_split_UCSD_Jax \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/comp-ancom_JAX_full_collapsedl7table_hsl.qzv




"
######## STOPPED CODE HERE ###############
######## CONTINUE DIVERSITY ANALYSES BELOW OR SKIP TO CST USING EXPORTED FILES ###############
######## BEGAN ESTROUS HEATMAP GENERATION in R ############
"


# other diversity analyses
#need to change file names accordingly
mkdir subset_GBS/Output
mkdir subset_GBS/Visualization
mkdir subset_GBS/exported

qiime feature-table summarize \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--m-sample-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--o-visualization subset_GBS/Visualization/Summary_HMBfinal_GBS_featuretable_filt.qzv
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/Summary_featuretableMM.qzv

%%bash -e
qiime taxa barplot \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_barplot
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/_barplot.qzv

qiime diversity alpha-rarefaction \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--p-max-depth 1000 \
--p-metrics 'observed_features' \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--p-min-depth 5 \
--p-steps 25 \
--p-iterations 20 \
--o-visualization subset_GBS/Visualization/JK_otu_max_gbs.qzv \

qiime diversity alpha-rarefaction \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--p-max-depth 1000 \
--p-metrics 'shannon' \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--p-min-depth 5 \
--p-steps 25 \
--p-iterations 20 \
--o-visualization subset_GBS/Visualization/JK_shannon_max_gbs.qzv \
#chose 300

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

%%%%beta%%%%
  
"
#MEAN RAREFIED TABLE
%%bash -e
qiime diversity beta \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix
    
qiime diversity beta-group-significance \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--m-metadata-column CST_char \
--p-method permanova \
--p-pairwise TRUE \
--p-permutations 999 \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_Beta_diversity_permanova
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/Beta_diversity_permanova.qzv
"

qiime diversity beta \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--p-metric braycurtis \
--p-pseudocount none \
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
qiime diversity pcoa-biplot \
--i-pcoa subset_GBS/Output/HMBfinal_GBS__PCoA_BCD.qza \
--i-features subset_GBS/OTU_filtdecontamfilt_GBS_relfreq.qza \
--o-biplot subset_GBS/Output/HMBfinal_GBS_pcoa_BCD_biplot 
#Saved PCoAResults % Properties('biplot') to: subset_GBS/Output/pcoa_BCD_biplot.qza

%%bash -e
qiime emperor biplot \
--i-biplot subset_GBS/Output/HMBfinal_GBS_pcoa_BCD_biplot.qza \
--m-sample-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--m-feature-metadata-file mergedstudy.gg2.tabletaxonomy.qza \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_pcoaBCD_biplot
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/pcoaBCD_biplot.qzv

%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--m-metadata-column CST_char \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization subset_GBS/Visualization/HMBfinal_GBS_Beta_diversity_permanova
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--m-metadata-column CST_char \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_GBS_Beta_diversity_permdisp

"
qiime diversity adonis \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--p-permutations 999 \
--p-formula "site_line+tp_days" \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_GBS_bc_linetp_adonis \
--verbose

qiime diversity filter-distance-matrix \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--p-where "site_line='BCM_HMB'" \
--o-filtered-distance-matrix subset_GBS/Output/HMBfinal_HMB_filtered-distance-matrix.qza


qiime diversity adonis \
--i-distance-matrix subset_GBS/Output/HMBfinal_HMB_filtered-distance-matrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset_HMb.txt \
--p-permutations 999 \
--p-formula "Uterine_clearance+tp_days" \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_GBS_HMBcleartp_adonis \
--verbose

qiime tools view Visualization/unrarefied_pair-wise/HMBfinal_GBS_Beta_diversity_adonis.qzv
"

%%%%longitudinal%%%%
  
%%bash -e
qiime diversity filter-distance-matrix \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--p-where "treatment='GBS'" \
--o-filtered-distance-matrix subset_GBS/Output/HMBfinal_GBS_filtered-distance-matrix.qza
#Saved DistanceMatrix to: filtered-distance-matrix.qza

%%bash -e
qiime longitudinal first-distances \
--i-distance-matrix subset_GBS/Output/HMBfinal_GBS_filtered-distance-matrix.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--p-state-column tp_days \
--p-individual-id-column mouse_origid \
--p-replicate-handling error \
--o-first-distances subset_GBS/Output/HMBfinal_GBStp_BC_FirstDistances 
#Saved SampleData[FirstDifferences] to: _BC_FirstDistances.qza

qiime feature-table filter-samples \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS_relfreq.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--p-where "[site_line]='BCM_HMB'" \
--o-filtered-table subset_GBS/HMB_GBSsubset_relfreqtable.qza

qiime longitudinal volatility \
--i-table subset_GBS/HMB_GBSsubset_relfreqtable.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset_noNAN.txt \
--m-metadata-file subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_shannon.qza \
--p-default-metric shannon_entropy \
--p-default-group-column CST_char \
--p-state-column tp_days \
--p-individual-id-column mouse_origid \
--o-visualization subset_GBS/Visualization/shannon_GBSvolatility.qzv

qiime longitudinal volatility \
--i-table subset_GBS/HMB_GBSsubset_relfreqtable.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset_noNAN.txt \
--m-metadata-file subset_GBS/Output/HMBfinal_GBS_HMBfinal_GBS_alpha_shannon.qza \
--p-default-metric swabGBS_CFU \
--p-default-group-column CST_char \
--p-state-column tp_days \
--p-individual-id-column mouse_origid \
--o-visualization subset_GBS/Visualization/swabCFU_GBSvolatility.qzv

qiime feature-table filter-samples \
--i-table subset_GBS/OTU_filtdecontamfilt_GBS.qza \
--m-metadata-file subset_GBS/samples-to-keep-GBSvolatility.txt \
--o-filtered-table subset_GBS/HMB_GBSvolatility_feattable.qza

qiime longitudinal feature-volatility \
--i-table subset_GBS/HMB_GBSvolatility_feattable.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset_noNAN.txt \
--p-state-column tp_days \
--p-individual-id-column mouse_origid \
--p-n-estimators 50 \
--output-dir subset_GBS/feature_volatility

qiime tools export \
--input-path subset_GBS/Output/HMBfinal_GBStp_BC_FirstDistances.qza \
--output-path subset_GBS/exported/HMBfinal_GBStp_BC_FirstDistances

qiime tools export \
--input-path subset_GBS/Output/HMBfinal_GBS_pcoa_BCD_biplot.qza \
--output-path subset_GBS/exported/HMBfinal_GBS_pcoa_BCD_biplot.qza

qiime tools export \
--input-path subset_GBS/Output/HMBfinal_GBS__BCDmatrix.qza \
--output-path subset_GBS/exported/HMBfinal_GBS__BCDmatrix.qza

qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--p-where "[study_topic]='HMB_GBSchallenge'" \
--o-filtered-table HMB_GBSuterine-table.qza

qiime feature-table filter-samples \
--i-table HMB_GBSuterine-table.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--p-where "[body_site]='vaginal tract'" \
--o-filtered-table HMB_GBSuterinee-VTtable.qza

#might need to run... add HMB before.qza
qiime feature-table filter-samples \
--i-table HMB_GBSuterine-table.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--p-where "[mouse_line]='HMB'" \
--o-filtered-table HMB_GBSuterine-VTtableHMB.qza

qiime taxa collapse \
--i-table HMB_GBSuterine-VTtableHMB.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-level 6 \
--o-collapsed-table HMB_GBSuterine-collapsedl6tableHMB.qza

qiime composition add-pseudocount \
--i-table HMB_GBSuterine-collapsedl6tableHMB.qza \
--o-composition-table comp-HMB_GBSuterine-collapsedl6tableHMB.qza

qiime composition ancom \
--i-table comp-HMB_GBSuterine-collapsedl6tableHMB.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column Uterine_clearance \
--o-visualization l6-ancom-Uterine_ColonizeClear.qzv

qiime taxa collapse \
--i-table HMB_GBSuterine-VTtableHMB.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-level 7 \
--o-collapsed-table HMB_GBSuterine-collapsedl7tableHMB.qza

qiime composition add-pseudocount \
--i-table HMB_GBSuterine-collapsedl7tableHMB.qza \
--o-composition-table comp-HMB_GBSuterine-collapsedl7tableHMB.qza

qiime composition ancom \
--i-table comp-HMB_GBSuterine-collapsedl7tableHMB.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column Uterine_clearance \
--o-visualization l7-ancom-Uterine_ColonizeClear.qzv

cd /Users/
qiime taxa collapse \
--i-table Proposals/Manuscriptste/Mine/2022_HMB_characterization/16S_analyses_compleHMB_GBSuterine-VTtableHMB.qza \
--i-taxonomy Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy.gg2.tabletaxonomy.qza \
--p-level 7 \
--o-collapsed-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMB_GBSuterine-collapsedl7tableHMB.qza

qiime composition add-pseudocount \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMB_GBSuterine-collapsedl7tableHMB.qza \
--o-composition-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/comp-HMB_GBSuterine-collapsedl7tableHMB.qza

qiime composition ancom \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/comp-HMB_GBSuterine-collapsedl7tableHMB.qza \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column Uterine_clearance \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/l7-ancom-Uterine_ColonizeClear.qzv

qiime composition ancom \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/comp-HMB_GBSuterine-collapsedl7tableHMB.qza \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/HMBpaper_metadata_GBSsubset_HMb.txt \
--m-metadata-column Uterine_clearance \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/l7-ancom-Uterine_ColonizeClear_HMb.qzv

qiime feature-table filter-samples \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMB_GBSuterine-table.qza \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMBpaper_metadata_FINAL_addseqs.txt \
--p-where "[mouse_line]='Conventional'" \
--o-filtered-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMB_GBSuterine-VTtableConv.qza

qiime taxa collapse \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/HMB_GBSuterine-VTtableConv.qza \
--i-taxonomy Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/mergedstudy.gg2.tabletaxonomy.qza \
--p-level 7 \
--o-collapsed-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/HMB_GBSuterine-collapsedl7tableconv.qza

qiime composition add-pseudocount \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/HMB_GBSuterine-collapsedl7tableconv.qza \
--o-composition-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/comp-HMB_GBSuterine-collapsedl7tableconv.qza

qiime composition ancom \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/comp-HMB_GBSuterine-collapsedl7tableconv.qza \
--m-metadata-file Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/HMBpaper_metadata_GBSsubset_conv.txt \
--m-metadata-column Uterine_clearance \
--o-visualization Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/l7-ancom-Uterine_ColonizeClear_conv.qzv





%%%%%%%%other diversity analyses - longitudinal_baseline %%%%%%%%
  #need to change file names accordingly
  
mkdir subset_longitudinal/Output
mkdir subset_longitudinal/Visualization
mkdir subset_longitudinal/exported

qiime feature-table summarize \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--m-sample-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--o-visualization subset_longitudinal/Visualization/Summary_HMBfinal_longitudinal_featuretable_filt.qzv
#Saved subset_longitudinal/Visualization to: subset_longitudinal/Visualization/Summary_featuretableMM.qzv

%%bash -e
qiime taxa barplot \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--o-visualization subset_longitudinal/Visualization/HMBfinal_longitudinal_barplot
#Saved subset_longitudinal/Visualization to: subset_longitudinal/Visualization/_barplot.qzv

qiime diversity alpha-rarefaction \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--p-max-depth 1000 \
--p-metrics 'observed_features' \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--p-min-depth 5 \
--p-steps 25 \
--p-iterations 50 \
--o-visualization subset_longitudinal/Visualization/JK_otu_max_longitudinal.qzv \

qiime diversity alpha-rarefaction \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--p-max-depth 1000 \
--p-metrics 'shannon' \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--p-min-depth 5 \
--p-steps 25 \
--p-iterations 50 \
--o-visualization subset_longitudinal/Visualization/JK_shannon_max_longitudinal.qzv \

%%%%alpha%%%%
  
  %%bash -e
qiime diversity alpha \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--p-metric observed_features \
--o-alpha-diversity subset_longitudinal/Output/HMBfinal_longitudinal_alpha_otus
#Saved SampleData[AlphaDiversity] to: subset_longitudinal/Output/HMBfinal_longitudinal_alpha_otus.qza
#no longer --p-metric observed_otus because updated qiime2

%%bash -e
qiime diversity alpha-correlation \
--i-alpha-diversity subset_longitudinal/Output/HMBfinal_longitudinal_alpha_otus.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--o-visualization subset_longitudinal/Visualization/HMBfinal_longitudinal_alpha_otus
#Saved subset_longitudinal/Visualization to: subset_longitudinal/Visualization/HMBfinal_longitudinal_alpha_otus.qzv

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_longitudinal/Output/HMBfinal_longitudinal_alpha_otus.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--o-visualization subset_longitudinal/Visualization/HMBfinal_longitudinal_alpha_otus_usethis
#Saved subset_longitudinal/Visualization to: subset_longitudinal/Visualization/HMBfinal_longitudinal_alpha_otus_usethis.qzv
qiime tools export \
--input-path subset_longitudinal/Output/HMBfinal_longitudinal_alpha_otus.qza \
--output-path subset_longitudinal/Output/HMBfinal_longitudinal_alpha_otus

%%bash -e
qiime diversity alpha \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--p-metric shannon \
--o-alpha-diversity subset_longitudinal/Output/HMBfinal_longitudinal_alpha_shannon
#Saved SampleData[AlphaDiversity] to: subset_longitudinal/Output/HMBfinal_longitudinal_alpha_shannon.qza

qiime diversity alpha-correlation \
--i-alpha-diversity subset_longitudinal/Output/HMBfinal_longitudinal_alpha_shannon.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--o-visualization subset_longitudinal/Visualization/HMBfinal_longitudinal_alpha_shannon

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_longitudinal/Output/HMBfinal_longitudinal_alpha_shannon.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--o-visualization subset_longitudinal/Visualization/HMBfinal_longitudinal_alpha_shannon
#Saved subset_longitudinal/Visualization to: subset_longitudinal/Visualization/HMBfinal_longitudinal_alpha_shannon.qzv


%%%%beta%%%%
" MEAN RAREFIED TABLE
%%bash -e
qiime diversity beta \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--p-metric braycurtis \
--p-pseudocount none \
--p-n-jobs 1 \
--o-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix
#Saved DistanceMatrix to: subset_longitudinal/Output/_BCDmatrix.qza
      
qiime diversity beta-group-significance \
--i-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-metadata-column cohort \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization subset_longitudinal/Visualization/HMBfinal_longitudinal_Beta_diversity_permanova
#Saved subset_longitudinal/Visualization to: subset_longitudinal/Visualization/Beta_diversity_permanova.qzv
      
qiime diversity beta-group-significance \
--i-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-metadata-column vivarium \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization subset_longitudinal/Visualization/HMBfinal_BCMroom_Beta_diversity_permanova
"

%%bash -e
qiime diversity beta \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix
#Saved DistanceMatrix to: subset_longitudinal/Output/_BCDmatrix.qza

%%bash -e
qiime diversity pcoa \
--i-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix.qza \
--p-number-of-dimensions 3 \
--o-pcoa subset_longitudinal/Output/HMBfinal_longitudinal_PCoA_BCD
#Saved PCoAResults to: subset_longitudinal/Output/_PCoA_BCD.qza

%%bash -e
qiime emperor plot \
--i-pcoa subset_longitudinal/Output/HMBfinal_longitudinal_PCoA_BCD.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--p-ignore-missing-samples False \
--o-visualization subset_longitudinal/Visualization/HMBfinal_longitudinal_PCoA_BCD
#Saved subset_longitudinal/Visualization to: subset_longitudinal/Visualization/_PCoA_BCD.qzv

%%bash -e
qiime feature-table relative-frequency \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--o-relative-frequency-table subset_longitudinal/Output/HMBfinal_longitudinal_relative
#Saved FeatureTable[RelativeFrequency] to: subset_longitudinal/Output/_relative.qza

%%bash -e
qiime diversity pcoa-biplot \
--i-pcoa subset_longitudinal/Output/HMBfinal_longitudinal_PCoA_BCD.qza \
--i-features subset_longitudinal/Output/HMBfinal_longitudinal_relative.qza \
--o-biplot subset_longitudinal/Output/HMBfinal_longitudinal_pcoa_BCD_biplot 
#Saved PCoAResults % Properties('biplot') to: subset_longitudinal/Output/pcoa_BCD_biplot.qza

%%bash -e
qiime emperor biplot \
--i-biplot subset_longitudinal/Output/HMBfinal_longitudinal_pcoa_BCD_biplot.qza \
--m-sample-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-feature-metadata-file mergedstudy.gg2.tabletaxonomy.qza \
--o-visualization subset_longitudinal/Visualization/HMBfinal_longitudinal_pcoaBCD_biplot
#Saved subset_longitudinal/Visualization to: subset_longitudinal/Visualization/pcoaBCD_biplot.qzv

%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-metadata-column cohort_number \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_longitudinalcombnum_Beta_diversity_permanova
#Saved subset_longitudinal/Visualization to: subset_longitudinal/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-metadata-column cohort_number \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_longitudinalcombnum_Beta_diversity_permdisp

qiime diversity beta-group-significance \
--i-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-metadata-column vivarium \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_BCMroom_Beta_diversity_permanova

qiime diversity beta-group-significance \
--i-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-metadata-column vivarium \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_BCMroom_Beta_diversity_permdisp

%%%%longitudinal%%%%  
  %%bash -e
qiime diversity filter-distance-matrix \
--i-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--p-where "mouse_line='HMB'" \
--o-filtered-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_filtered-distance-matrix.qza
#Saved DistanceMatrix to: filtered-distance-matrix.qza

%%bash -e
qiime longitudinal first-distances \
--i-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_filtered-distance-matrix.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--p-state-column tp_days \
--p-individual-id-column mouse_origid \
--p-replicate-handling error \
--o-first-distances subset_longitudinal/Output/HMBfinal_longitudinal_BC_FirstDistances 
#Saved SampleData[FirstDifferences] to: _BC_FirstDistances.qza

qiime tools export \
--input-path subset_longitudinal/Output/HMBfinal_longitudinal_BC_FirstDistances.qza \
--output-path subset_longitudinal/exported/HMBfinal_longitudinal_BC_FirstDistances

qiime tools export \
--input-path subset_longitudinal/Output/HMBfinal_longitudinal_pcoa_BCD_biplot.qza \
--output-path subset_longitudinal/exported/HMBfinal_longitudinal_pcoa_BCD_biplot.qza

qiime tools export \
--input-path subset_longitudinal/Output/HMBfinal_longitudinal_BCDmatrix.qza \
--output-path subset_longitudinal/exported/HMBfinal_longitudinal_BCDmatrix.qza


qiime diversity pcoa \
--i-distance-matrix Output/Unrarefied_wnu_matrices_longitudinal.qza \
--p-number-of-dimensions 3 \
--o-pcoa subset_longitudinal/Output/HMBfinal_longitudinal_PCoA_wnu
qiime emperor plot \
--i-pcoa subset_longitudinal/Output/HMBfinal_longitudinal_PCoA_wnu.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--p-ignore-missing-samples False \
--o-visualization subset_longitudinal/Visualization/HMBfinal_longitudinal_PCoA_wnu
%%bash -e
qiime diversity pcoa-biplot \
--i-pcoa subset_longitudinal/Output/HMBfinal_longitudinal_PCoA_wnu.qza \
--i-features subset_longitudinal/Output/HMBfinal_longitudinal_relative.qza \
--o-biplot subset_longitudinal/Output/HMBfinal_longitudinal_pcoa_wnu_biplot
%%bash -e
qiime emperor biplot \
--i-biplot subset_longitudinal/Output/HMBfinal_longitudinal_pcoa_wnu_biplot.qza \
--m-sample-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-feature-metadata-file mergedstudy.gg2.tabletaxonomy.qza \
--o-visualization subset_longitudinal/Visualization/HMBfinal_longitudinal_pcoawnu_biplot

qiime feature-table filter-samples \
--i-table Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--p-where "[body_site]='vaginal tract'" \
--o-filtered-table subset_longitudinal/HMB_longitudinal-table.qza

qiime longitudinal first-distances \
--i-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_filtered-distance-matrix.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset_noNAN.txt \
--p-state-column tp_consecswabs \
--p-individual-id-column mouse_throughstudies \
--p-replicate-handling error \
--o-first-distances subset_longitudinal/Output/HMBfinal_longitudinal_BC_FirstDistances_staging 

qiime tools export \
--input-path subset_longitudinal/Output/HMBfinal_longitudinal_BC_FirstDistances_staging.qza \
--output-path subset_longitudinal/exported/HMBfinal_longitudinal_BC_FirstDistances_staging

qiime longitudinal volatility \
--i-table subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset_noNAN.txt \
--m-metadata-file subset_longitudinal/Output/HMBfinal_longitudinal_alpha_shannon.qza \
--p-default-metric shannon_entropy \
--p-default-group-column cohort_number \
--p-state-column tp_consecswabs \
--p-individual-id-column mouse_throughstudies \
--o-visualization subset_longitudinal/Visualization/shannon_cohort_volatility.qzv

qiime longitudinal feature-volatility \
--i-table subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset_noNAN.txt \
--m-metadata-file subset_longitudinal/Output/HMBfinal_longitudinal_alpha_shannon.qza \
--m-metadata-file subset_longitudinal/Output/HMBfinal_longitudinal_alpha_otus.qza \
--p-state-column tp_consecswabs \
--p-individual-id-column mouse_throughstudies \
--p-n-estimators 100 \
--p-feature-count 20 \
--output-dir subset_longitudinal/visualization/feature_volatility_top20

qiime taxa collapse \
--i-table subset_longitudinal/HMB_longitudinal-table.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-level 6 \
--o-collapsed-table subset_longitudinal/HMB_longitudinal_room-collapsedl6tableHMB.qza

qiime composition add-pseudocount \
--i-table subset_longitudinal/HMB_longitudinal_room-collapsedl6tableHMB.qza \
--o-composition-table subset_longitudinal/comp-HMB_longitudinal_room-collapsedl6tableHMB.qza

qiime composition ancom \
--i-table subset_longitudinal/comp-HMB_longitudinal_room-collapsedl6tableHMB.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column vivarium \
--o-visualization subset_longitudinal/l6-ancom-room.qzv

qiime composition ancom \
--i-table subset_longitudinal/comp-HMB_longitudinal_room-collapsedl6tableHMB.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column tp_days \
--o-visualization subset_longitudinal/l6-ancom-days.qzv

qiime taxa collapse \
--i-table subset_longitudinal/HMB_longitudinal-table.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-level 7 \
--o-collapsed-table subset_longitudinal/HMB_longitudinal_room-collapsedl7tableHMB.qza

qiime composition add-pseudocount \
--i-table subset_longitudinal/HMB_longitudinal_room-collapsedl7tableHMB.qza \
--o-composition-table subset_longitudinal/comp-HMB_longitudinal_room-collapsedl7tableHMB.qza

qiime composition ancom \
--i-table subset_longitudinal/comp-HMB_longitudinal_room-collapsedl7tableHMB.qza \
--m-metadata-file HMBpaper_metadata_FINAL_addseqs.txt \
--m-metadata-column vivarium \
--o-visualization subset_longitudinal/l7-ancom-room.qzv

qiime composition ancom \
--i-table subset_longitudinal/comp-HMB_longitudinal_room-collapsedl7tableHMB.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset_noNAN.txt \
--m-metadata-column tp_consecswabs_cat \
--o-visualization subset_longitudinal/l7-ancom-duration.qzv

qiime composition ancom \
--i-table subset_longitudinal/comp-HMB_longitudinal_room-collapsedl7tableHMB.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset_noNAN.txt \
--m-metadata-column cohort \
--o-visualization subset_longitudinal/l7-ancom-cohort.qzv




%%%%%%%%other diversity analyses- Vivarium %%%%%%%%
  #need to change file names accordingly
  
mkdir subset_vivarium/Output
mkdir subset_vivarium/Visualization
mkdir subset_vivarium/exported

"
qiime feature-table summarize \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--m-sample-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--o-visualization subset_vivarium/Visualization/Summary_HMBfinal_vivarium_featuretable_filt.qzv
#Saved subset_vivarium/Visualization to: subset_vivarium/Visualization/Summary_featuretableMM.qzv
      
%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_vivarium/Output/HMBfinal_vivarium_BCDmatrix.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--m-metadata-column birthplace \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization subset_vivarium/Visualization/HMBfinal_vivarium_Beta_diversity_permanova
#Saved subset_vivarium/Visualization to: subset_vivarium/Visualization/Beta_diversity_permanova.qzv
"
%%bash -e
qiime taxa barplot \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--o-visualization subset_vivarium/Visualization/HMBfinal_vivarium_barplot
#Saved subset_vivarium/Visualization to: subset_vivarium/Visualization/_barplot.qzv

#run alpha rarefaction script with vivariumHMB since it encompasses all these and a few more

%%%%beta%%%%
  %%bash -e
qiime diversity beta \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_vivarium/Output/HMBfinal_vivarium_BCDmatrix
#Saved DistanceMatrix to: subset_vivarium/Output/_BCDmatrix.qza

%%bash -e
qiime diversity beta \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_vivarium/Output/HMBfinal_vivarium_BCDmatrix
#Saved DistanceMatrix to: subset_vivarium/Output/_BCDmatrix.qza

%%bash -e
qiime diversity pcoa \
--i-distance-matrix subset_vivarium/Output/HMBfinal_vivarium_BCDmatrix.qza \
--p-number-of-dimensions 3 \
--o-pcoa subset_vivarium/Output/HMBfinal_vivarium_PCoA_BCD
#Saved PCoAResults to: subset_vivarium/Output/_PCoA_BCD.qza

%%bash -e
qiime emperor plot \
--i-pcoa subset_vivarium/Output/HMBfinal_vivarium_PCoA_BCD.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--p-ignore-missing-samples False \
--o-visualization subset_vivarium/Visualization/HMBfinal_vivarium_PCoA_BCD
#Saved subset_vivarium/Visualization to: subset_vivarium/Visualization/_PCoA_BCD.qzv

%%bash -e
qiime feature-table relative-frequency \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--o-relative-frequency-table subset_vivarium/Output/HMBfinal_vivarium_relative
#Saved FeatureTable[RelativeFrequency] to: subset_vivarium/Output/_relative.qza

%%bash -e
qiime diversity pcoa-biplot \
--i-pcoa subset_vivarium/Output/HMBfinal_vivarium_PCoA_BCD.qza \
--i-features subset_vivarium/Output/HMBfinal_vivarium_relative.qza \
--o-biplot subset_vivarium/Output/HMBfinal_vivarium_pcoa_BCD_biplot 
#Saved PCoAResults % Properties('biplot') to: subset_vivarium/Output/pcoa_BCD_biplot.qza

%%bash -e
qiime emperor biplot \
--i-biplot subset_vivarium/Output/HMBfinal_vivarium_pcoa_BCD_biplot.qza \
--m-sample-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--m-feature-metadata-file mergedstudy.gg2.tabletaxonomy.qza \
--o-visualization subset_vivarium/Visualization/HMBfinal_vivarium_pcoaBCD_biplot
#Saved subset_vivarium/Visualization to: subset_vivarium/Visualization/pcoaBCD_biplot.qzv

%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_vivarium/Output/HMBfinal_vivarium_BCDmatrix.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--m-metadata-column birthplace_UCSD \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization subset_vivarium/Visualization/HMBfinal_vivarium_Beta_diversity_permanova_UCSD
#Saved subset_vivarium/Visualization to: subset_vivarium/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix subset_vivarium/Output/HMBfinal_vivarium_BCDmatrix.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--m-metadata-column birthplace_UCSD \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_vivarium_Beta_diversity_permdisp_UCSD
#Rotate through metadata-column to get combination where BCM is first, Jackson is first, and then UCSD is first to be compared to
#example outputs below where birthplace is original, starts with BCM. Then birthplace_JAX and birthplace_UCSD
#HMBfinal_vivarium_Beta_diversity_permdisp
#HMBfinal_vivarium_Beta_diversity_permdisp_JAX
#HMBfinal_vivarium_Beta_diversity_permdisp_UCSD




%%%%alpha%%%%
  
  %%bash -e
qiime diversity alpha \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--p-metric observed_features \
--o-alpha-diversity subset_vivarium/Output/HMBfinal_vivarium_alpha_otus
#Saved SampleData[AlphaDiversity] to: subset_vivarium/Output/HMBfinal_vivarium_alpha_otus.qza
#no longer --p-metric observed_otus because updated qiime2

%%bash -e
qiime diversity alpha-correlation \
--i-alpha-diversity subset_vivarium/Output/HMBfinal_vivarium_alpha_otus.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--o-visualization subset_vivarium/Visualization/HMBfinal_vivarium_alpha_otus
#Saved subset_vivarium/Visualization to: subset_vivarium/Visualization/HMBfinal_vivarium_alpha_otus.qzv

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_vivarium/Output/HMBfinal_vivarium_alpha_otus.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--o-visualization subset_vivarium/Visualization/HMBfinal_vivarium_alpha_otus_usethis
#Saved subset_vivarium/Visualization to: subset_vivarium/Visualization/HMBfinal_vivarium_alpha_otus_usethis.qzv
qiime tools export \
--input-path subset_vivarium/Output/HMBfinal_vivarium_alpha_otus.qza \
--output-path subset_vivarium/Output/HMBfinal_vivarium_alpha_otus

%%bash -e
qiime diversity alpha \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--p-metric shannon \
--o-alpha-diversity subset_vivarium/Output/HMBfinal_vivarium_alpha_shannon
#Saved SampleData[AlphaDiversity] to: subset_vivarium/Output/HMBfinal_vivarium_alpha_shannon.qza

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_vivarium/Output/HMBfinal_vivarium_alpha_shannon.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--o-visualization subset_vivarium/Visualization/HMBfinal_vivarium_alpha_shannon
#Saved subset_vivarium/Visualization to: subset_vivarium/Visualization/HMBfinal_vivarium_alpha_shannon.qzv

qiime taxa collapse \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-level 6 \
--o-collapsed-table subset_vivarium/HMB_vivarium-collapsedl6tableHMB.qza

qiime composition add-pseudocount \
--i-table subset_vivarium/HMB_vivarium-collapsedl6tableHMB.qza \
--o-composition-table subset_vivarium/comp-HMB_vivarium-collapsedl6tableHMB.qza

qiime composition ancom \
--i-table subset_vivarium/comp-HMB_vivarium-collapsedl6tableHMB.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--m-metadata-column birthplace \
--o-visualization subset_vivarium/l6-ancom-birthplace_.qzv

qiime composition ancom \
--i-table subset_vivarium/comp-HMB_vivarium-collapsedl6tableHMB.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--m-metadata-column cohort \
--o-visualization subset_vivarium/l6-ancom-bpcohort.qzv

qiime taxa collapse \
--i-table subset_vivarium/OTU_filtdecontamfilt_vivarium.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-level 7 \
--o-collapsed-table subset_vivarium/HMB_vivarium-collapsedl7tableHMB.qza

qiime composition add-pseudocount \
--i-table subset_vivarium/HMB_vivarium-collapsedl7tableHMB.qza \
--o-composition-table subset_vivarium/comp-HMB_vivarium-collapsedl7tableHMB.qza

qiime composition ancom \
--i-table subset_vivarium/comp-HMB_vivarium-collapsedl7tableHMB.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--m-metadata-column birthplace \
--o-visualization subset_vivarium/l7-ancom-birthplace.qzv





%%%%%%%%other diversity analyses- vivariumHMB %%%%%%%%
  #need to change file names accordingly
  
mkdir subset_vivariumHMB/Output
mkdir subset_vivariumHMB/Visualization
mkdir subset_vivariumHMB/exported

qiime feature-table summarize \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--m-sample-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--o-visualization subset_vivariumHMB/Visualization/Summary_HMBfinal_vivariumHMB_featuretable_filt.qzv
#Saved subset_vivariumHMB/Visualization to: subset_vivariumHMB/Visualization/Summary_featuretableMM.qzv

%%bash -e
qiime taxa barplot \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--o-visualization subset_vivariumHMB/Visualization/HMBfinal_vivariumHMB_barplot
#Saved subset_vivariumHMB/Visualization to: subset_vivariumHMB/Visualization/_barplot.qzv

qiime diversity alpha-rarefaction \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--p-max-depth 1000 \
--p-metrics 'observed_features' \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--p-min-depth 5 \
--p-steps 25 \
--p-iterations 50 \
--o-visualization subset_vivariumHMB/Visualization/JK_otu_max_vivHMB.qzv

qiime diversity alpha-rarefaction \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--p-max-depth 1000 \
--p-metrics 'shannon' \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--p-min-depth 5 \
--p-steps 25 \
--p-iterations 50 \
--o-visualization subset_vivariumHMB/Visualization/JK_shannon_max_vivHMb.qzv

%%%%alpha%%%%
  
  %%bash -e
qiime diversity alpha \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--p-metric observed_features \
--o-alpha-diversity subset_vivariumHMB/Output/HMBfinal_vivariumHMB_alpha_otus
#Saved SampleData[AlphaDiversity] to: subset_vivariumHMB/Output/HMBfinal_vivariumHMB_alpha_otus.qza
#no longer --p-metric observed_otus because updated qiime2

%%bash -e
qiime diversity alpha-correlation \
--i-alpha-diversity subset_vivariumHMB/Output/HMBfinal_vivariumHMB_alpha_otus.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--o-visualization subset_vivariumHMB/Visualization/HMBfinal_vivariumHMB_alpha_otus
#Saved subset_vivariumHMB/Visualization to: subset_vivariumHMB/Visualization/HMBfinal_vivariumHMB_alpha_otus.qzv

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_vivariumHMB/Output/HMBfinal_vivariumHMB_alpha_otus.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--o-visualization subset_vivariumHMB/Visualization/HMBfinal_vivariumHMB_alpha_otus_usethis
#Saved subset_vivariumHMB/Visualization to: subset_vivariumHMB/Visualization/HMBfinal_vivariumHMB_alpha_otus_usethis.qzv
qiime tools export \
--input-path subset_vivariumHMB/Output/HMBfinal_vivariumHMB_alpha_otus.qza \
--output-path subset_vivariumHMB/Output/HMBfinal_vivariumHMB_alpha_otus

%%bash -e
qiime diversity alpha \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--p-metric shannon \
--o-alpha-diversity subset_vivariumHMB/Output/HMBfinal_vivariumHMB_alpha_shannon
#Saved SampleData[AlphaDiversity] to: subset_vivariumHMB/Output/HMBfinal_vivariumHMB_alpha_shannon.qza

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_vivariumHMB/Output/HMBfinal_vivariumHMB_alpha_shannon.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--o-visualization subset_vivariumHMB/Visualization/HMBfinal_vivariumHMB_alpha_shannon
#Saved subset_vivariumHMB/Visualization to: subset_vivariumHMB/Visualization/HMBfinal_vivariumHMB_alpha_shannon.qzv

%%%%beta%%%%
"
%%bash -e
qiime diversity beta \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_vivariumHMB/Output/HMBfinal_vivariumHMB_BCDmatrix
#Saved DistanceMatrix to: subset_vivariumHMB/Output/_BCDmatrix.qza
      
%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_vivariumHMB/Output/HMBfinal_vivariumHMB_BCDmatrix.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column site_line \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization subset_vivariumHMB/Visualization/HMBfinal_vivariumHMB_Beta_diversity_permanova
#Saved subset_vivariumHMB/Visualization to: subset_vivariumHMB/Visualization/Beta_diversity_permanova.qzv
"

%%bash -e
qiime diversity beta \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_vivariumHMB/Output/HMBfinal_vivariumHMB_BCDmatrix
#Saved DistanceMatrix to: subset_vivariumHMB/Output/_BCDmatrix.qza

%%bash -e
qiime diversity pcoa \
--i-distance-matrix subset_vivariumHMB/Output/HMBfinal_vivariumHMB_BCDmatrix.qza \
--p-number-of-dimensions 3 \
--o-pcoa subset_vivariumHMB/Output/HMBfinal_vivariumHMB_PCoA_BCD
#Saved PCoAResults to: subset_vivariumHMB/Output/_PCoA_BCD.qza

%%bash -e
qiime emperor plot \
--i-pcoa subset_vivariumHMB/Output/HMBfinal_vivariumHMB_PCoA_BCD.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--p-ignore-missing-samples False \
--o-visualization subset_vivariumHMB/Visualization/HMBfinal_vivariumHMB_PCoA_BCD
#Saved subset_vivariumHMB/Visualization to: subset_vivariumHMB/Visualization/_PCoA_BCD.qzv

%%bash -e
qiime feature-table relative-frequency \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--o-relative-frequency-table subset_vivariumHMB/Output/HMBfinal_vivariumHMB_relative
#Saved FeatureTable[RelativeFrequency] to: subset_vivariumHMB/Output/_relative.qza

%%bash -e
qiime diversity pcoa-biplot \
--i-pcoa subset_vivariumHMB/Output/HMBfinal_vivariumHMB_PCoA_BCD.qza \
--i-features subset_vivariumHMB/Output/HMBfinal_vivariumHMB_relative.qza \
--o-biplot subset_vivariumHMB/Output/HMBfinal_vivariumHMB_pcoa_BCD_biplot 
#Saved PCoAResults % Properties('biplot') to: subset_vivariumHMB/Output/pcoa_BCD_biplot.qza

%%bash -e
qiime emperor biplot \
--i-biplot subset_vivariumHMB/Output/HMBfinal_vivariumHMB_pcoa_BCD_biplot.qza \
--m-sample-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-feature-metadata-file mergedstudy.gg2.tabletaxonomy.qza \
--o-visualization subset_vivariumHMB/Visualization/HMBfinal_vivariumHMB_pcoaBCD_biplot
#Saved subset_vivariumHMB/Visualization to: subset_vivariumHMB/Visualization/pcoaBCD_biplot.qzv

%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_vivariumHMB/Output/HMBfinal_vivariumHMB_BCDmatrix.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column site_line_HMB \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_vivariumHMB_Beta_diversity_permanova_HMB
#Saved subset_vivariumHMB/Visualization to: subset_vivariumHMB/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix subset_vivariumHMB/Output/HMBfinal_vivariumHMB_BCDmatrix.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column site_line_HMB \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_vivariumHMB_Beta_diversity_permdisp_HMB

qiime taxa collapse \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-level 6 \
--o-collapsed-table subset_vivariumHMB/HMB_vivariumHMB-collapsedl6tableHMB.qza

qiime composition add-pseudocount \
--i-table subset_vivariumHMB/HMB_vivariumHMB-collapsedl6tableHMB.qza \
--o-composition-table subset_vivariumHMB/comp-HMB_vivariumHMB-collapsedl6tableHMB.qza

qiime composition ancom \
--i-table subset_vivariumHMB/comp-HMB_vivariumHMB-collapsedl6tableHMB.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column site_line \
--o-visualization subset_vivariumHMB/l6-ancom-birthplace_.qzv

qiime taxa collapse \
--i-table subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-level 7 \
--o-collapsed-table subset_vivariumHMB/HMB_vivariumHMB-collapsedl7tableHMB.qza

qiime composition add-pseudocount \
--i-table subset_vivariumHMB/HMB_vivariumHMB-collapsedl7tableHMB.qza \
--o-composition-table subset_vivariumHMB/comp-HMB_vivariumHMB-collapsedl7tableHMB.qza

qiime composition ancom \
--i-table subset_vivariumHMB/comp-HMB_vivariumHMB-collapsedl7tableHMB.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column site_line \
--o-visualization subset_vivariumHMB/l7-ancom-birthplace.qzv

qiime composition ancom \
--i-table subset_vivariumHMB/comp-HMB_vivariumHMB-collapsedl7tableHMB.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column house_siteline \
--o-visualization subset_vivariumHMB/l7-ancom-birthplace.qzv



%%%%%%%%other diversity analyses- FPVS %%%%%%%%
  #need to change file names accordingly
  
mkdir subset_FPVS/Output
mkdir subset_FPVS/Visualization
mkdir subset_FPVS/exported

qiime feature-table summarize \
--i-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--m-sample-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--o-visualization subset_FPVS/Visualization/Summary_HMBfinal_FPVS_featuretable_filt.qzv
#Saved subset_FPVS/Visualization to: subset_FPVS/Visualization/Summary_featuretableMM.qzv

%%bash -e
qiime taxa barplot \
--i-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--o-visualization subset_FPVS/Visualization/HMBfinal_FPVS_barplot
#Saved subset_FPVS/Visualization to: subset_FPVS/Visualization/_barplot.qzv

qiime diversity alpha-rarefaction \
--i-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--p-max-depth 1000 \
--p-metrics 'observed_features' \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--p-min-depth 5 \
--p-steps 10 \
--p-iterations 50 \
--o-visualization subset_FPVS/Visualization/JK_otu_max_FPVS.qzv

qiime diversity alpha-rarefaction \
--i-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--p-max-depth 1000 \
--p-metrics 'shannon' \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--p-min-depth 5 \
--p-steps 10 \
--p-iterations 50 \
--o-visualization subset_FPVS/Visualization/JK_shannon_max_FPVS.qzv

%%%%alpha%%%%
  
  %%bash -e
qiime diversity alpha \
--i-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--p-metric observed_features \
--o-alpha-diversity subset_FPVS/Output/HMBfinal_FPVS_alpha_otus

%%bash -e
qiime diversity alpha-correlation \
--i-alpha-diversity subset_FPVS/Output/HMBfinal_FPVS_alpha_otus.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--o-visualization subset_FPVS/Visualization/HMBfinal_FPVS_alpha_otus
#Saved subset_FPVS/Visualization to: subset_FPVS/Visualization/HMBfinal_FPVS_alpha_otus.qzv

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_FPVS/Output/HMBfinal_FPVS_alpha_otus.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--o-visualization subset_FPVS/Visualization/HMBfinal_FPVS_alpha_otus_usethis
#Saved subset_FPVS/Visualization to: subset_FPVS/Visualization/HMBfinal_FPVS_alpha_otus_usethis.qzv
qiime tools export \
--input-path subset_FPVS/Output/HMBfinal_FPVS_alpha_otus.qza \
--output-path subset_FPVS/Output/HMBfinal_FPVS_alpha_otus

%%bash -e
qiime diversity alpha \
--i-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--p-metric shannon \
--o-alpha-diversity subset_FPVS/Output/HMBfinal_FPVS_alpha_shannon
#Saved SampleData[AlphaDiversity] to: subset_FPVS/Output/HMBfinal_FPVS_alpha_shannon.qza

%%bash -e
qiime diversity alpha-correlation \
--i-alpha-diversity subset_FPVS/Output/HMBfinal_FPVS_alpha_shannon.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--o-visualization subset_FPVS/Visualization/HMBfinal_FPVS_alpha_shannon

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_FPVS/Output/HMBfinal_FPVS_alpha_shannon.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--o-visualization subset_FPVS/Visualization/HMBfinal_FPVS_alpha_shannon
#Saved subset_FPVS/Visualization to: subset_FPVS/Visualization/HMBfinal_FPVS_alpha_shannon.qzv

%%%%beta%%%%
"
%%bash -e
qiime diversity beta \
--i-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_FPVS/Output/HMBfinal_FPVS_BCDmatrix
#Saved DistanceMatrix to: subset_FPVS/Output/_BCDmatrix.qza
      
%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_FPVS/Output/HMBfinal_FPVS_BCDmatrix.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--m-metadata-column body_site \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization subset_FPVS/Visualization/HMBfinal_FPVS_Beta_diversity_permanova
#Saved subset_FPVS/Visualization to: subset_FPVS/Visualization/Beta_diversity_permanova.qzv
"
#OUTPUTS {#f9e,5}
%%bash -e
qiime diversity beta \
--i-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_FPVS/Output/HMBfinal_FPVS_BCDmatrix
#Saved DistanceMatrix to: subset_FPVS/Output/_BCDmatrix.qza

%%bash -e
qiime diversity pcoa \
--i-distance-matrix subset_FPVS/Output/HMBfinal_FPVS_BCDmatrix.qza \
--p-number-of-dimensions 3 \
--o-pcoa subset_FPVS/Output/HMBfinal_FPVS_PCoA_BCD
#Saved PCoAResults to: subset_FPVS/Output/_PCoA_BCD.qza

%%bash -e
qiime emperor plot \
--i-pcoa subset_FPVS/Output/HMBfinal_FPVS_PCoA_BCD.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--p-ignore-missing-samples False \
--o-visualization subset_FPVS/Visualization/HMBfinal_FPVS_PCoA_BCD
#Saved subset_FPVS/Visualization to: subset_FPVS/Visualization/_PCoA_BCD.qzv

%%bash -e
qiime feature-table relative-frequency \
--i-table subset_FPVS/OTU_filtdecontamfilt_FPVS.qza \
--o-relative-frequency-table subset_FPVS/Output/HMBfinal_FPVS_relative
#Saved FeatureTable[RelativeFrequency] to: subset_FPVS/Output/_relative.qza

%%bash -e
qiime diversity pcoa-biplot \
--i-pcoa subset_FPVS/Output/HMBfinal_FPVS_PCoA_BCD.qza \
--i-features subset_FPVS/Output/HMBfinal_FPVS_relative.qza \
--o-biplot subset_FPVS/Output/HMBfinal_FPVS_pcoa_BCD_biplot 
#Saved PCoAResults % Properties('biplot') to: subset_FPVS/Output/pcoa_BCD_biplot.qza

%%bash -e
qiime emperor biplot \
--i-biplot subset_FPVS/Output/HMBfinal_FPVS_pcoa_BCD_biplot.qza \
--m-sample-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--m-feature-metadata-file mergedstudy.gg2.tabletaxonomy.qza \
--o-visualization subset_FPVS/Visualization/HMBfinal_FPVS_pcoaBCD_biplot
#Saved subset_FPVS/Visualization to: subset_FPVS/Visualization/pcoaBCD_biplot.qzv

%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_FPVS/Output/HMBfinal_FPVS_BCDmatrix.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--m-metadata-column body_site \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_FPVS_Beta_diversity_permanova
#Saved subset_FPVS/Visualization to: subset_FPVS/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix subset_FPVS/Output/HMBfinal_FPVS_BCDmatrix.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--m-metadata-column body_site \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_FPVS_Beta_diversity_permdisp





%%%%%%%%other diversity analyses- estrousstaging %%%%%%%%
  #need to change file names accordingly
  
mkdir subset_estrousstaging/Output
mkdir subset_estrousstaging/Visualization
mkdir subset_estrousstaging/exported

qiime feature-table summarize \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--m-sample-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--o-visualization subset_estrousstaging/Visualization/Summary_HMBfinal_estrousstaging_featuretable_filt.qzv
#Saved subset_estrousstaging/Visualization to: subset_estrousstaging/Visualization/Summary_featuretableMM.qzv

%%bash -e
qiime taxa barplot \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--o-visualization subset_estrousstaging/Visualization/HMBfinal_estrousstaging_barplot
#Saved subset_estrousstaging/Visualization to: subset_estrousstaging/Visualization/_barplot.qzv

qiime diversity alpha-rarefaction \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--p-max-depth 1000 \
--p-metrics 'observed_features' \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--p-min-depth 5 \
--p-steps 25 \
--p-iterations 50 \
--o-visualization subset_estrousstaging/Visualization/JK_otu_max_es.qzv

qiime diversity alpha-rarefaction \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--p-max-depth 1000 \
--p-metrics 'shannon' \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--p-min-depth 5 \
--p-steps 25 \
--p-iterations 50 \
--o-visualization subset_estrousstaging/Visualization/JK_shannon_max_es.qzv

%%%%alpha%%%%
  
  %%bash -e
qiime diversity alpha \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--p-metric observed_features \
--o-alpha-diversity subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_otus

%%bash -e
qiime diversity alpha-correlation \
--i-alpha-diversity subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_otus.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--o-visualization subset_estrousstaging/Visualization/HMBfinal_estrousstaging_alpha_otus

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_otus.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--o-visualization subset_estrousstaging/Visualization/HMBfinal_estrousstaging_alpha_otus_usethis
qiime tools export \
--input-path subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_otus.qza \
--output-path subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_otus

%%bash -e
qiime diversity alpha \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--p-metric shannon \
--o-alpha-diversity subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_shannon
#Saved SampleData[AlphaDiversity] to: subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_shannon.qza

%%bash -e
qiime diversity alpha-group-significance \
--i-alpha-diversity subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_shannon.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--o-visualization subset_estrousstaging/Visualization/HMBfinal_estrousstaging_alpha_shannon

%%%%beta%%%%
"
%%bash -e
qiime diversity beta \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_estrousstaging/Output/HMBfinal_estrousstaging_BCDmatrix
      
%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_estrousstaging/Output/HMBfinal_estrousstaging_BCDmatrix.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column estrous_stage \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization subset_estrousstaging/Visualization/HMBfinal_estrousstaging_Beta_diversity_permanova
"

%%bash -e
qiime diversity beta \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--p-metric braycurtis \
--p-pseudocount 1 \
--p-n-jobs 1 \
--o-distance-matrix subset_estrousstaging/Output/HMBfinal_estrousstaging_BCDmatrix

%%bash -e
qiime diversity pcoa \
--i-distance-matrix subset_estrousstaging/Output/HMBfinal_estrousstaging_BCDmatrix.qza \
--p-number-of-dimensions 3 \
--o-pcoa subset_estrousstaging/Output/HMBfinal_estrousstaging_PCoA_BCD

%%bash -e
qiime emperor plot \
--i-pcoa subset_estrousstaging/Output/HMBfinal_estrousstaging_PCoA_BCD.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--p-ignore-missing-samples False \
--o-visualization subset_estrousstaging/Visualization/HMBfinal_estrousstaging_PCoA_BCD

%%bash -e
qiime feature-table relative-frequency \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--o-relative-frequency-table subset_estrousstaging/Output/HMBfinal_estrousstaging_relative

%%bash -e
qiime diversity pcoa-biplot \
--i-pcoa subset_estrousstaging/Output/HMBfinal_estrousstaging_PCoA_BCD.qza \
--i-features subset_estrousstaging/Output/HMBfinal_estrousstaging_relative.qza \
--o-biplot subset_estrousstaging/Output/HMBfinal_estrousstaging_pcoa_BCD_biplot 

%%bash -e
qiime emperor biplot \
--i-biplot subset_estrousstaging/Output/HMBfinal_estrousstaging_pcoa_BCD_biplot.qza \
--m-sample-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-feature-metadata-file mergedstudy.gg2.tabletaxonomy.qza \
--o-visualization subset_estrousstaging/Visualization/HMBfinal_estrousstaging_pcoaBCD_biplot

%%bash -e
qiime diversity beta-group-significance \
--i-distance-matrix subset_estrousstaging/Output/HMBfinal_estrousstaging_BCDmatrix.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column estrous_stage_A_estrus \
--p-method permanova \
--p-pairwise TRUE \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_estrousstaging_Beta_diversity_permanova_ESTRUS

qiime diversity beta-group-significance \
--i-distance-matrix subset_estrousstaging/Output/HMBfinal_estrousstaging_BCDmatrix.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column estrous_stage_A_estrus \
--p-method permdisp \
--p-pairwise TRUE \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_estrousstaging_Beta_diversity_permdisp_ESTRUS

qiime diversity beta-group-significance \
--i-distance-matrix subset_estrousstaging/Output/HMBfinal_estrousstaging_BCDmatrix.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column CST_char \
--p-method permanova \
--p-pairwise TRUE \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_CST_Beta_diversity_permanova

qiime diversity beta-group-significance \
--i-distance-matrix subset_estrousstaging/Output/HMBfinal_estrousstaging_BCDmatrix.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column CST_char \
--p-method permdisp \
--p-pairwise TRUE \
--p-permutations 999 \
--o-visualization Visualization/unrarefied_pair-wise/HMBfinal_CST_Beta_diversity_permdisp
#qiime diversity filter-distance-matrix \
#--i-distance-matrix subset_estrousstaging/Output/HMBfinal_estrousstaging_BCDmatrix.qza \
#--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
#--p-where "mouse_line='HMB'" \
#--o-filtered-distance-matrix subset_longitudinal/Output/HMBfinal_longitudinal_filtered-distance-matrix.qza

%%bash -e
qiime longitudinal first-distances \
--i-distance-matrix subset_estrousstaging/Output/HMBfinal_estrousstaging_BCDmatrix.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--p-state-column tp_days \
--p-individual-id-column mouse_origid \
--p-replicate-handling error \
--o-first-distances subset_estrousstaging/Output/HMBfinal_estrous_BC_FirstDistances 
#Saved SampleData[FirstDifferences] to: _BC_FirstDistances.qza

qiime longitudinal volatility \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging_noNAN.txt \
--m-metadata-file subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_shannon.qza \
--m-metadata-file subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_otus.qza \
--p-default-metric shannon_entropy \
--p-default-group-column CST_char \
--p-state-column tp_days \
--p-individual-id-column mouse_origid \
--o-visualization subset_estrousstaging/Visualization/shannon_estvolatility.qzv

qiime longitudinal volatility \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous_relfreq.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging_noNAN.txt \
--m-metadata-file subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_shannon.qza \
--m-metadata-file subset_estrousstaging/Output/HMBfinal_estrousstaging_alpha_otus.qza \
--p-default-metric shannon_entropy \
--p-default-group-column CST_char \
--p-state-column tp_days \
--p-individual-id-column mouse_origid \
--o-visualization subset_estrousstaging/Visualization/CST_estvolatility.qzv

qiime longitudinal feature-volatility \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging_noNAN.txt \
--p-state-column tp_days \
--p-individual-id-column mouse_origid \
--p-n-estimators 20 \
--output-dir subset_estrousstaging/feature_volatility

qiime tools export \
--input-path subset_estrousstaging/Output/HMBfinal_estrous_BC_FirstDistances.qza \
--output-path subset_estrousstaging/exported/HMBfinal_estrous_BC_FirstDistances

qiime taxa collapse \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-level 6 \
--o-collapsed-table subset_estrousstaging/HMB_estrousstaging-collapsedl6tableHMB.qza

qiime composition add-pseudocount \
--i-table subset_estrousstaging/HMB_estrousstaging-collapsedl6tableHMB.qza \
--o-composition-table subset_estrousstaging/comp-HMB_estrousstaging-collapsedl6tableHMB.qza

qiime composition ancom \
--i-table subset_estrousstaging/comp-HMB_estrousstaging-collapsedl6tableHMB.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column estrous_stage \
--o-visualization subset_estrousstaging/l6-ancom-eststage_.qzv

qiime taxa collapse \
--i-table subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza \
--i-taxonomy mergedstudy.gg2.tabletaxonomy.qza \
--p-level 7 \
--o-collapsed-table subset_estrousstaging/HMB_estrousstaging-collapsedl7tableHMB.qza

qiime composition add-pseudocount \
--i-table subset_estrousstaging/HMB_estrousstaging-collapsedl7tableHMB.qza \
--o-composition-table subset_estrousstaging/comp-HMB_estrousstaging-collapsedl7tableHMB.qza

qiime composition ancom \
--i-table subset_estrousstaging/comp-HMB_estrousstaging-collapsedl7tableHMB.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column estrous_stage \
--o-visualization subset_estrousstaging/l7-ancom-eststage.qzv








#all beta
#gbs - don't need PCoA plot, thank heavens
# qiime diversity beta-rarefaction \
# --i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/subset_GBS/OTU_filtdecontamfilt_GBS.qza  \
# --i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
# --p-metric weighted_normalized_unifrac \
# --p-clustering-method upgma \
# --m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/HMBpaper_metadata_GBSsubset.txt \
# --p-sampling-depth 100 \
# --p-iterations 100 \
# --p-correlation-method spearman \
# --p-color-scheme BrBG_r \
# --o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/Visualization/JK_w_norm_unif_gbs.qzv \
# --output-dir /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/Visualization/Jackknife_beta_gbs
# 
# qiime diversity beta-rarefaction \
# --i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/subset_GBS/OTU_filtdecontamfilt_GBS.qza  \
# --i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
# --p-metric braycurtis \
# --p-clustering-method upgma \
# --m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/HMBpaper_metadata_GBSsubset.txt \
# --p-sampling-depth 100 \
# --p-iterations 100 \
# --p-correlation-method spearman \
# --p-color-scheme BrBG_r \
# --o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/Visualization/JK_w_norm_unif_gbs.qzv \
# --output-dir /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/Visualization/Jackknife_beta_gbs

#longitudinal
qiime diversity beta-rarefaction \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza  \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme RdGy_r \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/Visualization/JK_w_norm_unif_100long.qzv \
--output-dir /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/Visualization/Jackknife_betaWNU_100longitudinal

qiime diversity beta-rarefaction \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza  \
--p-metric braycurtis \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme RdGy_r \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/Visualization/JK_bc_100long.qzv \
--output-dir /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/Visualization/Jackknife_betaBC_100longitudinal
#first attempt ran depth of 40.
#second attempt ran depth 300 for 50 iterations

#vivarium - just run the vivarium HMB one and visualize the non-HMB subset for this figure, then show HMB as "overlay"
# qiime diversity beta-rarefaction \
# --i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivarium/OTU_filtdecontamfilt_vivarium.qza  \
# --i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
# --p-metric weighted_normalized_unifrac \
# --p-clustering-method upgma \
# --m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
# --p-sampling-depth 100 \
# --p-iterations 100 \
# --p-correlation-method spearman \
# --p-color-scheme RdYlBu_r \
# --o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivarium/Visualization/JK_w_norm_unif_vivairum.qzv
# 
# qiime diversity beta-rarefaction \
# --i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivarium/OTU_filtdecontamfilt_vivarium.qza  \
# --i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
# --p-metric braycurtis \
# --p-clustering-method upgma \
# --m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
# --p-sampling-depth 100 \
# --p-iterations 100 \
# --p-correlation-method spearman \
# --p-color-scheme RdYlBu_r \
# --o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivarium/Visualization/JK_bc_vivarium.qzv
#first attempt depth 100, iterations 20

#vivarium HMb
qiime diversity beta-rarefaction \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza  \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme RdYlBu_r \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivariumHMB/Visualization/JK_w_norm_unif_vivariumHMB.qzv

qiime diversity beta-rarefaction \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza  \
--p-metric braycurtis \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme RdYlBu_r \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivariumHMB/Visualization/JK_bc_vivariumHMB.qzv


#FPVS
qiime diversity beta-rarefaction \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_FPVS/OTU_filtdecontamfilt_FPVS.qza  \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_FPVS/HMBpaper_metadata_FPVS.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme PuOr_r \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_FPVS/Visualization/JK_w_norm_unif_FPVS.qzv \
--verbose
#running 08092023 at 9:07pm
#10:03 was on 8th iteration
#10:10 was on 10th iteration
#10:14 was on 12th iteration
#10:51 was on 26th iteration
#11:16 was on 40th iteration
#by 1:28ish all iterations were complete and code finished cleanup (new comand line)

qiime diversity beta-rarefaction \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_FPVS/OTU_filtdecontamfilt_FPVS.qza  \
--p-metric braycurtis \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_FPVS/HMBpaper_metadata_FPVS.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme PuOr_r \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_FPVS/Visualization/JK_bc_FPVS.qzv \
--verbose
#first attempt ran at 40 depth

#estrous staging
qiime diversity beta-rarefaction \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza  \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme BrBG \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/Visualization/JK_w_norm_unif_es.qzv

qiime diversity beta-rarefaction \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza  \
--p-metric braycurtis \
--p-clustering-method upgma \
--m-metadata-file /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--p-sampling-depth 100 \
--p-iterations 100 \
--p-correlation-method spearman \
--p-color-scheme BrBG \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/Visualization/JK_bc_es.qzv
#FPVS 080923 night, vivarium 081023 7am, longitudinal 081023 11am


"  
for i in {1..100}; do \


qiime diversity beta-phylogenetic \
--i-phylogeny Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--i-table Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Decontam_process/decontam_items_final/Output/filtered_OTU_filtdecontam.qza \
--p-metric weighted_normalized_unifrac \
--p-variance-adjusted TRUE \
--o-distance-matrix Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/
  

; done
"
cd ~
#best to make one giant matrix of entire feat table and then partition, ahould've added metaata column for each study with variables "yes" and "no"
qiime diversity beta-phylogenetic \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_GBS/OTU_filtdecontamfilt_GBS.qza  \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-variance-adjusted FALSE \
--p-bypass-tips FALSE \
--o-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_GBS

qiime diversity beta-phylogenetic \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/OTU_filtdecontamfilt_long_baseline.qza  \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-variance-adjusted FALSE \
--p-bypass-tips FALSE \
--o-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_longitudinal

qiime diversity beta-phylogenetic \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_vivariumHMB/OTU_filtdecontamfilt_vivariumHMB.qza  \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-variance-adjusted FALSE \
--p-bypass-tips FALSE \
--o-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_vivariumHMB

qiime diversity beta-phylogenetic \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_FPVS/OTU_filtdecontamfilt_FPVS.qza  \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-variance-adjusted FALSE \
--p-bypass-tips FALSE \
--o-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_FPVS

qiime diversity beta-phylogenetic \
--i-table /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/OTU_filtdecontamfilt_estrous.qza  \
--i-phylogeny /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/2022.10.phylogeny.asv.nwk.qza \
--p-metric weighted_normalized_unifrac \
--p-variance-adjusted FALSE \
--p-bypass-tips FALSE \
--o-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_estrous

#old diversity matrix with "TREU" variance adjustment
qiime tools export \
--input-path /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_GBS.qza \
--output-path  /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/Unrarefied_wnu_matrices_GBS_adj
#current v with FALSE
qiime tools export \
--input-path /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_GBS.qza \
--output-path  /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/exported/Unrarefied_wnu_matrices_GBS


#repeats of above tests but for wnu matrices
"GBS"
qiime diversity filter-distance-matrix \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_GBS.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSHMB.txt \
--p-where "FullDataSet='Yes'" \
--o-filtered-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_GBSHMB

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_GBS.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--m-metadata-column CST_char \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_GBS_wnupermanova
#Saved subset_GBS/Visualization to: subset_GBS/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_GBS.qza \
--m-metadata-file subset_GBS/HMBpaper_metadata_GBSsubset.txt \
--m-metadata-column CST_char \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_GBS_wnupermdisp

"#Longitudinal"
qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_longitudinal.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-metadata-column cohort_number \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_longCohortcombnum_wnupermanova
#Saved subset_longitudinal/Visualization to: subset_longitudinal/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_longitudinal.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-metadata-column cohort_number \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_longCohortcombnum_wnupermdisp

qiime tools export \
--input-path /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_longCohortcombnum_wnupermanova.qzv \
--output-path  /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/exported/HMBfinal_longCohortcombnum_wnupermanova
qiime tools export \
--input-path /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_longCohortcombnum_wnupermdisp.qzv \
--output-path  /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/exported/HMBfinal_longCohortcombnum_wnupermdisp

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_longitudinal.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-metadata-column vivarium \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_longVivarium_wnupermanova

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_longitudinal.qza \
--m-metadata-file subset_longitudinal/HMBpaper_metadata_longitudinalsubset.txt \
--m-metadata-column vivarium \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_longVivarium_wnupermdisp

"vivarium"
qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_vivarium_conv.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--m-metadata-column birthplace_UCSD \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_vivarium_birthplace_wnupermanova_UCSD
#Saved subset_vivarium/Visualization to: subset_vivarium/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_vivarium_conv.qza \
--m-metadata-file subset_vivarium/HMBpaper_metadata_vivariumsubset.txt \
--m-metadata-column birthplace_UCSD \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_vivarium_birthplace_wnupermdisp_UCSD

"vivariumHMB"
qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_vivariumHMB.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column site_line_HMB \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_vivariumHMB_birthplace_wnupermanova_HMB
#Saved subset_vivariumHMB/Visualization to: subset_vivariumHMB/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_vivariumHMB.qza \
--m-metadata-file subset_vivariumHMB/HMBpaper_metadata_vivariumHMB.txt \
--m-metadata-column site_line_HMB \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_vivariumHMB_birthplace_wnupermdisp_HMB
#uses site_line originally but added site_line_HMB  
"FPVS"
qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_FPVS.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--m-metadata-column body_site \
--p-method permanova \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_FPVS_wnupermanova
#Saved subset_FPVS/Visualization to: subset_FPVS/Visualization/Beta_diversity_permanova.qzv

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_FPVS.qza \
--m-metadata-file subset_FPVS/HMBpaper_metadata_FPVS.txt \
--m-metadata-column body_site \
--p-method permdisp \
--p-pairwise True \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_FPVS_wnupermdips

"estrous"
qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_estrous.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column estrous_stage \
--p-method permanova \
--p-pairwise TRUE \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_eststage_wnupermanova

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_estrous.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column estrous_stage \
--p-method permdisp \
--p-pairwise TRUE \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_eststage_wnupermdisp
#estrous_stage_A_proestrus

"CST"
qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_estrous.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column CST_char \
--p-method permanova \
--p-pairwise TRUE \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_CST_char_wnupermanova

qiime diversity beta-group-significance \
--i-distance-matrix /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Output/Unrarefied_wnu_matrices_estrous.qza \
--m-metadata-file subset_estrousstaging/HMBpaper_metadata_estrousstaging.txt \
--m-metadata-column CST_char \
--p-method permdisp \
--p-pairwise TRUE \
--p-permutations 999 \
--o-visualization /Users/marlydmejia/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/Visualization/Unrarefied_wnu_pair-wise/HMBfinal_CST_char_wnupermdisp

#------------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%% R studio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


######################CST-HEATMAP!!!
"Community State Types Baseline"
#perform in R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("factoextra")
library(factoextra) # clustering visualization
library(phyloseq) #OTU table handling
#package stats is automatically loaded in R, used for kmeans, wss, and heatmapping

setwd("~/Desktop")

%%%%%%%%%%%%%%% on estrous dataset %%%%%%%%%%%%%%%
#########################################################3
#use relative frequency file, consolidated by taxa and transposed to have sample_name column and taxa headings  
HMBfinalEstrous <- read.csv("Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/OTU_filtdecontamfilt_longitudinal_relfreq/HMB_baseline_taxahead_allincluded.csv", header=TRUE, sep = ",", row.names="sample_name")
meta_HMBfinalEstrousstages2 <- read.csv("Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_longitudinal/HMBpaper_metadata_estrous.csv", header=TRUE, sep = ",", row.names="sample_name")

setwd("~/Desktop/Patras_Lab/Proposals/Manuscripts/Mine/2022_HMB_characterization/16S_analyses_complete/subset_estrousstaging/CST")

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
rect.hclust(fitHMBfinalEs, k=3, border="lightblue")
rect.hclust(fitHMBfinalEs, k=4, border="red")
rect.hclust(fitHMBfinalEs, k=5, border="green")
rect.hclust(fitHMBfinalEs, k=6, border="Blue")
rect.hclust(fitHMBfinalEs, k=7, border="purple")
rect.hclust(fitHMBfinalEs, k=8, border="orange")

print(group5HMBfinalEs)   
write.csv(group5HMBfinalEs,"newCST5_HMBfinalEs_mapping.csv")
write.csv(group6HMBfinalEs,"newCST6_HMBfinalEs_mapping.csv")
write.csv(group7HMBfinalEs,"newCST7_HMBfinalEs_mapping.csv")

#printed and wrote CST 5 - CST 7 just to see who would be categorized differently
#Manually add first column heading "sample_name"
meta_HMBfinal7Estrousstudy <- read.csv("newCST7_HMBfinalEs_mapping.csv", header=TRUE, sep = ",", row.names="sample_name")
meta_HMBfinal6Estrousstudy <- read.csv("newCST6_HMBfinalEs_mapping.csv", header=TRUE, sep = ",", row.names="sample_name")
meta_HMBfinal5Estrousstudy <- read.csv("newCST5_HMBfinalEs_mapping.csv", header=TRUE, sep = ",", row.names="sample_name")

#below is a metadata file I made one column at a time, but I already loaded the one above in the same order, so I don't need to recreate a metadata sheet.... I don't think
#CST_HMBinitEsstudy_metadata <- read.csv("CST_metadata_HMBinitEsstudy.csv", header=TRUE, sep = ",", row.names="sample_name")
#This is to select the top taxa and then I rearrange/ consolidate to "Other", but I have already done this
topN = 15
HMBfinalEstrousT <- t(HMBfinalEstrous)
HMBfinalEstrous_taxa <- otu_table(HMBfinalEstrousT, taxa_are_rows = TRUE)
most_abundant_taxaHMBinitEsstudy = sort(taxa_sums(HMBfinalEstrous_taxa), TRUE)[1:topN]
print(most_abundant_taxaHMBinitEsstudy)
write.csv(most_abundant_taxaHMBinitEsstudy,"top15HMBinitEsst_taxa.csv")

most_abundant_taxa20HMBinitEsst <- sort(taxa_sums(HMBfinalEstrous_taxa), TRUE)[1:20]
print(most_abundant_taxa20HMBinitEsst)
write.csv(most_abundant_taxa20HMBinitEsst,"top20GHMBinitEsst_taxa.csv")
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

OTU_HMBinitEsre_ordername30 <- read.csv("RelFreq_30taxahead_HMBinitEsstudy.csv", header=TRUE, sep = ",", row.names="sample_name")
OTU_clustHMBfinalEs30_matrix <- as.matrix(OTU_HMBinitEsre_ordername30)

#for CST
my_group_HMBfinalEs <- as.numeric(as.factor(substr(meta_HMBfinal7Estrousstudy$CST,1,1)))
mycolHMBfinalEs <- c("seashell2","lightgoldenrod2","lightsteelblue4","lemonchiffon2","sienna3","dodgerblue4","skyblue4")
colSideHMBfinalEs <- mycolHMBfinalEs[my_group_HMBfinalEs]

hclust_rows <- as.dendrogram(hclust(dist(as.matrix(HMBfinalEstrous),method="euclidean"),method="ward.D"))  # Calculate hclust dendrograms

heatmap(OTU_clustHMBfinalEs30_matrix, cexCol = 0.8, scale= "none", RowSideColors=colSideHMBfinalEs, col=colMainHMBfinalEs, Colv= NA, Rowv=hclust_rows)
heatmap(as.matrix(HMBfinalEstrous), cexCol = 0.8, scale= "none", RowSideColors=colSideHMBfinalEs, col=colMainHMBfinalEs, Colv= NA, distfun = function(x) dist(x, method="euclidean"), hclustfun = function(x) hclust(x, method="ward.D"))

legend(x="right",legend)
# 
# my_group_6HMBfinalEs <- as.numeric(as.factor(substr(meta_6HMBfinalEstrousstudy$CST,1,1)))
# mycol6HMBfinalEs <- c("tomato","yellow","orange","lightseagreen","royalblue", "pink")
# colSide6HMBfinalEs <- mycol6HMBfinalEs[my_group_6HMBfinalEs]
# heatmap(OTU_clustHMBfinalEs_matrix, cexCol = 0.8, scale= "none", RowSideColors=colSide6HMBfinalEs, col=colMainHMBfinalEs, Colv= NA, distfun = function(x) dist(x, method="euclidean"), hclustfun = function(x) hclust(x, method="ward.D"))

#need files to be only subject ID and parameter, so 2 columns total. But check if (1,1) means something about array column length and variable column

#for estrous stage
my_group_HMBfinalEstrx <- as.numeric(as.factor(substr(meta_HMBfinalEstrousstages2$estrous_stage,1,1)))
esstcol_blues <- c("#FFFFFF", "#6f8da6", "#E8E1A6", "#254f82", "#c1d8d1") 
#white-ish, med blue, dull yellow, dark blue, seafoam
#blank, diestrus, estrus, metestrus, proestrus
colSideHMBfinalEstrxblues <- esstcol_blues[my_group_HMBfinalEstrx]
heatmap(OTU_clustHMBfinalEs30_matrix, cexCol = 0.8, scale= "none", RowSideColors=colSideHMBfinalEstrxblues, col=colMainHMBfinalEs,Colv= NA, Rowv=hclust_rows)

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

#image legend scale for body of heatmap
image(1:nrow(OTU_clustHMBfinalEs30_matrix), 1, as.matrix(1:nrow(OTU_clustHMBfinalEs30_matrix)), 
      col=colorRampPalette(brewer.pal(9, "BuPu"))(25),
      xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
