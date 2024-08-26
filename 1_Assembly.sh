############################################################################
######################## Step 1: Assembly ##################################
############################################################################

############################################################################
# Using Trimmomatic to clean Illumina reads

for i in Clat_larva_Ant1 Clat_larva_Ant2 Clat_larva_Ant3 #etc
do

	echo "trimming $i"
     trimmomatic PE -phred33 \
         -trimlog "$i.trimmomatic.log" \
         "raw_reads/${i}_1.fq.gz" "raw_reads/${i}_2.fq.gz" \
         "${i}_output_forward_paired.fq.gz" "${i}_output_forward_unpaired.fq.gz" \
         "${i}_output_reverse_paired.fq.gz" "${i}_output_reverse_unpaired.fq.gz" \
         HEADCROP:10 TRAILING:20 MINLEN:30 SLIDINGWINDOW:4:20
done

############################################################################
# Transcriptomic assembly (one per triplicate)

Trinity --seqType fq --samples_file sample_file.txt --max_memory 50G

############################################################################
# Running TransDecoder on each transcriptome to identify coding sequences

TransDecoder.LongOrfs -t Trinity_larvae_ant.fasta -m 50
TransDecoder.Predict -t Trinity_larvae_ant.fasta

############################################################################
# Using BUSCO to estimate the completeness of each transcriptome by organ

busco -i larva_antenna_Trinity_transdecoder_cds.fasta -l endopterygota_odb10 -m transcriptome -o larva_antenna_Trinity_transdecoder_cds.fasta.busco

############################################################################
# Performing statistics on each transcriptome assembly by organ

TrinityStats.pl larva_antenna_Trinity_transdecoder_cds.fasta
quast.py larva_antenna_Trinity_transdecoder_cds.fasta

############################################################################
# Concatenation of all transcriptomes in a file "cat_total_Trinity_transdecoder_cds.fasta"

cat *_Trinity_transdecoder_cds.fasta > cat_total_Trinity_transdecoder_cds.fasta

############################################################################
# Using cd-HIT for clustering and removing redundancies to create a global transcriptome assembly

cd-hit-est -i cat_total_Trinity_transdecoder_cds.fasta -o cat_total_Trinity_transdecoder_cds_CD-HIT.fasta -n 8 -c 0.9 -M 1000

############################################################################
# Removing very small contigs from cat_total_Trinity_transdecoder_cds_CD-HIT.fasta

python3 remove_small_contigs.py cat_total_Trinity_transdecoder_cds_CD-HIT.fasta

############################################################################
# Using BUSCO to estimate the completeness of the consensus

busco -i cat_total_Trinity_transdecoder_cds_CD-HIT.fasta_sup300 -l endopterygota_odb10 -m transcriptome -o cat_total_Trinity_transdecoder_cds_CD-HIT.fasta_sup300.busco

############################################################################
# Performing statistics on transcriptome assembly

TrinityStats.pl cat_total_Trinity_transdecoder_cds_CD-HIT.fasta_sup300
quast.py cat_total_Trinity_transdecoder_cds_CD-HIT.fasta_sup300

############################################################################
# Translating CDS into amino acid sequences

transeq cat_total_Trinity_transdecoder_cds_CD-HIT.fasta cat_total_Trinity_transdecoder_cds_CD-HIT.pep.fasta -sformat pearson

