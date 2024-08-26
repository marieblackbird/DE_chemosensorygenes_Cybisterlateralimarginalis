############################################################################
######################### Step 3: Quantification ###########################
############################################################################

############################################################################
# Using Kallisto to quantify the expression of each transcript from the consensus 
# transcriptomic assembly using Illumina reads (cleaned with Trimmomatic) per condition
############################################################################

# Creating the index
kallisto index -i cat_total_Trinity_transdecoder_cds_CD-HIT.fasta_sup300_final.idx cat_total_Trinity_transdecoder_cds_CD-HIT.fasta_sup300_final

# Creating one quantification file for each organ per stage per replicat

for i in Clat_larva_Ant1 Clat_larva_Ant2 Clat_larva_Ant3 #etc
do
  	kallisto quant -i cat_total_Trinity_transdecoder_cds_CD-HIT.fasta_sup300_final.idx -o $i.quantif \
$i.output_forward_paired.fq.gz $i.output_reverse_paired.fq.gz
done
