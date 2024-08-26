############################################################################
######################### Step 2: Annotation ###############################
############################################################################

############################################################################
# Remove asterisks as InterProScan does not support them

sed -i 's/\*//g' cat_larvae_Trinity_transdecoder_cds_CD-HIT.pep.fasta > cat_larvae_Trinity_transdecoder_cds_CD-HIT.wo_stop.pep.fasta

############################################################################
# Using interproscan to annotate chemosensory genes "de novo"

interproscan.sh -i ccat_larvae_Trinity_transdecoder_cds_CD-HIT.wo_stop.pep.fasta -f TSV --goterms --pathways

############################################################################
# Using BLAST to annotate chemosensory genes based on known coleopteran chemosensory genes

makeblastdb -in cat_total_Trinity_transdecoder_cds_CD-HIT.fasta -dbtype nucl
tblastn -query OR_dataset_all_Clat.fasta -out db_ClattotalCDHIT_query_OR_Clat -evalue 0.001 -db cat_total_Trinity_transdecoder_cds_CD-HIT.fasta -outfmt 6
tblastn -query IR_dataset_all_Clat.fasta -out db_ClattotalCDHIT_query_IR_Clat -evalue 0.001 -db cat_total_Trinity_transdecoder_cds_CD-HIT.fasta -outfmt 6
tblastn -query GR_dataset_all_Clat.fasta -out db_ClattotalCDHIT_query_GR_Clat -evalue 0.001 -db cat_total_Trinity_transdecoder_cds_CD-HIT.fasta -outfmt 6
tblastn -query OBP_dataset_all_Clat.fasta -out db_ClattotalCDHIT_query_OBP_Clat -evalue 0.001 -db cat_total_Trinity_transdecoder_cds_CD-HIT.fasta -outfmt 6

############################################################################
# Parsing BLAST output

python3 parse_blast.py db_ClattotalCDHIT_query_OR_Clat

############################################################################
# Manual checking
# Cleaning the consensus transcriptome (removing artefactual duplications after checking closely chemosensory genes) -> cat_total_Trinity_transdecoder_cds_CD-HIT.fasta_sup300_final
