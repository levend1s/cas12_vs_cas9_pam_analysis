#!/bin/sh
set -e

WORKING_DIR=/Users/jlevendis/crypto_cas12_analysis
OFFSET_DISTANCE=100
OUTPUT_DIR=${WORKING_DIR}/pam_analysis
mkdir -p ${OUTPUT_DIR}

# using rqc commit d0ca7f0b4b06b77351d859c3bd2cdbdd195cb34d
RQC_PATH=${WORKING_DIR}/rqc/rqc.py

# Cryptosporidium parvum BGF
# PROJECT_PREFIX=cpbgf
# ANNOTATION=${WORKING_DIR}/CpBGF_Repository/Genomes/CpBGF_genome_V1.1.gff
# GENOME=${WORKING_DIR}/CpBGF_Repository/Genomes/CpBGF_genome_v1.fasta
# # In CPBGF, a protein coding gene will look like:
# # CP141118.1      Genbank gene    3333    4291    .       -       .       ID=gene-cpbgf_1004000;Name=cpbgf_1004000;gbkey=Gene;gene_biotype=protein_coding;locus_tag=cpbgf_1004000
# # So the below command makes sure the third column is "mRNA" and then extracts the gene name from the attributes column ($9). Specifically, we split on "|" and ";" and take the FOURTH element which is the gene name
# awk '$3=="mRNA" {split ($9,x,/[|-]/); print x[4]}' ${ANNOTATION} > ${OUTPUT_DIR}/${PROJECT_PREFIX}_protein_coding_genes.tsv

# Toxoplasma gondii ME49
PROJECT_PREFIX=toxo_me49
ANNOTATION=${WORKING_DIR}/TgondiiME49/gff/data/ToxoDB-68_TgondiiME49.gff
GENOME=${WORKING_DIR}/TgondiiME49/fasta/data/ToxoDB-68_TgondiiME49_Genome.fasta
# CM002034.1      Genbank mRNA    18521   26575   .       +       .       ID=rna-gnl|WGS:ABPA|mrna.TGME49_292920;Parent=gene-TGME49_292920;Note=transcript TGME49_292920;gbkey=mRNA;locus_tag=TGME49_292920;orig_protein_id=gnl|WGS:ABPA|cds.TGME49_292920;orig_transcript_id=gnl|WGS:ABPA|mrna.TGME49_292920;product=heat shock protein 75%2C putative
# So the below command makes sure the third column is "mRNA" and then extracts the gene name from the attributes column ($9). Specifically, we split on "." and ";" and take the SECOND element which is the gene name.
awk '$3=="protein_coding_gene" {split ($9,x,/[=;]/); print x[2]}' ${ANNOTATION} > ${OUTPUT_DIR}/${PROJECT_PREFIX}_protein_coding_genes.tsv

# Plasmodium falciparum 3D7
# PROJECT_PREFIX=pf3d7
# ANNOTATION=${WORKING_DIR}/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff
# GENOME=${WORKING_DIR}/Pfalciparum3D7/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta
# # Pf3D7_14_v3     VEuPathDB       protein_coding_gene     2940178 2941916 .       +       .       ID=PF3D7_1472000;Name=ISY1;description=pre-mRNA-splicing factor ISY1%2C putative;ebi_biotype=protein_coding
# # So the below command makes sure the third column is "mRNA" and then extracts the gene name from the attributes column ($9). Specifically, we split on "=" and ";" and take the second element which is the gene name.
# awk '$3=="protein_coding_gene" {split ($9,x,/[=;]/); print x[2]}' ${ANNOTATION} > ${OUTPUT_DIR}/${PROJECT_PREFIX}_protein_coding_genes.tsv

# find start and and stop codons
python ${RQC_PATH} motif_finder -a ${ANNOTATION} -g ${GENOME} -o ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codons_all.tsv --feature_filter CDS -m "^ATG"
python ${RQC_PATH} motif_finder -a ${ANNOTATION} -g ${GENOME} -o ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codons_all.tsv --feature_filter CDS -m "(TA[AG]|TGA)$"
python ${RQC_PATH} motif_finder -a ${ANNOTATION} -g ${GENOME} -o ${OUTPUT_DIR}/${PROJECT_PREFIX}_NGG_all.tsv -m ".GG"
python ${RQC_PATH} motif_finder -a ${ANNOTATION} -g ${GENOME} -o ${OUTPUT_DIR}/${PROJECT_PREFIX}_TTTN_all.tsv -m "TTT."

# filter to only include one start and stop codon per gene
# because we used a feature filter, the ID field looks like this
# CP141118.1      4201    4204    ATc     0       -       cpbgf_1004000-CDS1
# so we split on "-" and take the first element which is the gene name, store it in a new column called ID
# This means we just see if the base genem name exists, not isoform specific checking
# However we do a version sort on the 7th column (the feature_ID) so that we keep the first CDS seen (for start codon) seen, or last CDS (for stop codon, we reverse sort it and take the first one seen) for genes with multiple CDS features (exons)
printf "contig\tstart\tend\tname\tscore\tstrand\tfeature_ID\tID\n" > ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codons.tsv
printf "contig\tstart\tend\tname\tscore\tstrand\tfeature_ID\tID\n" > ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codons.tsv
awk '{split ($7,x,/[.-]/); print $0"\t"x[1]}' ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codons_all.tsv | tail -n +2 | sort -V -k7 - | awk '!seen[$8]++' >> ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codons.tsv
awk '{split ($7,x,/[.-]/); print $0"\t"x[1]}' ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codons_all.tsv | tail -n +2 | sort -rV -k7 - | awk '!seen[$8]++' >> ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codons.tsv

# calculate offsets for NGG and TTTN PAMs relative to start and stop codons
python ${RQC_PATH} calculate_offsets -d ${OFFSET_DISTANCE} -r ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codons.tsv -o ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codon_offsets.tsv -i NGG ${OUTPUT_DIR}/${PROJECT_PREFIX}_NGG_all.tsv TTTN ${OUTPUT_DIR}/${PROJECT_PREFIX}_TTTN_all.tsv
python ${RQC_PATH} calculate_offsets -d ${OFFSET_DISTANCE} -r ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codons.tsv -o ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codon_offsets.tsv -i NGG ${OUTPUT_DIR}/${PROJECT_PREFIX}_NGG_all.tsv TTTN ${OUTPUT_DIR}/${PROJECT_PREFIX}_TTTN_all.tsv

# plot offsets
python ${RQC_PATH} plot_relative_distance -l "start codon" -d ${OFFSET_DISTANCE} -i ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codon_offsets.tsv -o ${OUTPUT_DIR}/${PROJECT_PREFIX}_${OFFSET_DISTANCE}_start.eps
python ${RQC_PATH} plot_relative_distance -l "stop codon" -d ${OFFSET_DISTANCE} -i ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codon_offsets.tsv -o ${OUTPUT_DIR}/${PROJECT_PREFIX}_${OFFSET_DISTANCE}_stop.eps

# generate summary statistics, missing genes IDS etc
printf "looking for missing start codons...\n"
while read line; do grep -q "$line" ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codons.tsv || echo $line; done < ${OUTPUT_DIR}/${PROJECT_PREFIX}_protein_coding_genes.tsv > ${OUTPUT_DIR}/${PROJECT_PREFIX}_missing_start_codon_ids.tsv
printf "looking for missing stop codons...\n"
while read line; do grep -q "$line" ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codons.tsv || echo $line; done < ${OUTPUT_DIR}/${PROJECT_PREFIX}_protein_coding_genes.tsv > ${OUTPUT_DIR}/${PROJECT_PREFIX}_missing_stop_codon_ids.tsv

printf "number of protein coding genes: "; wc -l < ${OUTPUT_DIR}/${PROJECT_PREFIX}_protein_coding_genes.tsv # no header
printf "number of genes with start codons: "; tail -n +2 ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codons.tsv | wc -l
printf "number of genes with stop codons: "; tail -n +2 ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codons.tsv | wc -l
printf "number of genes missing start codons: "; wc -l < ${OUTPUT_DIR}/${PROJECT_PREFIX}_missing_start_codon_ids.tsv # no header
printf "number of genes missing stop codons: "; wc -l < ${OUTPUT_DIR}/${PROJECT_PREFIX}_missing_stop_codon_ids.tsv # no header
printf "total NGG in genome: "; tail -n +2 ${OUTPUT_DIR}/${PROJECT_PREFIX}_NGG_all.tsv | wc -l
printf "total TTTN in genome: "; tail -n +2 ${OUTPUT_DIR}/${PROJECT_PREFIX}_TTTN_all.tsv | wc -l
printf "average NGG count around start codons: "; awk -F'\t' 'NR>1 {sum += $4; count++} END { if (NR > 0) print sum / count}' ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codon_offsets.tsv
printf "average TTTN count around start codons: "; awk -F'\t' 'NR>1 {sum += $6; count++} END { if (NR > 0) print sum / count}' ${OUTPUT_DIR}/${PROJECT_PREFIX}_start_codon_offsets.tsv
printf "average NGG count around stop codons: "; awk -F'\t' 'NR>1 {sum += $4; count++} END { if (NR > 0) print sum / count}' ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codon_offsets.tsv
printf "average TTTN count around stop codons: "; awk -F'\t' 'NR>1 {sum += $6; count++} END { if (NR > 0) print sum / count}' ${OUTPUT_DIR}/${PROJECT_PREFIX}_stop_codon_offsets.tsv