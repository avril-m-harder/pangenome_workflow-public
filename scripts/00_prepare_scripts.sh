#!/bin/bash
#
# prepare_scripts.sh
#

source /home/aharder_scratch_f13/pangenome_workflow/scripts/99_init_script_vars.sh 

# -----------------------------------------------------------------------------
# Set 01.sh Minigraph-Cactus options
# -----------------------------------------------------------------------------
## Chrom-level script directory
MC_CHROM_SCRIPT_DIR="${SCRIPT_DIR}/01_minigraph-cactus_chrom_scripts"
mkdir -p ${MC_CHROM_SCRIPT_DIR}

## Prep chrom-level scripts
for CHROM in ${chroms[@]}
do

	sed "s/sCHROM/${CHROM}/g" ${TEMPLATE_DIR}/01_minigraph-cactus_sCHROM.TEMPLATE.sh > \
	${MC_CHROM_SCRIPT_DIR}/01_minigraph-cactus_${CHROM}.sh
	
	sed -i "s/MCNTHREADS/${MC_NTHREADS}/g" \
	${MC_CHROM_SCRIPT_DIR}/01_minigraph-cactus_${CHROM}.sh
	
	sed -i "s/QUEUE01/${QUEUE_01}/g" \
	${MC_CHROM_SCRIPT_DIR}/01_minigraph-cactus_${CHROM}.sh
	
	sed -i "s/pangenome_workflow/${WDNAME}/g" \
	${MC_CHROM_SCRIPT_DIR}/01_minigraph-cactus_${CHROM}.sh

done

# -----------------------------------------------------------------------------
# Set 02.sh vg haplotype options
# -----------------------------------------------------------------------------
sed "s/VGHAPNTHREADS/${VG_HAP_NTHREADS}/g" \
${TEMPLATE_DIR}/02_vg_combine_and_haplos.TEMPLATE.sh > \
${SCRIPT_DIR}/02_vg_combine_and_haplos.sh
	
sed -i "s/QUEUE02/${QUEUE_02}/g" \
${SCRIPT_DIR}/02_vg_combine_and_haplos.sh
	
sed -i "s/pangenome_workflow/${WDNAME}/g" \
${SCRIPT_DIR}/02_vg_combine_and_haplos.sh

# -----------------------------------------------------------------------------
# Set 03.sh vg giraffe options
# -----------------------------------------------------------------------------
## Sample-level script directory
VG_SAMP_SCRIPT_DIR="${SCRIPT_DIR}/03_vg_giraffe-subgraph_samp_scripts"
mkdir -p ${VG_SAMP_SCRIPT_DIR}

## prepare correct scripts depending on whether FASTQs are interleaved or not
NCOLS=$(awk --field-separator="\t" "{ print NF }" ${FQ_LIST} | uniq)

if [ $NCOLS -eq 3 ]; then
	echo "separate FASTQ files"

	while read -a line
	do
		
		sed "s/SAMP_NAME/${line[0]}/g" \
		${TEMPLATE_DIR}/03_vg_giraffe-subgraph.TEMPLATE.sh > \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
		
		sed -i "s|FQ_FILE_1|${line[1]}|g" \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
		
		sed -i "s|FQ_FILE_2|${line[2]}|g" \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
	
		sed -i "s/QUEUE03/${QUEUE_03}/g" \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
	
		sed -i "s/VGGIRNTHREADS/${VG_GIR_NTHREADS}/g" \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
		
		sed -i "s/pangenome_workflow/${WDNAME}/g" \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
	
	done < ${FQ_LIST}


elif [ $NCOLS -eq 2 ]; then
	echo "interleaved"
	
	while read -a line
	do
		
		sed "s/SAMP_NAME/${line[0]}/g" \
		${TEMPLATE_DIR}/03_vg_giraffe-subgraph.TEMPLATE-INTERLEAVED.sh > \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
		
		sed -i "s|FQ_FILE|${line[1]}|g" \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
	
		sed -i "s/QUEUE03/${QUEUE_03}/g" \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
	
		sed -i "s/VGGIRNTHREADS/${VG_GIR_NTHREADS}/g" \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
		
		sed -i "s/pangenome_workflow/${WDNAME}/g" \
		${VG_SAMP_SCRIPT_DIR}/03_vg_giraffe-subgraph_${line[0]}.sh
	
	done < ${FQ_LIST}
	
else
	echo "incorrect number of columns in ${FQ_LIST}"
fi


# -----------------------------------------------------------------------------
# Set 04.sh vg call options
# -----------------------------------------------------------------------------
## Sample-level script directory
VG_CALL_SAMP_SCRIPT_DIR="${SCRIPT_DIR}/04_vg_call_SVs_samp_scripts"
mkdir -p ${VG_CALL_SAMP_SCRIPT_DIR}

while read -a line
do
		
	sed "s/SAMP_NAME/${line[0]}/g" \
		${TEMPLATE_DIR}/04_vg_call_SVs.TEMPLATE.sh > \
		${VG_CALL_SAMP_SCRIPT_DIR}/04_vg_call_SV_${line[0]}.sh
	
	sed -i "s/QUEUE04/${QUEUE_04}/g" \
	${VG_CALL_SAMP_SCRIPT_DIR}/04_vg_call_SV_${line[0]}.sh
	
	sed -i "s/VGCALLNTHREADS/${VG_CALL_NTHREADS}/g" \
	${VG_CALL_SAMP_SCRIPT_DIR}/04_vg_call_SV_${line[0]}.sh
	
	sed -i "s/pangenome_workflow/${WDNAME}/g" \
	${VG_CALL_SAMP_SCRIPT_DIR}/04_vg_call_SV_${line[0]}.sh
		
done < ${FQ_LIST}


# -----------------------------------------------------------------------------
# Set 05.sh varscan options
# -----------------------------------------------------------------------------
## Sample-level script directory
VARSCAN_SAMP_SCRIPT_DIR="${SCRIPT_DIR}/05_varscan_SNPs_indels_samp_scripts"
mkdir -p ${VARSCAN_SAMP_SCRIPT_DIR}

while read -a line
do
		
	sed "s/SAMP_NAME/${line[0]}/g" \
		${TEMPLATE_DIR}/05_varscan_SNPs_indels.TEMPLATE.sh > \
		${VARSCAN_SAMP_SCRIPT_DIR}/05_varscan_SNPs_indels_${line[0]}.sh
	
	sed -i "s/QUEUE05/${QUEUE_05}/g" \
	${VARSCAN_SAMP_SCRIPT_DIR}/05_varscan_SNPs_indels_${line[0]}.sh
	
	sed -i "s/VARSCANNTHREADS/${VG_CALL_NTHREADS}/g" \
	${VARSCAN_SAMP_SCRIPT_DIR}/05_varscan_SNPs_indels_${line[0]}.sh
	
	sed -i "s/pangenome_workflow/${WDNAME}/g" \
	${VARSCAN_SAMP_SCRIPT_DIR}/05_varscan_SNPs_indels_${line[0]}.sh
		
done < ${FQ_LIST}

# -----------------------------------------------------------------------------
# Set accessory script options
# -----------------------------------------------------------------------------
sed -i "s/pangenome_workflow/${WDNAME}/g" \
${SCRIPT_DIR}/PREP_01_create_single_chrom_fastas.sh

sed -i "s/pangenome_workflow/${WDNAME}/g" \
${SCRIPT_DIR}/accessory_functions/ACC_01_odgi_ad_hoc_plots.sh

sed -i "s/pangenome_workflow/${WDNAME}/g" \
${SCRIPT_DIR}/accessory_functions/ACC_02_clipped_seq_analysis.sh

sed -i "s/pangenome_workflow/${WDNAME}/g" \
${SCRIPT_DIR}/accessory_functions/ACC_03_getting_coordinates_from_node_IDs.sh
