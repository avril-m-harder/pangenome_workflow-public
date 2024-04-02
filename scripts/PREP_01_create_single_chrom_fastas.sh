#!/bin/bash
#$ -N sorghum_chrom_split
#$ -cwd
#$ -V
#$ -q all.q
#$ -pe smp 2
#$ -m ae
#$ -M #######@hudsonalpha.org

source ~/.bashrc
source /home/#######_scratch_f13/pangenome_workflow/scripts/99_init_script_vars.sh

TMP_DIR=`/bin/mktemp -d -p /mnt/data1/tmp`
cd ${TMP_DIR}
echo "${TMP_DIR}"

MAMBA="/home/#######/.conda/envs"


# -----------------------------------------------------------------------------
# Copy assemblies to REF_DIR -- EDIT IF NECESSARY
# -----------------------------------------------------------------------------

for FILE in $(find ${REF_DIR} -name "*mainGenome.fasta")
do
	cp ${FILE} ${REF_DIR}
done

# -----------------------------------------------------------------------------
# Create chroms-only FASTA files -- EDIT IF NECESSARY
# -----------------------------------------------------------------------------

mamba activate ${MAMBA}/bioinfo_tools

for f in $(ls ${REF_DIR}/*.mainGenome.fasta)
do
	samtools faidx ${f} > ${REF_DIR}/${f}.fai

	base=$(basename ${f} .mainGenome.fasta)
	for c in ${chroms[@]}
	do
		echo ${c} > tmp.list
		seqtk subseq ${f} tmp.list | bgzip > ${base}.${c}.fa.gz
	done
done

mamba deactivate

# -----------------------------------------------------------------------------
# Clean up temp dir
# -----------------------------------------------------------------------------

rsync -avuP ${TMP_DIR}/*.fa* ${REF_DIR}
cd ${OUTDIR}
rm -rf ${TMP_DIR}