#! /bin/bash
#$ -q "all.q"
#$ -pe smp 2
#$ -cwd
#$ -V
#$ -N odgi_bed_plots
#$ -m ae
#$ -M [email]

source ~/.bashrc
source /home/aharder_scratch_f13/pangenome_workflow/scripts/99_init_script_vars.sh

## create temp directory and echo path to .o file
TMP_DIR=`/bin/mktemp -d -p /mnt/data1/tmp`
cd ${TMP_DIR}
echo "${TMP_DIR}"

# -----------------------------------------------------------------------------
# odgi plot for ROIs
#
# format for input text file should be 1 ROI per line:
#
# chrom / path / start / end
#
# -----------------------------------------------------------------------------

mamba activate pg_tools

INDIR="${ODGI_OUT}"
OUTDIR="${ODGI_OUT}"

while read -a line
do

	CHROM=${line[0]}
	IN_BED="${CHROM}_${line[2]}-${line[3]}.bed"
	echo ${line[1]}$'\t'${line[2]}$'\t'${line[3]}$'\n' > ${IN_BED}
	base=$(basename ${IN_BED} .bed)
	IN_OG=$(find ${MC_OUT_DIR} -name "*${CHROM}.sorted.og")
	OUT_OG="${base}.og"
	SORT_OUT_OG="${base}.sorted.og"
	OUT_PNG="${base}.sorted.png"
	
	if [ ! -f ${OUTDIR}/${OUT_PNG} ]; then
		
		rsync -avuP ${IN_OG} .
		IN_OG=$(basename ${IN_OG})
		
		odgi extract \
			-i ${IN_OG} \
			-o ${OUT_OG} \
			-b ${IN_BED} \
			-c 0 \
			-E \
			-P

		odgi sort -i ${OUT_OG} -o ${SORT_OUT_OG}
		odgi viz -i ${SORT_OUT_OG} -o ${OUT_PNG} -s '#' -P 
		
	fi
done < ${INFO_DIR}/odgi_ranges.txt

# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------

rsync -avuP ./*.png ${OUTDIR}
cd ${WORKDIR}
rm -rf ${TMP_DIR}

mamba deactivate
