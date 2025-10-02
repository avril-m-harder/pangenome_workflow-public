#! /bin/bash
#$ -q "QUEUE03"
#$ -pe smp VGGIRNTHREADS
#$ -cwd
#$ -V
#$ -N vg_giraffe_03_SAMP_NAME
#$ -m ae
#$ -M [email]

source ~/.bashrc
source /home/aharder_scratch_f13/pangenome_workflow/scripts/99_init_script_vars.sh

## create temp directory and echo path to .o file
TMP_DIR=`/bin/mktemp -d -p /mnt/data1/tmp`
cd ${TMP_DIR}
echo "${TMP_DIR}"

# -----------------------------------------------------------------------------
# Activate env if necessary, set variable names
# -----------------------------------------------------------------------------
mamba activate ${MAMBA}/pg_tools
echo "activating pg_tools environment - " $(TZ=${ZONE} date) >> ${LOGFILE}
mamba list >> ${LOGFILE}

INDIR="${MC_OUT_DIR}"
OUTDIR="${GIR_OUT_DIR}"

if [ ${#REFARR[@]} -eq 1 ]; then
	GFA="mc_${TAXON}_all_chroms.gfa"
elif [ ${#REFARR[@]} -gt 1 ]; then
	GFA="mc_${TAXON}_all_chroms_ALL_PATHS_AS_REFS.gfa"
fi

XG="${GFA}.xg"
GBZ="${GFA}.gbz"
DIST="${GFA}.gbz.dist"
RI="${GFA}.gbz.ri"
MIN="${GFA}.gbz.min"

rsync -avuP ${INDIR}/${GFA}* .

SAMP="SAMP_NAME"
GAM="${SAMP}_allhapsgraph.gam"
BAM="${SAMP}_allhapsgraph.bam"
SORTGAM="${SAMP}_allhapsgraph.sorted.gam"
SORTBAM="${SAMP}_allhapsgraph.sorted.bam"
DEDUPBAM="${SAMP}_allhapsgraph.sorted.dedup.bam"
STATS="dup_metrics_${SAMP}_allhapsgraph_giraffeBAM.txt"

LOGFILE="${LOG_DIR}/03_giraffe_allhapsgraph_${SAMP}_$(date +"%Y_%m_%d_%I_%M_%p").log"
touch ${LOGFILE}

# -----------------------------------------------------------------------------
# Map short reads to graph with all haplotypes - for SV calling
# ----------------------------------------------------------------------------- 

echo "SAMP_NAME - mapping to allhapsgraph k=${KLEN} -> GAM - " $(TZ=${ZONE} date) >> ${LOGFILE}

if [ ! -f ${OUTDIR}/${SORTGAM} ]; then

	if [ -z "${FQ1}" ]; then

		rsync -LvP FQ_FILE_1 .
		rsync -LvP FQ_FILE_2 .

		FQ1=$(ls ./*R1*fastq*)
		FQ2=$(ls ./*R2*fastq*)

		if [[ ${FQ1} == *.bz2 ]]; then
			bzip2 -d ${FQ1}
			bzip2 -d ${FQ2}
			FQ1=$(basename ${FQ1} .bz2)
			FQ2=$(basename ${FQ2} .bz2)
		fi

	fi

	apptainer exec ${VG_IMAGE} \
	vg giraffe \
		-Z ${GBZ} \
		--fastq-in ${FQ1} \
		--fastq-in ${FQ2} \
		--interleaved \
		--max-multimaps 1 \
		-o gam \
		--sample SAMP_NAME_k${KLEN} \
		--threads ${VG_GIR_NTHREADS} \
		--progress > \
		${GAM}

	echo "SAMP_NAME - preparing sTM index for allhapsgraph GAM - " $(TZ=${ZONE} date) >> ${LOGFILE}

	apptainer exec ${VG_IMAGE} \
	vg gamsort \
		-i ${GAMINDEX} \
		--threads ${VG_GIR_NTHREADS} \
		${GAM} > \
		${SORTGAM}
	
	rsync -avuP ${GAMINDEX} ${STMVIZ_DIR}
	rsync -avuP ${GAMINDEX} ${OUTDIR}
	rsync -avuP ${SORTGAM} ${STMVIZ_DIR}
	rsync -avuP ${GAM} ${OUTDIR}
	rsync -avuP ${SORTGAM} ${OUTDIR}

	vg stats \
		-a ${SORTGAM} \
		${GBZ} > \
		${GIR_STATS_DIR}/${SORTGAM}_stats.txt
		
else

	rsync -avuP ${OUTDIR}/${SORTGAM} .

fi


# -----------------------------------------------------------------------------
# Map short reads to graph with all haplotypes --> BAM - for SNP calling
# ----------------------------------------------------------------------------- 

echo "SAMP_NAME - mapping to allhapsgraph k=${KLEN} -> BAM - " $(TZ=${ZONE} date) >> ${LOGFILE}


if [ ! -f ${OUTDIR}/${SORTBAM} ]; then

## replacing the below section with vg surject - keeping in case that doesn't work
# 	if [ -z "${FQ1}" ]; then
# 
# 		rsync -LvP FQ_FILE_1 .
# 		rsync -LvP FQ_FILE_2 .
# 
# 		FQ1=$(ls ./*R1*fastq*)
# 		FQ2=$(ls ./*R2*fastq*)
# 
# 		if [[ ${FQ1} == *.bz2 ]]; then
# 			bzip2 -d ${FQ1}
# 			bzip2 -d ${FQ2}
# 			FQ1=$(basename ${FQ1} .bz2)
# 			FQ2=$(basename ${FQ2} .bz2)
# 		fi
# 
# 	fi
# 
# 	apptainer exec ${VG_IMAGE} \
# 	vg giraffe \
# 		-Z ${GBZ} \
# 		--fastq-in ${FQ1} \
# 		--fastq-in ${FQ2} \
# 		--max-multimaps 1 \
# 		-o BAM \
# 		--sample SAMP_NAME_k${KLEN} \
# 		--threads ${VG_GIR_NTHREADS} \
# 		--progress > \
# 		${BAM}

	apptainer exec ${VG_IMAGE} \
	vg surject \
		--threads ${VG_GIR_NTHREADS} \
		-x ${XG} \
		--interleaved \
		--progress \
		--bam-output \
		${SORTGAM} > \
		${BAM}
	
	samtools sort \
		--threads $((${VG_GIR_NTHREADS}-1)) \
		${BAM} > \
		${SORTBAM}
	
	samtools index ${SORTBAM}

	rsync -avuP ${SORTBAM}* ${OUTDIR}
	
	samtools flagstat ${SORTBAM} > \
	${BAM_STATS_DIR}/${SORTBAM}_flagstat.txt
else

	rsync -avuP ${OUTDIR}/${SORTBAM} .

fi

	# -------------------------------------------------------------------------
	# Dedup BAM if it'll go
	# -------------------------------------------------------------------------

	echo "SAMP_NAME - dedup-ing BAM - " $(TZ=${ZONE} date) >> ${LOGFILE}
	
	picard MarkDuplicates \
		I=${SORTBAM} \
		O=${DEDUPBAM} \
		M=${STATS} \
		VALIDATION_STRINGENCY=LENIENT \
		ASSUME_SORTED=true \
		SORTING_COLLECTION_SIZE_RATIO=0.05 \
		REMOVE_DUPLICATES=true \
		TMP_DIR=${TMP_DIR} \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		MAX_RECORDS_IN_RAM=2750000 \
		READ_NAME_REGEX='[a-zA-Z0-9]+[:_][0-9]+[:_][a-zA-Z0-9]+[:_][0-9]+[:_]([0-9]+)[:_]([0-9]+)[:_]([0-9]+)' \
		CREATE_INDEX=true
		
	rsync -avuP ${STATS} ${BAM_STATS_DIR}
	rsync -avuP ${DEDUPBAM} ${OUTDIR}
	
	
echo "SAMP_NAME - complete - " $(TZ=${ZONE} date) >> ${LOGFILE}
	

# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------

rsync -avuP *allhapsgraph* ${OUTDIR}
cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
