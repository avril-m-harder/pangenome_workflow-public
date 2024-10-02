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

INDIR="${MC_OUT_DIR}"
OUTDIR="${GIR_OUT_DIR}"

GFA="mc_${TAXON}_all_chroms.gfa"
GBZ="mc_${TAXON}_all_chroms.gfa.gbz"
DIST="mc_${TAXON}_all_chroms.gfa.gbz.dist"
RI="mc_${TAXON}_all_chroms.gfa.gbz.ri"
MIN="mc_${TAXON}_all_chroms.gfa.gbz.min"
HAPL="mc_${TAXON}_all_chroms.gfa.gbz.k${KLEN}.hapl"

rsync -avuP ${INDIR}/${GFA}* .

# -----------------------------------------------------------------------------
# Sample haplotypes based on kmer counts
# -----------------------------------------------------------------------------

SAMP="SAMP_NAME_k${KLEN}"
SUBGBZ="mc_${TAXON}_all_chroms.gfa.${SAMP}.subgraph.gbz"
GAM="${SAMP}_subgraph.gam"
BAM="${SAMP}_subgraph.bam"
SORTGAM="${SAMP}_subgraph.sorted.gam"
GAMINDEX="${SAMP}_subgraph.sorted.gam.gai"
SORTBAM="${SAMP}_subgraph.sorted.bam"
DEDUPBAM="${SAMP}_subgraph.sorted.dedup.bam"
STATS="dup_metrics_${SAMP}_subgraph_giraffeBAM.txt"
KFF="${SAMP}.kff"

LOGFILE="${LOG_DIR}/03_giraffe_${SAMP}_$(date +"%Y_%m_%d_%I_%M_%p").log"
touch ${LOGFILE}

## count k-mers of various lengths in sample reads
if [ ! -f ${KMC_OUT}/${KFF} ]; then

	rsync -LvP FQ_FILE .

	FQ=$(ls ./*fastq*)

	if [[ ${FQ} == *.bz2 ]]; then
		bzip2 -d ${FQ}
		FQ=$(basename ${FQ} .bz2)
	fi

	## if using hammer, ~7.8 Gb memory/CPU? 192 CPU, 1.5 Tb
	echo "SAMP_NAME - counting ${KLEN}-mers - " $(TZ=${ZONE} date) >> ${LOGFILE}
	kmc \
		-k${KLEN} \
		-m200 \
		-okff \
		-t${VG_GIR_NTHREADS} \
		${FQ} \
		${SAMP} \
		.

	rsync -avuP ${KFF} ${KMC_OUT}
	
else
	
	rsync -avuP ${KMC_OUT}/${KFF} .

fi


## sample haplotypes to build subgraph (include ref path only if there's only 1)
echo "SAMP_NAME - haplotype sampling k=${KLEN} to build subgraph - " $(TZ=${ZONE} date) >> ${LOGFILE}

if [ ! -f ${OUTDIR}/${SUBGBZ} ]; then

	if [ ${#REFARR[@]} -eq 1 ]; then

		vg haplotypes \
			--threads ${VG_GIR_NTHREADS} \
			--include-reference \
			--haplotype-input ${HAPL} \
			--kmer-input ${KFF} \
			--kmer-length ${KLEN} \
			--num-haplotypes ${NHAPLOS} \
			--verbosity 2 \
			--subchain-length ${CHAINLEN} \
			-g ${SUBGBZ} \
			${GBZ}
			
	elif [ ${#REFARR[@]} -gt 1 ]; then
	
		vg haplotypes \
			--threads ${VG_GIR_NTHREADS} \
			--haplotype-input ${HAPL} \
			--kmer-input ${KFF} \
			--kmer-length ${KLEN} \
			--num-haplotypes ${NHAPLOS} \
			--verbosity 2 \
			--subchain-length ${CHAINLEN} \
			-g ${SUBGBZ} \
			${GBZ}
			
	fi
		
	echo "SAMP_NAME - preparing sTM index for subgraph GBZ - " $(TZ=${ZONE} date) >> ${LOGFILE}
	${STM_DIR}/prepare_vg.sh ${SUBGBZ}
	rsync -avuP ${SUBGBZ}* ${STMVIZ_DIR}
	rsync -avuP ${SUBGBZ}* ${OUTDIR}
else
	rsync -avuP ${OUTDIR}/${SUBGBZ}* .
fi



# -----------------------------------------------------------------------------
# Map short reads to haplotype-sampled graph file - for SV calling
# ----------------------------------------------------------------------------- 

echo "SAMP_NAME - mapping to subgraph k=${KLEN} -> GAM - " $(TZ=${ZONE} date) >> ${LOGFILE}

if [ ! -f ${OUTDIR}/${SORTGAM} ]; then

	if [ -z "${FQ}" ]; then
	
		rsync -LvP FQ_FILE .

		FQ=$(ls ./*fastq*)

		if [[ ${FQ} == *.bz2 ]]; then
			bzip2 -d ${FQ}
			FQ=$(basename ${FQ} .bz2)
		fi

	fi
	
	vg giraffe \
		-Z ${SUBGBZ} \
		--fastq-in ${FQ} \
		--interleaved \
		--max-multimaps 1 \
		-o gam \
		--sample SAMP_NAME_k${KLEN} \
		--threads ${VG_GIR_NTHREADS} \
		--progress > \
		${GAM}

	echo "SAMP_NAME - preparing sTM index for subgraph GAM - " $(TZ=${ZONE} date) >> ${LOGFILE}
	
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
		${SUBGBZ} > \
		${GIR_STATS_DIR}/${SORTGAM}_stats.txt
		
else

	rsync -avuP ${OUTDIR}/${SORTGAM} .

fi


# -----------------------------------------------------------------------------
# Map short reads to haplotype-sampled graph file --> BAM - for SNP calling
# ----------------------------------------------------------------------------- 

echo "SAMP_NAME - mapping to subgraph k=${KLEN} -> BAM - " $(TZ=${ZONE} date) >> ${LOGFILE}


if [ ! -f ${OUTDIR}/${SORTBAM} ]; then

	if [ -z "${FQ}" ]; then
	
		rsync -LvP FQ_FILE .

		FQ=$(ls ./*fastq*)

		if [[ ${FQ} == *.bz2 ]]; then
			bzip2 -d ${FQ}
			FQ=$(basename ${FQ} .bz2)
		fi

	fi

	vg giraffe \
		-Z ${SUBGBZ} \
		--fastq-in ${FQ} \
		--interleaved \
		--max-multimaps 1 \
		-o BAM \
		--sample SAMP_NAME_k${KLEN} \
		--threads ${VG_GIR_NTHREADS} \
		--progress > \
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

mamba deactivate

	# -------------------------------------------------------------------------
	# Dedup BAM if it'll go
	# -------------------------------------------------------------------------

	echo "SAMP_NAME - dedup-ing BAM - " $(TZ=${ZONE} date) >> ${LOGFILE}
	
	mamba activate ${MAMBA}/bioinfo_tools
	
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

rsync -avuP *subgraph* ${OUTDIR}
cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
