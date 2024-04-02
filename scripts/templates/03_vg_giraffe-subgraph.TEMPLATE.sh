#! /bin/bash
#$ -q "QUEUE03"
#$ -pe smp VGGIRNTHREADS
#$ -cwd
#$ -V
#$ -N vg_giraffe_03_SAMP_NAME
#$ -m ae
#$ -M #######@hudsonalpha.org

source ~/.bashrc
source /home/#######_scratch_f13/pangenome_workflow/scripts/99_init_script_vars.sh

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
DIST="mc_${TAXON}_all_chroms.gfa.dist"
RI="mc_${TAXON}_all_chroms.gfa.ri"
MIN="mc_${TAXON}_all_chroms.gfa.min"
HAPL="mc_${TAXON}_all_chroms.gfa.k${KLEN}.hapl"

rsync -avuP ${INDIR}/${GFA}* .

# -----------------------------------------------------------------------------
# Sample haplotypes based on kmer counts
# -----------------------------------------------------------------------------

SAMP="SAMP_NAME_k${KLEN}"
SUBGBZ="mc_${TAXON}_all_chroms.gfa.${SAMP}.subgraph.gbz"
GAM="${SAMP}_subgraph.gam"
BAM="${SAMP}_subgraph.bam"
SORTGAM="${SAMP}_subgraph.sorted.gam"
SORTBAM="${SAMP}_subgraph.sorted.bam"
KFF="${SAMP}.kff"

LOGFILE="${LOG_DIR}/03_giraffe_${SAMP}_$(date +"%Y_%m_%d_%I_%M_%p").log"
touch ${LOGFILE}

## count k-mers of various lengths in sample reads
rsync -avuP ${READ_DIR}/FQ_FILE_1 .
rsync -avuP ${READ_DIR}/FQ_FILE_2 .

echo "FQ_FILE_1" >> files.lst
echo "FQ_FILE_2" >> files.lst


if [ ! -f ${KMC_OUT}/${KFF} ]; then

	## if using hammer, ~7.8 Gb memory/CPU? 192 CPU, 1.5 Tb
	echo "SAMP_NAME - counting ${KLEN}-mers - " $(date -u) >> ${LOGFILE}
	apptainer exec \
		-H $(pwd) \
		${KMCDOCK} \
		kmc \
		-k${KLEN} \
		-m200 \
		-okff \
		-t${VG_GIR_NTHREADS} \
		@files.lst \
		${SAMP} \
		.

	rsync -avuP ${KFF} ${KMC_OUT}
	
else
	
	rsync -avuP ${KMC_OUT}/${KFF} .

fi


## sample haplotypes to build subgraph
echo "SAMP_NAME - haplotype sampling k=${KLEN} to build subgraph - " $(date -u) >> ${LOGFILE}

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
		
	echo "SAMP_NAME - preparing sTM index for subgraph GBZ - " $(date -u) >> ${LOGFILE}
	${STM_DIR}/prepare_vg.sh ${SUBGBZ}
	rsync -avuP ${SUBGBZ}* ${STMVIZ_DIR}
	rsync -avuP ${SUBGBZ}* ${OUTDIR}
else
	rsync -avuP ${OUTDIR}/${SUBGBZ}* .
fi



# -----------------------------------------------------------------------------
# Map short reads to haplotype-sampled graph file - for SV calling
# ----------------------------------------------------------------------------- 

echo "SAMP_NAME - mapping to subgraph k=${KLEN} -> GAM - " $(date -u) >> ${LOGFILE}

if [ ! -f ${OUTDIR}/${SORTGAM} ]; then

	vg giraffe \
		-Z ${SUBGBZ} \
		--fastq-in FQ_FILE_1 \
		--fastq-in FQ_FILE_2 \
		--max-multimaps 1 \
		-o gam \
		--sample SAMP_NAME_k${KLEN} \
		--threads ${VG_GIR_NTHREADS} \
		--progress > \
		${GAM}

	echo "SAMP_NAME - preparing sTM index for subgraph GAM - " $(date -u) >> ${LOGFILE}
	${STM_DIR}/prepare_gam.sh ${GAM}
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

echo "SAMP_NAME - mapping to subgraph k=${KLEN} -> BAM - " $(date -u) >> ${LOGFILE}


if [ ! -f ${OUTDIR}/${SORTBAM} ]; then
	vg giraffe \
		-Z ${SUBGBZ} \
		--fastq-in FQ_FILE_1 \
		--fastq-in FQ_FILE_2 \
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

fi

	
echo "SAMP_NAME - complete - " $(date -u) >> ${LOGFILE}
	

# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------

rsync -avuP *subgraph* ${OUTDIR}
cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
