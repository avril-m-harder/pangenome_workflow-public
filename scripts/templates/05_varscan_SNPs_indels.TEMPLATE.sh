#! /bin/bash
#$ -q "QUEUE05"
#$ -pe smp VARSCANNTHREADS
#$ -cwd
#$ -V
#$ -N varscan_05_SAMP_NAME
#$ -m ae
#$ -M [email]

source ~/.bashrc
source /home/######/pangenome_workflow/scripts/99_init_script_vars.sh

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

INDIR="${GIR_OUT_DIR}"
OUTDIR="${VARSCAN_OUT_DIR}"

if [ ${#REFARR[@]} -eq 1 ]; then
	GFA="mc_${TAXON}_all_chroms.gfa"
elif [ ${#REFARR[@]} -gt 1 ]; then
	GFA="mc_${TAXON}_all_chroms_ALL_PATHS_AS_REFS.gfa"
fi
REF_PATH_FA_BGZIP="$(basename ${GFA} .gfa)_${PRIMREF}_path.fa.gz"

SAMP="SAMP_NAME_k${KLEN}"
SORTBAM="${SAMP}_subgraph.sorted.bam"
DEDUPBAM="${SAMP}_subgraph.sorted.dedup.bam"
MPILEUP="${SAMP}_subgraph.giraffeBAM.mpileup"
TMP_VCF="${SAMP}_subgraph.giraffeBAM.tmp.vcf.gz"
FINAL_VCF="${SAMP}_subgraph.giraffeBAM.vcf.gz"

if [ ! -f ${INDIR}/${DEDUPBAM} ]; then
	rsync -avuP ${INDIR}/${SORTBAM} .
	INBAM=${SORTBAM}
else
	rsync -avuP ${INDIR}/${DEDUPBAM} .
	INBAM=${DEDUPBAM}
fi
FILTBAM=$(basename ${INBAM} .bam).filt.bam

LOGFILE="${LOG_DIR}/05_varscan_${SAMP}_$(date +"%Y_%m_%d_%I_%M_%p").log"
touch ${LOGFILE}

rsync -avuP ${MC_OUT_DIR}/${REF_PATH_FA_BGZIP} .

# -----------------------------------------------------------------------------
# Filter BAM file
# -----------------------------------------------------------------------------

samtools view \
	--min-MQ 15 \
	-@ ${VARSCAN_NTHREADS} \
	-b \
	-o ${FILTBAM} \
	${INBAM}
	
rsync -avuP ${FILTBAM} ${INDIR}


# -----------------------------------------------------------------------------
# Prep reference and header info
# -----------------------------------------------------------------------------

echo "SAMP_NAME - prepping reference and VCF header info - " $(TZ=${ZONE} date) >> ${LOGFILE}

samtools faidx ${REF_PATH_FA_BGZIP}

cat ${REF_PATH_FA_BGZIP}.fai | awk '{print "##contig=<ID="$1",length="$2">"}' > header_mid.tmp

cat ${INFO_DIR}/header_top.tmp header_mid.tmp ${INFO_DIR}/header_bottom.tmp > header

refbase=$(basename ${REF_PATH_FA_BGZIP})
sed -i "s/sREF/${refbase}/g" header

sed "s/ssample/${SAMP}/g" header | \
bgzip -c -@ ${VARSCAN_NTHREADS} > ${SAMP}.header.gz


# -----------------------------------------------------------------------------
# Calculate mean read depth and lower/upper depth limits for keeping variants
# -----------------------------------------------------------------------------

samtools coverage ${FILTBAM} > tmp.coverage.txt

AVG_DP=$(grep -v "^#" tmp.coverage.txt | awk '{ sum1+=($3*$7); } { sum2+=$3; } END{ print sum1/sum2;}')
LO_LIM=$(awk -vdp=$AVG_DP -vlo=$LO_PROP 'BEGIN{printf "%.0f" ,dp * lo}')
if (( $(echo "$LO_LIM $MIN" | awk '{print ($1 < $2)}') )); then
	LO_LIM=$MIN
fi
UP_LIM=$(awk -vdp=$AVG_DP -vlo=$UP_PROP 'BEGIN{printf "%.0f" ,dp * lo}')

DEPTHFILT_VCF="${SAMP}_subgraph.giraffeBAM.d${LO_LIM}-${UP_LIM}.vcf.gz"
D8_VCF="${SAMP}_subgraph.giraffeBAM.d8.vcf.gz"


# -----------------------------------------------------------------------------
# Generate pileup and call variants with varscan
# -----------------------------------------------------------------------------

echo "SAMP_NAME - generating mpileup - " $(TZ=${ZONE} date) >> ${LOGFILE}
samtools mpileup \
	-d 500 \
	-Q 20 \
	-f ${REF_PATH_FA_BGZIP} \
	-o ${MPILEUP} \
	${FILTBAM}

echo "SAMP_NAME - calling variants - " $(TZ=${ZONE} date) >> ${LOGFILE}
varscan mpileup2cns \
	${MPILEUP} \
	--min-coverage 1 \
	--min-reads2 0 \
	--output-vcf 1 | \
	grep -v "^#" | \
	awk 'FS=OFS="\t" {split($10,s,":"); if(s[1] != "./.") print $1,$2,$3,$4,$5,$6,$7,".","GT:RD:AD", s[1]":"s[5]":"s[6]}' | \
	bgzip -c -@ ${VARSCAN_NTHREADS} > ${TMP_VCF}

cat ${SAMP}.header.gz ${TMP_VCF} > ${FINAL_VCF}

## filter to keep sites with acceptable read depths
filt="(FMT/RD + FMT/AD)>=${LO_LIM} & (FMT/RD + FMT/AD)<=${UP_LIM}"
bcftools filter -i "$filt" ${FINAL_VCF} | bgzip -c -@ ${VARSCAN_NTHREADS} > ${DEPTHFILT_VCF}

## filter to keep sites with >= 8 reads
LO_LIM=8
filt="(FMT/RD + FMT/AD)>=${LO_LIM}"
bcftools filter -i "$filt" ${FINAL_VCF} | bgzip -c -@ ${VARSCAN_NTHREADS} > ${D8_VCF}

# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------

echo "SAMP_NAME - complete - " $(TZ=${ZONE} date) >> ${LOGFILE}

rsync -avuP ${FINAL_VCF} ${OUTDIR}
rsync -avuP ${DEPTHFILT_VCF} ${OUTDIR}
rsync -avuP ${D8_VCF} ${OUTDIR}
cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
