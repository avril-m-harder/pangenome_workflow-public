#! /bin/bash
#$ -q "QUEUE05"
#$ -pe smp VARSCANNTHREADS
#$ -cwd
#$ -V
#$ -N varscan_05_SAMP_NAME
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

mamba activate ${MAMBA}/bioinfo_tools

INDIR="${GIR_OUT_DIR}"
OUTDIR="${VARSCAN_OUT_DIR}"

REF_PATH_FA_BGZIP="mc_${TAXON}_all_chroms_${PRIMREF}_path.fa.gz"
SAMP="SAMP_NAME_k${KLEN}"
SORTBAM="${SAMP}_subgraph.sorted.bam"
MPILEUP="${SAMP}_subgraph.giraffeBAM.mpileup"
TMP_VCF="${SAMP}_subgraph.giraffeBAM.tmp.vcf.gz"
FINAL_VCF="${SAMP}_subgraph.giraffeBAM.vcf.gz"
D8_VCF="${SAMP}_subgraph.giraffeBAM.d8.vcf.gz"

LOGFILE="${LOG_DIR}/05_varscan_${SAMP}_$(date +"%Y_%m_%d_%I_%M_%p").log"
touch ${LOGFILE}

rsync -avuP ${MC_OUT_DIR}/${REF_PATH_FA_BGZIP} .
rsync -avuP ${INDIR}/${SORTBAM} .

# -----------------------------------------------------------------------------
# Prep reference and header info
# -----------------------------------------------------------------------------

echo "SAMP_NAME - prepping reference and VCF header info - " $(date -u) >> ${LOGFILE}

samtools faidx ${REF_PATH_FA_BGZIP}

cat ${REF_PATH_FA_BGZIP}.fai | awk '{print "##contig=<ID="$1",length="$2">"}' > header_mid.tmp

cat ${INFO_DIR}/header_top.tmp header_mid.tmp ${INFO_DIR}/header_bottom.tmp > header

refbase=$(basename ${REF_PATH_FA_BGZIP})
sed -i "s/sREF/${refbase}/g" header

sed "s/ssample/${SAMP}/g" header | \
bgzip -c -@ ${VARSCAN_NTHREADS} > ${SAMP}.header.gz


# -----------------------------------------------------------------------------
# Generate pileup and call variants with varscan
# -----------------------------------------------------------------------------

echo "SAMP_NAME - generating mpileup - " $(date -u) >> ${LOGFILE}
samtools mpileup \
	-d 500 \
	-Q 20 \
	-f ${REF_PATH_FA_BGZIP} \
	-o ${MPILEUP} \
	${SORTBAM}

echo "SAMP_NAME - calling variants - " $(date -u) >> ${LOGFILE}
varscan mpileup2cns \
	${MPILEUP} \
	--min-coverage 1 \
	--min-reads2 0 \
	--output-vcf 1 | \
	grep -v "^#" | \
	awk 'FS=OFS="\t" {split($10,s,":"); if(s[1] != "./.") print $1,$2,$3,$4,$5,$6,$7,".","GT:RD:AD", s[1]":"s[5]":"s[6]}' | \
	bgzip -c -@ ${VARSCAN_NTHREADS} > ${TMP_VCF}

cat ${SAMP}.header.gz ${TMP_VCF} > ${FINAL_VCF}

## filter to keep sites with read depth >= 8
gzip -dck ${FINAL_VCF} | \
awk '/^#/ {print} OFS="\t" {split($10,a,":"); split(a[2],b,","); $10 = a[1]; $9="GT"; if(b[1] + b[2] >= 8) print $0}' | \
bgzip -c ${VARSCAN_NTHREADS} > ${D8_VCF}


# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------

echo "SAMP_NAME - complete - " $(date -u) >> ${LOGFILE}

rsync -avuP ${FINAL_VCF} ${OUTDIR}
rsync -avuP ${D8_VCF} ${OUTDIR}
cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
