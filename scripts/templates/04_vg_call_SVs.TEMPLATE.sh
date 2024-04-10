#! /bin/bash
#$ -q "QUEUE04"
#$ -pe smp VGCALLNTHREADS
#$ -cwd
#$ -V
#$ -N vg_call_04_SAMP_NAME
#$ -m ae
#$ -M aharder@hudsonalpha.org

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

INDIR="${GIR_OUT_DIR}"
OUTDIR="${CALL_OUT_DIR}"

SAMP="SAMP_NAME_k${KLEN}"
SUBXG="mc_${TAXON}_all_chroms.gfa.${SAMP}.subgraph.gbz.xg"
SORTGAM="${SAMP}_subgraph.sorted.gam"
FILTGAM="${SAMP}_subgraph.sorted.filtered.gam"
OUTVCF="${SAMP}_mc_subgraph.vcf"
FILTVCF="${SAMP}_mc_subgraph_PASS.vcf"
NORMVCF="${SAMP}_mc_subgraph_PASS.norm.vcf"
SVVCF="${SAMP}_mc_subgraph_PASS.norm.SVsONLY.vcf.gz"

LOGFILE="${LOG_DIR}/04_vg_call_${SAMP}_$(date +"%Y_%m_%d_%I_%M_%p").log"
touch ${LOGFILE}

# -----------------------------------------------------------------------------
# Compute read support (pack) and call variants (call)
# -----------------------------------------------------------------------------

rsync -avuP ${INDIR}/${SUBXG} .
rsync -avuP ${INDIR}/${SORTGAM} .

echo "(1/4) running vg filter - " $(date -u) >> ${LOGFILE}
vg filter \
	${SORTGAM} \
	--threads ${VG_CALL_NTHREADS} \
	-r 0.90 \
	-fu \
	-m 1 \
	-q 15 \
	-D 999 \
	-x ${SUBXG} > \
	${FILTGAM}
	
vg stats \
	-a ${FILTGAM} \
	${SUBGBZ} > \
	${GIR_STATS_DIR}/${FILTGAM}_stats.txt

echo "(2/4) running vg pack - " $(date -u) >> ${LOGFILE}
vg pack \
	-x ${SUBXG} \
	-g ${FILTGAM} \
	-Q 5 \
	--threads ${VG_CALL_NTHREADS} \
	-o tmp.pack

echo "(3/4) running vg call - " $(date -u) >> ${LOGFILE}
## -a = genotype every snarl, including reference calls (gVCF)
vg call \
	${SUBXG} \
	-a \
	--sample ${SAMP} \
	--threads ${VG_CALL_NTHREADS} \
	-k tmp.pack > \
	${OUTVCF}
	
rm tmp.pack

## edit weirdly formatted chrom names output by Minigraph-Cactus in header
sed -i "s/${PRIMREF}#0#//g" ${OUTVCF}

# -----------------------------------------------------------------------------
# Normalize and filter to keep passing variants >= 50 bp in length
# -----------------------------------------------------------------------------

echo "(4/4) running VCF filtering - " $(date -u) >> ${LOGFILE}
bcftools view \
	-f PASS \
	${OUTVCF} > \
	${FILTVCF}
	
bcftools norm \
	-m-any \
	-O z \
	-o ${NORMVCF} \
	${FILTVCF}
	
grep "^#" ${NORMVCF} > tmp1.header1
grep -v "^#CHROM" tmp1.header1 > tmp2.header1
grep "^#CHROM" ${NORMVCF} > tmp.header2
cat tmp2.header1 tmp.header2 > full.header

grep -v "^#" ${NORMVCF} | \
awk '{if (sqrt(((length($4)-length($5))^2))<50) print $0}' > \
tmp.${SAMP}.vcf

cat full.header tmp.${SAMP}.vcf | \
gzip > \
${SVVCF}

rm tmp.${SAMP}.vcf
rm ${FILTVCF}
rm ${NORMVCF}

gzip *.vcf
rsync -avuP *.vcf.gz ${OUTDIR}

echo "SAMP_NAME - complete - " $(date -u) >> ${LOGFILE}
	

# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------

cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
