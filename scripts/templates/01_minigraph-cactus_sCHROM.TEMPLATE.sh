#! /bin/bash
#$ -q "QUEUE01"
#$ -pe smp MCNTHREADS
#$ -cwd
#$ -V
#$ -N cactus_sCHROM
#$ -m ae
#$ -M aharder@hudsonalpha.org


source ~/.bashrc
source /home/aharder_scratch_f13/pangenome_workflow/scripts/99_init_script_vars.sh

## create temp directory and echo path to .o file
TMP_DIR=`/bin/mktemp -d -p /mnt/data1/tmp`
cd ${TMP_DIR}
echo "${TMP_DIR}"


# -----------------------------------------------------------------------------
# Declare variables and prep files
# -----------------------------------------------------------------------------
RUN_ID="${TAXON}_sCHROM"

readonly JOBSTORE_IMAGE=jobStore_${RUN_ID}.img
readonly CACTUS_SCRATCH=${TMP_DIR}/cactus-${RUN_ID}
OUTDIR="${MC_OUT_DIR}/${RUN_ID}"
mkdir -p ${OUTDIR}

LOGFILE="${LOG_DIR}/01_minigraph-cactus_sCHROM_$(date +"%Y_%m_%d_%I_%M_%p").log"
touch ${LOGFILE}

echo "sCHROM - prepping for construction - " $(date -u) >> ${LOGFILE}

## not using the restart option unless really need to -- would work with tmp system
# restart=''

# creating the local working directories
mkdir -p -m 777 ${CACTUS_SCRATCH}/upper ${CACTUS_SCRATCH}/work
truncate -s 300M "${JOBSTORE_IMAGE}"
apptainer exec ${CACTUS_IMAGE} mkfs.ext3 -d ${CACTUS_SCRATCH} "${JOBSTORE_IMAGE}"

mkdir -m 700 -p ${CACTUS_SCRATCH}/tmp
mkdir cactus_wd

# copy input files or cactus can't find them
cp ${SEQFILE} ./sCHROM_seqfile.txt
sed -i 's/xCHROMx/sCHROM/g' sCHROM_seqfile.txt

cp ${CPFILE} ./sCHROM_file_copy.txt
sed -i 's/xCHROMx/sCHROM/g' sCHROM_file_copy.txt

while read -a line
do
	cp ${line[0]} .
done < sCHROM_file_copy.txt

# -----------------------------------------------------------------------------
# Run Minigraph-Cactus pangenome
# -----------------------------------------------------------------------------

echo "sCHROM - running minigraph-cactus graph construction - " $(date -u) >> ${LOGFILE}
apptainer exec --cleanenv \
  --overlay ${JOBSTORE_IMAGE} \
  --bind ${CACTUS_SCRATCH}/tmp:/tmp \
  --env PYTHONNOUSERSITE=1 \
  ${CACTUS_IMAGE} \
  cactus-pangenome \
  	--workDir=cactus_wd \
  	--binariesMode local \
  	--reference ${REFSAMP} \
  	--chop \
  	--vcf clip \
  	--giraffe clip \
  	--gfa clip \
  	--gbz full clip \
  	--xg clip \
  	--chrom-vg clip \
  	--chrom-og full clip \
  	--outDir ${RUN_ID} \
  	--outName ${RUN_ID} \
  	js \
  	sCHROM_seqfile.txt

echo "sCHROM - graph construction complete - " $(date -u) >> ${LOGFILE}

# -----------------------------------------------------------------------------
# Viz with odgi
# -----------------------------------------------------------------------------

mamba activate ${MAMBA}/pg_tools

cut -f 1 sCHROM_seqfile.txt > samp.list

echo "sCHROM - starting odgi - " $(date -u) >> ${LOGFILE}
OG=$(find . -name "*.og")
for OGFN in ${OG}
do

	base=$(basename ${OGFN} .og)
	SORTOG="${base}.sorted.og"
	PNG="${TAXON}_${base}.sorted.png"
	INV_PNG="${TAXON}_${base}.sorted.inversions.png"
	COL_PNG="${TAXON}_${base}.sorted.collapsed.png"

	odgi sort \
		-i ${OGFN} \
		--threads ${MC_NTHREADS} \
		-P \
		-O \
		-Y \
		-o ${SORTOG}

	odgi viz \
		-i ${SORTOG} \
		--threads ${MC_NTHREADS} \
		-o ${PNG}
		
	odgi viz \
		-i ${SORTOG} \
		--threads ${MC_NTHREADS} \
		-o ${INV_PNG} \
		-z
	
	odgi viz \
		-i ${SORTOG} \
		--threads ${MC_NTHREADS} \
		-o ${COL_PNG} \
		-M samp.list

done

echo "sCHROM - odgi complete - " $(date -u) >> ${LOGFILE}
mamba deactivate

# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------
rm *.list
rm *.fa
rm *.fa.gz
rm *.img
rm sCHROM_seqfile.txt
rm sCHROM_file_copy.txt
rsync -avuP *.png ${ODGI_OUT}
rm *.png
rsync -avuP $TMP_DIR/* ${OUTDIR}

cd ${OUTDIR}
rm -rf ${TMP_DIR}