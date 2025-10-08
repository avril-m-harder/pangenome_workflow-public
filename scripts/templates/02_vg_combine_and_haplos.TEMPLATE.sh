#! /bin/bash
#$ -q "QUEUE02"
#$ -pe smp VGHAPNTHREADS
#$ -cwd
#$ -V
#$ -N vg_combine_02
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

INDIR="${MC_OUT_DIR}"
OUTDIR="${GIR_OUT_DIR}"
LOGFILE="${LOG_DIR}/02_haplos_$(TZ=${ZONE} date +"%Y_%m_%d_%I_%M_%p").log"

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
HAPL="${GFA}.gbz.k${KLEN}.hapl"
REF_PATH_FA_BGZIP="$(basename ${GFA} .gfa)_${PRIMREF}_path.fa.gz"


mamba activate ${MAMBA}/pg_tools
touch ${LOGFILE}
echo "activating pg_tools environment - " $(TZ=${ZONE} date) >> ${LOGFILE}
mamba list >> ${LOGFILE}

# -----------------------------------------------------------------------------
# Collect chrom-level VG files and merge
# -----------------------------------------------------------------------------

if [ ! -f ${INDIR}/${GFA} ]; then

	echo "merging chrom-level GFAs - " $(TZ=${ZONE} date) >> ${LOGFILE}

	for CHROM in ${chroms[@]}
	do
		
		rsync -avuP $(find ${INDIR} -name "${TAXON}_${CHROM}.gfa.gz") .
		gunzip ${TAXON}_${CHROM}.gfa.gz

	done

	ALLFNS=$( echo $(for CHROM in ${chroms[@]}; do echo "${TAXON}_${CHROM}.gfa"; done) )

	apptainer exec ${VG_IMAGE} \
	vg combine \
		${ALLFNS} > \
		${GFA}

	rsync -avuP ${GFA} ${INDIR}
	
	rm ${TAXON}*.gfa
	
else

	rsync -avuP ${INDIR}/${GFA} .

fi


# -----------------------------------------------------------------------------
# Build GBZ from GFA and extract reference path for downstream apps
# 
# https://github.com/vgteam/vg/wiki/Haplotype-Sampling
# -----------------------------------------------------------------------------

echo "building GBZ from GFA - " $(TZ=${ZONE} date) >> ${LOGFILE}
if ! [[ -f ${INDIR}/${GBZ} ]] || ! [[ -f ${INDIR}/${XG} ]]; then


		## writes GBZ from GFA
		if [ ! -f ${INDIR}/${GBZ} ]; then
		apptainer exec ${VG_IMAGE} \
		vg gbwt \
			--gbz-format \
			--num-threads ${VG_HAP_NTHREADS} \
			-g ${GBZ} \
			-G ${GFA}
		fi
		
		## writes path metadata to text file
		if [ ! -f ${GRAPH_STATS_DIR}/${GBZ}.pathstats.txt ]; then
		apptainer exec ${VG_IMAGE} \
		vg paths \
			-x ${GBZ} \
			--threads ${VG_HAP_NTHREADS} \
			-M > \
			${GRAPH_STATS_DIR}/${GBZ}.pathstats.txt
		fi

		## writes reference path as FASTA
		if [ ! -f ${INDIR}/${REF_PATH_FA_BGZIP} ]; then
		apptainer exec ${VG_IMAGE} \
		vg paths \
			-x ${GBZ} \
			--extract-fasta \
			--threads ${VG_HAP_NTHREADS} \
			--sample ${PRIMREF} | \
			bgzip > \
			${REF_PATH_FA_BGZIP}
		fi

		## writes XG from GFA
		if [ ! -f ${INDIR}/${XG} ]; then
		apptainer exec ${VG_IMAGE} \
		vg convert \
			-g ${GFA} \
			--threads ${VG_HAP_NTHREADS} \
			-x > \
			${XG}
		fi
		
		rsync -avuP ${XG} ${INDIR}	
		rsync -avuP ${XG} ${STMVIZ_DIR}	
		rsync -avuP ${GBZ}* ${INDIR}
		rsync -avuP ${REF_PATH_FA_BGZIP} ${INDIR}

else
	
	rsync -avuP ${INDIR}/${GBZ}* .

fi

# -----------------------------------------------------------------------------
# Generate indices and haplotype file for full graph
# -----------------------------------------------------------------------------

if [ ! -f ${INDIR}/${HAPL} ]; then

	echo "${line[0]} - building giraffe indices - " $(TZ=${ZONE} date) >> ${LOGFILE}
	
	## writes distance index needed for writing haplotype index
	if [ ! -f ${INDIR}/${DIST} ]; then				
	apptainer exec ${VG_IMAGE} \
	vg index \
		-t ${VG_HAP_NTHREADS} \
		--dist-name ${DIST} \
		--no-nested-distance \
		${GBZ}
	fi

	## writes r-index needed for writing haplotype index
	if [ ! -f ${INDIR}/${RI} ]; then				
	apptainer exec ${VG_IMAGE} \
	vg gbwt \
		--num-threads ${VG_HAP_NTHREADS} \
		--gbz-input ${GBZ} \
		--r-index ${RI}
	fi
		
	## writes graph haplotype info for subgraph construction
	apptainer exec ${VG_IMAGE} \
	vg haplotypes \
		--distance-index ${DIST} \
		--r-index ${RI} \
		--threads ${VG_HAP_NTHREADS} \
		--include-reference \
		--haplotype-output ${HAPL} \
		--kmer-length ${KLEN} \
		--subchain-length ${CHAINLEN} \
		--verbosity 2 \
		${GBZ}
	
	rsync -avuP ${GBZ}* ${INDIR}

fi

echo "${line[0]} - complete - " $(TZ=${ZONE} date) >> ${LOGFILE}

# -----------------------------------------------------------------------------
# Collecting summary stats + making some plots
# -----------------------------------------------------------------------------

echo "running summary scripts - " $(TZ=${ZONE} date) >> ${LOGFILE}

## Coverage of each sample X each chromosome in the graph
cat $(find ${INDIR} -name "minigraph.split.log") > \
${GRAPH_STATS_DIR}/minigraph_split_logs.txt

## horizontal bar/line plots of proportion retained per sample, 1 plot per chromosome
Rscript --vanilla \
	${SCRIPT_DIR}/R01_chrom_coverage_stats.R \
	${GRAPH_STATS_DIR}/minigraph_split_logs.txt \
	${GRAPH_STATS_PLOTS}/${TAXON}_chrom_coverage_props.pdf \
	${TAXON}
	
## dot plots of proportion retained by chromosome and by sample (2 plots)
## get coordinates of sequences kept as subpaths for each input sample (haplotype)
grep -e '^W' ${GFA} | \
cut -f2-6 | \
awk '{ print $1 "\t" $3 "\t" $4 "\t" $5 }' | \
awk '{ print $2 "\t" $3 "\t" $4 >$1"_subpath_included_coords.bed" }'
for SAMP in $(cut -f1 ${INFO_DIR}/cactus_CHROM_seqfile.txt)
do
	SAMP1=$(sed 's/[^[:alnum:]]/*/g' <<< "$SAMP")
	if [ -f ${REF_DIR}/*${SAMP1}*chr_only*.fai ]; then
		FAI=$(ls ${REF_DIR}/*${SAMP}*chr_only*.fai)
	else
		FAI=$(ls ${REF_DIR}/*${SAMP1}*.fai)
	fi
	cut -f1,2 ${FAI} > ${SAMP}.genome
done
Rscript --vanilla \
	${SCRIPT_DIR}/R03_plotting_clipped_by_chrom_and_samp.R \
	${TAXON}
	
## Panacus plots: can be customized by editing ${INFO_DIR}/panacus_report.yaml
BASE=$(basename ${GFA} .gfa)

panacus report \
	${INFO_DIR}/panacus_report.yaml > \
	${GRAPH_STATS_PLOTS}/${BASE}_panacus_report.html

## Plotting coverage of samples x chromosomes by representation in the 
## final (clipped) graph - chromosome diagram format

grep -e '^W' ${GFA} | \
cut -f2-6 | \
awk '{ print $1 "#" $2 "#" $3 ":" $4 "-" $5 }' > \
${BASE}.paths.txt

Rscript --vanilla \
	${SCRIPT_DIR}/R02_chrom_coverage_graphics.R \
	./${BASE}.paths.txt \
	${GRAPH_STATS_PLOTS}/${BASE}_chrom_coverage_graphic.pdf

echo "complete - " $(TZ=${ZONE} date) >> ${LOGFILE}

# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------

cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
