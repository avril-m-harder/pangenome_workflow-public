#! /bin/bash
#$ -q "QUEUE02"
#$ -pe smp VGHAPNTHREADS
#$ -cwd
#$ -V
#$ -N vg_combine_02
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

INDIR="${MC_OUT_DIR}"
OUTDIR="${GIR_OUT_DIR}"
LOGFILE="${LOG_DIR}/02_haplos_$(date +"%Y_%m_%d_%I_%M_%p").log"

if [ ${#REFARR[@]} -eq 1 ]; then
	GFA="mc_${TAXON}_all_chroms.gfa"
elif [ ${#REFARR[@]} -gt 1 ]; then
	GFA="mc_${TAXON}_all_chroms_ALL_PATHS_AS_REFS.gfa"
fi

if [ ${#REFARR[@]} -eq 1 ]; then
	GBZ="mc_${TAXON}_all_chroms.gfa.gbz"
elif [ ${#REFARR[@]} -gt 1 ]; then
	TMP_GBZ="mc_${TAXON}_all_chroms_ALL_PATHS_AS_REFS.gfa.gbz"
	GBZ="mc_${TAXON}_all_chroms.gfa.gbz"
fi

DIST="mc_${TAXON}_all_chroms.gfa.gbz.dist"
RI="mc_${TAXON}_all_chroms.gfa.gbz.ri"
MIN="mc_${TAXON}_all_chroms.gfa.gbz.min"
HAPL="mc_${TAXON}_all_chroms.gfa.gbz.k${KLEN}.hapl"
REF_PATH_FA_BGZIP="mc_${TAXON}_all_chroms_${PRIMREF}_path.fa.gz"


mamba activate ${MAMBA}/pg_tools
touch ${LOGFILE}

# -----------------------------------------------------------------------------
# Collect chrom-level VG files and merge
# -----------------------------------------------------------------------------

if [ ! -f ${INDIR}/${GFA} ]; then

	echo "merging chrom-level GFAs - " $(date -u) >> ${LOGFILE}

	for CHROM in ${chroms[@]}
	do
		
		rsync -avuP $(find ${INDIR} -name "${TAXON}_${CHROM}.gfa.gz") .
		gunzip ${TAXON}_${CHROM}.gfa.gz

	done

	ALLFNS=$( echo $(for CHROM in ${chroms[@]}; do echo "${TAXON}_${CHROM}.gfa"; done) )

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
# https://github.com/vgteam/vg/issues/3303
# -----------------------------------------------------------------------------

echo "building GBZ from GFA - " $(date -u) >> ${LOGFILE}
if [ ! -f ${INDIR}/${GBZ} ]; then

	if [ ${#REFARR[@]} -eq 1 ]; then

		vg gbwt \
			--gbz-format \
			-g ${GBZ} \
			-G ${GFA}

		vg paths \
			-x ${GBZ} \
			-M > \
			${GRAPH_STATS_DIR}/${GBZ}.pathstats.txt

		vg paths \
			-x ${GBZ} \
			--extract-fasta \
			--sample ${PRIMREF} |
			bgzip > \
			${REF_PATH_FA_BGZIP}


		${STM_DIR}/prepare_vg.sh ${GBZ}
		rsync -avuP ${GBZ}* ${INDIR}
		rsync -avuP ${GBZ}* ${STMVIZ_DIR}
		rsync -avuP ${REF_PATH_FA} ${INDIR}

	elif [ ${#REFARR[@]} -gt 1 ]; then

		vg gbwt \
			--gbz-format \
			-g ${TMP_GBZ} \
			-G ${GFA}

		vg gbwt \
			-Z \
			--set-reference ${PRIMREF} \
			--gbz-format \
			-g ${GBZ} \
			${TMP_GBZ}

		vg paths \
			-x ${GBZ} \
			-M > \
			${GRAPH_STATS_DIR}/${GBZ}.pathstats.txt

		${STM_DIR}/prepare_vg.sh ${GBZ}
		${STM_DIR}/prepare_vg.sh ${TMP_GBZ}
		rsync -avuP ${GBZ}* ${INDIR}
		rsync -avuP ${GBZ}* ${STMVIZ_DIR}
		rsync -avuP ${TMP_GBZ}* ${INDIR}
		rsync -avuP ${TMP_GBZ}* ${STMVIZ_DIR}
		rsync -avuP ${REF_PATH_FA_BGZIP} ${INDIR}
	
	fi

else
	
	rsync -avuP ${INDIR}/${GBZ} .

fi

# -----------------------------------------------------------------------------
# Generate indices and haplotype file for full graph
# -----------------------------------------------------------------------------

if [ ! -f ${INDIR}/${HAPL} ]; then

	echo "${line[0]} - building giraffe indices - " $(date -u) >> ${LOGFILE}
	vg index \
		--dist-name ${DIST} \
		${GBZ}

	vg gbwt \
		--gbz-input ${GBZ} \
		--r-index ${RI}

	vg minimizer \
		-d ${DIST} \
		-t ${VG_HAP_NTHREADS} \
		-o ${MIN} \
		${GBZ}
		
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

echo "${line[0]} - complete - " $(date -u) >> ${LOGFILE}

# -----------------------------------------------------------------------------
# Collecting summary stats
# -----------------------------------------------------------------------------

echo "running summary scripts - " $(date -u) >> ${LOGFILE}

## Coverage of each sample X each chromosome in the graph
cat $(find ${INDIR} -name "minigraph.split.log") > \
${GRAPH_STATS_DIR}/minigraph_split_logs.txt

Rscript --vanilla \
	${SCRIPT_DIR}/R01_chrom_coverage_stats.R \
	${GRAPH_STATS_DIR}/minigraph_split_logs.txt \
	${GRAPH_STATS_PLOTS}/${TAXON}_chrom_coverage_plots.pdf
	
## Panacus plots: node coverage and graph growth with added samples
base=$(basename ${GFA} .gfa)

grep -e '^W' ${GFA} | \
cut -f2-6 | \
awk '{ print $1 "#" $2 "#" $3 ":" $4 "-" $5 }' > \
${base}.paths.txt

# grep -ve '${REFSAMP}' ${base}.paths.txt > \
# ${base}.paths.haplotypes.txt
cp ${base}.paths.txt ${base}.paths.haplotypes.txt

if [ ! -f ${GRAPH_STATS_PLOTS}/panacus_${base}.histgrowth.node.pdf ]; then
	RUST_LOG=info panacus histgrowth \
		-t${VG_HAP_NTHREADS} \
		-l 1,2,1,1,1 \
		-q 0,0,1,0.5,0.1 \
		-S \
		-a \
		-s ${base}.paths.haplotypes.txt \
		-c all \
		${GFA} > \
		${base}.histgrowth.node.tsv

	panacus-visualize \
		-e ${base}.histgrowth.node.tsv \
		> ${GRAPH_STATS_PLOTS}/panacus_${base}.histgrowth.node.pdf
fi

if [ ! ${GRAPH_STATS_PLOTS}/panacus_${base}.histgrowth.html ]; then
	RUST_LOG=info panacus histgrowth \
		-t${VG_HAP_NTHREADS} \
		-l 1,2,1,1,1 \
		-q 0,0,1,0.5,0.1 \
		-S \
		-s ${base}.paths.haplotypes.txt \
		-c all \
		-a \
		-o html \
		${GFA} > \
		${GRAPH_STATS_PLOTS}/panacus_${base}.histgrowth.html
fi

## Plotting coverage of samples x chromosomes by representation in the 
## final (clipped) graph
Rscript --vanilla \
	${SCRIPT_DIR}/R02_chrom_coverage_graphics.R \
	./${base}.paths.haplotypes.txt \
	${REF_DIR} \
	${GRAPH_STATS_PLOTS}/${base}_chrom_coverage_graphic.pdf

echo "complete - " $(date -u) >> ${LOGFILE}

# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------

cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
