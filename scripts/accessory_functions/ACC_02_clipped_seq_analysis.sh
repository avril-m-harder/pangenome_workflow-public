#! /bin/bash
#$ -q "all.q"
#$ -pe smp 2
#$ -cwd
#$ -V
#$ -N clipped_ACC02
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
OUTDIR="${CLIP_OUT}"

GFA="mc_${TAXON}_all_chroms.gfa"

mamba activate ${MAMBA}/pg_tools

# -----------------------------------------------------------------------------
# Collecting summary stats
# -----------------------------------------------------------------------------

rsync -avuP ${INDIR}/${GFA} .

rsync -avuP ${INFO_DIR}/clipped_analysis_files.list .

## get coordinates of sequences kept as subpaths for each input sample (haplotype)
grep -e '^W' ${GFA} | \
cut -f2-6 | \
awk '{ print $1 "\t" $3 "\t" $4 "\t" $5 }' | \
awk '{ print $2 "\t" $3 "\t" $4 >$1"_subpath_included_coords.bed" }'

## take complement to find coordinates of excluded (clipped) sequences and intersect excluded sequences with GFF files
while read -a line
do

	SAMP=${line[0]}
	REP_GFF=${line[1]}
	GEN_GFF=${line[2]}
	EXCL_BED="${SAMP}_subpath_excluded_coords.bed"
		
	if [ -f ${REF_DIR}/*${SAMP}*chr_only*.fai ]; then
		FAI=$(ls ${REF_DIR}/*${SAMP}*chr_only*.fai)
	else
		FAI=$(ls ${REF_DIR}/*${SAMP}*.fai)
	fi

	cut -f1,2 ${FAI} > ${SAMP}.genome
	
	bedtools sort \
		-g ${SAMP}.genome \
		-i ${SAMP}_subpath_included_coords.bed > \
		${SAMP}_subpath_included_coords.sorted.bed
	
	bedtools complement \
		-i ${SAMP}_subpath_included_coords.sorted.bed \
		-g ${SAMP}.genome > \
		${EXCL_BED}
		
	rsync -avuP ${REF_DIR}/${REP_GFF} .
	rsync -avuP ${REF_DIR}/${GEN_GFF} .
	
	## get repeat intersects
	BASE=$(basename ${REP_GFF} .gff3.gz)
	REP_GFF_CHR="${BASE}.CHR_ONLY.gff3.gz"

	gzip -cd ${REP_GFF} | \
	grep "^Chr" | gzip -c > \
	${REP_GFF_CHR}

	bedtools intersect \
		-a ${EXCL_BED} \
		-b ${REP_GFF_CHR} \
		-wao > \
		${SAMP}_clips_repeats_intersects.txt
	
	## get gene intersects
	BASE=$(basename ${GEN_GFF} .gff3.gz)
	GEN_GFF_CHR="${BASE}.CHR_ONLY.gff3.gz"

	gzip -cd ${GEN_GFF} | \
	grep "^Chr" | gzip -c > \
	${GEN_GFF_CHR}
		
	bedtools intersect \
		-a ${EXCL_BED} \
		-b ${GEN_GFF_CHR} \
		-wao > \
		${SAMP}_clips_genes_intersects.txt		
	
done < clipped_analysis_files.list

Rscript --vanilla \
	${ACC_FUNC_DIR}/ACC_R02_collapsing_and_plotting_overlaps.R \
	${TAXON}


# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------

mkdir -p ${OUTDIR}/bed_files
rsync -avuP ./*.bed ${OUTDIR}/bed_files
rsync -avuP ./*.pdf ${OUTDIR}
rsync -avuP ./*.pdf ${GRAPH_STATS_PLOTS}
rsync -avuP ./*clipped_sequence_summary.txt ${OUTDIR}

cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
