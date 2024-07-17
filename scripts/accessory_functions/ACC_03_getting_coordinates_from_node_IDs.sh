#! /bin/bash
#$ -q "all.q"
#$ -pe smp 2
#$ -cwd
#$ -V
#$ -N node_coords_ACC03
#$ -m ae
#$ -M aharder@hudsonalpha.org

##—<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>—<>—<>—<>—##
##                                                                                     ##
##  Takes input file (*_breakpoint_nodes.txt, described in first section below) and    ##
##  determines haplotype-specific coordinates for either the first or last position of ##
##  a known node of interest. Coordinates are with respect to the original sequence    ##
##  input to the graph.																   ##
##  																				   ##
##  Output file format is a tsv (*_breakpoint_coords.txt) with the columns:			   ##
##  1)  node annotation = can be whatever desciptive name you set for each node	       ##
##		>>  if pulling sequences (below), include *_start and *_end with match 		   ##
##			prefixes to indicate bounds for ROIs									   ##
##  2)  sample = ... sample/haplotype name                                             ##
##  3)  node ID = node ID in the graph (integer)                                       ##
##  4)  chromosome																	   ##
##  5)  start/end coordinate = coordinate of start or end of node of interest in	   ##
##  	path of interest (start or end is indicated in input file)                     ##
##  6) start coordinate of subpath traversing node of interest						   ##
##  																				   ##
##  Can also auto extract sequences corresponding to ROIs if indicated by PULL_SEQS    ##
##  setting in first section below.													   ##
##                                                                                     ##
##—<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>——<>—<>—<>—<>—##

source ~/.bashrc
source /home/aharder_scratch_f13/pangenome_workflow/scripts/99_init_script_vars.sh

## create temp directory and echo path to .o file
TMP_DIR=`/bin/mktemp -d -p /mnt/data1/tmp`
cd ${TMP_DIR}
echo "${TMP_DIR}"

# -----------------------------------------------------------------------------
# Set which set of breakpoint nodes you want coordinates for
# -----------------------------------------------------------------------------

## node list file: breakpoint ID / sample / node ID / start or end
## 1) breakpoint ID can be any text useful for describing the ROI
## 2) sample can be a single sample name or "all"
##		>> all = all paths traversing the node of interest will be identified and processed
## 3) node ID can be identified by examining a graph VCF or (more easily) checking ROIs in
## 	  sequenceTubeMap and manually identifying nodes of interest
##
## 4) start or end colum:
##		>> 0 = want coordinate of first position in node
##		>> 1 = want coordinate of last position in node

NODELIST="${INFO_DIR}/example_breakpoint_nodes.txt"
base=$(basename ${NODELIST} _breakpoint_nodes.txt)

# -----------------------------------------------------------------------------
# Set variable names, activate env (automated)
# -----------------------------------------------------------------------------

LOGFILE="${LOG_DIR}/ACC_03_${base}_breakpoints_$(date +"%Y_%m_%d_%I_%M_%p").log"
INDIR="${MC_OUT_DIR}"
OUTDIR="${NCOORDS_OUT}"

if [ ${#REFARR[@]} -eq 1 ]; then
	GFA="mc_${TAXON}_all_chroms.gfa"
elif [ ${#REFARR[@]} -gt 1 ]; then
	GFA="mc_${TAXON}_all_chroms_ALL_PATHS_AS_REFS.gfa"
fi

XG="${GFA}".xg

mamba activate ${MAMBA}/pg_tools
touch ${LOGFILE}

rsync -avuP ${INDIR}/${GFA} .
rsync -avuP ${INDIR}/${XG} .

# -----------------------------------------------------------------------------
# Collect breakpoint information (automated)
# -----------------------------------------------------------------------------
## for each path x node of interest combo:
while read -a line
do
	NODEANNOT=${line[0]}
	HAP=${line[1]}
	NODEID=${line[2]}
	S_OR_E=${line[3]}
	
	echo "starting search for node annot = ${NODEANNOT} - " $(date -u) >> ${LOGFILE}
	
	if [ ${HAP} == "all" ]; then
		
		## identify which haplotypes traverse that node
		declare -a HAPS=( $(grep "^W" ${GFA} | grep "[^0-9]${NODEID}[^0-9]" | awk '{print $2}') )
		echo ">>> ${#HAPS[@]} haps found for node ${NODEID} - " $(date -u) >> ${LOGFILE}
		
		for HAP in ${HAPS[@]}
		do
		
		## then find the subpath traversing the node		
		SUBP=$(grep "^W" ${GFA} | grep ${HAP} | grep "[^0-9]${NODEID}[^0-9]" | cut -f2,3,4,5,6)
		echo -e "${SUBP}"'\t'"${NODEID}"'\t'"${NODEANNOT}"'\t'"${S_OR_E}" >> subpaths.txt
		
		done
	
	else
	
		## find the subpath traversing the node
		SUBP=$(grep "^W" ${GFA} | grep ${HAP} | grep "[^0-9]${NODEID}[^0-9]" | cut -f2,3,4,5,6)
		echo -e "${SUBP}"'\t'"${NODEID}"'\t'"${NODEANNOT}"'\t'"${S_OR_E}" >> subpaths.txt	
		
	fi
	
	echo ">>> node annot = ${NODEANNOT} search complete - " $(date -u) >> ${LOGFILE}
	
done < ${NODELIST}

echo "starting node coordinate pull - " $(date -u) >> ${LOGFILE}
while read -a line
do
	SAMP=${line[0]}
	CHROM=${line[2]}
	STARTCOORD=${line[3]}
	NODEID=${line[5]}
	NODEANNOT=${line[6]}
	S_OR_E=${line[7]}

	if [ ${STARTCOORD} == 0 ]; then
		PNAME="${line[0]}#${line[1]}#${line[2]}"
	else
		PNAME="${line[0]}#${line[1]}#${line[2]}[${line[3]}-${line[4]}]"
	fi

	echo ${PNAME}
	echo ${NODEID}

	FINDOUT=$(vg find \
		-x ${XG} \
		-n ${NODEID} \
		-P ${PNAME})
	
	echo -e "${NODEANNOT}"'\t'"${SAMP}"'\t'"${FINDOUT}"'\t'"${CHROM}"'\t'"${STARTCOORD}"'\t'"${S_OR_E}" >> node_coords.txt	
		
done < subpaths.txt
echo ">>> node coordinate pull complete - " $(date -u) >> ${LOGFILE}

echo "starting breakpoint coordinate calculation - " $(date -u) >> ${LOGFILE}
while read -a line
do
	if [[ "${line[5]}" -eq 1 ]]; then
		NODEANNOT=${line[0]}
		SAMP=${line[1]}
		NODEID=${line[2]}
		NODESTART=${line[3]}
		CHROM=${line[4]}
		SUBPATHSTART=${line[5]}
		NODELEN=$(grep "^S" ${GFA} | grep ${NODEID} | awk '{ print length($3)}')
		ENDCOORD=$((${NODELEN} + ${NODESTART} + ${SUBPATHSTART}))
		
		echo -e "${NODEANNOT}"'\t'"${SAMP}"'\t'"${NODEID}"'\t'"${CHROM}"'\t'"${ENDCOORD}"'\t'"${line[5]}" >> ${base}_breakpoint_coords.txt
		
	else
		NODEANNOT=${line[0]}
		SAMP=${line[1]}
		NODEID=${line[2]}
		NODESTART=${line[3]}
		CHROM=${line[4]}
		SUBPATHSTART=${line[5]}
		STARTCOORD=$((${NODESTART} + ${SUBPATHSTART} + 1))
		
		echo -e "${NODEANNOT}"'\t'"${SAMP}"'\t'"${NODEID}"'\t'"${CHROM}"'\t'"${STARTCOORD}"'\t'"${line[5]}" >> ${base}_breakpoint_coords.txt
	
	fi

done < node_coords.txt
echo ">>> breakpoint coordinate calculation complete - " $(date -u) >> ${LOGFILE}

# -----------------------------------------------------------------------------
# Pull sequences if option is set (automated)
# -----------------------------------------------------------------------------

if [ ${PULL_SEQS} -eq 1 ]; then

	OUTFAS=${base}_ROI_sequences.fa

	while read -a line
	do

		if [[ ${line[0]} == *"start"*  ]]; then
			
			START=${line[0]}
			BASE=$(basename ${START} _start)
			END=$(echo "${BASE}_end")
			SAMP=${line[1]}
			CHROM=${line[3]}
			
			## subtract 1 to go from base-1 to base-0 for BED input to seqtk
			S_COORD=$((${line[4]}-1))
			S_SUBPATH=${line[5]}
			
			E_COORD=$(grep ${END} ${base}_breakpoint_coords.txt | \
			grep ${SAMP} | grep ${CHROM} | awk '{print $5}')
			E_SUBPATH=$(grep ${END} ${base}_breakpoint_coords.txt | \
			grep ${SAMP} | grep ${CHROM} | awk '{print $6}')
		
		
			INFAS=$(ls ${REF_DIR}/*${SAMP}*${CHROM}*)
			echo -e "${CHROM}\t${S_COORD}\t${E_COORD}" > tmp.bed
			seqtk subseq ${INFAS} tmp.bed > tmp.fas
		
			if [[ "${S_SUBPATH}" != "${E_SUBPATH}" ]]; then
				sed -i "s/>${CHROM}/>${SAMP}|${CHROM}|${BASE}|CONTAINS_CLIPPED/g" tmp.fas
			else
				sed -i "s/>${CHROM}/>${SAMP}|${CHROM}|${BASE}/g" tmp.fas
			fi
			
			cat tmp.fas >> ${OUTFAS}
		
			rm tmp.bed
			rm tmp.fas

		fi

	done < ${base}_breakpoint_coords.txt
	
	rsync -avuP ${OUTFAS} ${OUTDIR}

fi

# -----------------------------------------------------------------------------
# Clean up tmp dir (automated)
# -----------------------------------------------------------------------------

rsync -avuP ${base}_breakpoint_coords.txt ${OUTDIR}

cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
