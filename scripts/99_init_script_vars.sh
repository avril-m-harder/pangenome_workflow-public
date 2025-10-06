#!/bin/bash
#
# initialize_script_variables.sh
#
# This script initializes variables that are used in multiple scripts in the
# pangenome graph workflow. It should be included in all scripts for each step of
# the workflow. Put the following line in the scripts:
#
#     source /home/######/pangenome_workflow/scripts/99_init_script_vars.sh
#

## Set main project directory
WDNAME="pangenome_workflow"
WORKDIR="/home/######/${WDNAME}"

## Set taxon name (can be anything useful for ID)
TAXON="camelinaTEST"

## Can set time zone for timestamps in log files
ZONE="America/Chicago"

# -----------------------------------------------------------------------------
# Create input and output directories
# -----------------------------------------------------------------------------
INFO_DIR="${WORKDIR}/info"
SCRIPT_DIR="${WORKDIR}/scripts"
TEMPLATE_DIR="${SCRIPT_DIR}/templates"
ACC_FUNC_DIR="${SCRIPT_DIR}/accessory_functions"
REF_DIR="${WORKDIR}/ref_data"
STATS_DIR="${WORKDIR}/stats"

## analysis-specific directories
MC_OUT_DIR="${WORKDIR}/mc_output"
GIR_OUT_DIR="${WORKDIR}/vg_giraffe_output"
CALL_OUT_DIR="${WORKDIR}/vg_call_output"
VARSCAN_OUT_DIR="${WORKDIR}/varscan_output"
GIR_STATS_DIR="${STATS_DIR}/vg_giraffe_gam_stats"
BAM_STATS_DIR="${STATS_DIR}/vg_bam_stats"
GRAPH_STATS_DIR="${STATS_DIR}/graph_stats_info"
GRAPH_STATS_PLOTS="${GRAPH_STATS_DIR}/plots"
DEDUP_STATS_DIR="${STATS_DIR}/dedup_bam_stats"
READ_DIR="${WORKDIR}/prepped_reads"
KMC_OUT="${WORKDIR}/kmc_output"
LOG_DIR="${WORKDIR}/logfiles"
ODGI_OUT="${WORKDIR}/odgi_output"
CLIP_OUT="${WORKDIR}/clipped_sequence_analyses"
NCOORDS_OUT="${WORKDIR}/node_coord_search_output"
PANKMER_OUT="${WORKDIR}/pankmer_output"
PGSCUNK_OUT="${WORKDIR}/pgscunk_output"


mkdir -p ${REF_DIR} ${STATS_DIR} ${MC_OUT_DIR} ${GIR_OUT_DIR} ${GIR_STATS_DIR} \
${READ_DIR} ${KMC_OUT} ${LOG_DIR} ${ODGI_OUT} ${GRAPH_STATS_DIR} ${GRAPH_STATS_PLOTS} \
${BAM_STATS_DIR} ${CALL_OUT_DIR} ${VARSCAN_OUT_DIR} ${DEDUP_STATS_DIR} ${CLIP_OUT} \
${NCOORDS_OUT} ${PANKMER_OUT} ${PGSCUNK_OUT}

## program directories
STM_DIR="/home/aharder/bin/sequenceTubeMap/scripts"
STMVIZ_DIR="/home/######/seqTubeMap_input_data"
MAMBA="/home/######/.conda/envs"


# -----------------------------------------------------------------------------
# Set 01.sh Minigraph-Cactus options
# -----------------------------------------------------------------------------

## Cactus version and number of threads
readonly CACTUS_IMAGE=docker://quay.io/comparative-genomics-toolkit/cactus:v3.0.0

## Files with sample info
readonly SEQFILE="${INFO_DIR}/cactus_CHROM_seqfile.txt"
readonly CPFILE="${INFO_DIR}/cactus_CHROM_file_copy.txt"

## Set queue and number of threads
QUEUE_01=hammer.q
MC_NTHREADS=32

## Set reference sample -- 2 options here: --------------

## >>> Option 1: Set all input haplotypes as reference paths in the full graph. The first
##				 haplotype/sequence in the list will be the primary reference included in
##				 all subgraphs, but only the number of haplotypes set below for 03.sh will
##				 be included as additional paths in subgraphs.
REFSAMP=$(awk '{print $1}' ${SEQFILE})
declare -a REFARR=( $(awk '{print $1}' ${SEQFILE}) )
PRIMREF=${REFARR[0]}

## >>> Option 2: Set a single input haplotype as a reference path in the full graph (this
##				 will also be the primary reference included in all subgraphs).
##				 Gotta match something in the first column of the SEQFILE.
# REFARR=("CsativaPrytzh")
# REFSAMP=${REFARR[0]}
# PRIMREF=${REFARR[0]}
## ------------------------------------------------------

## Define chrom names
declare -a chroms=(Chr01 Chr02)

# -----------------------------------------------------------------------------
# Set 02.sh vg haplotype options
# -----------------------------------------------------------------------------

## Set version of vg to use (docker updated with every release)
readonly VG_IMAGE=docker://quay.io/vgteam/vg:v1.68.0

## Set number of threads and queue and edit script header
VG_HAP_NTHREADS=32
QUEUE_02=all.q

## set target block size for haplotype sampling (default: 10000);
## applies to 03.sh, too.
CHAINLEN=10000

## Set k-mer length for subgraph construction; applies to 03.sh too
KLEN=29

# -----------------------------------------------------------------------------
# Set 03.sh vg giraffe options
# -----------------------------------------------------------------------------

## Set number of threads and queue and edit script header
VG_GIR_NTHREADS=32
QUEUE_03=all.q

## Set number of haplotypes to include in each subgraph
NHAPLOS=8

## Set list of FASTQ files for samples (must be .gz, not .bz2)
FQ_LIST="${INFO_DIR}/fq_list.txt"


# -----------------------------------------------------------------------------
# Set 04.sh vg call options
# -----------------------------------------------------------------------------

## Set number of threads and queue and edit script header
VG_CALL_NTHREADS=32
QUEUE_04=all.q


# -----------------------------------------------------------------------------
# Set 05.sh varscan options
# -----------------------------------------------------------------------------
VARSCAN_NTHREADS=32
QUEUE_05=all.q

## Minimum acceptable read depth for keeping a variant call
MIN=8

## Minimum acceptable read depth: LO_PROP * mean depth of initial BAM
## Will use higher of this value of MIN (e.g., if LO_LIM = 6, value used will be 8)
LO_PROP=0.25

## Maximum acceptable read depth: UP_PROP * mean depth of initial BAM
UP_PROP=2
