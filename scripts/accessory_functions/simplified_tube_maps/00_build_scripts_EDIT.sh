#! /bin/bash

# -----------------------------------------------------------------------------
# Define variables for ROI
# -----------------------------------------------------------------------------

## uses minigraph-cactus naming convention
TAXON="sorghum31"
CHROM="Chr01"
ROI="LGS1"
REF="BTx623"
START=75097270
END=75132500
MIN_L=5000
MIN_KB=$((MIN_L/1000))
REF_SEQ="/home/aharder_scratch_f13/sorghum_pangenome/ref_data/Sorghum_bicolor_var_BTx623.${CHROM}.fa"

## name files according to ROI or coordinates
# REGION="${TAXON}_${CHROM}_${START}-${END}"
REGION="${TAXON}_${CHROM}_${ROI}"

# -----------------------------------------------------------------------------
# Auto mkdir and write scripts to it
# -----------------------------------------------------------------------------

source ../../99_init_script_vars.sh

SET="${TAXON}_${CHROM}_${START}-${END}-vcfw-${MIN_KB}kb"
OUTDIR="${WORKDIR}/simplified_tube_maps/${SET}"
mkdir -p ${OUTDIR}

## prep 01.sh
cp 01_cluster_TEMPLATE.sh ${OUTDIR}/01_cluster-${REGION}.sh
sed -i "s|sTAXON|${TAXON}|g" ${OUTDIR}/01_cluster-${REGION}.sh
sed -i "s|sCHROM|${CHROM}|g" ${OUTDIR}/01_cluster-${REGION}.sh
sed -i "s|sSTART|${START}|g" ${OUTDIR}/01_cluster-${REGION}.sh
sed -i "s|sEND|${END}|g" ${OUTDIR}/01_cluster-${REGION}.sh
sed -i "s|sREGION|${REGION}|g" ${OUTDIR}/01_cluster-${REGION}.sh

## prep 02.sh
cp 02_local_TEMPLATE.sh ${OUTDIR}/02_local-${REGION}.sh
sed -i "s|sREGION|${REGION}|g" ${OUTDIR}/02_local-${REGION}.sh
sed -i "s|sOUTDIR|${OUTDIR}|g" ${OUTDIR}/02_local-${REGION}.sh
sed -i "s|sMIN_L|${MIN_L}|g" ${OUTDIR}/02_local-${REGION}.sh
sed -i "s|sREF|${REF}|g" ${OUTDIR}/02_local-${REGION}.sh
sed -i "s|sMIN_KB|${MIN_KB}|g" ${OUTDIR}/02_local-${REGION}.sh

## prep 03.sh
cp 03_cluster_TEMPLATE.sh ${OUTDIR}/03_cluster-${REGION}.sh
sed -i "s|sSET|${SET}|g" ${OUTDIR}/03_cluster-${REGION}.sh
sed -i "s|sREF_SEQ|${REF_SEQ}|g" ${OUTDIR}/03_cluster-${REGION}.sh
