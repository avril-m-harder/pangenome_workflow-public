#! /bin/bash

source ../../scripts/99_init_script_vars.sh
mamba activate ${MAMBA}/pg_tools

rsync -avuP ${MC_OUT_DIR}/sTAXON_sCHROM/sTAXON_sCHROM/sTAXON_sCHROM.vcf.gz ${OUTDIR}

tabix sTAXON_sCHROM.vcf.gz

bcftools view \
	-r sCHROM:sSTART-sEND \
	-a \
	-m2 \
	-O z \
	-o sREGION.vcf.gz \
	sTAXON_sCHROM.vcf.gz

mamba deactivate
