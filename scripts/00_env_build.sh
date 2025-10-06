#!/bin/bash
#
# prepare_scripts.sh
#

source /home/######/pangenome_workflow/scripts/99_init_script_vars.sh 

# -----------------------------------------------------------------------------
# Create conda env for all tools needed in the MC graph workflow.
# (other than cactus + vg --> docker and pg-scunk which is picky about dependencies
# and needs its own environment, below)
# -----------------------------------------------------------------------------

mamba create -p ${MAMBA}/pg_tools \
	kmc \
	odgi \
	panacus \
	samtools \
	picard \
	bcftools \
	varscan

# Using mamba 
mamba create -p ${MAMBA}/pg-scunk \
	bioconda::kmc=3.2.4 \
	bioconda::samtools=1.21 \
	bioconda::vg=1.65.0 \
	conda-forge::zlib=1.3.1 \
	conda-forge::r-base \
	bioconda::bwa=0.7.19 \
	bioconda::bedtools=2.31.1 
