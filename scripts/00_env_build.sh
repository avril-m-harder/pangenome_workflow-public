#!/bin/bash
#
# prepare_scripts.sh
#

source /home/aharder_scratch_f13/pangenome_workflow/scripts/99_init_script_vars.sh 

# -----------------------------------------------------------------------------
# Create conda env for all tools needed in the MC graph workflow.
# (other than cactus + vg --> docker)
# -----------------------------------------------------------------------------

mamba create -p ${MAMBA}/pg_tools \
	kmc \
	odgi \
	panacus \
	samtools \
	picard \
	bcftools \
	varscan
