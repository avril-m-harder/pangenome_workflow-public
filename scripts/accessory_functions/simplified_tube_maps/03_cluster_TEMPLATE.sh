#! /bin/bash

source ../../scripts/99_init_script_vars.sh
mamba activate ${MAMBA}/pg_tools

for VARS in sSET.norm sSET.norm.simplified sSET.simplified.multi
do

	if [[ -f ${VARS}.vcf.gz ]]; then

		tabix -f ${VARS}.vcf.gz

		vg construct \
			-m 1024 \
			-a \
			-S \
			-r sREF_SEQ \
			-v ${VARS}.vcf.gz > \
			${VARS}.vg
	
		vg gbwt \
			-v ${VARS}.vcf.gz \
			-x ${VARS}.vg \
			-g ${VARS}.gbz \
			--gbz-format

		vg convert -f ${VARS}.gbz > ${VARS}.gfa
		sed -i "1s/.*/H	VN:Z:1.1	RS:Z:BTx623/" ${VARS}.gfa
		vg convert -x ${VARS}.gbz > ${VARS}.xg
	
	fi

done
	
rsync -avuP *.xg /home/aharder_scratch_f13/seqTubeMap_input_data/

mamba deactivate
