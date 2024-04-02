# pangenome_workflow

## Starting ingredients:
```
.
|____info
| |____fq_list.txt
| |____cactus_CHROM_file_copy.txt
| |____cactus_CHROM_seqfile.txt
| |____header_top.tmp
| |____header_bottom.tmp
|
|____scripts
| |____PREP_01_create_single_chrom_fastas.sh
| |____00_prepare_scripts.sh
| |____99_init_script_vars.sh
| |____R01_chrom_coverage_stats.R
| |____R02_chrom_coverage_graphics.R
| |____templates
|   |____01_minigraph-cactus_sCHROM.TEMPLATE.sh
|   |____02_vg_combine_and_haplos.TEMPLATE.sh
|   |____03_vg_giraffe-subgraph.TEMPLATE.sh
|   |____03_vg_giraffe-subgraph.TEMPLATE-INTERLEAVED.sh
|   |____04_vg_call_SVs.TEMPLATE.sh
|   |____05_varscan_SNPs_indels.TEMPLATE.sh
|
|____ref_data
| |____[all input reference FASTAs]
|
|____prepped_reads
  |____[all gzipped (not bgzipped) input FASTQs, interleaved or separate]

```


order of operations:
1. add reference FASTA and read FASTQ files to appropriate dirs
2. update /info/fq_list.txt and /info/cactus*.txt files
3. update 00.sh and 99.sh as necessary
4. edit and run PREP_01.sh as necessary
5. run 00.sh
6. qsub each 01.sh script, run rest of numbered bash jobs in order as previous steps finish
