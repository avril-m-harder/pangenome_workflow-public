# pangenome_workflow

### Starting ingredients:
```
.
|____info
| |____fq_list.txt
| |____cactus_CHROM_file_copy.txt
| |____cactus_CHROM_seqfile.txt
| |____header_top.tmp
| |____header_bottom.tmp
| |____panacus_report_TEMPLATE.yaml
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
|   |____03_vg_giraffe-allhapsgraph.TEMPLATE-INTERLEAVED.sh
|   |____03_vg_giraffe-allhapsgraph.TEMPLATE.sh
|   |____03_vg_giraffe-subgraph.TEMPLATE-INTERLEAVED.sh
|   |____03_vg_giraffe-subgraph.TEMPLATE.sh
|   |____04_vg_call-allhapsgraph.SVs.TEMPLATE.sh
|   |____04_vg_call-subgraph.SVs.TEMPLATE.sh
|   |____05_varscan_SNPs_indels.TEMPLATE.sh
|
|____ref_data
| |____[all input reference FASTAs]
|
|____prepped_reads
  |____[all gzipped (not bgzipped) input FASTQs, interleaved or separate]

```


### Order of operations:
1. add reference FASTA and read FASTQ files (if aligning reads to graph) to appropriate dirs
2. update /info/fq_list.txt (if aligning reads) and /info/cactus*.txt (necessary for graph building) files
3. update 00.sh and 99.sh as necessary
   * in addition to conda env creation, [PG-SCUnK](https://github.com/cumtr/PG-SCUnK/) needs to be installed
4. edit and run PREP_01.sh as necessary to split genome FASTAs into chromosome FASTAs
5. run 00.sh
6. qsub each 01.sh script, run rest of numbered bash jobs in order as previous steps finish


