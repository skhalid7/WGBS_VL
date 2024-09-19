# WGBS_VL
Example command:
```bash
singularity exec --bind /gpfs/path/to/reads:/mnt --bind /additional/path/:/genome \
WGBS_pipeline.sif \
python /mnt/processing/WGBS_script_launcher.py \
--set1_reads /mnt/path/to/reads_1.fq.gz \
--set2_reads /mnt/path/to/reads_2.fq.gz \
--output_dir /mnt/output/folder \
--ref /genome \
--filename sample_name\
--adapter_sequence1 ATGCA \
--adapter_sequence2 AGATC
````


## Output Files
sample_name_bismark_bt2_pe.sorted.RG.mkdups.resort.bismark.cov.gz: 
#Chr Pos_start Pos_end methylated% methylated_read un_methylated_reads

sample_name_bismark_bt2_pe.sorted.RG.mkdups.resort.bismark.bedGraph.gz: (UCSC filename)
#Chr Pos_start Pos_end %methylated

sample_name_bismark_bt2_pe.sorted.RG.mkdups.resort.M-bias.txt:
Per read methylation contexts, header is self explanatory
