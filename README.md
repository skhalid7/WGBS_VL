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
1. sample_name_bismark_bt2_pe.sorted.RG.mkdups.resort.bismark.cov.gz:<br>
   Chr  Pos_start  Pos_end  methylated%  methylated_read  un_methylated_reads
2. sample_name_bismark_bt2_pe.sorted.RG.mkdups.resort.bismark.bedGraph.gz: (UCSC format)<br>
  Chr  Pos_start  Pos_end  %methylated
3. sample_name_bismark_bt2_pe.sorted.RG.mkdups.resort.M-bias.txt:<br>
  Per read methylation contexts, header is self explanatory
