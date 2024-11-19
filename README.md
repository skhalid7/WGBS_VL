# WGBS_VL
# Example command for WGBS commands:
Python script launches two scripts step 1 and step 2
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
--step [1/2]
````
## Output Files
1. sample_name_bismark_bt2_pe.sorted.RG.mkdups.resort.bismark.cov.gz:<br>
   Chr  Pos_start  Pos_end  methylated%  methylated_read  un_methylated_reads
2. sample_name_bismark_bt2_pe.sorted.RG.mkdups.resort.bismark.bedGraph.gz: (UCSC format)<br>
  Chr  Pos_start  Pos_end  %methylated
3. sample_name_bismark_bt2_pe.sorted.RG.mkdups.resort.M-bias.txt:<br>
  Per read methylation contexts, header is self explanatory


# Example command for WGS processing files:
```bash
#!/bin/bash

# Input arguments
set1_reads=$1
set2_reads=$2
sample_name=$3
RG=$4 #Readgroup file
adapter1=$5
adapter2=$6
ref_file=$7 #reference file
LB=$8 #lab name for AddReadGroups command
PU=$9 #Processing Unit for AddReadGroups command
CN=${10} #center name for AddReadGroups command
output_dir=${11} 

singularity exec --bind /scratch:/mnt --bind /gpfs/shared_drive:/mnt2 \
WGBS_pipeline.sif \
bash complete_WGS_processing_pipeline_step1.sh \
$set1_reads $set2_reads $sample_name $RG $adapter1 $adapter2 $ref_file $LB $PU $CN $output_dir

#!/bin/bash

sample_name=$1
path_2_files=$2

singularity exec --bind /scratch:/mnt --bind /gpfs/shared_drive:/mnt2 \
WGBS_pipeline.sif \
bash complete_WGS_processing_pipeline_step2.sh \
$sample_name $path_2_files
```
## Output Files
1. sample_name.RG.sorted.mkdups.bam
2. Coverage and flagstat file
