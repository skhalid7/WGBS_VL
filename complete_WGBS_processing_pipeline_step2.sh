#!/bin/bash

# Input arguments for the second group
output_dir=$1
ref=$2
filename=$3
path_2_picard=$4 #optional
path_2_gatk=$5 #optional

cores=$(awk '/^processor/{n+=1}END{print n - 1}' /proc/cpuinfo) #use all available cores -1
ram_gb=$(free | grep Mem | awk '{print int($7 * 0.75 / 1000000)}') #get available RAM

# Step 5: AddReadGroups
echo "Adding Read Groups..."
java -Xmx${ram_gb}g -jar $path_2_picard AddOrReplaceReadGroups \
        I=${output_dir}/${filename}_bismark_bt2_pe.sorted.bam \
        O=${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.bam \
        LB=Veeramah PL=COMPLETE PU=BGI SM=${filename} \
        CN=BGI RGID=${filename}_1 VALIDATION_STRINGENCY=SILENT

# Step 6: Mark duplicates
echo "Marking duplicates..."
java -Xmx${ram_gb}g -jar $path_2_picard MarkDuplicates \
  I=${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.bam \
  O=${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.bam \
  AS=TRUE \
  M=${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.metrics \
  REMOVE_DUPLICATES=TRUE \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=TRUE

# Step 7: Re-sort using readgroups
echo "Re-sorting and indexing for Bismark..."
samtools sort -@$cores -n -o ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.resort.bam ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.bam

# Step 8: Generate flagstat report
echo "Generating flagstat report..."
samtools flagstat -@$cores ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.resort.bam > ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.resort.flagstat

# Step 9: Run bismark_methylation_extractor
echo "Extracting methylation stats..."
bismark_methylation_extractor -p --comprehensive --bedGraph --counts --cytosine_report --genome_folder $ref --multicore $cores ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.resort.bam -o $output_dir

# Step 10: Run bismark_methylation_extractor --mbias-only
echo "Running M-bias extraction..."
bismark_methylation_extractor -p --mbias_only ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.resort.bam --multicore $cores -o $output_dir
