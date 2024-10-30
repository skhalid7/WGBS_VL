#!/bin/bash

# Input arguments
sample_name=$1
path_2_files=$2

# Determine available cores and RAM, WARNING if this is being run on a shared cluster this should be manually configured!!
cores=$(awk '/^processor/{n+=1}END{print n - 1}' /proc/cpuinfo) #use all available cores -1
ram_gb=$(free | grep Mem | awk '{print int($7 * 0.75 / 1000000)}') #get available RAM

echo $path_2_files
echo $sample_name

ls ${path_2_files}/*${sample_name}*RG.bam

echo 'Merginging bams...'
samtools merge -@ $cores ${path_2_files}/${sample_name}.RG.merged.bam ${path_2_files}/*${sample_name}_*RG.bam

echo 'Merge and index...'
samtools sort -@ $cores -o ${path_2_files}/${sample_name}.RG.merged.sort.bam ${path_2_files}/${sample_name}.RG.merged.bam
samtools index -@ $cores ${path_2_files}/${sample_name}.RG.merged.sort.bam

echo 'Marking duplicates...'
java -Xmx${ram_gb}g -jar /opt/picard/picard.jar MarkDuplicates I=${path_2_files}/${sample_name}.RG.merged.sort.bam \
	O=${path_2_files}/${sample_name}.RG.merged.sort.mkdups.bam \
	AS=TRUE \
	M=${path_2_files}/${sample_name}.RG.merged.sort.mkdups.metrics \
	REMOVE_DUPLICATES=FALSE \
	VALIDATION_STRINGENCY=LENIENT \
	CREATE_INDEX=TRUE

echo 'Create stat files'
samtools flagstat -@ $cores ${path_2_files}/${sample_name}.RG.merged.sort.mkdups.bam > ${path_2_files}/${sample_name}.RG.merged.sort.mkdups.flagstat
samtools coverage -q20 -Q20 ${path_2_files}/${sample_name}.RG.merged.sort.mkdups.bam > ${path_2_files}/${sample_name}.RG.merged.sort.mkdups.coverage.stats


rm ${path_2_files}/${sample_name}.RG.merged.bam
rm ${path_2_files}/${sample_name}.RG.merged.sort.bam
