#!/bin/bash

# Input arguments
set1_reads=$1
set2_reads=$2
output_dir=$3
ref=$4
adapter_sequence1=$5
adapter_sequence2=$6
filename=$7 #optional
path_2_picard=$8 #optional
path_2_gatk=$9 #optional

# Determine available cores and RAM, WARNING if this is being run on a shared cluster this should be manually configured!!
cores=$(awk '/^processor/{n+=1}END{print n - 1}' /proc/cpuinfo) #use all available cores -1
ram_gb=$(free | grep Mem | awk '{print int($7 * 0.75 / 1000000)}') #get available RAM

# Step 1: Run trim_galore on paired-end reads
echo "Running trim_galore..."
trim_galore --output_dir $output_dir --paired --cores $cores --quality 20 --fastqc --adapter $adapter_sequence1 --adapter2 $adapter_sequence2 $set1_reads $set2_reads --basename ${filename}_trimmed

# Define the output from trim_galore
trimmed_r1="${output_dir}/${filename}_trimmed_val_1.fq.gz"
trimmed_r2="${output_dir}/${filename}_trimmed_val_2.fq.gz"

# Step 2: Run bismark using the trimmed reads
echo "Running bismark..."
cores_bismark=$(echo $cores | awk '{print int($0 / 3)}') #bismark needs some extra cores because it runs bowtie on 4 cores. If you specify 4 cores it will run on 8 as per its documentation
bismark -o $output_dir $ref -1 $trimmed_r1 -2 $trimmed_r2 --parallel $cores_bismark --temp_dir $output_dir

# Step 3: Sorting the BAM file
echo "Sorting BAM file..."
samtools sort -@$cores -o ${output_dir}/${filename}_bismark_bt2_pe.sorted.bam ${output_dir}/${filename}_trimmed_val_1_bismark_bt2_pe.bam #note the change in basename at this point!

# Step 4: Indexing the BAM file
echo "Indexing BAM file..."
samtools index -@$cores ${output_dir}/${filename}_bismark_bt2_pe.sorted.bam

# Step 5: AddReadGroups
echo "Adding Read Groups... "
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
echo "Re-sorting and indexing for Bismark... "
samtools sort -@$cores -n -o ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.resort.bam ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.bam #indexing this file will throw an error becuase it isn't sorted on position anymore

# Step 8: Generate flagstat report
echo "Generating flagstat report..."
samtools flagstat -@$cores ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.resort.bam > ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.resort.flagstat

#Step 9: Run bismark_methylation_extractor
echo "Extracting methylation sats"
bismark_methylation_extractor -p --comprehensive --bedGraph --counts --cytosine_report --genome_folder $ref --multicore $cores ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.resort.bam -o $output_dir

#Step 10: Run bismark_methylation_extractor --mbias-only
bismark_methylation_extractor -p --mbias_only ${output_dir}/${filename}_bismark_bt2_pe.sorted.RG.mkdups.resort.bam --multicore $cores -o $output_dir

#Step 11
#Run GATK DepthOfCoverage
