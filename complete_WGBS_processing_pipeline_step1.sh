#!/bin/bash

# Input arguments
set1_reads=$1
set2_reads=$2
output_dir=$3
ref=$4
adapter_sequence1=$5
adapter_sequence2=$6
filename=$7 #optional

# Determine available cores and RAM
cores=$(awk '/^processor/{n+=1}END{print n - 1}' /proc/cpuinfo) #use all available cores -1

# Step 1: Run trim_galore on paired-end reads
echo "Running trim_galore..."
trim_galore --output_dir $output_dir --paired --cores $cores --quality 20 --fastqc --adapter $adapter_sequence1 --adapter2 $adapter_sequence2 $set1_reads $set2_reads --basename ${filename}_trimmed

# Define the output from trim_galore
trimmed_r1="${output_dir}/${filename}_trimmed_val_1.fq.gz"
trimmed_r2="${output_dir}/${filename}_trimmed_val_2.fq.gz"

# Step 2: Run bismark using the trimmed reads
echo "Running bismark..."
cores_bismark=$(echo $cores | awk '{print int($0 / 3)}') #bismark needs some extra cores because it runs bowtie on 4 cores.
bismark -o $output_dir $ref -1 $trimmed_r1 -2 $trimmed_r2 --parallel $cores_bismark --temp_dir $output_dir

# Step 3: Sorting the BAM file
echo "Sorting BAM file..."
samtools sort -@$cores -o ${output_dir}/${filename}_bismark_bt2_pe.sorted.bam ${output_dir}/${filename}_trimmed_val_1_bismark_bt2_pe.bam

# Step 4: Indexing the BAM file
echo "Indexing BAM file..."
samtools index -@$cores ${output_dir}/${filename}_bismark_bt2_pe.sorted.bam

