#!/bin/bash

# Input arguments
set1_reads=$1
set2_reads=$2
sample_name=$3
RG=$4
adapter1=$5
adapter2=$6
ref_file=$7
LB=$8
PU=$9
CN=${10}
output_dir=${11}

echo $CN
echo $output_dir

# Determine available cores and RAM, WARNING if this is being run on a shared cluster this should be manually configured!!
cores=$(awk '/^processor/{n+=1}END{print n - 1}' /proc/cpuinfo) #use all available cores -1
ram_gb=$(free | grep Mem | awk '{print int($7 * 0.75 / 1000000)}') #get available RAM

# Step 1: Run trim_galore on paired-end reads
echo "Running trim_galore..."
echo $output_dir
echo $cores
echo $adapter1
echo $adapter2
echo $set1_reads
echo $set2_reads
echo $RG

trim_galore --output_dir $output_dir --paired --cores $cores --quality 20 --fastqc --adapter $adapter1 --adapter2 $adapter2 $set1_reads $set2_reads --basename ${sample_name}_${RG}_trimmed

# Define the output from trim_galore
trimmed_r1="${output_dir}/${sample_name}_${RG}_trimmed_val_1.fq.gz"
trimmed_r2="${output_dir}/${sample_name}_${RG}_trimmed_val_2.fq.gz"

#Step 2: Run BWA-MEM and samtools to align and sort reads
bwa mem -t $cores $ref_file $trimmed_r1 $trimmed_r2 | samtools view -b -h -@ $cores | samtools sort -@ $cores -O BAM -o ${output_dir}/${sample_name}_${RG}.bam

# Step 3: AddReadGroups
echo "Adding Read Groups..."
java -Xmx${ram_gb}g -jar /opt/picard/picard.jar AddOrReplaceReadGroups \
        I=${output_dir}/${sample_name}_${RG}.bam \
        O=${output_dir}/${sample_name}_${RG}.RG.bam \
        LB=$LB PL=COMPLETE PU=$PU SM=${sample_name} \
        CN=BGI RGID=$RG VALIDATION_STRINGENCY=SILENT


#removing initial bam file without read groups
rm ${output_dir}/${sample_name}_${RG}.bam
