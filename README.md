# WGBS_VL
Example command:
```bash
singularity exec --bind /gpfs/path/to/reads:/mnt --bind /additional/path/:/genome WGBS_pipeline.sif python /mnt/processing/WGBS_script_launcher.py --set1_reads /mnt/processing/GEPE-PCHA-12_1_mini_test.fq.gz --set2_reads /mnt/processing/GEPE-PCHA-12_2_mini_test.fq.gz --output_dir /mnt/output_2 --ref /genome/Rachael/P_papua_genome/ --filename GEPE-PCHA-12 --adapter_sequence1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence2 AGATCGGAAGCGTCGTGTAGGGAAAGAGTGT
````
