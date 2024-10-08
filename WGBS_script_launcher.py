#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Run the complete WGBS processing pipeline.")
    
    # Add flag-based arguments
    parser.add_argument('--set1_reads', required=True, help="Path to the first set of paired-end reads.")
    parser.add_argument('--set2_reads', required=True, help="Path to the second set of paired-end reads.")
    parser.add_argument('--output_dir', required=True, help="Directory where output will be saved.")
    parser.add_argument('--ref', required=True, help="Reference genome directory for Bismark.")
    parser.add_argument('--adapter_sequence1', required=True, help="Adapter sequence 1 for Trim Galore.")
    parser.add_argument('--adapter_sequence2', required=True, help="Adapter sequence 2 for Trim Galore.")
    
    # Add optional arguments with default values
    parser.add_argument('--filename', default=None, help="Base filename for outputs (default: based on set1_reads).")
    parser.add_argument('--path_2_picard', default="/opt/picard/picard.jar", help="Path to Picard jar file (default: '/opt/picard/picard.jar').")
    parser.add_argument('--path_2_gatk', default="/opt/gatk4/gatk-package-4.6.0.0-local.jar", help="Path to GATK4 jar file (default: '/opt/gatk4/gatk-package-4.6.0.0-local.jar').")
    parser.add_argument('--path_2_WGBS_script', default=None, help="Path to WGBS bash script (default: Assumes it is in the same path as this python script")
    parser.add_argument('--step', choices=['step1', 'step2', 'complete'], required=True, help="Specify which part of the pipeline to run: 'step1' for Steps TrimGalore + Bismark Alignment, 'step2' for MarkDuplicates + Bismark methylation extractor", default = 'complete')

    # Parse the arguments
    args = parser.parse_args()

    # Set default filename if not provided, derived from set1_reads by removing .fq.gz
    if args.filename is None:
        args.filename = os.path.basename(args.set1_reads).replace('.fq.gz', '')

    #Set defaul WGBS script path if not provided
    if args.path_2_WGBS_script is None:
        args.path_2_WGBS_script = "/".join(sys.argv[0].split("/")[0:-1])

    # Expand the path to Picard if needed
    picard_path = os.path.expanduser(args.path_2_picard)
    gatk_path = os.path.expanduser(args.path_2_gatk)

    # Build the bash command for step 1 and step 2
    bash_command_step1 = [
            "bash", "{}/complete_WGBS_processing_pipeline_step1.sh".format(args.path_2_WGBS_script),
            args.set1_reads, args.set2_reads, args.output_dir,
            args.ref, args.adapter_sequence1, args.adapter_sequence2,
            args.filename
    ]
    
    bash_command_step2 = [
            "bash", "{}/complete_WGBS_processing_pipeline_step2.sh".format(args.path_2_WGBS_script),
            args.output_dir, args.ref, args.filename, picard_path, gatk_path
    ]

    # Run the command
    try:
        if args.step == 'step1' or args.step == 'complete':
            subprocess.run(bash_command_step1, check=True)
            print("Step1 completed successfully.")

        if args.step == 'complete' or args.step == 'step2':
            subprocess.run(bash_command_step2, check=True)
            print("Step2 complete successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running the pipeline: {e}")

if __name__ == "__main__":
    main()

