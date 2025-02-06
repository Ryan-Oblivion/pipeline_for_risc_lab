#!/bin/env bash

#SBATCH --mem=40GB
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --job-name=nextflow_chip

source $HOME/.bashrc_rj_test.sh

conda activate nextflow_two

# --PE paramater for single end reads
# --SE parameter for pair end reads
# when useing SE do --single_end_reads and give the path to your single end reads with a glob pattern if you have other reads in that directory you dont want
# if you have an adapter sequence use --ada_seq, then specify the sequence with --adapter_seq_str which will take the string sequence
nextflow run nasa_pipeline.nf -profile 'nasa_pipeline ' \
-resume \
--SE \
--single_end_reads '/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz' \
--ada_seq --adapter_seq_str 'AGATCGGAAGAGC'