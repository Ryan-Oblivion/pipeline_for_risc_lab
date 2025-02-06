// first i will activate dsl2
nextflow.enable.dsl=2

// NOTE: i would prefer if we could just keep all the reads in one directory and just use glob pattern to grab what we want

// first process is to download genomes that will be used or just get the path to the genome
//process download_genomes{}

// this is the process for if the adapter sequence is known


process fastp_SE_adapter_known {
    // using this conda yml file that I created. Nextflow will now make its own conda environment from the dependencies found within.
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastp_rj_env.yml'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    publishDir './fastp_qc_single_end', mode: 'copy', pattern:'*_fp_filt.fastq'
    publishDir './fastp_qc_single_end/html_reports', mode: 'copy', pattern:'*.html'

    input:
    // input names dont have to be the exact same as what is seen in the workflow section of this script
    // as long and they are in the same order, you will have the correct input
    path(fastq_files) // the files
    val(fastq_names) // the names of the files as a value input channel
    val(adapter_seq)


    output:
    path("${out_name}"), emit: filtered_fastqs
    path("${fastq_names}*.html"), emit: fastp_html_reports


    script:

    // getting the output name
    out_name = "${fastq_names}_fp_filt.fastq"


    """
    #!/usr/bin/env bash

    # ASK HERA ABOUT THE ADAPTERS USED AND FIND THEM FROM THE ILLUMINA KIT

    # remember this process is for the single end reads so only one read in and one filtered qc read out
    # adapter trimming is enabled by defualt for single end but use --detect_adapter_for_pe to enable it in paired end data
    # you can specify the adapter sequence if you have them. look at documentation here: https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp
    
    # --adapter_sequence: string sequence that represents the adapter used. in the illumina NEBNext libraries adapter read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA. Adapter read 2 is AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    # --dedup enables the deduplication to drop the duplicated reads
    # --dup_calc_accuracy: defualt is 3 but can go from 1~6. higher levels use more memory
    # --trim_poly_g: novaSeq data polyG can happen in read tails since G means no signal in the illumina two-color system. fastp can detect and trim them.
    # --qualified_quality_phred: the quality value that a base is qualified defualt is 15
    # --unqualified_percent_limit: how many percent of bases are allowed to be unqualified defualt is 40
    # --umi: not using this yet have to use --umi_loc also and i dont know that, fastp can extract the unique molecular identifiers and append them to the first part of the read names, so the umi's will be present in SAM/BAM records If the UMI is in the reads, then it will be shifted from read so that the read will become shorter. If the UMI is in the index, it will be kept.
    # --overrepresentation_analysis and --overrepresentation_sampling see fastp documentation
    # --html: this is how the report file will be made

    # their fastq2bam single end pipeline used trim_galore tool for trimming adapters.
    # I want to implement a way to choose to let fastp do this by default or use a known sequence that the user can input.
    ###### can add or remove more options when needed #####
    
    fastp \
    --in1 "${fastq_files}" \
    --out1 "${out_name}" \
    --adapter_sequence "${adapter_seq}" \
    --dedup \
    --dup_calc_accuracy 5 \
    --trim_poly_g \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --overrepresentation_analysis \
    --overrepresentation_sampling 20 \
    --html "${fastq_names}_fastp.html"






    """

}







// This next process will run some qc to look at the fastq files and trim adapters from the single end reads
process fastp_SE {
    // using this conda yml file that I created. Nextflow will now make its own conda environment from the dependencies found within.
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastp_rj_env.yml'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    publishDir './fastp_qc_single_end', mode: 'copy', pattern:'*_fp_filt.fastq'
    publishDir './fastp_qc_single_end/html_reports', mode: 'copy', pattern:'*.html'

    input:
    // input names dont have to be the exact same as what is seen in the workflow section of this script
    // as long and they are in the same order, you will have the correct input
    path(fastq_files) // the files
    val(fastq_names) // the names of the files as a value input channel



    output:
    path("${out_name}"), emit: filtered_fastqs
    path("${fastq_names}*.html"), emit: fastp_html_reports


    script:

    // getting the output name
    out_name = "${fastq_names}_fp_filt.fastq"


    """
    #!/usr/bin/env bash

    # ASK HERA ABOUT THE ADAPTERS USED AND FIND THEM FROM THE ILLUMINA KIT

    # remember this process is for the single end reads so only one read in and one filtered qc read out
    # adapter trimming is enabled by defualt for single end but use --detect_adapter_for_pe to enable it in paired end data
    # you can specify the adapter sequence if you have them. look at documentation here: https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp
    
    # --adapter_sequence: string sequence that represents the adapter used. in the illumina NEBNext libraries adapter read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA. Adapter read 2 is AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    # --dedup enables the deduplication to drop the duplicated reads
    # --dup_calc_accuracy: defualt is 3 but can go from 1~6. higher levels use more memory
    # --trim_poly_g: novaSeq data polyG can happen in read tails since G means no signal in the illumina two-color system. fastp can detect and trim them.
    # --qualified_quality_phred: the quality value that a base is qualified defualt is 15
    # --unqualified_percent_limit: how many percent of bases are allowed to be unqualified defualt is 40
    # --umi: not using this yet have to use --umi_loc also and i dont know that, fastp can extract the unique molecular identifiers and append them to the first part of the read names, so the umi's will be present in SAM/BAM records If the UMI is in the reads, then it will be shifted from read so that the read will become shorter. If the UMI is in the index, it will be kept.
    # --overrepresentation_analysis and --overrepresentation_sampling see fastp documentation
    # --html: this is how the report file will be made

    # their fastq2bam single end pipeline used trim_galore tool for trimming adapters.
    # I want to implement a way to choose to let fastp do this by default or use a known sequence that the user can input.
    ###### can add or remove more options when needed #####
    
    fastp \
    --in1 "${fastq_files}" \
    --out1 "${out_name}" \
    --dedup \
    --dup_calc_accuracy 5 \
    --trim_poly_g \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --overrepresentation_analysis \
    --overrepresentation_sampling 20 \
    --html "${fastq_names}_fastp.html"






    """

}

// this next process is for fastqc tool
// dont forget about multiqc m
process fastqc_SE {
    // using the conda environment 
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastqc_rj_env.yml'
    publishDir './fastqc_htmls', mode: 'copy', pattern: '*.html'


    input:

    path(fastq_filt_files)
    val(fastq_filt_names)


    output:
    path("*.html"), emit: fastqc_htmls
    path("*.zip"), emit: fastqc_zip_files

    script:
    out_name = "${fastq_filt_files}"

    """
    #!/usr/bin/env bash

    # I need to add the adapter sequences later, ask hera.

    fastqc "${fastq_filt_files}"
    



    """
}

process multiqc_SE {
    // this yml file doesnt work
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/multiqc_rj_env.yml'

    //conda '/ru-auth/local/home/rjohnson/miniconda3/envs/multiqc_rj'
    
    publishDir './multiQC_collection', mode: 'copy', pattern: '*.html'

    input:
    path(fastp_filt_html)



    output:

    path("*.html"), emit: multiqc_html_collection

    script:


    """
    #!/usr/bin/env bash

    # I think i jsut have to pass the html files to multiqc for it to complie them
    multiqc . \
    --interactive \
    --profile-runtime \
    --title "Single-end QC"




    """

}


// Creating two processes that will index the reference genome

process bwa_index_genome {
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    publishDir './genome_index_bwa', mode: 'copy', pattern: '*'

    input:
    path(ref_genome)
    


    output:

    path("*"), emit: genome_index_files


    script:
    // not getting the basename because bwa expects the index files will have the exact file name of the genome but with an .{ext} on the end. 
    // example genome.fa, genome.fa.amb, genome.fa.ann
    // not genome.fa, genome.amb, genome.ann
    //genome_file_name = "${ref_genome.baseName}"

    """
    #!/usr/bin/env bash

    ############### parameters used ###############
    # -p: a string representing the prefix of the output database [same as the db filename] so i think just the base name of the genome file
    # -a: a string, choosing the algorithm to construct the BWT index. either is or bwtsw. I'll use bwtsw


    bwa index \
    -p "${ref_genome}" \
    -a bwtsw \
    "${ref_genome}"


    """

}

// creating a process that will align the reads to the genome. i will take in the reference genome, the index files, the filtered fastq's and their names

process bwa_align_SE {
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    publishDir './bwa_outputs_singleEnd_SAM', mode: 'copy', pattern: '*.sam'
    publishDir './sai_alignment_files', mode: 'copy', pattern: '*.sai'


    input:
    path(ref_genome)
    path(genome_index_files)
    path(fastq_filt_files)
    val(fastq_filt_names)


    output:

    path("*.sam"), emit: sam_se_files
    path("*.sai"), emit: sai_align_files


    script:
    sai_output_file = "${fastq_filt_names}_out.sai"
    sam_name = "${fastq_filt_names}.sam"

    """
    #!/usr/bin/env bash

    ############# Parameters used ############
    # first i need to get the sai file by using bwa aln. This gives the SA coordinates of the input reads
    # -t (nThrds) number of threads for multi threading mode. defualt is 1

    # using bwa samse : this will generate alignments in the SAM format given single-end reads
    # the two parameters may be used in the future.
    # -n: takes an integer. max number of alignments to output in the XA tag for reads paired properly
    # -r: takes a string. specify the read group
    #
    ##########################################

    ls .

    bwa aln \
    -t 20 \
    "${ref_genome}"  \
    "${fastq_filt_files}" \
    > "${sai_output_file}"


    bwa samse \
    "${ref_genome}" \
    "${sai_output_file}" \
    "${fastq_filt_files}" \
    > "${sam_name}"


    """
}

process samtools_sort {
    // using the conda yml file for samtools
    // it doesnt work
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/samtools_rj_env.yml'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools_rj'

    publishDir './sorted_bam_files', mode: 'copy', pattern: '*_sorted.bam'
    publishDir './indexed_bam_files', mode: 'copy', pattern: '*.{bai, csi}'

    input:
    path(sam_files)


    output:

    path("*_sorted.bam"), emit: sorted_bams
    //tuple path("*.{bai,csi}"), emit: indexed_bams
    //tuple path("*.bai"), path("*.csi"), emit: indexed_bams
    path("*.bai"), emit: indexed_bams


    script:

    // i will start using baseName inside the process since its easier to keep track of different names an uses less inputs into a process
    out_bam = "${sam_files.baseName}_sorted.bam"



    """
    #!/usr/bin/env bash

    ################# samtools parameters used ################
    # for samtools sort
    # -o : takes a file. it writes the final sorted output to file rather than standard output
    # -O : write the final output as sam, bam, or cram

    # now for samtools index to get index files
    # -b, --bai: create a bai index; this version of samtools does not support --bai --csi just use -b -c
    # -c, --csi: create a csi index
    # -o, --output: write the output index to a file specified  only when one alignment file is being indexed

    ###########################################################

    samtools sort \
    -o "${out_bam}" \
    -O bam \
    "${sam_files}"

    # so i will use the out_bam for input to samtools index since it has to be coordinate sorted
    # I will not make an out file name since I hope samtools index will just add the prefix
    # i was not able to use both -b and -c in the same samtools index call. i can just write another samtools index with -c instead if i want that index also.

    samtools index \
    -b \
    "${out_bam}" 






    """
}

workflow {

    // this is the end seq alignment steps first


    // i will use a path already in the hpc as the defualt human genome but the user can change the genome by using -human_genome parameter and putting the path to a new genome in the command line when calling nextflow run
    params.human_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/genome/Sequence/Bowtie2Index/genome.fa')

    // putting the human genome in a channel
    // keeping the human genome in a value channel so i can have other processes run more than once.
    human_genome_ch = Channel.value(params.human_genome)



    // i want to add an if then logic to the pipeline so i know which type of reads are comming in paired end or single end

    //if ( params.PE )
         // lets get the channel for the reads first

      //   fastp_PE()
         
         //align_PE_reads(human_genome_ch)
    
    
    
    
    if ( params.SE )

        // lets get the channel for the single end reads first
        // only use the single end read 1 data from the end seq which are already stored here: /rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs

    

         params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
         se_reads_files = Channel.fromPath(params.single_end_reads)
         
         // now let's get the basename of the single end reads
         // removing both the .gz and the .fastq
         // I would normally use file.baseName here but it had the .gz and the .fastq
         se_reads_files.flatten()
                        .map{ file -> file.name.replace('.fastq.gz','')}
                        .set{se_reads_name}
         
         // let's view both the files and the names to make sure they match in order
         //se_reads_files.view()
         //se_reads_name.view()
         // this is where i send both the input file and their corresponding basenames to the fastp_SE process
         

         // if the adapter sequence is known then input it as a string if not dont use the parameter
        if ( params.ada_seq ) {

            params.adapter_seq_str = 'AGATCGGAAGAGC' // this is just a place holder value for the adapter sequence
            adapter_ch = Channel.value(params.adapter_seq_str)

            fastp_SE_adapter_known(se_reads_files.take(3), se_reads_name.take(3), adapter_ch) // will have to make a new process for if the adapter sequence is known

            fastq_filts = fastp_SE_adapter_known.out.filtered_fastqs
            //fastp_SE.out.view()
            fastq_filts.map{file -> file.baseName}
                        .set{fastq_filts_name}

            // now getting the html files since i think fastqc combines them into one, that might be multiqc
            fastp_filt_html = fastp_SE_adapter_known.out.fastp_html_reports

        }    
        else {

            fastp_SE(se_reads_files.take(3), se_reads_name.take(3))

                // take all of the filtered fastq files and put them in a channel name
            // since the fastq files might be in a different order, if i need to get their base names I will have to do it from this new channel below
            fastq_filts = fastp_SE.out.filtered_fastqs
            //fastp_SE.out.view()
            fastq_filts.map{file -> file.baseName}
                        .set{fastq_filts_name}

            // now getting the html files since i think fastqc combines them into one, that might be multiqc
            fastp_filt_html = fastp_SE.out.fastp_html_reports



        }
        

         //fastp_SE(se_reads_files.take(3), se_reads_name.take(3)) // REMEMBER TO REMOVE THIS TESTING FEATURE WHERE IT WILL ONLY TAKE THE FIRST 3

         

         
         //fastp_SE.out.view()

         //fastq_filts.view()
         //fastp_filt_html.view()
        //fastp_filt_html.collect().view()


         // now creating a fastqc process
        fastqc_SE(fastq_filts, fastq_filts_name)

        fastqc_html_files = fastqc_SE.out.fastqc_htmls
        fastqc_zips = fastqc_SE.out.fastqc_zip_files

        // now using multiqc to combine all of the se zip files. Multiqc takes the zip files generated by fastqc and puts them in a single html file
        multiqc_SE(fastqc_zips.collect())

        // first have a seprate process that indexes the reference genome using bwa or bwa mem so this part doesnt have to be done again and will be cached
        bwa_index_genome(human_genome_ch)

        // collecting the genome index files from the last process 
        // not sure if i should keep track of the order the files are in first
        // it looks like they are in the same order that they appeared in the published dir using ll
        //bwa_index_genome.out.genome_index_files.view()
        genome_index_files_ch = bwa_index_genome.out.genome_index_files

        // Now I need to pass the human genome file to the process to index the genome file. Also I will add the filtered fastq files from fastp into this process that will be aligned to the genome. 
        // each run of this only takes 20-30 min to run but since the hpc only is allowing 2-3 to run at one time it takes 3 hours
        bwa_align_SE(human_genome_ch, genome_index_files_ch, fastq_filts, fastq_filts_name )

        //bwa_align_SE.out.sam_se_files.view()

        // making a channel for the sam files generated
        sam_files = bwa_align_SE.out.sam_se_files

        // now I want to take any sam files generated by the bwa and use samtools to order them and convert them into bam files
        // I will hopefully be able to do this outside of the if else statement so the sam file from both conditions can be passed to the same samtools process
        // since i am just testing the pipeline i should find a way to do this on only a few files (about 3-4)
        
        samtools_sort(sam_files.take(3)) // using take 3 should only take the first 3 files from the sam_files channel which should have 64 sam files. This is just for production and testing. will remove when running pipeline for real.

        samtools_sort.out.sorted_bams.view()
        samtools_sort.out.indexed_bams.view()
}