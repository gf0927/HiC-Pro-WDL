version 1.0

workflow hic_docker {
    
    input {
        #Setting
        String docker_image

        String BOWTIE2_GLOBAL_OPTIONS
        String BOWTIE2_LOCAL_OPTIONS

        String N_cpu
        String N_mem

        String cutsite

        String OUTPUT_DATA_PATH

        #PATH
        String HiC_PATH


        #inputs
        String Sample_Name
        String REFERENCE_GENOME
        String Bowtie2Index_PATH

        File Sample_r1
        File Sample_r2
        File BED_FILE


        File Bowtie_idx_1
        File Bowtie_idx_2
        File Bowtie_idx_3
        File Bowtie_idx_4
        File Bowtie_idx_re1
        File Bowtie_idx_re2
        File Bowtie_makesh

    }

    # call unzip_bowtie2index {
    #     input:
    #         bwt_index = Bowtie2Index_PATH,
    #         docker_image =docker_image,
    #         REFERENCE_GENOME = REFERENCE_GENOME
    # }

    call bowtie_global_mapping {
        input:
            bowtie_idx_1 = Bowtie_idx_1,
            bowtie_idx_2 = Bowtie_idx_2,
            bowtie_idx_3 = Bowtie_idx_3,
            bowtie_idx_4 = Bowtie_idx_4,
            bowtie_idx_re1 = Bowtie_idx_re1,
            bowtie_idx_re2 = Bowtie_idx_re2,
            bowtie_makesh = Bowtie_makesh,
            sample_name  = Sample_Name,
            REFERENCE_GENOME = REFERENCE_GENOME,
            BOWTIE2_GLOBAL_OPTIONS = BOWTIE2_GLOBAL_OPTIONS,
            odir = OUTPUT_DATA_PATH,
            BOWTIE2_IDX = Bowtie2Index_PATH,
            SAMPLE_R1 = Sample_r1,
            SAMPLE_R2 = Sample_r2,
            bwt_cpu = N_cpu,
            docker_image = docker_image,
            task_mem = N_mem
    }

    call bowtie_local_trimming {
        input:
            gol_unmap_fastq1=bowtie_global_mapping.gol_unmap,
            gol_unmap_fastq2=bowtie_global_mapping.gol_unmap2,
            HiC_PATH=HiC_PATH,
            cutsite=cutsite,
            sample_name = Sample_Name,
            REFERENCE_GENOME=REFERENCE_GENOME,
            docker_image = docker_image,
            task_cpu = N_cpu,
            task_mem = N_mem
    }

    call bowtie_local_mapping{
        input:
            unmap_trimmed_r1=bowtie_local_trimming.unmap_trimmed_r1,
            unmap_trimmed_r2=bowtie_local_trimming.unmap_trimmed_r2,
            bowtie_idx_1 = Bowtie_idx_1,
            bowtie_idx_2 = Bowtie_idx_2,
            bowtie_idx_3 = Bowtie_idx_3,
            bowtie_idx_4 = Bowtie_idx_4,
            bowtie_idx_re1 = Bowtie_idx_re1,
            bowtie_idx_re2 = Bowtie_idx_re2,
            bowtie_makesh = Bowtie_makesh,
            sample_name  = Sample_Name,
            task_cpu = N_cpu,
            REFERENCE_GENOME = REFERENCE_GENOME,
            sample_name = Sample_Name,
            BOWTIE2_LOCAL_OPTIONS = BOWTIE2_LOCAL_OPTIONS,
            docker_image = docker_image,
            task_mem = N_mem
    }

    call mapping_combine {
        input:
            global_mapped_r1 = bowtie_global_mapping.gol_mapped,
            global_mapped_r2 = bowtie_global_mapping.gol_mapped2,
            local_mapped_r1 = bowtie_local_mapping.local_mapped_r1,
            local_mapped_r2 = bowtie_local_mapping.local_mapped_r2,
            sample_name = Sample_Name,
            REFERENCE_GENOME = REFERENCE_GENOME,
            docker_image = docker_image,
            task_cpu = N_cpu,
            task_mem = N_mem
    }
    call merge_pairs {
        input:
            merged_r1 = mapping_combine.bwt2merged_r1,
            merged_r2 = mapping_combine.bwt2merged_r2,
            REFERENCE_GENOME = REFERENCE_GENOME,
            sample_name = Sample_Name,
            HiC_PATH = HiC_PATH,
            docker_image = docker_image,
            task_cpu = N_cpu,
            task_mem = N_mem
    }

    call mapped_hic_fragments {
        input:
            bwt2pairs = merge_pairs.bwt2pairs,
            HiC_PATH=HiC_PATH,
            BED_FILE = BED_FILE,
            docker_image = docker_image,
            task_cpu = N_cpu,
            task_mem = N_mem
    }


}


# task unzip_bowtie2index {
#     input {
#         File bwt_index = bwt_index
#         String docker_image
#         String REFERENCE_GENOME
#     }
#     runtime {
#         docker: docker_image
#     }
#     command {
#         mkdir ./input
#         unzip ${bwt_index} -d ./input
#     }
#     output {
#         File bowtie2idx = "./input"
#         File bowtie2idx_makefile = "./input/make_${REFERENCE_GENOME}.sh"
#     }
# }

task bowtie_global_mapping {
    input {
        String sample_name
        String REFERENCE_GENOME
        String BOWTIE2_GLOBAL_OPTIONS
        String odir
        String bwt_cpu
        String BOWTIE2_IDX
        File SAMPLE_R1
        File SAMPLE_R2
        String docker_image
        String task_mem
        File bowtie_idx_1
        File bowtie_idx_2
        File bowtie_idx_3
        File bowtie_idx_4
        File bowtie_idx_re1
        File bowtie_idx_re2
        File bowtie_makesh

    }
    String ldir = odir+"/logs"
    String bowtie_idx_path = "../inputs/1446856593/"+REFERENCE_GENOME
    #String odir = odir+"/bowtie_results/bwt2_global/dixon_2M"
    command {
        mkdir ./output
        mkdir ./output/logs
        date > ./output/logs/time.log
        echo "##HiC-Pro mapping" > ./output/logs/${sample_name}_r1_bowtie2.log
        echo "##HiC-Pro mapping" > ./output/logs/${sample_name}_r2_bowtie2.log
        bowtie2 ${BOWTIE2_GLOBAL_OPTIONS} --un ./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap.fastq --rg-id BMG --rg SM:${sample_name}_r1 -p ${bwt_cpu} -x ${bowtie_idx_path} -U ${SAMPLE_R1} 2>> ./output/logs/${sample_name}_r1_bowtie2.log | samtools view -F 4 -bS - > ./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.bam
        bowtie2 ${BOWTIE2_GLOBAL_OPTIONS} --un ./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap.fastq --rg-id BMG --rg SM:${sample_name}_r2 -p ${bwt_cpu} -x ${bowtie_idx_path} -U ${SAMPLE_R2} 2>> ./output/logs/${sample_name}_r2_bowtie2.log | samtools view -F 4 -bS - > ./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.bam
        rm -r ${REFERENCE_GENOME}
        date >> ./output/logs/time.log
    }
    runtime {
        docker: docker_image
        cpu : bwt_cpu
        memory : task_mem + "GB"
    }
    output {
        File gol_unmap = "./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap.fastq"
        File gol_mapped = "./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.bam"
        File gol_log = "./output/logs/${sample_name}_r1_bowtie2.log"
        File gol_unmap2 = "./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap.fastq"
        File gol_mapped2 = "./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.bam"
        File gol_log2 = "./output/logs/${sample_name}_r2_bowtie2.log"
    }
}

task bowtie_local_trimming {
    input {
        File gol_unmap_fastq1
        File gol_unmap_fastq2
        String HiC_PATH
        String cutsite
        String sample_name
        String REFERENCE_GENOME
        String docker_image
        Int task_cpu
        Int task_mem
    }
    runtime {
        docker :docker_image
        cpu : task_cpu
        memory : task_mem + "GB"
    }
    command {
        mkdir output
        mkdir output/logs
        date > ./output/logs/time.log
        ${HiC_PATH}/scripts/cutsite_trimming --fastq ${gol_unmap_fastq1} --cutsite ${cutsite} --out ./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap_trimmed.fastq > ./output/logs/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap_readsTrimming.log 2>&1
        ${HiC_PATH}/scripts/cutsite_trimming --fastq ${gol_unmap_fastq2} --cutsite ${cutsite} --out ./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap_trimmed.fastq > ./output/logs/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap_readsTrimming.log 2>&1
        date >> ./output/logs/time.log
    }
    output {
        File unmap_trimmed_r1="./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap_trimmed.fastq"
        File unmap_trimmed_r2="./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap_trimmed.fastq"
        File trim_log1="./output/logs/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap_readsTrimming.log"
        File trim_log2="./output/logs/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap_readsTrimming.log"
    }
}

task bowtie_local_mapping {
    input {
        File unmap_trimmed_r1
        File unmap_trimmed_r2
        File bowtie_idx_1
        File bowtie_idx_2
        File bowtie_idx_3
        File bowtie_idx_4
        File bowtie_idx_re1
        File bowtie_idx_re2
        File bowtie_makesh


        String BOWTIE2_LOCAL_OPTIONS

        String REFERENCE_GENOME
        String sample_name

        String docker_image

        Int task_cpu
        Int task_mem

    }
    runtime {
        docker :docker_image
        cpu : task_cpu
        memory :task_mem +"GB"
    }
    String bowtie_idx_path = "../inputs/1446856593/"+REFERENCE_GENOME
    command {
        mkdir ./output
        mkdir ./output/logs
        date > ./output/logs/time.log
        bowtie2 ${BOWTIE2_LOCAL_OPTIONS} --rg-id BML --rg SM:${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap -p ${task_cpu} -x ${bowtie_idx_path} -U ${unmap_trimmed_r1} 2>> ./output/logs/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap_bowtie2.log | samtools view -bS - > ./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap_bwt2loc.bam
        bowtie2 ${BOWTIE2_LOCAL_OPTIONS} --rg-id BML --rg SM:${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap -p ${task_cpu} -x ${bowtie_idx_path} -U ${unmap_trimmed_r2} 2>> ./output/logs/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap_bowtie2.log | samtools view -bS - > ./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap_bwt2loc.bam
        rm -r ${REFERENCE_GENOME}
        date >> ./output/logs/time.log
    }

     output {
        File local_mapped_r1 = "./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap_bwt2loc.bam"
        File local_mapped_r2 = "./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap_bwt2loc.bam"
        File local_mapped_r1_log = "./output/logs/${sample_name}_r1_${REFERENCE_GENOME}.bwt2glob.unmap_bowtie2.log"
        File local_mapped_r2_log = "./output/logs/${sample_name}_r2_${REFERENCE_GENOME}.bwt2glob.unmap_bowtie2.log"
     }
}

task mapping_combine {
    input {
        File global_mapped_r1
        File global_mapped_r2
        File local_mapped_r1
        File local_mapped_r2
        String sample_name
        String REFERENCE_GENOME
        String docker_image
        Int task_cpu
        Int task_mem
    }
    runtime {
        docker : docker_image
        cpu :task_cpu
        memory :task_mem + "GB"
    }
    command {
        mkdir ./output
        mkdir ./output/logs
        mkdir ./tmp
        date > ./output/logs/time.log
        samtools merge -@ 20 -n -f ./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2merged.bam ${global_mapped_r1} ${local_mapped_r1}
        samtools merge -@ 20 -n -f ./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2merged.bam ${global_mapped_r2} ${local_mapped_r2}
        samtools sort -@ 20 -n -T tmp/${sample_name}_r1_${REFERENCE_GENOME} -o ./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2merged.sorted.bam ./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2merged.bam
        samtools sort -@ 20 -n -T tmp/${sample_name}_r2_${REFERENCE_GENOME} -o ./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2merged.sorted.bam ./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2merged.bam
        mv ./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2merged.sorted.bam ./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2merged.bam
        mv ./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2merged.sorted.bam ./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2merged.bam
        date >> ./output/logs/time.log
    }

    output {
        File bwt2merged_r1 = "./output/${sample_name}_r1_${REFERENCE_GENOME}.bwt2merged.bam"
        File bwt2merged_r2 = "./output/${sample_name}_r2_${REFERENCE_GENOME}.bwt2merged.bam"
    }
}

task merge_pairs {
    input {
        File merged_r1
        File merged_r2
        String HiC_PATH
        String sample_name
        String REFERENCE_GENOME
        String docker_image
        Int task_cpu
        Int task_mem
    }
    runtime {
        docker : docker_image
        cpu : task_cpu
        memory : task_mem +"GB"
    }
    command {
        mkdir ./output
        mkdir ./output/logs
        date > ./output/logs/time.log
        python ${HiC_PATH}/scripts/mergeSAM.py -q 0 -t -v -f ${merged_r1} -r ${merged_r2} -o ./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.bam
        date >> ./output/logs/time.log
    }
    output {
        File bwt2pairs = "./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.bam"
    }
}

task mapped_hic_fragments {
    input {
        File bwt2pairs
        String HiC_PATH
        File BED_FILE
        String docker_image
        Int task_cpu
        Int task_mem
    }
    runtime {
        docker : docker_image
        cpu : task_cpu
        memory : task_mem + "GB"
    }
    command {
        mkdir ./output
        mkdir ./output/logs
        date > ./output/logs/time.log
        python ${HiC_PATH}/scripts/mapped_2hic_fragments.py -v -S -t 100 -m 100000 -s 100 -l 600 -a -f ${BED_FILE} -r ${bwt2pairs} -o ./output
        date >> ./output/logs/time.log
    }
    output {
        File bam_validpairs = "./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs_interaction.bam"
        File validPairs =  "./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.validPairs"
        File DEPairs = "./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.DEPairs"
        File DumpedPairs = "./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.DumpPairs"
        File FiltPairs = "./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.FiltPairs"
        File REPairs = "./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.REPairs"
        File SinglePairs ="./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.SinglePairs"
        File SCPairs = "./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.SCPairs"
        File RSstat = "./output/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.RSstat"
    }

}
