version 1.0
workflow HiCPro {
    #input sample files
    input {

    File input_sample1_r1
    File input_sample1_r2



    String sample_1_r1_1
    String sample_1_r2_2

    String sample_name

    String BED_FILE_NAME

    Int N_cpu

    #PATH
    String Sample_PATH
    String BOWTIE2_PATH
    String BOWTIE2_IDX_PATH
    String RAWDATA_PATH
    String SAMTOOLS_PATH
    String HiC_PATH
    String PYTHON_PATH
    String BED_FILE_PATH

    String OUTPUT_DATA_PATH

    #OPTIONS
    String BOWTIE2_GLOBAL_OPTIONS
    String BOWTIE2_LOCAL_OPTIONS
    String cutsite

    #   File Config_SYS
    }
    #call system
    call preparation{
       input:
          odir = OUTPUT_DATA_PATH
    }

    call bowtie_global {
        input:
            BOWTIE2_GLOBAL_OPTIONS = BOWTIE2_GLOBAL_OPTIONS,
            sample_1_r1 = sample_1_r1_1,
            sample_1_r2 = sample_1_r2_2,
            REFERENCE_GENOME = sample_name,
            bwt_cpu = N_cpu,
            odir = OUTPUT_DATA_PATH,
            BOWTIE2_IDX = BOWTIE2_IDX_PATH,
            BOWTIE2_PATH = BOWTIE2_PATH,
            SAMTOOLS_PATH =SAMTOOLS_PATH,
            RAWDATA_PATH = RAWDATA_PATH,
            SAMPLE1_R1 = input_sample1_r1,
            SAMPLE1_R2 = input_sample1_r2


    }

    call bowtie_local_trimming {
        input:
            gol_unmap_fastq1=bowtie_global.gol_unmap,
            gol_unmap_fastq2=bowtie_global.gol_unmap2,
            HiC_PATH=HiC_PATH,
            cutsite=cutsite,
            sample_1_r1 =sample_1_r1_1,
            sample_1_r2 =sample_1_r2_2,
            REFERENCE_GENOME=sample_name
    }

    call bowtie_local_mapping{
        input:
            unmap_trimmed_r1=bowtie_local_trimming.unmap_trimmed_r1,
            unmap_trimmed_r2=bowtie_local_trimming.unmap_trimmed_r2,
            N_cpu = N_cpu,
            REFERENCE_GENOME = sample_name,
            sample_1_r1 = sample_1_r1_1,
            sample_1_r2 = sample_1_r2_2,
            BOWTIE2_LOCAL_OPTIONS = BOWTIE2_LOCAL_OPTIONS,
            Bowtie2Index = BOWTIE2_IDX_PATH,
            BOWTIE2_PATH =BOWTIE2_PATH,
            SAMTOOLS_PATH = SAMTOOLS_PATH,
            odir = OUTPUT_DATA_PATH
    }
    
    call mapping_combine {
        input:
            SAMTOOLS_PATH = SAMTOOLS_PATH,
            odir = OUTPUT_DATA_PATH,
            global_mapped_r1 = bowtie_global.gol_mapped,
            global_mapped_r2 = bowtie_global.gol_mapped2,
            local_mapped_r1 = bowtie_local_mapping.local_mapped_r1,
            local_mapped_r2 = bowtie_local_mapping.local_mapped_r2,
            sample_1_r1 = sample_1_r1_1,
            sample_1_r2 = sample_1_r2_2,
            REFERENCE_GENOME = sample_name
    }
    call merge_pairs {
        input:
            merged_r1 = mapping_combine.bwt2merged_r1,
            merged_r2 = mapping_combine.bwt2merged_r2,
            REFERENCE_GENOME = sample_name,
            sample_name = "SRR400264_00",
            odir = OUTPUT_DATA_PATH,
            PYTHON_PATH = PYTHON_PATH,
            HiC_PATH = HiC_PATH
    }

    call mapped_hic_fragments {
        input:
            bwt2pairs = merge_pairs.bwt2pairs,
            PYTHON_PATH =PYTHON_PATH,
            HiC_PATH=HiC_PATH,
            odir = OUTPUT_DATA_PATH,
            BED_FILE_NAME = BED_FILE_NAME,
            BED_FILE_PATH = BED_FILE_PATH
    }
}


task preparation {
    input {
        String odir
    }
    String LOGS_DIR = "logs"
    command {
        mkdir -p ${odir}/${LOGS_DIR}
    }
}

task bowtie_global {
    input {
        String sample_1_r1
        String sample_1_r2
        String REFERENCE_GENOME
        String BOWTIE2_GLOBAL_OPTIONS
        String odir
        Int bwt_cpu
        String BOWTIE2_IDX
        String BOWTIE2_PATH
        String SAMTOOLS_PATH
        String RAWDATA_PATH
        File SAMPLE1_R1
        File SAMPLE1_R2

    }
    String ldir = odir+"/logs/dixon2M"
    #String odir = odir+"/bowtie_results/bwt2_global/dixon_2M"
    command {
        echo "##HiC-Pro mapping" > ${ldir}/${sample_1_r1}_bowtie2.log
        echo "##HiC-Pro mapping" > ${ldir}/${sample_1_r2}_bowtie2.log
        ${BOWTIE2_PATH} ${BOWTIE2_GLOBAL_OPTIONS} --un ${odir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap.fastq --rg-id BMG --rg SM:${sample_1_r1} -p ${bwt_cpu} -x ${BOWTIE2_IDX}/hg19 -U ${SAMPLE1_R1} 2>> ${ldir}/${sample_1_r1}_bowtie2.log | ${SAMTOOLS_PATH} view -F 4 -bS - > ${odir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.bam
        ${BOWTIE2_PATH} ${BOWTIE2_GLOBAL_OPTIONS} --un ${odir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap.fastq --rg-id BMG --rg SM:${sample_1_r2} -p ${bwt_cpu} -x ${BOWTIE2_IDX}/hg19 -U ${SAMPLE1_R2} 2>> ${ldir}/${sample_1_r2}_bowtie2.log | ${SAMTOOLS_PATH} view -F 4 -bS - > ${odir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.bam
    }

    output {
        File gol_unmap = "${odir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap.fastq"
        File gol_mapped = "${odir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.bam"
        File gol_log = "${ldir}/${sample_1_r1}_bowtie2.log"
        File gol_unmap2 = "${odir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap.fastq"
        File gol_mapped2 = "${odir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.bam"
        File gol_log2 = "${ldir}/${sample_1_r2}_bowtie2.log"
    }
}

task bowtie_local_trimming {
    input {
        File gol_unmap_fastq1
        File gol_unmap_fastq2
        String HiC_PATH
        String cutsite
        String sample_1_r1
        String sample_1_r2
        String REFERENCE_GENOME
    }
    command {
        ${HiC_PATH}/scripts/cutsite_trimming --fastq ${gol_unmap_fastq1} --cutsite ${cutsite} --out ${HiC_PATH}/output/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap_trimmed.fastq > ${HiC_PATH}/output/logs/dixon2M/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap_readsTrimming.log 2>&1
        ${HiC_PATH}/scripts/cutsite_trimming --fastq ${gol_unmap_fastq2} --cutsite ${cutsite} --out ${HiC_PATH}/output/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap_trimmed.fastq > ${HiC_PATH}/output/logs/dixon2M/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap_readsTrimming.log 2>&1

    }
    output {
        File unmap_trimmed_r1="${HiC_PATH}/output/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap_trimmed.fastq"
        File unmap_trimmed_r2="${HiC_PATH}/output/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap_trimmed.fastq"
        File trim_log1="${HiC_PATH}/output/logs/dixon2M/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap_readsTrimming.log"
        File trim_log2="${HiC_PATH}/output/logs/dixon2M/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap_readsTrimming.log"
    }
}

task bowtie_local_mapping {
    input {
        File unmap_trimmed_r1
        File unmap_trimmed_r2
        String BOWTIE2_PATH
        String BOWTIE2_LOCAL_OPTIONS
        String Bowtie2Index
        Int N_cpu
        String REFERENCE_GENOME
        String sample_1_r1
        String sample_1_r2
        String SAMTOOLS_PATH
        String odir
        

    }
    String ldir = odir+"/logs/dixon2M"
    command {
        ${BOWTIE2_PATH} ${BOWTIE2_LOCAL_OPTIONS} --rg-id BML --rg SM:${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap -p ${N_cpu} -x ${Bowtie2Index}/${REFERENCE_GENOME} -U ${unmap_trimmed_r1} 2>> ${ldir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap_bowtie2.log | ${SAMTOOLS_PATH} view -bS - > ${odir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap_bwt2loc.bam
        ${BOWTIE2_PATH} ${BOWTIE2_LOCAL_OPTIONS} --rg-id BML --rg SM:${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap -p ${N_cpu} -x ${Bowtie2Index}/${REFERENCE_GENOME} -U ${unmap_trimmed_r2} 2>> ${ldir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap_bowtie2.log | ${SAMTOOLS_PATH} view -bS - > ${odir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap_bwt2loc.bam

    }

     output {
        File local_mapped_r1 = "${odir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap_bwt2loc.bam"
        File local_mapped_r2 = "${odir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap_bwt2loc.bam"
        File local_mapped_r1_log = "${ldir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2glob.unmap_bowtie2.log"
        File local_mapped_r2_log = "${ldir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2glob.unmap_bowtie2.log"
     }
}

task mapping_combine {
    input {
        String SAMTOOLS_PATH
        String odir
        File global_mapped_r1
        File global_mapped_r2
        File local_mapped_r1
        File local_mapped_r2
        String sample_1_r1
        String sample_1_r2
        String REFERENCE_GENOME
    }
    
    command {
        ${SAMTOOLS_PATH} merge -@ 20 -n -f ${odir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2merged.bam ${global_mapped_r1} ${local_mapped_r1}
        ${SAMTOOLS_PATH} merge -@ 20 -n -f ${odir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2merged.bam ${global_mapped_r2} ${local_mapped_r2}
        ${SAMTOOLS_PATH} sort -@ 20 -n -T tmp/${sample_1_r1}_${REFERENCE_GENOME} -o ${odir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2merged.sorted.bam ${odir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2merged.bam
        ${SAMTOOLS_PATH} sort -@ 20 -n -T tmp/${sample_1_r2}_${REFERENCE_GENOME} -o ${odir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2merged.sorted.bam ${odir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2merged.bam
        mv ${odir}/${sample_1_r1}_${REFERENCE_GENOME}.bwt2merged.sorted.bam ${sample_1_r1}_${REFERENCE_GENOME}.bwt2merged.bam
        mv ${odir}/${sample_1_r2}_${REFERENCE_GENOME}.bwt2merged.sorted.bam ${sample_1_r2}_${REFERENCE_GENOME}.bwt2merged.bam

    }

    output {
        File bwt2merged_r1 = "${sample_1_r1}_${REFERENCE_GENOME}.bwt2merged.bam"
        File bwt2merged_r2 = "${sample_1_r2}_${REFERENCE_GENOME}.bwt2merged.bam"
    }
}

task merge_pairs {
    input {
        File merged_r1
        File merged_r2
        String PYTHON_PATH
        String HiC_PATH
        String odir 
        String sample_name
        String REFERENCE_GENOME
    }
    command {
        ${PYTHON_PATH} ${HiC_PATH}/scripts/mergeSAM.py -q 0 -t -v -f merged_r1 -r ${merged_r2} -o ${odir}/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.bam

    }
    output {
        File bwt2pairs = "${odir}/${sample_name}_${REFERENCE_GENOME}.bwt2pairs.bam"
    }
}

task mapped_hic_fragments {
    input {
        File bwt2pairs
        String PYTHON_PATH
        String odir
        String HiC_PATH
        String BED_FILE_NAME
        String BED_FILE_PATH
    }
    command {
        ${PYTHON_PATH} ${HiC_PATH}/scripts/mapped_2hic_fragments.py -v -S -t 100 -m 100000 -s 100 -l 600 -a -f ${BED_FILE_PATH}/${BED_FILE_NAME} -r ${bwt2pairs} -o ${odir}/result
        
    }
}