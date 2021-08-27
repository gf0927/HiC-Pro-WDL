version 1.0

workflow Hic_Docker {
    input {
        #Setting
        String dockerImage
        String bowtie2GlobalOptions
        String bowtie2LocalOptions
        String nCPU
        String nMem
        String cutsite
        Int binSize

        #PATH
        String hicPath

        #inputs
        String sampleName
        String referenceGenome
        File bowtieIndexPath
        File sampleR1
        File sampleR2
        File bedFile
        File tableFile
    }

    call Bowtie_Global_Mapping {
        input:
            sampleName  = sampleName,
            refGenome = referenceGenome,
            bowtie2GlobalOpts = bowtie2GlobalOptions,
            bowtieIndexPath = bowtieIndexPath,
            sampleR1 = sampleR1,
            sampleR2 = sampleR2,
            dockerCPU = nCPU,
            dockerImage = dockerImage,
            dockerMEM = nMem
    }

    call Bowtie_Local_Trimming {
        input:
            globalUnmappedFastq1=Bowtie_Global_Mapping.gol_unmap,
            globalUnmappedFastq2=Bowtie_Global_Mapping.gol_unmap2,
            hicPath=hicPath,
            cutsite=cutsite,
            sampleName = sampleName,
            refGenome=referenceGenome,
            dockerImage = dockerImage,
            dockerCPU = nCPU,
            dockerMEM = nMem
    }

    call Bowtie_Local_Mapping{
        input:
            unmapTrimmedR1 = Bowtie_Local_Trimming.unmapTrimmedR1,
            unmapTrimmedR2 = Bowtie_Local_Trimming.unmapTrimmedR2,
            bowtieIndexPath = bowtieIndexPath,
            sampleName  = sampleName,
            dockerCPU = nCPU,
            refGenome = referenceGenome,
            bowtie2LocalOptions = bowtie2LocalOptions,
            dockerImage = dockerImage,
            dockerMEM = nMem
    }

    call Mapping_Combine {
        input:
            globalMappedR1 = Bowtie_Global_Mapping.gol_mapped,
            globalMappedR2 = Bowtie_Global_Mapping.gol_mapped2,
            localMappedR1 = Bowtie_Local_Mapping.localMappedR1,
            localMappedR2 = Bowtie_Local_Mapping.localMappedR2,
            sampleName = sampleName,
            refGenome = referenceGenome,
            dockerImage = dockerImage,
            dockerCPU = nCPU,
            dockerMEM = nMem
    }

    call Mapping_Stat {
        input:
            mappedGlob_R1 = Bowtie_Global_Mapping.gol_mapped,
            mappedGlob_R2 = Bowtie_Global_Mapping.gol_mapped2,
            mappedLoc_R1 = Bowtie_Local_Mapping.localMappedR1,
            mappedLoc_R2 = Bowtie_Local_Mapping.localMappedR2,
            mappedMerged_R1 = Mapping_Combine.bwt2merged_R1,
            mappedMerged_R2 = Mapping_Combine.bwt2merged_R2,
            dockerImage = dockerImage,
            dockerCPU = nCPU,
            dockerMEM = nMem,
            sampleName = sampleName,
            refGenome = referenceGenome
    }

    call Merge_Pairs {
        input:
            merged_R1 = Mapping_Combine.bwt2merged_R1,
            merged_R2 = Mapping_Combine.bwt2merged_R2,
            refGenome = referenceGenome,
            sampleName = sampleName,
            hicPath = hicPath,
            dockerImage = dockerImage,
            dockerCPU = nCPU,
            dockerMEM = nMem
    }

    call Mapped_Hic_Fragments {
        input:
            bowtieMappedPairs = Merge_Pairs.bowtieMappedPairs,
            hicPath=hicPath,
            bedFile = bedFile,
            dockerImage = dockerImage,
            dockerCPU = nCPU,
            dockerMEM = nMem,
            sampleName = sampleName,
            refGenome = referenceGenome
    }

    call Build_Matrix {
        input:
            allValidPairs = Mapped_Hic_Fragments.validPairs,
            tableFile = tableFile,
            hicPath = hicPath,
            dockerImage = dockerImage,
            dockerCPU = nCPU,
            dockerMEM = nMem,
            binSize = binSize,
            sampleName = sampleName

    }
}

task Bowtie_Global_Mapping {
    input {
        String sampleName
        String refGenome
        String bowtie2GlobalOpts
        String dockerCPU
        File bowtieIndexPath
        File sampleR1
        File sampleR2
        String dockerImage
        String dockerMEM
    }

    command {
        set -e -o pipefail
        mkdir ./output
        mkdir ./output/logs
        date > ./output/logs/time.log
        echo "##HiC-Pro mapping" > ./output/logs/~{sampleName}_r1_bowtie2.log
        echo "##HiC-Pro mapping" > ./output/logs/~{sampleName}_r2_bowtie2.log
        bowtie2 ~{bowtie2GlobalOpts} \
            --un ./output/~{sampleName}_r1_~{refGenome}.bwt2glob.unmap.fastq \
            --rg-id BMG \
            --rg SM:~{sampleName}_r1 \
            -p ~{dockerCPU} \
            -x ~{bowtieIndexPath}/~{refGenome} \
            -U ~{sampleR1} \
            2>> ./output/logs/~{sampleName}_r1_bowtie2.log | \
        samtools view -F 4 -bS - > ./output/~{sampleName}_r1_~{refGenome}.bwt2glob.bam
        bowtie2 ~{bowtie2GlobalOpts} \
            --un ./output/~{sampleName}_r2_~{refGenome}.bwt2glob.unmap.fastq \
            --rg-id BMG \
            --rg SM:~{sampleName}_r2 \
            -p ~{dockerCPU} \
            -x ~{bowtieIndexPath}/~{refGenome} \
            -U ~{sampleR2} \
            2>> ./output/logs/~{sampleName}_r2_bowtie2.log | \
        samtools view -F 4 -bS - > ./output/~{sampleName}_r2_~{refGenome}.bwt2glob.bam
        date >> ./output/logs/time.log
    }

    runtime {
        docker: dockerImage
        cpu : dockerCPU
        memory : dockerMEM + "GB"
    }

    output {
        File gol_unmap = "./output/~{sampleName}_r1_~{refGenome}.bwt2glob.unmap.fastq"
        File gol_mapped = "./output/~{sampleName}_r1_~{refGenome}.bwt2glob.bam"
        File gol_log = "./output/logs/~{sampleName}_r1_bowtie2.log"
        File gol_unmap2 = "./output/~{sampleName}_r2_~{refGenome}.bwt2glob.unmap.fastq"
        File gol_mapped2 = "./output/~{sampleName}_r2_~{refGenome}.bwt2glob.bam"
        File gol_log2 = "./output/logs/~{sampleName}_r2_bowtie2.log"
    }
}

task Bowtie_Local_Trimming {
    input {
        File globalUnmappedFastq1
        File globalUnmappedFastq2
        String hicPath
        String cutsite
        String sampleName
        String refGenome
        String dockerImage
        String dockerCPU
        String dockerMEM
    }

    runtime {
        docker : dockerImage
        cpu : dockerCPU
        memory : dockerMEM + "GB"
    }

    command {
        set -e -o pipefail
        mkdir output
        mkdir output/logs
        date > ./output/logs/time.log
        ~{hicPath}/scripts/cutsite_trimming \
            --fastq ~{globalUnmappedFastq1} \
            --cutsite ~{cutsite} \
            --out ./output/~{sampleName}_r1_~{refGenome}.bwt2glob.unmap_trimmed.fastq > ./output/logs/~{sampleName}_r1_~{refGenome}.bwt2glob.unmap_readsTrimming.log 2>&1
        ~{hicPath}/scripts/cutsite_trimming \
            --fastq ~{globalUnmappedFastq2} --cutsite ~{cutsite} \
            --out ./output/~{sampleName}_r2_~{refGenome}.bwt2glob.unmap_trimmed.fastq > ./output/logs/~{sampleName}_r2_~{refGenome}.bwt2glob.unmap_readsTrimming.log 2>&1
        date >> ./output/logs/time.log
    }

    output {
        File unmapTrimmedR1="./output/~{sampleName}_r1_~{refGenome}.bwt2glob.unmap_trimmed.fastq"
        File unmapTrimmedR2="./output/~{sampleName}_r2_~{refGenome}.bwt2glob.unmap_trimmed.fastq"
        File trim_log1="./output/logs/~{sampleName}_r1_~{refGenome}.bwt2glob.unmap_readsTrimming.log"
        File trim_log2="./output/logs/~{sampleName}_r2_~{refGenome}.bwt2glob.unmap_readsTrimming.log"
    }
}

task Bowtie_Local_Mapping {
    input {
        File unmapTrimmedR1
        File unmapTrimmedR2
        File bowtieIndexPath
        String bowtie2LocalOptions
        String refGenome
        String sampleName
        String dockerImage
        String dockerCPU
        String dockerMEM
    }

    runtime {
        docker : dockerImage
        cpu : dockerCPU
        memory : dockerMEM +"GB"
    }

    command {
        set -e -o pipefail
        mkdir ./output
        mkdir ./output/logs
        date > ./output/logs/time.log
        bowtie2 ~{bowtie2LocalOptions} \
            --rg-id BML \
            --rg SM:~{sampleName}_r1_~{refGenome}.bwt2glob.unmap \
            -p ~{dockerCPU} \
            -x ~{bowtieIndexPath}/~{refGenome} \
            -U ~{unmapTrimmedR1} \
            2>> ./output/logs/~{sampleName}_r1_~{refGenome}.bwt2glob.unmap_bowtie2.log | \
        samtools view -bS - > ./output/~{sampleName}_r1_~{refGenome}.bwt2glob.unmap_bwt2loc.bam
        bowtie2 ~{bowtie2LocalOptions} \
            --rg-id BML \
            --rg SM:~{sampleName}_r2_~{refGenome}.bwt2glob.unmap \
            -p ~{dockerCPU} \
            -x ~{bowtieIndexPath}/~{refGenome} \
            -U ~{unmapTrimmedR2} \
            2>> ./output/logs/~{sampleName}_r2_~{refGenome}.bwt2glob.unmap_bowtie2.log | \
        samtools view -bS - > ./output/~{sampleName}_r2_~{refGenome}.bwt2glob.unmap_bwt2loc.bam
        date >> ./output/logs/time.log
    }

    output {
        File localMappedR1 = "./output/~{sampleName}_r1_~{refGenome}.bwt2glob.unmap_bwt2loc.bam"
        File localMappedR2 = "./output/~{sampleName}_r2_~{refGenome}.bwt2glob.unmap_bwt2loc.bam"
        File localMappedR1_log = "./output/logs/~{sampleName}_r1_~{refGenome}.bwt2glob.unmap_bowtie2.log"
        File localMappedR2_log = "./output/logs/~{sampleName}_r2_~{refGenome}.bwt2glob.unmap_bowtie2.log"
     }
}

task Mapping_Combine {
    input {
        File globalMappedR1
        File globalMappedR2
        File localMappedR1
        File localMappedR2
        String sampleName
        String refGenome
        String dockerImage
        String dockerCPU
        String dockerMEM
    }

    runtime {
        docker : dockerImage
        cpu :dockerCPU
        memory :dockerMEM + "GB"
    }

    command {
        set -e -o pipefail
        mkdir ./output
        mkdir ./output/logs
        mkdir ./tmp
        date > ./output/logs/time.log
        samtools merge \
            -@ ~{dockerCPU} \
            -n \
            -f ./output/~{sampleName}_r1_~{refGenome}.bwt2merged.bam ~{globalMappedR1} ~{localMappedR1}
        samtools merge \
            -@ ~{dockerCPU} \
            -n \
            -f ./output/~{sampleName}_r2_~{refGenome}.bwt2merged.bam ~{globalMappedR2} ~{localMappedR2}
        samtools sort \
            -@ ~{dockerCPU} \
            -n \
            -T tmp/~{sampleName}_r1_~{refGenome} \
            -o ./output/~{sampleName}_r1_~{refGenome}.bwt2merged.sorted.bam ./output/~{sampleName}_r1_~{refGenome}.bwt2merged.bam
        samtools sort \
            -@ ~{dockerCPU} \
            -n \
            -T tmp/~{sampleName}_r2_~{refGenome} \
            -o ./output/~{sampleName}_r2_~{refGenome}.bwt2merged.sorted.bam ./output/~{sampleName}_r2_~{refGenome}.bwt2merged.bam
        mv ./output/~{sampleName}_r1_~{refGenome}.bwt2merged.sorted.bam ./output/~{sampleName}_r1_~{refGenome}.bwt2merged.bam
        mv ./output/~{sampleName}_r2_~{refGenome}.bwt2merged.sorted.bam ./output/~{sampleName}_r2_~{refGenome}.bwt2merged.bam
        date >> ./output/logs/time.log
    }

    output {
        File bwt2merged_R1 = "./output/~{sampleName}_r1_~{refGenome}.bwt2merged.bam"
        File bwt2merged_R2 = "./output/~{sampleName}_r2_~{refGenome}.bwt2merged.bam"
    }
}

task Mapping_Stat {
    input {
        File mappedGlob_R1
        File mappedGlob_R2
        File mappedLoc_R1
        File mappedLoc_R2
        File mappedMerged_R1
        File mappedMerged_R2
        String dockerImage
        String dockerCPU
        String dockerMEM
        String sampleName
        String refGenome
    }

    runtime {
        docker : dockerImage
        cpu : dockerCPU
        memory : dockerMEM + "GB"
    }
    
    command {
        set -e -o pipefail
        mkdir ./output
        echo "##~{sampleName}_R1_~{refGenome}.mapstat" > ./output/~{sampleName}_R1_~{refGenome}.mapstat
        echo -n -e "total\t" >> ./output/~{sampleName}_R1_~{refGenome}.mapstat
        samtools view -c ~{mappedMerged_R1} >> ./output~{sampleName}_R1_~{refGenome}.mapstat
        echo -n -e "mapped\t" >> ./output/~{sampleName}_R1_~{refGenome}.mapstat
        samtools view -c -F 4 ~{mappedMerged_R1} >> ./output/~{sampleName}_R1_~{refGenome}.mapstat
        echo -n -e "global\t" >> ./output~{sampleName}_R1_~{refGenome}.mapstat
        samtools view -c -F 4 ~{mappedGlob_R1} >> ./output/~{sampleName}_R1_~{refGenome}.mapstat
        echo -n -e "local\t" >> ./output/~{sampleName}_R1_~{refGenome}.mapstat
        samtools view -c -F 4 ~{mappedLoc_R1} >> ./output/~{sampleName}_R1_~{refGenome}.mapstat

        echo "##~{sampleName}_R2_~{refGenome}.mapstat" > ./output/~{sampleName}_R2_~{refGenome}.mapstat
        echo -n -e "total\t" >> ./output/~{sampleName}_R2_~{refGenome}.mapstat
        samtools view -c ~{mappedMerged_R2} >> ./output/~{sampleName}_R2_~{refGenome}.mapstat
        echo -n -e "mapped\t" >> ./output/~{sampleName}_R2_~{refGenome}.mapstat
        samtools view -c -F 4 ~{mappedMerged_R2} >> ./output/~{sampleName}_R2_~{refGenome}.mapstat
        echo -n -e "global\t" >> ./output/~{sampleName}_R2_~{refGenome}.mapstat
        samtools view -c -F 4 ~{mappedGlob_R2} >> ./output/~{sampleName}_R2_~{refGenome}.mapstat
        echo -n -e "local\t" >> ./output/~{sampleName}_R2_~{refGenome}.mapstat
        samtools view -c -F 4 ~{mappedLoc_R2} >> ./output/~{sampleName}_R2_~{refGenome}.mapstat

    }

    output {
        File R1_mapstat = "./output/~{sampleName}_R1_~{refGenome}.mapstat"
        File R2_mapstat = "./output/~{sampleName}_R2_~{refGenome}.mapstat"
    }
}


task Merge_Pairs {
    input {
        File merged_R1
        File merged_R2
        String hicPath
        String sampleName
        String refGenome
        String dockerImage
        Int dockerCPU
        Int dockerMEM
    }

    runtime {
        docker : dockerImage
        cpu : dockerCPU
        memory : dockerMEM +"GB"
    }

    command {
        set -e -o pipefail
        mkdir ./output
        mkdir ./output/logs
        date > ./output/logs/time.log
        python3 ~{hicPath}/scripts/mergeSAM.py \
            -q 0 \
            -t \
            -v \
            -f ~{merged_R1} \
            -r ~{merged_R2} \
            -o ./output/~{sampleName}_~{refGenome}.bwt2pairs.bam
        date >> ./output/logs/time.log
    }

    output {
        File bowtieMappedPairs = "./output/~{sampleName}_~{refGenome}.bwt2pairs.bam"
    }
}

task Mapped_Hic_Fragments {
    input {
        File bowtieMappedPairs
        String hicPath
        File bedFile
        String dockerImage
        String dockerCPU
        String dockerMEM
        String sampleName
        String refGenome
    }

    runtime {
        docker : dockerImage
        cpu : dockerCPU
        memory : dockerMEM + "GB"
    }

    command {
        set -e -o pipefail
        mkdir ./output
        mkdir ./output/logs
        date > ./output/logs/time.log
        python3 ~{hicPath}/scripts/mapped_2hic_fragments.py \
            -v \
            -S \
            -t 100 \
            -m 100000 \
            -s 100 \
            -l 600 \
            -a \
            -f ~{bedFile} \
            -r ~{bowtieMappedPairs} \
            -o ./output
        LANG=en; sort \
            -T tmp \
            -k2,2V -k3,3n -k5,5V -k6,6n \
            -o ./output/~{sampleName}_~{refGenome}.bwt2pairs.validPairs \
            ./output/~{sampleName}_~{refGenome}.bwt2pairs.validPairs
        date >> ./output/logs/time.log
    }

    output {
        File bam_validpairs = "./output/~{sampleName}_~{refGenome}.bwt2pairs_interaction.bam"
        File validPairs =  "./output/~{sampleName}_~{refGenome}.bwt2pairs.validPairs"
        File DEPairs = "./output/~{sampleName}_~{refGenome}.bwt2pairs.DEPairs"
        File DumpedPairs = "./output/~{sampleName}_~{refGenome}.bwt2pairs.DumpPairs"
        File FiltPairs = "./output/~{sampleName}_~{refGenome}.bwt2pairs.FiltPairs"
        File REPairs = "./output/~{sampleName}_~{refGenome}.bwt2pairs.REPairs"
        File SinglePairs ="./output/~{sampleName}_~{refGenome}.bwt2pairs.SinglePairs"
        File SCPairs = "./output/~{sampleName}_~{refGenome}.bwt2pairs.SCPairs"
        File RSstat = "./output/~{sampleName}_~{refGenome}.bwt2pairs.RSstat"
    }
}

task Build_Matrix {
    input {
        File allValidPairs
        File tableFile
        String hicPath
        String dockerImage
        String dockerCPU
        String dockerMEM
        String sampleName
        Int binSize

    }

    runtime {
        docker : dockerImage
        cpu : dockerCPU
        memory : dockerMEM + "GB"
    }

    command {
        set -e -o pipefail
        mkdir ./output
        cat ~{allValidPairs} | ~{hicPath}/scripts/build_matrix \
            --matrix-format upper \
            --binsize ~{binSize} \
            --chrsizes ~{tableFile} \
            --ifile /dev/stdin \
            --oprefix ./output/~{sampleName}_~{binSize}
    }

    output {
        File matrix = "./output/~{sampleName}_~{binSize}.matrix"
        File abs_bed = "./output/~{sampleName}_~{binSize}_abs.bed"
    }
}