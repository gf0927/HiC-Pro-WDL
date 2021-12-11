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
        String tagR1
        String tagR2

        #PATH
        String hicPath

        #inputs
        String sampleName
        String resName
        String referenceGenome
        File bowtieIndexPath
        Array[Array[File]] inputSamples
        Array[Array[String]] inputSettings

        File bedFile
        File tableFile
        File fakeFile
    }

    scatter (idx in range(length(inputSamples))) {
        call Bowtie_Global_Mapping {
            input:
                sampleName = inputSettings[idx][0],
                refGenome = referenceGenome,
                bowtie2GlobalOpts = bowtie2GlobalOptions,
                bowtieIndexPath = bowtieIndexPath,
                sampleR1 = inputSamples[idx][0],
                sampleR2 = inputSamples[idx][1],
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
                sampleName = inputSettings[idx][0],
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
                sampleName  = inputSettings[idx][0],
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
                sampleName = inputSettings[idx][0],
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
                sampleName = inputSettings[idx][0],
                refGenome = referenceGenome
        }

        call Merge_Pairs {
            input:
                merged_R1 = Mapping_Combine.bwt2merged_R1,
                merged_R2 = Mapping_Combine.bwt2merged_R2,
                refGenome = referenceGenome,
                sampleName = inputSettings[idx][0],
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
                sampleName = inputSettings[idx][0],
                refGenome = referenceGenome
        }

        call Merge_Valid_Interaction {
            input:
                dockerImage = dockerImage,
                dockerCPU = nCPU,
                dockerMEM = nMem,
                resName = resName,
                validPairs = Mapped_Hic_Fragments.validPairs,
                fakeFile = fakeFile
        }

        call Build_Matrix {
            input:
                allValidPairs = Merge_Valid_Interaction.allValidPairs,
                tableFile = tableFile,
                hicPath = hicPath,
                dockerImage = dockerImage,
                dockerCPU = nCPU,
                dockerMEM = nMem,
                binSize = binSize,
                sampleName = inputSettings[idx][0]
        }

        call Making_Plot {
            input:
                sampleName = inputSettings[idx][0],
                refGenome=referenceGenome,
                dockerCPU = nCPU,
                dockerMEM = nMem,
                dockerImage = dockerImage,
                tagR1 = tagR1,
                tagR2 = tagR2,
                hicPath = hicPath,
                rmMulti = "1",
                rmSingle = "1",
                mappingStatR1 = Mapping_Stat.R1_mapstat,
                mappingStatR2 = Mapping_Stat.R2_mapstat,
                pairStat = Merge_Pairs.pairStat,
                rsStat = Mapped_Hic_Fragments.RSstat,
                mergeStat = Merge_Valid_Interaction.mergeStat,
                validPairs = Merge_Valid_Interaction.allValidPairs
        }

        call Ice_Normalization {
            input:
                dockerImage = dockerImage,
                dockerCPU = nCPU,
                dockerMEM = nMem,
                hicPath = hicPath,
                sampleName = inputSettings[idx][0],
                rawMatrix = Build_Matrix.matrix
        }
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
        mkdir -p ./output
        mkdir -p ./output/logs
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
        File timer = "./output/logs/time.log"
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
        mkdir -p ./output
        mkdir -p ./output/logs
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
        File timer = "./output/logs/time.log"
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
        mkdir -p ./output
        mkdir -p ./output/logs
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
        File timer = "./output/logs/time.log"
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
        mkdir -p ./output
        mkdir -p ./output/logs
        mkdir -p ./tmp
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
        File timer = "./output/logs/time.log"
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
        mkdir -p ./output
        mkdir -p ./output/logs
        date > ./output/logs/time.log
        echo "##~{sampleName}_R1_~{refGenome}.mapstat" > ./output/~{sampleName}_R1_~{refGenome}.mapstat
        echo -n -e "total\t" >> ./output/~{sampleName}_R1_~{refGenome}.mapstat
        samtools view -c ~{mappedMerged_R1} >> ./output/~{sampleName}_R1_~{refGenome}.mapstat
        echo -n -e "mapped\t" >> ./output/~{sampleName}_R1_~{refGenome}.mapstat
        samtools view -c -F 4 ~{mappedMerged_R1} >> ./output/~{sampleName}_R1_~{refGenome}.mapstat
        echo -n -e "global\t" >> ./output/~{sampleName}_R1_~{refGenome}.mapstat
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

        date >> ./output/logs/time.log

    }

    output {
        File R1_mapstat = "./output/~{sampleName}_R1_~{refGenome}.mapstat"
        File R2_mapstat = "./output/~{sampleName}_R2_~{refGenome}.mapstat"
        File timer = "./output/logs/time.log"
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
        mkdir -p ./output
        mkdir -p ./output/logs
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
        File pairStat = "./output/~{sampleName}_~{refGenome}.bwt2pairs.pairstat"
        File timer = "./output/logs/time.log"
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

    command <<<
        set -e -o pipefail
        mkdir -p ./output
        mkdir -p ./output/logs
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
    >>>

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
        File timer = "./output/logs/time.log"
    }
}

task Merge_Valid_Interaction {
    input {
        File validPairs
        File fakeFile
        String dockerImage
        String dockerCPU
        String dockerMEM
        String resName
    }

    runtime {
        docker : dockerImage
        cpu : dockerCPU
        memory : dockerMEM + "GB"
    }

    command <<<
        set -e -o pipefail

        cat ~{validPairs}/* > sample.bwt2pairs.validPairs

        mkdir -p ./output
        mkdir -p ./tmp
        mkdir -p ./output/logs
        date > ./output/logs/time.log
        cp sample.bwt2pairs.validPairs ./tmp
        LANG=en; sort -T ./tmp -S 50% -k2,2V -k3,3n -k5,5V -k6,6n -m ./tmp/*.validPairs | \
        awk -F"\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=$2 || c2!=$5 || s1!=$3 || s2!=$6){print;c1=$2;c2=$5;s1=$3;s2=$6}' > ./output/~{resName}.allValidPairs
        echo -e -n "valid_interaction\t" > ./output/~{resName}_allValidPairs.mergestat
        cat sample.bwt2pairs.validPairs | wc -l >> ./output/~{resName}_allValidPairs.mergestat
        echo -e -n "valid_interaction_rmdup\t" >> ./output/~{resName}_allValidPairs.mergestat
        cat ./output/*.allValidPairs | wc -l >> ./output/~{resName}_allValidPairs.mergestat
        awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} \
            $2 == $5 \
            {cis=cis+1; d=$6>$3?$6-$3:$3-$6; \
            if (d<=20000){sr=sr+1}else{lr=lr+1}} \
            $2!=$5{trans=trans+1}\
            END{print "trans_interaction\t"trans"\ncis_interaction\t"cis"\ncis_shortRange\t"sr"\ncis_longRange\t"lr}' \
            ./output/~{resName}.allValidPairs >> ./output/~{resName}_allValidPairs.mergestat
        date >> ./output/logs/time.log
    >>>

    output {
        File allValidPairs = "./output/~{resName}.allValidPairs"
        File mergeStat = "./output/~{resName}_allValidPairs.mergestat"
        File timer = "./output/logs/time.log"
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
        mkdir -p ./output
        mkdir -p ./output/logs
        date > ./output/logs/time.log
        cat ~{allValidPairs} | ~{hicPath}/scripts/build_matrix \
            --matrix-format upper \
            --binsize ~{binSize} \
            --chrsizes ~{tableFile} \
            --ifile /dev/stdin \
            --oprefix ./output/~{sampleName}_~{binSize}
        date >> ./output/logs/time.log
    }

    output {
        File matrix = "./output/~{sampleName}_~{binSize}.matrix"
        File abs_bed = "./output/~{sampleName}_~{binSize}_abs.bed"
        File timer = "./output/logs/time.log"
    }
}

task Making_Plot {
    input {
        String sampleName
        String dockerCPU
        String dockerMEM
        String tagR1
        String tagR2
        String hicPath
        String dockerImage
        String rmMulti #是否給出多端比對的結果，預設為1，不給
        String rmSingle #是否給出單端比對的結果，預設為1，不給
        String refGenome

        File mappingStatR1
        File mappingStatR2
        File pairStat
        File rsStat
        File mergeStat
        File validPairs
    }
    
    command <<<
        set -e -o pipefail

        mkdir -p ./~{sampleName}
        mkdir -p ./temDir
        mkdir -p ./output
        mkdir -p ./output/logs
        date > ./output/logs/time.log
        cp -r ~{mappingStatR1} ./~{sampleName}
        cp -r ~{mappingStatR2} ./~{sampleName}
        cp -r ~{pairStat} ./~{sampleName}
        cp -r ~{rsStat} ./~{sampleName}
        python ~{hicPath}/scripts/merge_statfiles.py -d ./~{sampleName}/~{sampleName}_R1_~{refGenome}.mapstat -p "*_R1*.mapstat" -v > ./~{sampleName}/~{sampleName}.mmapStatR1
        python ~{hicPath}/scripts/merge_statfiles.py -d ./~{sampleName}/~{sampleName}_R2_~{refGenome}.mapstat -p "*_R2*.mapstat" -v > ./~{sampleName}/~{sampleName}.mmapStatR2
        python ~{hicPath}/scripts/merge_statfiles.py -d ./~{sampleName}/~{sampleName}_~{refGenome}.bwt2pairs.pairstat -p "*.pairstat" -v > ./~{sampleName}/~{sampleName}.mPairStat
        python ~{hicPath}/scripts/merge_statfiles.py -d ./~{sampleName}/~{sampleName}_~{refGenome}.bwt2pairs.RSstat -p "*.RSstat" -v > ./~{sampleName}/~{sampleName}.mRSstat

        cp ./~{sampleName}/~{sampleName}.mmapStatR1 ./temDir/~{sampleName}_R1_~{refGenome}.mapstat
        cp ./~{sampleName}/~{sampleName}.mmapStatR2 ./temDir/~{sampleName}_R2_~{refGenome}.mapstat
        cp ./~{sampleName}/~{sampleName}.mPairStat ./temDir/~{sampleName}_~{refGenome}.bwt2pairs.pairstat
        cp ./temDir/~{sampleName}_~{refGenome}.bwt2pairs.pairstat ./output/~{sampleName}.pairstat
        cp ./~{sampleName}/~{sampleName}.mRSstat ./temDir/~{sampleName}_~{refGenome}.bwt2pairs.RSstat
        cp -r ~{mergeStat} ./temDir/~{sampleName}.mergestat
        cp -r ~{validPairs} ./temDir/~{sampleName}.validPairs
        R CMD BATCH \
            --no-save \
            --no-restore \
            "--args picDir='./output' bwtDir='./temDir' sampleName='~{sampleName}' r1tag='~{tagR1}' r2tag='~{tagR2}' " \
            ~{hicPath}/scripts/plot_mapping_portion.R \
            ./output/logs/plot_mapping_portion.Rout

        R CMD BATCH --no-save --no-restore \
            "--args picDir='./output' bwtDir='./temDir' sampleName='~{sampleName}' rmMulti='~{rmMulti}' rmSingle='~{rmSingle}' " \
            ~{hicPath}/scripts/plot_pairing_portion.R \
            ./output/logs/plot_pairing_portion.Rout

        R CMD BATCH --no-save --no-restore \
            "--args picDir='./output' hicDir='./temDir' sampleName='~{sampleName}'" \
            ~{hicPath}/scripts/plot_hic_fragment.R  \
            ./output/logs/plot_hic_fragment.Rout

        R CMD BATCH --no-save --no-restore \
            "--args picDir='./output' hicDir='./temDir' statsDir='./temDir' sampleName='~{sampleName}'" \
            ~{hicPath}/scripts/plot_hic_contacts.R \
            ./output/logs/plot_hic_contacts.Rout
        date >> ./output/logs/time.log
    >>>

    runtime {
        docker : dockerImage
        cpu : dockerCPU
        memory : dockerMEM + "GB"
    }

    output {
        File plotHiCContactRanges = "./output/plotHiCContactRanges_~{sampleName}.pdf"
        File plotHiCFragment = "./output/plotHiCFragment_~{sampleName}.pdf"
        File plotHiCFragmentSize = "./output/plotHiCFragmentSize_~{sampleName}.pdf"
        File plotMapping = "./output/plotMapping_~{sampleName}.pdf"
        File plotMappingPairing = "./output/plotMappingPairing_~{sampleName}.pdf"
        File pairstat = "./output/~{sampleName}.pairstat"
        File timer = "./output/logs/time.log"
    }

}

task Ice_Normalization {
    input {
        String dockerImage
        String dockerCPU
        String dockerMEM
        String hicPath
        String sampleName
        File rawMatrix
    }

    command {
        set -e -o pipefail
        mkdir -p ./output
        mkdir -p ./output/logs
        date > ./output/logs/time.log
        python ~{hicPath}/scripts/ice \
            --results_filename ./output/~{sampleName}_20000_iced.matrix\
            --filter_low_counts_perc 0.02 \
            --filter_high_counts_perc 0 \
            --max_iter 100 \
            --eps 0.1 \
            --remove-all-zeros-loci \
            --output-bias 1 \
            --verbose 1 \
            ~{rawMatrix}
        date >> ./output/logs/time.log
    }

    runtime {
        docker : dockerImage
        cpu : dockerCPU
        memory : dockerMEM + "GB"
    }

    output {
        File icedMatrix = "./output/~{sampleName}_20000_iced.matrix"
        File timer = "./output/logs/time.log"
    }
}
