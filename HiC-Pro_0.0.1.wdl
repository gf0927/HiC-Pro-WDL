version 1.0

workflow hic_docker {

    input {
        #Setting
        String docker_image        
        String bowtie2_global_options
        String bowtie2_local_options        
        String n_cpu
        String n_mem        
        String cutsite        
        #PATH
        String hic_path
        #inputs
        String sample_name
        String reference_genome
        File bowtie_idx_path
        File sample_r1
        File sample_r2
        File bed_file
    }    
    call bowtie_global_mapping {
        input:
            sample_name  = sample_name,
            ref_genome = reference_genome,
            bowtie2_global_opts = bowtie2_global_options,
            bowtie_idx_path = bowtie_idx_path,
            sample_r1 = sample_r1,
            sample_r2 = sample_r2,
            bwt_cpu = n_cpu,
            docker_image = docker_image,
            task_mem = n_mem
    }    
    call bowtie_local_trimming {
        input:
            gol_unmap_fastq1=bowtie_global_mapping.gol_unmap,
            gol_unmap_fastq2=bowtie_global_mapping.gol_unmap2,
            hic_path=hic_path,
            cutsite=cutsite,
            sample_name = sample_name,
            ref_genome=reference_genome,
            docker_image = docker_image,
            task_cpu = n_cpu,
            task_mem = n_mem
    }    
    call bowtie_local_mapping{
        input:
            unmap_trimmed_r1 = bowtie_local_trimming.unmap_trimmed_r1,
            unmap_trimmed_r2 = bowtie_local_trimming.unmap_trimmed_r2,
            bowtie_idx_path = bowtie_idx_path,
            sample_name  = sample_name,
            task_cpu = n_cpu,
            ref_genome = reference_genome,
            bowtie2_local_options = bowtie2_local_options,
            docker_image = docker_image,
            task_mem = n_mem
    }    
    call mapping_combine {
        input:
            global_mapped_r1 = bowtie_global_mapping.gol_mapped,
            global_mapped_r2 = bowtie_global_mapping.gol_mapped2,
            local_mapped_r1 = bowtie_local_mapping.local_mapped_r1,
            local_mapped_r2 = bowtie_local_mapping.local_mapped_r2,
            sample_name = sample_name,
            ref_genome = reference_genome,
            docker_image = docker_image,
            task_cpu = n_cpu,
            task_mem = n_mem
    }    
    call merge_pairs {
        input:
            merged_r1 = mapping_combine.bwt2merged_r1,
            merged_r2 = mapping_combine.bwt2merged_r2,
            ref_genome = reference_genome,
            sample_name = sample_name,
            hic_path = hic_path,
            docker_image = docker_image,
            task_cpu = n_cpu,
            task_mem = n_mem
    }    
    call mapped_hic_fragments {
        input:
            bwt2pairs = merge_pairs.bwt2pairs,
            hic_path=hic_path,
            bed_file = bed_file,
            docker_image = docker_image,
            task_cpu = n_cpu,
            task_mem = n_mem,
            sample_name = sample_name,
            ref_genome = reference_genome
    }
}

task bowtie_global_mapping {
    input {
        String sample_name
        String ref_genome
        String bowtie2_global_opts
        String bwt_cpu
        File bowtie_idx_path
        File sample_r1
        File sample_r2
        String docker_image
        String task_mem
    }
    command {
        mkdir ./output
        mkdir ./output/logs
        date > ./output/logs/time.log
        echo "##HiC-Pro mapping" > ./output/logs/~{sample_name}_r1_bowtie2.log
        echo "##HiC-Pro mapping" > ./output/logs/~{sample_name}_r2_bowtie2.log
        bowtie2 ~{bowtie2_global_opts} \
            --un ./output/~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap.fastq \
            --rg-id BMG \
            --rg SM:~{sample_name}_r1 \
            -p ~{bwt_cpu} \
            -x ~{bowtie_idx_path}/~{ref_genome} \
            -U ~{sample_r1} \
            2>> ./output/logs/~{sample_name}_r1_bowtie2.log | \
        samtools view -F 4 -bS - > ./output/~{sample_name}_r1_~{ref_genome}.bwt2glob.bam        
        bowtie2 ~{bowtie2_global_opts} \
            --un ./output/~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap.fastq \
            --rg-id BMG \
            --rg SM:~{sample_name}_r2 \
            -p ~{bwt_cpu} \
            -x ~{bowtie_idx_path}/~{ref_genome} \
            -U ~{sample_r2} \
            2>> ./output/logs/~{sample_name}_r2_bowtie2.log | \
        samtools view -F 4 -bS - > ./output/~{sample_name}_r2_~{ref_genome}.bwt2glob.bam        
        rm -r ~{ref_genome}
        date >> ./output/logs/time.log
    }
    runtime {
        docker: docker_image
        cpu : bwt_cpu
        memory : task_mem + "GB"
    }
    output {
        File gol_unmap = "./output/~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap.fastq"
        File gol_mapped = "./output/~{sample_name}_r1_~{ref_genome}.bwt2glob.bam"
        File gol_log = "./output/logs/~{sample_name}_r1_bowtie2.log"
        File gol_unmap2 = "./output/~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap.fastq"
        File gol_mapped2 = "./output/~{sample_name}_r2_~{ref_genome}.bwt2glob.bam"
        File gol_log2 = "./output/logs/~{sample_name}_r2_bowtie2.log"
    }
}
task bowtie_local_trimming {
    input {
        File gol_unmap_fastq1
        File gol_unmap_fastq2
        String hic_path
        String cutsite
        String sample_name
        String ref_genome
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
        ~{hic_path}/scripts/cutsite_trimming \
            --fastq ~{gol_unmap_fastq1} \
            --cutsite ~{cutsite} \
            --out ./output/~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap_trimmed.fastq > ./output/logs/~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap_readsTrimming.log 2>&1
        ~{hic_path}/scripts/cutsite_trimming \
            --fastq ~{gol_unmap_fastq2} --cutsite ~{cutsite} \
            --out ./output/~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap_trimmed.fastq > ./output/logs/~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap_readsTrimming.log 2>&1
        date >> ./output/logs/time.log
    }
    output {
        File unmap_trimmed_r1="./output/~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap_trimmed.fastq"
        File unmap_trimmed_r2="./output/~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap_trimmed.fastq"
        File trim_log1="./output/logs/~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap_readsTrimming.log"
        File trim_log2="./output/logs/~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap_readsTrimming.log"
    }
}
task bowtie_local_mapping {
    input {
        File unmap_trimmed_r1
        File unmap_trimmed_r2
        File bowtie_idx_path        
        String bowtie2_local_options        
        String ref_genome
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
    command {
        mkdir ./output
        mkdir ./output/logs
        date > ./output/logs/time.log
        bowtie2 ~{bowtie2_local_options} \
            --rg-id BML \
            --rg SM:~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap \
            -p ~{task_cpu} \
            -x ~{bowtie_idx_path}/~{ref_genome} \
            -U ~{unmap_trimmed_r1} \
            2>> ./output/logs/~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap_bowtie2.log | \
        samtools view -bS - > ./output/~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap_bwt2loc.bam
        bowtie2 ~{bowtie2_local_options} \
            --rg-id BML \
            --rg SM:~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap \
            -p ~{task_cpu} \
            -x ~{bowtie_idx_path}/~{ref_genome} \
            -U ~{unmap_trimmed_r2} \
            2>> ./output/logs/~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap_bowtie2.log | \
        samtools view -bS - > ./output/~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap_bwt2loc.bam
        rm -r ~{ref_genome}
        date >> ./output/logs/time.log
    }     
    output {
        File local_mapped_r1 = "./output/~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap_bwt2loc.bam"
        File local_mapped_r2 = "./output/~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap_bwt2loc.bam"
        File local_mapped_r1_log = "./output/logs/~{sample_name}_r1_~{ref_genome}.bwt2glob.unmap_bowtie2.log"
        File local_mapped_r2_log = "./output/logs/~{sample_name}_r2_~{ref_genome}.bwt2glob.unmap_bowtie2.log"
     }
}
task mapping_combine {
    input {
        File global_mapped_r1
        File global_mapped_r2
        File local_mapped_r1
        File local_mapped_r2
        String sample_name
        String ref_genome
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
        samtools merge \
            -@ 20 \
            -n \
            -f ./output/~{sample_name}_r1_~{ref_genome}.bwt2merged.bam ~{global_mapped_r1} ~{local_mapped_r1}
        samtools merge \
            -@ 20 \
            -n \
            -f ./output/~{sample_name}_r2_~{ref_genome}.bwt2merged.bam ~{global_mapped_r2} ~{local_mapped_r2}
        samtools sort \
            -@ 20 \
            -n \
            -T tmp/~{sample_name}_r1_~{ref_genome} \
            -o ./output/~{sample_name}_r1_~{ref_genome}.bwt2merged.sorted.bam ./output/~{sample_name}_r1_~{ref_genome}.bwt2merged.bam
        samtools sort \
            -@ 20 \
            -n \
            -T tmp/~{sample_name}_r2_~{ref_genome} \
            -o ./output/~{sample_name}_r2_~{ref_genome}.bwt2merged.sorted.bam ./output/~{sample_name}_r2_~{ref_genome}.bwt2merged.bam
        mv ./output/~{sample_name}_r1_~{ref_genome}.bwt2merged.sorted.bam ./output/~{sample_name}_r1_~{ref_genome}.bwt2merged.bam
        mv ./output/~{sample_name}_r2_~{ref_genome}.bwt2merged.sorted.bam ./output/~{sample_name}_r2_~{ref_genome}.bwt2merged.bam
        date >> ./output/logs/time.log
    }    
    output {
        File bwt2merged_r1 = "./output/~{sample_name}_r1_~{ref_genome}.bwt2merged.bam"
        File bwt2merged_r2 = "./output/~{sample_name}_r2_~{ref_genome}.bwt2merged.bam"
    }
}
task merge_pairs {
    input {
        File merged_r1
        File merged_r2
        String hic_path
        String sample_name
        String ref_genome
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
        python3 ~{hic_path}/scripts/mergeSAM.py \
            -q 0 \
            -t \
            -v \
            -f ~{merged_r1} \
            -r ~{merged_r2} \
            -o ./output/~{sample_name}_~{ref_genome}.bwt2pairs.bam
        date >> ./output/logs/time.log
    }
    output {
        File bwt2pairs = "./output/~{sample_name}_~{ref_genome}.bwt2pairs.bam"
    }
}
task mapped_hic_fragments {
    input {
        File bwt2pairs
        String hic_path
        File bed_file
        String docker_image
        Int task_cpu
        Int task_mem
        String sample_name
        String ref_genome
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
        python3 ~{hic_path}/scripts/mapped_2hic_fragments.py \
            -v \
            -S \
            -t 100 \
            -m 100000 \
            -s 100 \
            -l 600 \
            -a \
            -f ~{bed_file} \
            -r ~{bwt2pairs} \
            -o ./output
        date >> ./output/logs/time.log
    }
    output {
        File bam_validpairs = "./output/~{sample_name}_~{ref_genome}.bwt2pairs_interaction.bam"
        File validPairs =  "./output/~{sample_name}_~{ref_genome}.bwt2pairs.validPairs"
        File DEPairs = "./output/~{sample_name}_~{ref_genome}.bwt2pairs.DEPairs"
        File DumpedPairs = "./output/~{sample_name}_~{ref_genome}.bwt2pairs.DumpPairs"
        File FiltPairs = "./output/~{sample_name}_~{ref_genome}.bwt2pairs.FiltPairs"
        File REPairs = "./output/~{sample_name}_~{ref_genome}.bwt2pairs.REPairs"
        File SinglePairs ="./output/~{sample_name}_~{ref_genome}.bwt2pairs.SinglePairs"
        File SCPairs = "./output/~{sample_name}_~{ref_genome}.bwt2pairs.SCPairs"
        File RSstat = "./output/~{sample_name}_~{ref_genome}.bwt2pairs.RSstat"

    }    
}
