HiC-WDL
===
## 主文件：hic-docker.wdl  

### bowtie_global_mapping : 用bowtie2跑mapping  
    main_input : bowtie2index.tar檔( tar cvf打包), sample_r1.fastq, sample_r2.fastq  
    main_output: unmap:.unmap.fastq, mapped:.bam  

### bowtie_local_trimming：將global沒map到的用酶切點位切一次  
    main_input: cutsite(酶切點位), global_mapping.unmap.fastq  
    main_output: .trimmed.fastq  

### bowtie_local_mapping:將切完後的基因再map一次  
    main_input: bowtie2index.tar, .trimmed.fastq  
    main_output: local_mapped.bam  

### mapping_combine:將global與local map到的檔案合一  
    main_input: global_mapped.bam, local_mapped.bam (r1,r2都會給)  
    main_output: r1_merged.bam, r2_merged.bam  

### merged_pairs: 用script中的mergeSAM.py將r1,r2 pairs對齊merge起來  
    main_input: r1_merged.bam, r2_merged.bam  
    main_output: merged_pairs.bam  

### mapped_hic_fragment:用mapped_2hic_fragments.py，根據bed file去對mapped出來的結果  
    main_input: .bed, merged_pairs.bam  
    main_output: .validPairs, SinglePairs等其他七項判定檔

