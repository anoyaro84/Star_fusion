

rule fastqc_fastq:
    input:
        PATH_FASTQ+"/trimmed/{sample}.fastq.gz"
    output:
        html=PATH_QC+"/{sample}_fastqc.html",
        zip=PATH_QC+"/{sample}_fastqc.zip"
    params:
        config["fastqc"]["extra"]
    wrapper:
        "0.17.0/bio/fastqc"


rule cutadapt:
    input:
        PATH_FASTQ+"/{sample}.fastq.gz"
    output:
        fastq=PATH_FASTQ+"/trimmed/{sample}.fastq.gz",
        qc=PATH_FASTQ+"/trimmed/{sample}.qc.txt"
    params:
        config["cutadapt"]["extra"]
    log:
        PATH_LOG + "/{sample}_cutadapt.log"
    wrapper:
        "0.17.2/bio/cutadapt/se"

rule star:
    input:
        PATH_FASTQ+"/trimmed/{sample}.fastq.gz"
    params:
        index = INDEX,
        thread = config['Ncores'],
        prefix = PATH_BAM+"/{sample}/"
    threads:
        config['Ncores']
    log:
        main=PATH_LOG + "/{sample}.star_log",
        final=PATH_LOG + "/{sample}.star_log_final"
    output:
        bam=PATH_BAM+"/{sample}.bam",
        junction=PATH_BAM+"/junctions/{sample}.out.junction"
    shell:
        """
        mkdir -p {params.prefix}
        STAR --genomeDir {params.index} --readFilesIn {input}\
             --runThreadN {params.thread} --limitBAMsortRAM 31532137230\
             --outSAMtype BAM SortedByCoordinate --twopassMode Basic\
             --outReadsUnmapped None --chimSegmentMin 5 \
             --chimJunctionOverhangMin 5 --alignSJDBoverhangMin 5 \
             --alignMatesGapMax 200000 --alignIntronMax 200000 \
             --chimSegmentReadGapMax parameter 3 --alignSJstitchMismatchNmax 5 -1 5 5\
             --readFilesCommand zcat\
             --outFileNamePrefix {params.prefix}
        mv {params.prefix}/Log.out {log.main}
        mv {params.prefix}/Log.final.out {log.final}
        mv {params.prefix}/Chimeric.out.junction {output.junction}
        mv {params.prefix}/Aligned.sortedByCoord.out.bam {output.bam}
        rm -rf {params.prefix}
        """

rule star_fusion:
    input:
        PATH_BAM+"/junctions/{sample}.out.junction"
    params:
        output_path=PATH_FUSION+"/{sample}_output/",
        reflib = REFLIB
    output:
        PATH_FUSION+"/{sample}.fusion_candidate"
    log:
        PATH_LOG + "/{sample}_star_fusion.log"
    shell:
        """
        STAR-Fusion --genome_lib_dir {params.reflib} -J {input}\
                    --output_dir {params.output_path} &> {log}
        cp {params.output_path}/star-fusion.fusion_candidates.final.abridged {output}
        """
