from os.path import join
import os
import pandas as pd
from scripts.Load import samplesheet
import json

snakedir = os.getcwd()
configfile: 'config.yaml'

print('Snakemake running dir:', snakedir)
print('Pipeline working dir:', config['workdir'])
print('#'*100,'\n')

print("Config setted to:")
print(json.dumps(config, indent=4, sort_keys=True), '\n')
print('#'*100,'\n')

sampledic, libdic, rundic = samplesheet(config['samplesheet'])
print("Feeding from samplesheet:")
print('sampledic:\n', json.dumps(sampledic, indent=4, sort_keys=True), '\n')
print('libdic:\n', json.dumps(libdic, indent=4, sort_keys=True), '\n')
print('rundic:\n', json.dumps(rundic, indent=4, sort_keys=True), '\n')
print('#'*100,'\n')

workdir: config['workdir']
itvfastp = ['%.4d'%(itv) for itv in range(1, int(config['parameter']['itvfastp'])+1)]
itv50 = ['%.4d'%(itv) for itv in range(1, 50+1)]
itv17 = ['%.4d'%(itv) for itv in range(1, 17+1)]
itvdv = ['%.5d'%(itv) for itv in range(int(config['parameter']['deepvariant_shards'])+1)]

## Global wildcards
#SAMPLES, = glob_wildcards(sampledic.keys()).sample


rule all:
    input:
        expand("02.Alignment/LevelPU/{run}/{run}.sort.md.metrics", run=rundic.keys()),
        expand("02.Alignment/LevelLB/{lib}/{lib}.sort.md.metrics", lib=libdic.keys()),
        expand("02.Alignment/chrM/{sample}/MD5.txt", sample=sampledic.keys()),
        expand("02.Alignment/BQSR/{sample}/{sample}.BQSR.cram", sample=sampledic.keys()),
        expand("02.Alignment/BQSR/{sample}/Metrics/Summary.ok", sample=sampledic.keys()),
        
        # BamQC
        lambda wildcards: [f"02.Alignment/BQSR/{sample}/Metrics/Stats.{sample}.SM.txt" for sample in sampledic],
        lambda wildcards: [f"02.Alignment/BQSR/{sample}/Metrics/Stats.{sample}.LB_{lib}.txt" for sample in sampledic for lib in sampledic[sample]],
        lambda wildcards: [f"02.Alignment/BQSR/{sample}/Metrics/Stats.{sample}.RG_{run}.txt" for sample in sampledic for lib in sampledic[sample] for run in libdic[lib]],
        lambda wildcards: [f"02.Alignment/BQSR/{sample}/Mosdepth/{sample}.mosdepth.summary.txt" for sample in sampledic],
        
        # Germline
        expand("03.Germline.deepvariant/{sample}/{sample}.deepvariant.vcf.gz", sample=sampledic.keys()),
        expand("03.Germline.freebayes/{sample}/{sample}.freebayes.vcf.gz", sample=sampledic.keys()),
        expand("03.Germline.GATK/{sample}/{sample}.vcf.gz", sample=sampledic.keys()),
        expand("03.Germline.Strelka/{sample}/results/variants/variants.vcf.gz", sample=sampledic.keys()),
        expand("03.Germline.Manta/{sample}/results/variants/candidateSV.vcf.gz", sample=sampledic.keys()),
        expand("03.Germline.Tiddit/{sample}/{sample}.tiddit.vcf.gz", sample=sampledic.keys()),
        expand("04.Germline.Gridss/{sample}/{sample}.gridss.vcf", sample=sampledic.keys()),

rule QC:
    input:
        reads1=lambda wildcards: rundic[wildcards.run]['Read1'],
        reads2=lambda wildcards: rundic[wildcards.run]['Read2'],
    output:
        reads1out=temp(expand("01.CleanData/{{run}}/{itv}.{{run}}.R1.cln.fq.gz", itv=itvfastp)),
        reads2out=temp(expand("01.CleanData/{{run}}/{itv}.{{run}}.R2.cln.fq.gz", itv=itvfastp)),
        htmlout="01.CleanData/{run}/{run}.QC.html",
        jsonout="01.CleanData/{run}/{run}.QC.json",
    params:
        reads1out="01.CleanData/{run}/{run}.R1.cln.fq.gz",
        reads2out="01.CleanData/{run}/{run}.R2.cln.fq.gz",
    log: 
        out = snakedir + "/logs/A1.QC/{run}.o",
        err = snakedir + "/logs/A1.QC/{run}.e",
    threads: 12
    resources: 
        mem = '12g',
        extra = ' --gres=lscratch:10 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][fastp]} "
        "fastp"
        " -i {input.reads1}"
        " -I {input.reads2}"
        " -o {params.reads1out}"
        " -O {params.reads2out}"
        " -h {output.htmlout}"
        " -j {output.jsonout}" 
        " -s {config[parameter][itvfastp]}"
        " -w {threads}"
        " >> {log.out} 2>> {log.err}"

        
rule Bwa_mem:
    input:
        read1="01.CleanData/{run}/{itv}.{run}.R1.cln.fq.gz",
        read2="01.CleanData/{run}/{itv}.{run}.R2.cln.fq.gz",
    output:
        bam=temp("02.Alignment/LevelPU/{run}/{itv}.{run}.bam"),
    log:
        out = snakedir+"/logs/B1.bwa/{run}.{itv}.o",
        err = snakedir+"/logs/B1.bwa/{run}.{itv}.e"

    threads: 36
    resources: 
        mem = '36g',
        extra = ' --gres=lscratch:200 ',
    params:
        bwa=lambda wildcards:' -K 100000000 -v 3 -Y -R "@RG\\tID:{}\\tLB:{}\\tPL:illumina\\tPU:{}\\tSM:{}" '.format(
                                wildcards.run, rundic[wildcards.run]['LB'], rundic[wildcards.run]['PU'],  rundic[wildcards.run]['SM']),
        samblaster="",
        sambambasort=" --tmpdir /lscratch/$SLURM_JOB_ID "
    run:
        bwacpu = threads-4,
        shell(
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][bwa]} "
        "bwa mem -t {bwacpu}"
        " {params.bwa}"
        " {config[references][bwaidx]}"
        " {input.read1} {input.read2} 2>> {log.err} |"
        " {config[singularity]} {config[simg][samtools]}"
        " samtools view"
        " -@ 8"
        " -O BAM"
        " -o {output.bam}"
        " >> {log.out} 2>> {log.err}"
        )

        
rule PuMarkdup:
    input:
        bam = lambda wildcards: expand("02.Alignment/LevelPU/{run}/{itv}.{run}.bam", itv=itvfastp, run=wildcards.run),# run=libdic[wildcards.lib]),
    output:
        bam = temp("02.Alignment/LevelPU/{run}/{run}.sort.md.bam"),
        #bai = temp("02.Alignment/LevelPU/{run}/{run}.sort.md.bam.bai"),
        #sbi = temp("02.Alignment/LevelPU/{run}/{run}.sort.md.bam.sbi"),
    log:
        out = snakedir+"/logs/B2.pu_mkdup/{run}.o",
        err = snakedir+"/logs/B2.pu_mkdup/{run}.e",
    threads:  56
    resources:
        mem = '80g',
        extra = ' --gres=lscratch:800 ',
    run:
        inputs = " ".join("-I {}".format(in_) for in_ in input.bam),
        shell(
            "module load singularity > {log.out} 2> {log.err}\n"
            "{config[singularity]} {config[simg][gatk]} "
            "gatk"
            " --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms60G -Xmx60G\""
            " MarkDuplicatesSpark"
            " {inputs}"
            " -O {output.bam}"
            " --optical-duplicate-pixel-distance 2500"
            " -R {config[references][fasta]}"
            " --spark-master local[{threads}]"
            " --tmp-dir /lscratch/$SLURM_JOB_ID"
            " --treat-unsorted-as-querygroup-ordered true"
            " >> {log.out} 2>> {log.err}"
        )            

rule PuMarkdupMetrics:
    input:
        bam= "02.Alignment/LevelPU/{run}/{run}.sort.md.bam",
        #bai = "02.Alignment/LevelPU/{run}/{run}.sort.md.bam.bai",
        #sbi = "02.Alignment/LevelPU/{run}/{run}.sort.md.bam.sbi",
    output:
        metrics = "02.Alignment/LevelPU/{run}/{run}.sort.md.metrics",
    log:
        out = snakedir+"/logs/B3.pu_mkdup_metrics/{run}.o",
        err = snakedir+"/logs/B3.pu_mkdup_metrics/{run}.e",
    threads:  4
    resources:
        mem = '8g',
        extra = ' --gres=lscratch:20 ',
    run:
        shell(
            "module load singularity > {log.out} 2> {log.err}\n"
            "{config[singularity]} {config[simg][gatk]} "
            "gatk"
            " CollectDuplicateMetrics"
            " -I {input.bam}"
            " -M {output.metrics}"
            " --CREATE_MD5_FILE true"
            " >> {log.out} 2>> {log.err}"
        )

rule LbMarkdup:
    input:
        bam = lambda wildcards: expand("02.Alignment/LevelPU/{run}/{run}.sort.md.bam", run=libdic[wildcards.lib]),
    output:
        bam = temp("02.Alignment/LevelLB/{lib}/{lib}.sort.md.bam"),
        #bai = temp("02.Alignment/LevelLB/{lib}/{lib}.sort.md.bam.bai"),
        #sbi = temp("02.Alignment/LevelLB/{lib}/{lib}.sort.md.bam.sbi"),
    log:
        out = snakedir+"/logs/B4.lb_mkdup/{lib}.o",
        err = snakedir+"/logs/B4.lb_mkdup/{lib}.e",
    threads:  56
    resources:
        mem = '80g',
        extra = ' --gres=lscratch:800 ',
    run:
        inputs = " ".join("-I {}".format(in_) for in_ in input.bam),
        shell(
            "module load singularity > {log.out} 2> {log.err}\n"
            "{config[singularity]} {config[simg][gatk]} "
            "gatk"
            " --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms60G -Xmx60G\""
            " MarkDuplicatesSpark"
            " {inputs}"
            " -O {output.bam}"
            " --optical-duplicate-pixel-distance 2500"
            " -R {config[references][fasta]}"
            " --spark-master local[{threads}]"
            " --tmp-dir /lscratch/$SLURM_JOB_ID"
            " --allow-multiple-sort-orders-in-input true"
            " >> {log.out} 2>> {log.err}"
        )

rule LbMarkdupMetrics:
    input:
        bam= "02.Alignment/LevelLB/{lib}/{lib}.sort.md.bam",
    output:
        metrics = "02.Alignment/LevelLB/{lib}/{lib}.sort.md.metrics",
    log:
        out = snakedir+"/logs/B5.lb_mkdup_metrics/{lib}.o",
        err = snakedir+"/logs/B5.lb_mkdup_metrics/{lib}.e",
    threads:  4
    resources:
        mem = '8g',
        extra = ' --gres=lscratch:20 ',
    run:
        shell(
            "module load singularity > {log.out} 2> {log.err}\n"
            "{config[singularity]} {config[simg][gatk]} "
            "gatk"
            " CollectDuplicateMetrics"
            " -I {input.bam}"
            " -M {output.metrics}"
            " --CREATE_MD5_FILE true"
            " >> {log.out} 2>> {log.err}"
        )

rule SmMerge:
    input:
        bam = lambda wildcards: ["02.Alignment/LevelLB/{}/{}.sort.md.bam".format(lib, lib) for lib in sampledic[wildcards.sample]],
    output:
        bam = temp("02.Alignment/LevelSM/{sample}/{sample}.sort.md.bam"),
        #bai = temp("02.Alignment/LevelSM/{sample}/{sample}.sort.md.bam.bai"),
    log:
        out = snakedir+"/logs/B6.sm_merge/{sample}.o",
        err = snakedir+"/logs/B6.sm_merge/{sample}.e",
    threads:  16
    resources:
        mem = '32g',
        extra = ' --gres=lscratch:10 ',
    run:
        if len(input.bam) > 1:
            inputs = " ".join(input.bam),
            shell(
                "module load singularity > {log.out} 2> {log.err}\n"
                "{config[singularity]} {config[simg][sambamba]} "
                "sambamba merge"
                " -t {threads}"
                " {output.bam}"
                " {inputs} "
                " >> {log.out} 2>> {log.err}"
            )
        else:
            shell(
                "module load singularity > {log.out} 2> {log.err}\n"
                "ln {input.bam[0]} {output.bam}"
                " >> {log.out} 2>> {log.err}\n"
                "{config[singularity]} {config[simg][samtools]} "
                " samtools index"
                " -@ {threads}"
                " {output.bam}"
                " >> {log.out} 2>> {log.err}"
            )

            
rule Mitobam:
    input:
        bam="02.Alignment/LevelSM/{sample}/{sample}.sort.md.bam",
        #bai="02.Alignment/LevelSM/{sample}/{sample}.sort.md.bam.bai",
    output:
        bam="02.Alignment/chrM/{sample}/{sample}.bam",
        md5="02.Alignment/chrM/{sample}/MD5.txt",
    log:
        out = snakedir+"/logs/B7.Mitobam/{sample}.o",
        err = snakedir+"/logs/B7.Mitobam/{sample}.e",
    threads:  8
    resources:
        mem = '16g',
        extra = ' --gres=lscratch:10 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        " {config[singularity]} {config[simg][samtools]}"
        " samtools view"
        " -@ 8"
        " -O BAM"
        " -o {output.bam}"
        " {input.bam}"
        " chrM"
        " >> {log.out} 2>> {log.err}\n"
        " {config[singularity]} {config[simg][samtools]}"
        " samtools index"
        " -@ 8"
        " {output.bam}"
        " >> {log.out} 2>> {log.err}\n"
        " cd `dirname {output.bam}` && md5sum `basename {output.bam}`* > `basename {output.md5}` 2>>{log.err}"
        
        
##############################################################################################

rule BQSR:
    input:
        bam="02.Alignment/LevelSM/{sample}/{sample}.sort.md.bam",
        #bai="02.Alignment/LevelSM/{sample}/{sample}.sort.md.bam.bai",
    params:
        itvbed="/data/yuk5/pipeline/wgs_germline/ref/hg38_chr_intervals/itv_{itv}.bed"
    output:
        metrics=temp("02.Alignment/BQSR/{sample}/{sample}.{itv}.BQSR.metrics"),
        bam=temp("02.Alignment/BQSR/{sample}/{sample}.{itv}.BQSR.bam"),
        bai=temp("02.Alignment/BQSR/{sample}/{sample}.{itv}.BQSR.bai"),
    log:
        out = snakedir+"/logs/B8.BQSR/{sample}.{itv}.o",
        err = snakedir+"/logs/B8.BQSR/{sample}.{itv}.e",
    threads:  2
    resources:
        mem  = '8g',
        extra = ' --gres=lscratch:400 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        " {config[singularity]} {config[simg][gatk]} "
        "gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms4G -Xmx4G -XX:ParallelGCThreads=2\" BaseRecalibrator "
        " -R {config[references][fasta]}"
        " -I {input.bam}"
        " -O {output.metrics}"
        " --use-original-qualities"
        " --known-sites {config[references][gatk_dbsnp]}"
        " --known-sites {config[references][gatk_1000g]}"
        " --known-sites {config[references][gatk_indel]}"
        " --intervals   {params.itvbed}"
        " >> {log.out} 2>> {log.err}\n"
        " {config[singularity]} {config[simg][gatk]} "
        "gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms4G -Xmx4G -XX:ParallelGCThreads=2\" ApplyBQSR"
        " --add-output-sam-program-record"
        " -R {config[references][fasta]}"
        " -I {input.bam}"
        " -O {output.bam}"
        " --use-original-qualities"
        " -bqsr {output.metrics}"
        " --static-quantized-quals 10"
        " --static-quantized-quals 20"
        " --static-quantized-quals 30"
        " -L {params.itvbed}"
        " >> {log.out} 2>> {log.err}\n"

        
rule BQSR_mergeM:
    input:
        expand("02.Alignment/BQSR/{{sample}}/{{sample}}.{itv}.BQSR.metrics", itv=itv17),
    output:
        metrics="02.Alignment/BQSR/{sample}/Metrics/BQSR.metrics",
    log:
        out = snakedir+"/logs/B9a.BQSR_mergereport/{sample}.o",
        err = snakedir+"/logs/B9a.BQSR_mergereport/{sample}.e",
    threads:  2
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:40 ',
    run:
        inputs = " ".join("-I {}".format(in_) for in_ in input),
        shell(
            "module load singularity > {log.out} 2> {log.err}\n"
            " {config[singularity]} {config[simg][gatk]} "
            " gatk --java-options \"-Xms3000m\" "
            " GatherBQSRReports"
            " {inputs}"
            " -O {output.metrics}"
            " >> {log.out} 2>> {log.err}\n"
            )
            


rule BQSR_mergeB:
    input:
        expand("02.Alignment/BQSR/{{sample}}/{{sample}}.{itv}.BQSR.bam", itv=itv17),
    output:
        bam="02.Alignment/BQSR/{sample}/{sample}.BQSR.bam",
    log:
        out = snakedir+"/logs/B9b.BQSR_mergebam/{sample}.o",
        err = snakedir+"/logs/B9b.BQSR_mergebam/{sample}.e",
    threads:  16
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:40 ',
    run:
        inputs = " ".join("{}".format(in_) for in_ in input),
        shell(
            "module load singularity > {log.out} 2> {log.err}\n"
            " {config[singularity]} {config[simg][sambamba]} "
            " sambamba merge "
            " -t {threads}"
            " {output.bam}"
            " {inputs}"
            " >> {log.out} 2>> {log.err}\n"
            )


rule BQSR_to_CRAM:
    input:
        bam="02.Alignment/BQSR/{sample}/{sample}.BQSR.bam",
    output:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
        crai="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram.crai",
        md5="02.Alignment/BQSR/{sample}/MD5.txt",
    log:
        out = snakedir+"/logs/B10.BQSR_to_CRAM/{sample}.o",
        err = snakedir+"/logs/B10.BQSR_to_CRAM/{sample}.e",
    threads:  12
    resources:
        mem  = '24g',
        extra = ' --gres=lscratch:100 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][samtools]} "
        "samtools"
        "    view"
        "    -@ {threads}"
        "    -T {config[references][fasta]}"
        "    {input.bam}"
        "    -O CRAM"
        "    -o {output.cram}"
        " >> {log.out} 2>> {log.err}\n"
        "{config[singularity]} {config[simg][samtools]} "
        "samtools"
        "    index"
        "    -@ {threads}"
        "    {output.cram}"
        " >> {log.out} 2>> {log.err}\n"
        "cd `dirname {output.cram}` \n"
        "md5sum $(basename {output.cram})* > $(basename {output.md5})"
        " >> {log.out} 2>> {log.err}\n"
        

rule BQSR_CollectQualityYieldMetrics:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        status=temp("02.Alignment/BQSR/{sample}/Metrics/QualityYield.metrics.ok"),
    params:
        prefix="02.Alignment/BQSR/{sample}/Metrics/QualityYield.metrics",
    log:
        out = snakedir+"/logs/B9d.BQSR_CollectQualityYieldMetrics/{sample}.o",
        err = snakedir+"/logs/B9d.BQSR_CollectQualityYieldMetrics/{sample}.e",
    threads:  2
    resources:
        mem  = '4g',
        extra = ' --gres=lscratch:40 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        " {config[singularity]} {config[simg][gatk]} "
        " gatk --java-options \"-Xms2000m -Xmx3000m\" "
        " CollectQualityYieldMetrics"
        " -I {input.cram}"
        " -O {params.prefix}"
        " --CREATE_MD5_FILE true"
        " -R {config[references][fasta]}"
        " >> {log.out} 2>> {log.err}\n"
        "touch {output.status}"
        " >> {log.out} 2>> {log.err}\n"
        
rule BQSR_CollectWgsMetrics:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        status=temp("02.Alignment/BQSR/{sample}/Metrics/Wgs.metrics.ok"),
    params:
        prefix="02.Alignment/BQSR/{sample}/Metrics/Wgs.metrics",
    log:
        out = snakedir+"/logs/B9e.BQSR_CollectWgsMetrics/{sample}.o",
        err = snakedir+"/logs/B9e.BQSR_CollectWgsMetrics/{sample}.e",
    threads:  2
    resources:
        mem  = '4g',
        extra = ' --gres=lscratch:40 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        " {config[singularity]} {config[simg][gatk]} "
        " gatk --java-options \"-Xms2000m -Xmx3000m\" "
        " CollectWgsMetrics"
        " -I {input.cram}"
        " -O {params.prefix}"
        " -R {config[references][fasta]}"
        " --INCLUDE_BQ_HISTOGRAM true"
        " --INTERVALS {config[references][wgscovitv]}"
        " --VALIDATION_STRINGENCY SILENT"
        " --USE_FAST_ALGORITHM true"
        " --CREATE_MD5_FILE true"
        " >> {log.out} 2>> {log.err}\n"
        "touch {output.status}"
        " >> {log.out} 2>> {log.err}\n"
        
rule BQSR_CollectAllReadsMultipleMetrics:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        status=temp("02.Alignment/BQSR/{sample}/Metrics/AllReadsMultiple.ok"),
    params:
        prefix="02.Alignment/BQSR/{sample}/Metrics/AllReadsMultiple",
    log:
        out = snakedir+"/logs/B9f.BQSR_CollectAllReadsMultipleMetrics/{sample}.o",
        err = snakedir+"/logs/B9f.BQSR_CollectAllReadsMultipleMetrics/{sample}.e",
    threads:  2
    resources:
        mem  = '8g',
        extra = ' --gres=lscratch:40 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        " {config[singularity]} {config[simg][gatk]} "
        " gatk --java-options \"-Xms5000m -Xmx6500m\" "
        " CollectMultipleMetrics"
        " -I {input.cram}"
        " -O {params.prefix}"
        " --ASSUME_SORTED true"
        " --CREATE_MD5_FILE true"
        " -R {config[references][fasta]}"
        " --VALIDATION_STRINGENCY LENIENT"
        " --PROGRAM null "
        " --PROGRAM CollectBaseDistributionByCycle "
        " --PROGRAM CollectInsertSizeMetrics "
        " --PROGRAM MeanQualityByCycle "
        " --PROGRAM QualityScoreDistribution "
        " --METRIC_ACCUMULATION_LEVEL null "
        " --METRIC_ACCUMULATION_LEVEL ALL_READS "
        " >> {log.out} 2>> {log.err}\n"
        "touch {output.status}"
        " >> {log.out} 2>> {log.err}\n"
        
rule BQSR_CollectReadGroupsMultipleMetrics:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        status=temp("02.Alignment/BQSR/{sample}/Metrics/ReadGroupsMultiple.ok"),
    params:
        prefix="02.Alignment/BQSR/{sample}/Metrics/ReadGroupsMultiple",
    log:
        out = snakedir+"/logs/B9g.BQSR_CollectReadGroupsMultipleMetrics/{sample}.o",
        err = snakedir+"/logs/B9g.BQSR_CollectReadGroupsMultipleMetrics/{sample}.e",
    threads:  2
    resources:
        mem  = '8g',
        extra = ' --gres=lscratch:40 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        " {config[singularity]} {config[simg][gatk]} "
        " gatk --java-options \"-Xms5000m -Xmx6500m\" "
        " CollectMultipleMetrics"
        " -I {input.cram}"
        " -O {params.prefix}"
        " --ASSUME_SORTED true"
        " --CREATE_MD5_FILE true"
        " -R {config[references][fasta]}"
        " --VALIDATION_STRINGENCY LENIENT"
        " --PROGRAM null "
        " --PROGRAM CollectBaseDistributionByCycle "
        " --PROGRAM CollectInsertSizeMetrics "
        " --PROGRAM CollectAlignmentSummaryMetrics "
        " --PROGRAM QualityScoreDistribution "
        " --METRIC_ACCUMULATION_LEVEL null "
        " --METRIC_ACCUMULATION_LEVEL READ_GROUP "
        " >> {log.out} 2>> {log.err}\n"
        "touch {output.status}"
        " >> {log.out} 2>> {log.err}\n"
        
rule BQSR_CollectAggregationMetrics:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        status=temp("02.Alignment/BQSR/{sample}/Metrics/SM_LB_Aggregation.ok"),
    params:
        prefix="02.Alignment/BQSR/{sample}/Metrics/SM_LB_Aggregation",
    log:
        out = snakedir+"/logs/B9h.BQSR_CollectAggregationMetrics/{sample}.o",
        err = snakedir+"/logs/B9h.BQSR_CollectAggregationMetrics/{sample}.e",
    threads:  2
    resources:
        mem  = '8g',
        extra = ' --gres=lscratch:40 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        " {config[singularity]} {config[simg][gatk]} "
        " gatk --java-options \"-Xms5000m -Xmx6500m\" "
        " CollectMultipleMetrics"
        " -I {input.cram}"
        " -O {params.prefix}"
        " --ASSUME_SORTED true"
        " --CREATE_MD5_FILE true"
        " -R {config[references][fasta]}"
        " --PROGRAM null "
        " --PROGRAM CollectAlignmentSummaryMetrics "
        " --PROGRAM CollectInsertSizeMetrics "
        " --PROGRAM CollectSequencingArtifactMetrics "
        " --PROGRAM QualityScoreDistribution "
        " --METRIC_ACCUMULATION_LEVEL null "
        " --METRIC_ACCUMULATION_LEVEL SAMPLE "
        " --METRIC_ACCUMULATION_LEVEL LIBRARY "
        " >> {log.out} 2>> {log.err}\n"
        "touch {output.status}"
        " >> {log.out} 2>> {log.err}\n"
        
rule BQSR_SummaryMetrics:
    input:
        status1="02.Alignment/BQSR/{sample}/Metrics/QualityYield.metrics.ok",
        status2="02.Alignment/BQSR/{sample}/Metrics/Wgs.metrics.ok",
        status3="02.Alignment/BQSR/{sample}/Metrics/AllReadsMultiple.ok",
        status4="02.Alignment/BQSR/{sample}/Metrics/ReadGroupsMultiple.ok",
        status5="02.Alignment/BQSR/{sample}/Metrics/SM_LB_Aggregation.ok",
        status6="02.Alignment/BQSR/{sample}/Metrics/SM_LB_Aggregation.ok",
    output:
        status="02.Alignment/BQSR/{sample}/Metrics/Summary.ok",
    log:
        out = snakedir+"/logs/B9z.BQSR_SummaryMetrics/{sample}.o",
        err = snakedir+"/logs/B9z.BQSR_SummaryMetrics/{sample}.e",
    threads:  2
    resources:
        mem  = '8g',
        extra = ' --gres=lscratch:40 ',
    shell:
        "touch {output.status}"

        
        
rule CRAM_TO_STAT_RG:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        stats="02.Alignment/BQSR/{sample}/Metrics/Stats.{sample}.RG_{run}.txt",
    log:
        out = snakedir+"/logs/B11.CRAM_TO_STAT_RG/{sample}.{run}.o",
        err = snakedir+"/logs/B11.CRAM_TO_STAT_RG/{sample}.{run}.e",
    threads:  12
    resources:
        mem  = '24g',
        extra = ' --gres=lscratch:100 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][samtools]} "
        "samtools"
        "    view"
        "    -@ 8"
        "    -b"
        "    -T {config[references][fasta]}"
        "    -r {wildcards.run}"
        "    {input.cram} | "
        "{config[singularity]} {config[simg][samtools]} "
        "samtools"
        "    stats"
        "    -@ 4"
        "    --reference {config[references][fasta]}"
        "    -"
        "    > {output.stats} 2>> {log.err}\n"
        
        
rule CRAM_TO_STAT_LB:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        stats="02.Alignment/BQSR/{sample}/Metrics/Stats.{sample}.LB_{lib}.txt",
    log:
        out = snakedir+"/logs/B12.CRAM_TO_STAT_LB/{sample}.{lib}.o",
        err = snakedir+"/logs/B12.CRAM_TO_STAT_LB/{sample}.{lib}.e",
    threads:  12
    resources:
        mem  = '24g',
        extra = ' --gres=lscratch:100 ',
    run:
        ids = ' '.join([' -r %s'%rg for rg in libdic[wildcards.lib]])
        shell(
            "module load singularity > {log.out} 2> {log.err}\n"
            "{config[singularity]} {config[simg][samtools]} "
            "samtools"
            "    view"
            "    -@ 8"
            "    -b"
            "    -T {config[references][fasta]}"
            "    {ids}"
            "    {input.cram} | "
            "{config[singularity]} {config[simg][samtools]} "
            "samtools"
            "    stats"
            "    -@ 4"
            "    --reference {config[references][fasta]}"
            "    -"
            "    > {output.stats} 2>> {log.err}\n"
        )
        
        
rule CRAM_TO_STAT_SM:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        stats="02.Alignment/BQSR/{sample}/Metrics/Stats.{sample}.SM.txt",
    log:
        out = snakedir+"/logs/B13.CRAM_TO_STAT_SM/{sample}.o",
        err = snakedir+"/logs/B13.CRAM_TO_STAT_SM/{sample}.e",
    threads:  12
    resources:
        mem  = '24g',
        extra = ' --gres=lscratch:100 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][samtools]} "
        "samtools"
        "    stats"
        "    -@ {threads}"
        "    --reference {config[references][fasta]}"
        "    {input.cram}"
        "    > {output.stats} 2>> {log.err}\n"
        
rule CRAM_TO_MOSDEPTH:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        stats="02.Alignment/BQSR/{sample}/Mosdepth/{sample}.mosdepth.summary.txt",
    params:
        prefix="02.Alignment/BQSR/{sample}/Mosdepth/{sample}"
    log:
        out = snakedir+"/logs/B14.CRAM_TO_MOSDEPTH/{sample}.o",
        err = snakedir+"/logs/B14.CRAM_TO_MOSDEPTH/{sample}.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:100 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "export MOSDEPTH_Q0=NO_COVERAGE \n"
        "export MOSDEPTH_Q1=LOW_COVERAGE \n"
        "export MOSDEPTH_Q2=CALLABLE \n"
        "export MOSDEPTH_Q3=HIGH_COVERAGE \n"
        "export MOSDEPTH_Q4=HIGH_COVERAGE_300 \n"
        "export MOSDEPTH_Q5=HIGH_COVERAGE_1000 \n"
        "{config[singularity]} {config[simg][mosdepth]} "
        "mosdepth"
        "    --threads {threads}"
        "    --fasta {config[references][fasta]}"
        "    --quantize 0:1:5:150:300:1000:"
        "    {params.prefix}"
        "    {input.cram}"
        " >> {log.out} 2>> {log.err}\n"
        
##################################################################
############################ Germline ############################
##################################################################

############## Germline: Haplotype Caller ##############

### HaplotypeCaller by intervals ###
rule HCitv:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        gvcf=temp("03.Germline.GATK/{sample}/itvs/{sample}.{itv}.vcf.gz"),
        bam=temp("03.Germline.GATK/{sample}/itvs/{sample}.{itv}.vcf.bam"),
    log:
        out = snakedir+"/logs/C1.HCitv/{sample}.{itv}.o",
        err = snakedir+"/logs/C1.HCitv/{sample}.{itv}.e",
    threads:  4
    resources:
        mem  = '24g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][gatk]} "
        "gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=2\" "
        "    HaplotypeCaller" 
        "    --reference {config[references][fasta]}" 
        "    --input {input.cram} "
        "    --output {output.gvcf} "
        "    --dbsnp {config[references][gatk_dbsnp]}"
        "    --intervals {config[references][hcitv]}/temp_{wildcards.itv}_of_50/scattered.interval_list"
        "    --bam-output {output.bam}"
        "    --tmp-dir /lscratch/$SLURM_JOB_ID/"
        " >> {log.out} 2>> {log.err}\n"
        

### Merge HaplotypeCaller intervals ###
rule HCmerge:
    input:
        vcf=expand("03.Germline.GATK/{{sample}}/itvs/{{sample}}.{itv}.vcf.gz",itv=itv50),
        bam=expand("03.Germline.GATK/{{sample}}/itvs/{{sample}}.{itv}.vcf.bam",itv=itv50),
    params:
        tmpbam="/lscratch/$SLURM_JOB_ID/{sample}.vcf.bamout.bam"
    output:
        gvcf="03.Germline.GATK/{sample}/{sample}.vcf.gz",
        gbam="03.Germline.GATK/{sample}/{sample}.vcf.bam",
        gmet="03.Germline.GATK/{sample}/{sample}.vcf.gz.variant_calling_summary_metrics",
    log:
        out = snakedir+"/logs/C2.HCmerge/{sample}.o",
        err = snakedir+"/logs/C2.HCmerge/{sample}.e",
    threads:  4
    resources:
        mem  = '24g',
        extra = ' --gres=lscratch:200 ',
    run:
        inputvcfs = " ".join("-I {}".format(in_) for in_ in input.vcf),
        inputbams = " ".join(" {}".format(in_) for in_ in input.bam),
        shell(
            "module load singularity > {log.out} 2> {log.err}\n"
            " {config[singularity]} {config[simg][gatk]} "
            "gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=2\" "
            "    MergeVcfs "
            "    -O {output.gvcf}"
            "    {inputvcfs}"
            " >> {log.out} 2>> {log.err}\n"
            "{config[singularity]} {config[simg][sambamba]} "
            "sambamba merge"
            "    -t {threads} "
            "    {params.tmpbam} "
            "    {inputbams}"
            " >> {log.out} 2>> {log.err}\n"
            "{config[singularity]} {config[simg][sambamba]} "
            "sambamba sort -t {threads}"
            "    --tmpdir /lscratch/$SLURM_JOB_ID "
            "    -o {output.gbam} "
            "    {params.tmpbam} "
            " >> {log.out} 2>> {log.err}\n"
            " {config[singularity]} {config[simg][gatk]} "
            "gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=2\" "
            "    CollectVariantCallingMetrics"
            "    -I {output.gvcf}"
            "    -O {output.gvcf}"
            "    --DBSNP {config[references][gatk_dbsnp]}"
            "    -SD {config[references][fadict]}"
            "    --GVCF_INPUT"
            "    --THREAD_COUNT {threads}"
            "    -TI {config[references][wgscallingitv]}"
            " >> {log.out} 2>> {log.err}\n"
        )
            
            
            
############## Germline: Freebayes ##############
            
rule Freebayes:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        vgz="03.Germline.freebayes/{sample}/{sample}.freebayes.vcf.gz",
    log:
        out = snakedir+"/logs/C3.freebayes/{sample}.o",
        err = snakedir+"/logs/C3.freebayes/{sample}.e",
    threads:  32
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:100 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][freebayes]} "
        "freebayes-parallel <("
        "{config[singularity]} {config[simg][freebayes]} "
        "    fasta_generate_regions.py"
        "    {config[references][fasta]}.fai"
        "    100000)"
        "    {threads}"
        "    -f {config[references][fasta]}"
        "    {input.cram} 2>> {log.err}|"
        "{config[singularity]} {config[simg][samtools]} "
        "bgzip > {output.vgz} 2>> {log.err}\n"
        "{config[singularity]} {config[simg][samtools]} "
        "tabix"
        "    -p vcf"
        "    {output.vgz}"
        " >> {log.out} 2>> {log.err}\n"
        
        
        
############## Germline: Deep Variant ##############

# DeepVariant Step 1
rule DeepVariant1:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        dv1a=temp( "03.Germline.deepvariant/{sample}/temp/gvcf.tfrecord-{itv}-of-%s.gz"%(itvdv[-1]) ),
        dv1b=temp( "03.Germline.deepvariant/{sample}/temp/make_examples.tfrecord-{itv}-of-%s.gz"%(itvdv[-1]) ),
        dv1c=temp( "03.Germline.deepvariant/{sample}/temp/make_examples.tfrecord-{itv}-of-%s.gz.example_info.json"%(itvdv[-1]) ),
    params:
        dv1a="03.Germline.deepvariant/{sample}/temp/gvcf.tfrecord@%d.gz"%(config['parameter']['deepvariant_shards']),
        dv1b="03.Germline.deepvariant/{sample}/temp/make_examples.tfrecord@%d.gz"%(config['parameter']['deepvariant_shards']),
    log:
        out = snakedir+"/logs/C4a.DeepVariant1/{sample}.{itv}.o",
        err = snakedir+"/logs/C4a.DeepVariant1/{sample}.{itv}.e",
    threads:  2
    resources:
        mem  = '4g',
        extra = ' --gres=lscratch:4',
    run:
        print(wildcards)
        itvi = int(wildcards.itv)
        shell(
            "module load singularity > {log.out} 2> {log.err}\n"
            "{config[singularity]} {config[simg][deepvariant]} "
            "/opt/deepvariant/bin/make_examples"
            "    --mode calling"
            "    --ref {config[references][fasta]}"
            "    --reads {input.cram}"
            "    --gvcf {params.dv1a}"
            "    --examples {params.dv1b}"
            "    --channels insert_size"
            "    --task {itvi}"
            " >> {log.out} 2>> {log.err}\n"
        )

# DeepVariant Step 2
rule DeepVariant2:
    input:
        dv1a = lambda wildcards: ["03.Germline.deepvariant/{sample}/temp/gvcf.tfrecord-{itv}-of-{all}.gz".format(sample=wildcards.sample, itv=i, all=itvdv[-1]) for i in itvdv[:-1]],
        dv1b = lambda wildcards: ["03.Germline.deepvariant/{sample}/temp/make_examples.tfrecord-{itv}-of-{all}.gz".format(sample=wildcards.sample, itv=i, all=itvdv[-1]) for i in itvdv[:-1]],
        dv1c = lambda wildcards: ["03.Germline.deepvariant/{sample}/temp/make_examples.tfrecord-{itv}-of-{all}.gz.example_info.json".format(sample=wildcards.sample, itv=i, all=itvdv[-1]) for i in itvdv[:-1]],
    output:
        call = temp( "03.Germline.deepvariant/{sample}/temp/call_variants_output.tfrecord.gz"),
    params:
        dv1b = "03.Germline.deepvariant/{sample}/temp/make_examples.tfrecord@%d.gz"%(config['parameter']['deepvariant_shards']),
    log:
        out = snakedir+"/logs/C4b.DeepVariant2/{sample}.o",
        err = snakedir+"/logs/C4b.DeepVariant2/{sample}.e",
    threads:  12
    resources:
        mem  = '24g',
        extra = ' --gres=lscratch:40,gpu:a100:1',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} --nv {config[simg][deepvariant]} "
        "/opt/deepvariant/bin/call_variants"
        "    --outfile {output.call}"
        "    --examples {params.dv1b}"
        "    --checkpoint /opt/models/wgs/model.ckpt"
        " >> {log.out} 2>> {log.err}\n"

# DeepVariant Step 3
rule DeepVariant3:
    input:
        dv1a = lambda wildcards: ["03.Germline.deepvariant/{sample}/temp/gvcf.tfrecord-{itv}-of-{all}.gz".format(sample=wildcards.sample, itv=i, all=itvdv[-1]) for i in itvdv[:-1]],
        call = "03.Germline.deepvariant/{sample}/temp/call_variants_output.tfrecord.gz",
    output:
        vcf  = "03.Germline.deepvariant/{sample}/{sample}.deepvariant.vcf.gz",
        gvcf = "03.Germline.deepvariant/{sample}/{sample}.deepvariant.gvcf.gz",
    params:
        vcf  = "03.Germline.deepvariant/{sample}/{sample}.deepvariant.vcf",
        gvcf = "03.Germline.deepvariant/{sample}/{sample}.deepvariant.gvcf",
        dv1a="03.Germline.deepvariant/{sample}/temp/gvcf.tfrecord@%d.gz"%(config['parameter']['deepvariant_shards']),
    log:
        out = snakedir+"/logs/C4c.DeepVariant3/{sample}.o",
        err = snakedir+"/logs/C4c.DeepVariant3/{sample}.e",
    threads:  4
    resources:
        mem  = '48g',
        extra = ' --gres=lscratch:40',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][deepvariant]} "
        "/opt/deepvariant/bin/postprocess_variants"
        "    --ref {config[references][fasta]}"
        "    --infile {input.call}"
        "    --nonvariant_site_tfrecord_path {params.dv1a}"
        "    --outfile {params.vcf}"
        "    --gvcf_outfile {params.gvcf}"
        " >> {log.out} 2>> {log.err}\n"
        
        "{config[singularity]} {config[simg][samtools]} "
        "bgzip {params.vcf}"
        " >> {log.out} 2>> {log.err}\n"
        
        "{config[singularity]} {config[simg][samtools]} "
        "bgzip {params.gvcf}"
        " >> {log.out} 2>> {log.err}\n"
        
        "{config[singularity]} {config[simg][samtools]} "
        "tabix -p vcf {output.vcf}"
        " >> {log.out} 2>> {log.err}\n"
        
        "{config[singularity]} {config[simg][samtools]} "
        "tabix -p vcf  {output.gvcf}"
        " >> {log.out} 2>> {log.err}\n"

        
############## Germline: Strelka ##############

rule STRELKA:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    params:
        dir="03.Germline.Strelka/{sample}",
    output:
        vgz="03.Germline.Strelka/{sample}/results/variants/variants.vcf.gz",
    log:
        out = snakedir+"/logs/C5.STRELKA/{sample}.o",
        err = snakedir+"/logs/C5.STRELKA/{sample}.e",
    threads:  24
    resources:
        mem  = '48g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][strelka]} "
        "configureStrelkaGermlineWorkflow.py"
        "    --bam {input.cram}"
        "    --referenceFasta {config[references][fasta]}"
        "    --runDir {params.dir}"
        "    >> {log.out} 2>> {log.err}\n"
        "cd {params.dir} \n"
        "{config[singularity]} {config[simg][strelka]} "
        "./runWorkflow.py"
        "    -m local"
        "    -j {threads}"
        "    >> {log.out} 2>> {log.err}\n"
        "    rm -rf workspace"
        "    >> {log.out} 2>> {log.err}\n"


rule MANTA:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    params:
        dir="03.Germline.Manta/{sample}",
    output:
        vgz="03.Germline.Manta/{sample}/results/variants/candidateSV.vcf.gz",
    log:
        out = snakedir+"/logs/C6.MANTA/{sample}.o",
        err = snakedir+"/logs/C6.MANTA/{sample}.e",
    threads:  24
    resources:
        mem  = '48g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][manta]} "
        "configManta.py"
        "    --bam {input.cram}"
        "    --reference {config[references][fasta]}"
        "    --runDir {params.dir}"
        "    >> {log.out} 2>> {log.err}\n"
        "cd {params.dir} \n"
        "{config[singularity]} {config[simg][manta]} "
        "./runWorkflow.py"
        "    -m local"
        "    -j {threads}"
        "    >> {log.out} 2>> {log.err}\n"
        "    rm -rf workspace"
        "    >> {log.out} 2>> {log.err}\n"
        
rule TIDDIT:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    params:
        prefix="03.Germline.Tiddit/{sample}/{sample}.tiddit",
    output:
        vgz="03.Germline.Tiddit/{sample}/{sample}.tiddit.vcf.gz",
    log:
        out = snakedir+"/logs/C7.TIDDIT/{sample}.o",
        err = snakedir+"/logs/C7.TIDDIT/{sample}.e",
    threads:  24
    resources:
        mem  = '48g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][tiddit]} "
        "tiddit"
        "    --sv"
        "    --threads {threads}"
        "    --bam {input.cram}"
        "    --ref {config[references][fasta]}"
        "    -o {params.prefix}"
        "    >> {log.out} 2>> {log.err}\n"
        "{config[singularity]} {config[simg][samtools]} "
        "bgzip {params.prefix}.vcf"
        "    >> {log.out} 2>> {log.err}\n"
        "{config[singularity]} {config[simg][samtools]} "
        "tabix -p vcf {output.vgz}"
        "    >> {log.out} 2>> {log.err}\n"

        
rule GRIDSS_preprocess:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    params:
        workspace="04.Germline.Gridss/{sample}/_gridss",
    output:
        ok="04.Germline.Gridss/{sample}/_gridss/preprocess.ok",
    log:
        out = snakedir+"/logs/C8.GRIDSS/pre.{sample}.o",
        err = snakedir+"/logs/C8.GRIDSS/pre.{sample}.e",
    threads:  8
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][gridss]} "
        "gridss "
        "  -s preprocess "
        "  -r {config[references][gridssfa]}"
        "  -t {threads} "
        "  -w {params.workspace} "
        "  {input.cram}"
        "  -b {config[references][blacklist]}"
        "  >> {log.out} 2>> {log.err}\n"
        "touch {output.ok}"
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.workspace}"
        "  >> {log.out} 2>> {log.err}\n"
        
rule GRIDSS_assemble:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
        ok="04.Germline.Gridss/{sample}/_gridss/preprocess.ok",
    params:
        workspace="04.Germline.Gridss/{sample}/_gridss",
        bam="04.Germline.Gridss/{sample}/{sample}.assemble.bam",
    output:
        ok ="04.Germline.Gridss/{sample}/_gridss/assemble_{shard}.ok",
    log:
        out = snakedir+"/logs/C8.GRIDSS/assemble.{sample}.{shard}.o",
        err = snakedir+"/logs/C8.GRIDSS/assemble.{sample}.{shard}.e",
    threads:  8
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][gridss]} "
        "gridss "
        "  -s assemble "
        "  -r {config[references][gridssfa]}"
        "  -t {threads}"
        "  -a {params.bam}"
        "  --jobnodes {config[parameter][gridss_shards]}"
        "  --jobindex {wildcards.shard}"
        "  {input.cram}"
        "  -w {params.workspace}"
        "  >> {log.out} 2>> {log.err}\n"
        "touch {output.ok}"
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.workspace}"
        "  >> {log.out} 2>> {log.err}\n"

rule GRIDSS_call:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
        ok = lambda wildcards: ["04.Germline.Gridss/{sample}/_gridss/assemble_{shard}.ok".format(sample=wildcards.sample, shard=i) for i in range(config['parameter']['gridss_shards'])],
    params:
        workspace="04.Germline.Gridss/{sample}/_gridss",
    output:
        bam="04.Germline.Gridss/{sample}/{sample}.assemble.bam",
        vcf="04.Germline.Gridss/{sample}/{sample}.gridss.vcf",
    log:
        out = snakedir+"/logs/C8.GRIDSS/call.{sample}.o",
        err = snakedir+"/logs/C8.GRIDSS/call.{sample}.e",
    threads:  8
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][gridss]} "
        "gridss "
        "  -s assemble,call "
        "  -r {config[references][gridssfa]}"
        "  -t {threads}"
        "  -a {output.bam}"
        "  -o {output.vcf}"
        "  -w {params.workspace}"
        "  {input.cram}"
        "  >> {log.out} 2>> {log.err}\n"

rule MELT_Ins:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    params:
        prefix="06.MEI.MELT/{sample}/ME_Insert",
    output:
        vcf="06.MEI.MELT/{sample}/ME_Insert/LINE1.final_comp.vcf",
    log:
        out = snakedir+"/logs/E1.MELT/ins.{sample}.o",
        err = snakedir+"/logs/E1.MELT/ins.{sample}.e",
    threads:  30
    resources:
        mem  = '72g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][melt]} "
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Single "
        "  -a "
        "  -c {threads} "
        "  -h {config[references][fasta]} "
        "  -bamfile {input.cram} "
        "  -n /opt/MELT/add_bed_files/Hg38/Hg38.genes.bed "
        "  -t /opt/MELT/me_refs/transposon_file_list.txt "
        "  -w {params.prefix}"
        "  >> {log.out} 2>> {log.err}\n"

rule MELT_Del1:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    params:
        prefix1="06.MEI.MELT/{sample}/ME_Deletion_LINE1",
    output:
        vcf1="06.MEI.MELT/{sample}/ME_Deletion_LINE1/DEL.final_comp.vcf",
    log:
        out = snakedir+"/logs/E1.MELT/del_line1.{sample}.o",
        err = snakedir+"/logs/E1.MELT/del_line1.{sample}.e",
    threads: 30
    resources:
        mem  = '72g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][melt]} "
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Deletion-Genotype"
        "  -bamfile {input.cram} "
        "  -w {params.prefix}/workspace"
        "  -bed /opt/MELT/add_bed_files/Hg38/LINE1.deletion.bed "
        "  -h {config[references][fasta]} "
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.prefix1}/workspace/*.del.tsv > /lscratch/$SLURM_JOB_ID/del_LINE1.output.list"
        "  2>> {log.err}\n"
        "{config[singularity]} {config[simg][melt]} "
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Deletion-Genotype"
        "  -mergelist /lscratch/$SLURM_JOB_ID/del_LINE1.output.list"
        "  -bed /opt/MELT/add_bed_files/Hg38/LINE1.deletion.bed"
        "  -h {config[references][fasta]} "
        "  -o params.prefix1"
        "  >> {log.out} 2>> {log.err}\n"
        
rule MELT_Del2:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    params:
        prefix="06.MEI.MELT/{sample}/ME_Deletion_AluY",
    output:
        vcf1="06.MEI.MELT/{sample}/ME_Deletion_AluY/DEL.final_comp.vcf",
    log:
        out = snakedir+"/logs/E1.MELT/del_aluy.{sample}.o",
        err = snakedir+"/logs/E1.MELT/del_aluy.{sample}.e",
    threads: 30
    resources:
        mem  = '72g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][melt]} "
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Deletion-Genotype"
        "  -bamfile {input.cram} "
        "  -w {params.prefix}/workspace"
        "  -bed /opt/MELT/add_bed_files/Hg38/AluY.deletion.bed "
        "  -h {config[references][fasta]} "
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.prefix1}/workspace/*.del.tsv > /lscratch/$SLURM_JOB_ID/del_AluY.output.list"
        "  2>> {log.err}\n"
        "{config[singularity]} {config[simg][melt]} "
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Deletion-Genotype"
        "  -mergelist /lscratch/$SLURM_JOB_ID/del_AluY.output.list"
        "  -bed /opt/MELT/add_bed_files/Hg38/AluY.deletion.bed"
        "  -h {config[references][fasta]} "
        "  -o params.prefix1"
        "  >> {log.out} 2>> {log.err}\n"
        
        
rule MSISENSORPRO:
    input:
        cram="02.Alignment/BQSR/{sample}/{sample}.BQSR.cram",
    output:
        prefix="06.MSI.MSISENSORPRO/{sample}/{sample}.msisensorpro",
    log:
        out = snakedir+"/logs/E1.MELT/del_aluy.{sample}.o",
        err = snakedir+"/logs/E1.MELT/del_aluy.{sample}.e",
    threads: 8
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:50 ',
    shell:
        "module load singularity > {log.out} 2> {log.err}\n"
        "{config[singularity]} {config[simg][msisensorpro]} "
        "msisensor-pro "
        "  pro"
        "  -d {config[references][msipon]} "
        "  -g {config[references][fasta]} "
        "  -t {input.cram} "
        "  -o {output.prefix}"
        "  -b {threads}"
        "  >> {log.out} 2>> {log.err}\n"
        