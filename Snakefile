from os.path import join
import os
import pandas as pd
from scripts.Load import samplesheet

snakedir = os.getcwd()
print(snakedir)
configfile: 'config.yaml'
print(config)
sampledic, libdic, rundic = samplesheet(config['samplesheet'])
workdir: config['workdir']
itv4 = ['%.4d'%(itv) for itv in range(1, int(config['references']['interval'])+1)]
itv50 = ['%.4d'%(itv) for itv in range(1, 50+1)]

## Global wildcards
#SAMPLES, = glob_wildcards(sampledic.keys()).sample

rule all:
    input:
        expand("02.Alignment/Level3/{sample}/{sample}.sort.md.bam", sample = sampledic.keys()),
        expand("02.Alignment/chrM/{sample}/{sample}.bam", sample = sampledic.keys()),
        expand("02.Alignment/Level3/{sample}/{sample}.BQSR.metrics", sample = sampledic.keys()),
        expand("02.Alignment/Level3/{sample}/{sample}.BQSR.bam.flagstat", sample = sampledic.keys()),
        expand("02.Alignment/Level3/{sample}/Callable/{sample}.callable.bed", sample = sampledic.keys()),
        
        expand("03.Germline.chrM/{sample}/{sample}.vcf", sample = sampledic.keys()), # mtoolbox
        expand("03.Germline.GATK/{sample}/{sample}.g.vcf.gz", sample = sampledic.keys()), #HaplotypeCaller
        expand("03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf.gz.anno/Merge.Anno.matrix.gz", sample=sampledic.keys()), #Lofreq
        expand("03.Germline.freebayes/{sample}/{sample}.freebayes.flt.vcf.gz.anno/Merge.Anno.matrix.gz", sample=sampledic.keys()), #(freebayes)
        expand("03.Germline.Strelka/{sample}/results/variants/variants.vcf.gz.anno/Merge.Anno.matrix.gz", sample=sampledic.keys()), #(strelka)
        
        
        expand("04.SV.Delly/{sample}/{sample}.delly.vcf.gz", sample=sampledic.keys()), #(Delly)
        expand("04.SV.Gridss/{sample}/{sample}.gridss.vcf.gz", sample=sampledic.keys()), #(Gridss)
        
        expand("02.Alignment/WindowDepth5K/{sample}/{sample}.bin.5K.depth.gz", sample=sampledic.keys()), #(WindowDepth)
        expand("05.CNV.Canvas/SingleSample/{sample}/CNV.vcf.gz", sample=sampledic.keys()), #(Canvas)
        
        expand("06.MEI.Scramble/{sample}/{sample}.mei.txt", sample=sampledic.keys()), #(Scramble)

        

rule QC:
    input:
        reads1=lambda wildcards: rundic[wildcards.run]['Read1'],
        reads2=lambda wildcards: rundic[wildcards.run]['Read2'],
    output:
        reads1out=temp(expand("01.CleanData/{{run}}/{itv}.{{run}}.R1.cln.fq.gz", itv=itv4)),
        reads2out=temp(expand("01.CleanData/{{run}}/{itv}.{{run}}.R2.cln.fq.gz", itv=itv4)),
        htmlout="01.CleanData/{run}/{run}.QC.html",
        jsonout="01.CleanData/{run}/{run}.QC.json",
    params:
        reads1out="01.CleanData/{run}/{run}.R1.cln.fq.gz",
        reads2out="01.CleanData/{run}/{run}.R2.cln.fq.gz",
    log: 
        out = snakedir + "/logs/A1.QC/{run}.o",
        err = snakedir + "/logs/A1.QC/{run}.e",
    threads: 10
    resources: 
        mem = '4g',
        extra = ' --gres=lscratch:10 '
    shell:
        "module load fastp &&"
        "/data/yuk5/script/fastp "
        " -i {input.reads1}" 
        " -I {input.reads2}" 
        " -o {params.reads1out}" 
        " -O {params.reads2out}" 
        " -h {output.htmlout}" 
        " -j {output.jsonout}" 
        " -s {config[references][interval]}"
        " -w {threads} > {log.out} 2> {log.err}"
        
rule Bwa_mem:
    input:
        read1="01.CleanData/{run}/{itv}.{run}.R1.cln.fq.gz",
        read2="01.CleanData/{run}/{itv}.{run}.R2.cln.fq.gz",
    output:
        bam=temp("02.Alignment/Level1/{run}/{itv}.{run}.sort.bam"),
        bai=temp("02.Alignment/Level1/{run}/{itv}.{run}.sort.bai"),
    log:
        out = snakedir+"/logs/B1.bwa/{run}.{itv}.o",
        err = snakedir+"/logs/B1.bwa/{run}.{itv}.e"
    params:
        bwa=lambda wildcards:' -K 10000000 -R "@RG\\tID:{}\\tLB:{}\\tPL:illumina\\tPU:{}\\tSM:{}" '.format(
                                wildcards.run, rundic[wildcards.run]['LB'],  wildcards.run,  rundic[wildcards.run]['SM']),
        samblaster="",
        sambambasort=" --tmpdir /lscratch/$SLURM_JOB_ID "
    threads: 32
    resources: 
        mem = '32g',
        extra = ' --gres=lscratch:200 '
    shell:
        "module load {config[modules][bwa]} {config[modules][gatk]} {config[modules][sambamba]} \n"
        "bwa mem -t {threads} "
        " {params.bwa} "
        " {config[references][bwaidx]} "
        " {input.read1} {input.read2} 2> {log.err} |"
        "gatk SortSam --java-options -Xmx30g "
        "  --MAX_RECORDS_IN_RAM 5000000 "
        "  -I /dev/stdin "
        "  --SORT_ORDER coordinate "
        "  --TMP_DIR /lscratch/$SLURM_JOB_ID/ "
        "  -O /dev/stdout  2>> {log.err} |"
        "gatk SetNmMdAndUqTags --java-options -Xmx30g "
        "  --INPUT /dev/stdin "
        "  --OUTPUT {output.bam} "
        "  --REFERENCE_SEQUENCE {config[references][fasta]}"
        "  --CREATE_INDEX true >> {log.out} 2>> {log.err}"
        
rule Markdup:
    input:
        bam = lambda wildcards: expand("02.Alignment/Level1/{run}/{itv}.{run}.sort.bam", itv=itv4, run=libdic[wildcards.lib]),
        bai = lambda wildcards: expand("02.Alignment/Level1/{run}/{itv}.{run}.sort.bai", itv=itv4, run=libdic[wildcards.lib]),
    output:
        bam    =temp("02.Alignment/Level2/{lib}/{lib}.sort.md.bam"),
        bai    =temp("02.Alignment/Level2/{lib}/{lib}.sort.md.bam.bai"),
        metrics="02.Alignment/Level2/{lib}/{lib}.sort.md.metrics",
    log:
        out = snakedir+"/logs/B2.mkdup/{lib}.o",
        err = snakedir+"/logs/B2.mkdup/{lib}.e",
    threads:  4
    resources:
        mem = '32g',
        extra = ' --gres=lscratch:200 ',
    run:
        inputs = " ".join("-I {}".format(in_) for in_ in input.bam),
        shell(
            """
            module load {config[modules][gatk]} {config[modules][sambamba]}
            gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms30G -Xmx30G -XX:ParallelGCThreads=2" MarkDuplicates  \
              {inputs} \
              -O {output.bam} \
              -M {output.metrics} \
              --MAX_RECORDS_IN_RAM 2000000 \
              --VALIDATION_STRINGENCY SILENT \
              --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 > {log.out} 2> {log.err}
            sambamba index -t {threads} {output.bam} >> {log.out} 2>> {log.err}
            """)
        
rule Merge_level3:
    input:
        bam = lambda wildcards: ["02.Alignment/Level2/{}/{}.sort.md.bam".format(lib, lib) for lib in sampledic[wildcards.sample]],
        bai = lambda wildcards: ["02.Alignment/Level2/{}/{}.sort.md.bam.bai".format(lib, lib) for lib in sampledic[wildcards.sample]],
    output:
        bam=temp("02.Alignment/Level3/{sample}/{sample}.sort.md.bam"),
        bai=temp("02.Alignment/Level3/{sample}/{sample}.sort.md.bam.bai"),
    log:
        out = snakedir+"/logs/B3.merge_level3/{sample}.o",
        err = snakedir+"/logs/B3.merge_level3/{sample}.e",
    threads:  8
    resources:
        mem = '16g',
        extra = ' --gres=lscratch:10 ',
    run:
        if len(input.bam) > 1:
            inputs = " ".join(input.bam),
            shell("""
            module load {config[modules][sambamba]}
            sambamba merge \
                -t {threads} \
                {output.bam} \
                {inputs} > {log.out} 2> {log.err}
            """)
        else:
            shell("""
            module load {config[modules][sambamba]}
            mv {input.bam[0]} {output.bam} > {log.out} 2> {log.err}
            sambamba index \
                -t {threads} \
                {output.bam} >>{log.out} 2>> {log.err}
            """)
        
rule chrMbam:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.sort.md.bam",
        bai="02.Alignment/Level3/{sample}/{sample}.sort.md.bam.bai",
    output:
        mbam=protected("02.Alignment/chrM/{sample}/{sample}.bam"),
    log:
        out = snakedir+"/logs/B4.chrMbam/{sample}.o",
        err = snakedir+"/logs/B4.chrMbam/{sample}.e",
    threads:  8
    resources:
        mem = '16g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        module load {config[modules][samtools]} {config[modules][sambamba]}
        samtools view -Sb -h {input.bam} chrM > {output.mbam} 2>> {log.err}
        sambamba index -t {threads} {output.mbam} >>{log.out} 2>> {log.err}
        """

rule BQSR:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.sort.md.bam",
        bai="02.Alignment/Level3/{sample}/{sample}.sort.md.bam.bai",
    params:
        itvbed="/data/yuk5/pipeline/wgs_germline/ref/hg38_chr_intervals/itv_{itv}.bed"
    output:
        metrics=temp("02.Alignment/Level3/{sample}/{sample}.{itv}.BQSR.metrics"),
        bam=temp("02.Alignment/Level3/{sample}/{sample}.{itv}.BQSR.bam"),
        bai=temp("02.Alignment/Level3/{sample}/{sample}.{itv}.BQSR.bai"),
    log:
        out = snakedir+"/logs/B5.BQSR/{sample}.{itv}.o",
        err = snakedir+"/logs/B5.BQSR/{sample}.{itv}.e",
    threads:  2
    resources:
        mem  = '8g',
        extra = ' --gres=lscratch:400 ',
    shell:
        """
        module load {config[modules][gatk]}
        gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms4G -Xmx4G -XX:ParallelGCThreads=2" BaseRecalibrator \
            -R {config[references][fasta]} \
            -I {input.bam} \
            -O {output.metrics} \
            --use-original-qualities \
            --known-sites {config[references][gatk_dbsnp]} \
            --known-sites {config[references][gatk_1000g]} \
            --known-sites {config[references][gatk_indel]} \
            --intervals   {params.itvbed} >> {log.out} 2>> {log.err}
        gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyBQSR \
            --add-output-sam-program-record \
            -R {config[references][fasta]} \
            -I {input.bam} \
            -O {output.bam} \
            --use-original-qualities \
            -bqsr {output.metrics} \
             --static-quantized-quals 10 \
             --static-quantized-quals 20 \
             --static-quantized-quals 30 \
            -L {params.itvbed} >> {log.out} 2>> {log.err}
        """
        
rule BQSR_mergeM:
    input:
        expand("02.Alignment/Level3/{{sample}}/{{sample}}.{itv}.BQSR.metrics", itv=itv4),
    output:
        metrics="02.Alignment/Level3/{sample}/{sample}.BQSR.metrics",
    log:
        out = snakedir+"/logs/B5.BQSR/{sample}.M.o",
        err = snakedir+"/logs/B5.BQSR/{sample}.M.e",
    threads:  2
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:40 ',
    run:
        inputs = " ".join("-I {}".format(in_) for in_ in input),
        shell("""
        module load {config[modules][gatk]} 
        gatk --java-options "-Xms3000m" \
          GatherBQSRReports \
          {inputs} \
          -O {output.metrics} >> {log.out} 2>> {log.err}
        """)

rule BQSR_mergeB:
    input:
        expand("02.Alignment/Level3/{{sample}}/{{sample}}.{itv}.BQSR.bam", itv=itv4),
    output:
        bam=protected("02.Alignment/Level3/{sample}/{sample}.BQSR.bam"),
    log:
        out = snakedir+"/logs/B5.BQSR/{sample}.B.o",
        err = snakedir+"/logs/B5.BQSR/{sample}.B.e",
    threads:  16
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:40 ',
    run:
        inputs = " ".join("{}".format(in_) for in_ in input),
        shell(
        """
        module load {config[modules][sambamba]}
        sambamba merge \
          -t {threads} \
          {output.bam} \
          {inputs} >> {log.out} 2>> {log.err}
        """
        )


rule stat_BQSR:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    params:
        stat1="02.Alignment/Level3/{sample}/stat/{sample}.stat1",
        stat2="02.Alignment/Level3/{sample}/stat/{sample}.stat2",
    output:
        statdir=directory("02.Alignment/Level3/{sample}/stat"),
        flags="02.Alignment/Level3/{sample}/{sample}.BQSR.bam.flagstat",
    log:
        out = snakedir+"/logs/B6.BQSRstat/{sample}.o",
        err = snakedir+"/logs/B6.BQSRstat/{sample}.e",
    threads:  4
    resources:
        mem  = '8g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        module load {config[modules][picard]} {config[modules][sambamba]} {config[modules][samtools]}
        mkdir -p {output.statdir}
        java -Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR CollectMultipleMetrics \
            INPUT={input.bam} \
            REFERENCE_SEQUENCE={config[references][fasta]} \
            OUTPUT={params.stat1} \
            ASSUME_SORTED=true \
            PROGRAM="null" \
            PROGRAM="CollectAlignmentSummaryMetrics" \
            PROGRAM="CollectInsertSizeMetrics" \
            PROGRAM="CollectSequencingArtifactMetrics" \
            PROGRAM="CollectGcBiasMetrics" \
            PROGRAM="QualityScoreDistribution" \
            METRIC_ACCUMULATION_LEVEL="null" \
            METRIC_ACCUMULATION_LEVEL="SAMPLE" \
            METRIC_ACCUMULATION_LEVEL="LIBRARY" >> {log.out} 2>> {log.err}
        java -Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR CollectMultipleMetrics \
            INPUT={input.bam} \
            REFERENCE_SEQUENCE={config[references][fasta]} \
            OUTPUT={params.stat2} \
            ASSUME_SORTED=true \
            PROGRAM="null" \
            PROGRAM="CollectAlignmentSummaryMetrics" \
            PROGRAM="CollectGcBiasMetrics" \
            METRIC_ACCUMULATION_LEVEL="null" \
            METRIC_ACCUMULATION_LEVEL="READ_GROUP"  >> {log.out} 2>> {log.err}
        sambamba flagstat -t {threads} {input.bam} > {output.flags} 2>> {log.err}
        """ 
        
rule CallableLoci:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        bed="02.Alignment/Level3/{sample}/Callable/{sample}.callable.bed",
    params:
        bed="02.Alignment/Level3/{sample}/Callable/{sample}.bed",
        sum="02.Alignment/Level3/{sample}/Callable/{sample}.summary",
    log:
        out = snakedir+"/logs/B8.Callable/{sample}.o",
        err = snakedir+"/logs/B8.Callable/{sample}.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        java -Xms8G -Xmx8G -XX:ParallelGCThreads=2 -jar {config[bins][gatk3]}  -T CallableLoci \
          -R {config[references][fasta]} \
          -L {config[references][wgscallingitv]} \
          -I {input.bam} \
          --maxDepth 1000 \
          --minBaseQuality 10 \
          --minMappingQuality 10 \
          --minDepth 5 \
          --minDepthForLowMAPQ 20 \
          --summary {params.sum} \
          -o {params.bed}  > {log.out} 2> {log.err}
        cat {params.bed} 2>> {log.err}|grep CALLABLE > {output.bed} 2>> {log.err}
        gzip {params.bed} 2>> {log.err}
        """

rule Mtoolbox:
    input:
        bam="02.Alignment/chrM/{sample}/{sample}.bam",
    output:
        config="03.Germline.chrM/{sample}/config",
        path=directory("03.Germline.chrM/{sample}"),
        vcf="03.Germline.chrM/{sample}/{sample}.vcf",
    log:
        out = snakedir+"/logs/C6.mtoolbox/{sample}.o",
        err = snakedir+"/logs/C6.mtoolbox/{sample}.e",
    threads:  4
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        module load {config[modules][python27]}
        mkdir -p {output.path}/OUT_{wildcards.sample}
        set +eu
        source /home/yuk5/yuk5/app/anaconda3/etc/profile.d/conda.sh >> {log.out} 2>> {log.err}
        conda activate mtoolbox >> {log.out} 2>> {log.err}
        source /home/yuk5/yuk5/app/anaconda3/envs/mtoolbox/MToolBox-1.2.1/setenv.sh >> {log.out} 2>> {log.err}
        set -eu
        {snakedir}/scripts/write_mtoolbox_config.py {input.bam} {output.path} {wildcards.sample} > {output.config} 2>> {log.err}
        MToolBox.sh -i {output.config} >> {log.out} 2>> {log.err}
        rm -rf {output.path}/tmp
        """

### HaplotypeCaller by intervals ###
rule HCitv:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    params:
        itv="/data/yuk5/pipeline/wgs_germline/ref/gatk_bundle_hg38_v0/scattered_calling_intervals/temp_{itv}_of_50/scattered.interval_list"
    output:
        gvcf=temp("03.Germline.GATK/{sample}/itvs/{sample}.{itv}.g.vcf.gz"),
        bam=temp("03.Germline.GATK/{sample}/itvs/{sample}.{itv}.g.vcf.bam"),
    log:
        out = snakedir+"/logs/C1.HCitv/{sample}.{itv}.o",
        err = snakedir+"/logs/C1.HCitv/{sample}.{itv}.e",
    threads:  4
    resources:
        mem  = '24g',
        extra = ' --gres=lscratch:50 ',
    shell:
        """
        module load {config[modules][gatk]} 
        gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=2" \
          HaplotypeCaller \
          -R {config[references][fasta]} \
          -O {output.gvcf} \
          -I {input.bam} \
          -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
          -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
          -ERC GVCF \
          -bamout {output.bam} \
          -L {params.itv}
        """

### Merge HaplotypeCaller intervals ###
rule HCmerge:
    input:
        vcf=expand("03.Germline.GATK/{{sample}}/itvs/{{sample}}.{itv}.g.vcf.gz",itv=itv50),
        bam=expand("03.Germline.GATK/{{sample}}/itvs/{{sample}}.{itv}.g.vcf.bam",itv=itv50),
    params:
        tmpbam="/lscratch/$SLURM_JOB_ID/{sample}.g.vcf.bamout.bam"
    output:
        gvcf="03.Germline.GATK/{sample}/{sample}.g.vcf.gz",
        gbam="03.Germline.GATK/{sample}/{sample}.g.vcf.bam",
        gmet="03.Germline.GATK/{sample}/{sample}.g.vcf.gz.variant_calling_summary_metrics",
    log:
        out = snakedir+"/logs/C2.HCmerge/{sample}.o",
        err = snakedir+"/logs/C2.HCmerge/{sample}.e",
    threads:  4
    resources:
        mem  = '24g',
        extra = ' --gres=lscratch:50 ',
    run:
        inputvcfs = " ".join("-I {}".format(in_) for in_ in input.vcf),
        inputbams = " ".join(" {}".format(in_) for in_ in input.bam),
        shell('''
        module load {config[modules][gatk]}  {config[modules][sambamba]}
        gatk --java-options " -Xms20G -Xmx20G  -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -XX:ParallelGCThreads=2" \
          MergeVcfs \
          -O {output.gvcf} \
          {inputvcfs} > {log.out} 2> {log.err}
        sambamba merge \
          -t {threads} \
          {params.tmpbam} \
          {inputbams} >> {log.out} 2>> {log.err}
        sambamba sort -t {threads} \
          -o {output.gbam} \
          {params.tmpbam} >> {log.out} 2>> {log.err}
        gatk --java-options " -Xms20G -Xmx20G  -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -XX:ParallelGCThreads=2" \
          CollectVariantCallingMetrics \
          -I {output.gvcf} \
          -O {output.gvcf} \
          --DBSNP {config[references][gatk_dbsnp]} \
          -SD {config[references][fadict]} \
          --GVCF_INPUT \
          --THREAD_COUNT {threads} \
          -TI {config[references][wgscallingitv]} >> {log.out} 2>> {log.err}
        ''')

# Lofreq 
rule Lofreqindelbam:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        bam=temp("02.Alignment/Lofreq/{sample}/{sample}.li.bam"),
        bai=temp("02.Alignment/Lofreq/{sample}/{sample}.li.bam.bai"),
    log:
        out = snakedir+"/logs/B7.lofreqindelbam/{sample}.o",
        err = snakedir+"/logs/B7.lofreqindelbam/{sample}.e",
    threads:  4
    resources:
        mem  = '8g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        module load {config[modules][lofreq]} {config[modules][sambamba]}
        lofreq indelqual \
          -f {config[references][fasta]} \
          --dindel \
          -o {output.bam} \
          {input.bam} >> {log.out} 2>> {log.err}
        sambamba index -t {threads} {output.bam}  >> {log.out} 2>> {log.err}
        """
rule Lofreq:
    input:
        bam="02.Alignment/Lofreq/{sample}/{sample}.li.bam",
        bai="02.Alignment/Lofreq/{sample}/{sample}.li.bam.bai",
    params:
        vcf="03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf",
    output:
        vgz="03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf.gz",
        tbi="03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf.gz.tbi",
    log:
        out = snakedir+"/logs/C3.lofreq/{sample}.o",
        err = snakedir+"/logs/C3.lofreq/{sample}.e",
    threads:  32
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:50 ',
    shell:
        """
        module load {config[modules][lofreq]} {config[modules][samtools]}
        lofreq call-parallel \
          --pp-threads {threads} \
          -f {config[references][fasta]} \
          -o {params.vcf} \
          --call-indels \
          {input.bam}  >> {log.out} 2>> {log.err}
        bgzip {params.vcf} >> {log.out} 2>> {log.err}
        tabix -p vcf {output.vgz}  >> {log.out} 2>> {log.err}
        """
rule anno_Lofreq:
    input:
        vgz="03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf.gz",
    output:
        mgz="03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf.gz.anno/Merge.Anno.matrix.gz",
    log:
        out = snakedir+"/logs/C3.anno_lofreq/{sample}.o",
        err = snakedir+"/logs/C3.anno_lofreq/{sample}.e",
    threads:  16
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:50 ',
    shell:
        """
        module load {config[modules][vcflib]} {config[modules][samtools]}
        /data/yuk5/pipeline/vcf_annotation/vcf_annotation.v1.wgs.py \
          {input.vgz} \
          {input.vgz}.anno {threads} n > {log.out} 2> {log.err}
        """

rule Freebayes:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        vgz="03.Germline.freebayes/{sample}/{sample}.freebayes.vcf.gz",
    log:
        out = snakedir+"/logs/C4.freebayes/{sample}.o",
        err = snakedir+"/logs/C4.freebayes/{sample}.e",
    threads:  32
    resources:
        mem  = '128g',
        extra = ' --gres=lscratch:100 ',
    shell:
        """
        module load {config[modules][freebayes]} {config[modules][samtools]}
        freebayes-parallel \
          <(fasta_generate_regions.py \
          {config[references][fasta]}.fai \
          100000) \
          32 \
          -f {config[references][fasta]} \
          {input.bam} 2>> {log.err}|\
        bgzip > {output.vgz}
        tabix -f -p vcf {output.vgz} >> {log.out} 2>> {log.err}
        """
        
rule anno_Freebayes:
    input:
        vgz="03.Germline.freebayes/{sample}/{sample}.freebayes.vcf.gz",
    output:
        fltvgz="03.Germline.freebayes/{sample}/{sample}.freebayes.flt.vcf.gz",
        mgz=   "03.Germline.freebayes/{sample}/{sample}.freebayes.flt.vcf.gz.anno/Merge.Anno.matrix.gz",
    log:
        out = snakedir+"/logs/C4.anno_freebayes/{sample}.o",
        err = snakedir+"/logs/C4.anno_freebayes/{sample}.e",
    threads:  16
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:50 ',
    shell:
        """
        module load {config[modules][vcflib]} {config[modules][samtools]}
        vcffilter \
          -f "QUAL > 20 & DP > 8 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
          {input.vgz} 2>{log.err} |bgzip > {output.fltvgz}
        tabix -f -p vcf {output.fltvgz} >> {log.out} 2>> {log.err}
        /data/yuk5/pipeline/vcf_annotation/vcf_annotation.v1.wgs.py \
          {output.fltvgz} \
          {output.fltvgz}.anno {threads} n >> {log.out} 2>> {log.err}
        """

rule Strelka:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    params:
        dir="03.Germline.Strelka/{sample}",
    output:
        vgz="03.Germline.Strelka/{sample}/results/variants/variants.vcf.gz",
    log:
        out = snakedir+"/logs/C5.strelka/{sample}.o",
        err = snakedir+"/logs/C5.strelka/{sample}.e",
    threads:  24
    resources:
        mem  = '48g',
        extra = ' --gres=lscratch:50 ',
    shell:
        """
        module load {config[modules][strelka]}
        configureStrelkaGermlineWorkflow.py \
          --bam {input.bam} \
          --referenceFasta {config[references][fasta]} \
          --runDir {params.dir} > {log.out} 2> {log.err}
        cd {params.dir}
        ./runWorkflow.py -m local -j {threads}
        rm -rf workspace
        """
rule anno_Strelka:
    input:
        vgz="03.Germline.Strelka/{sample}/results/variants/variants.vcf.gz",
    output:
        mgz="03.Germline.Strelka/{sample}/results/variants/variants.vcf.gz.anno/Merge.Anno.matrix.gz",
    log:
        out = snakedir+"/logs/C5.anno_strelka/{sample}.o",
        err = snakedir+"/logs/C5.anno_strelka/{sample}.e",
    threads:  16
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:50 ',
    shell:
        """
        /data/yuk5/pipeline/vcf_annotation/vcf_annotation.v1.wgs.py \
          {input.vgz} \
          {input.vgz}.anno {threads} n >> {log.out} 2>> {log.err}
        """

rule Gridss:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        vgz="04.SV.Gridss/{sample}/{sample}.gridss.vcf.gz",
    log:
        out = snakedir+"/logs/D01.Gridss/{sample}.o",
        err = snakedir+"/logs/D01.Gridss/{sample}.e",
    threads:  16
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:200 ',
    shell:
        """
        module load {config[modules][gridss]}
        cp {config[references][fasta]} /lscratch/${{SLURM_JOB_ID}}/ref.fa 
        gridss -r /lscratch/${{SLURM_JOB_ID}}/ref.fa \
          -t {threads} \
          -b {config[references][blacklist_encode]} \
          -o {output.vgz} \
          -a /lscratch/${{SLURM_JOB_ID}}/assembly.bam \
          -w /lscratch/${{SLURM_JOB_ID}} \
          {input.bam} > {log.out} 2> {log.err}
        """

rule Delly:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        sv=  "04.SV.Delly/{sample}/{sample}.delly.sv.bcf",
        svgz="04.SV.Delly/{sample}/{sample}.delly.sv.vcf.gz",
        cnv= "05.CNV.Delly/{sample}/{sample}.delly.cnv.bcf",
        cvgz="05.CNV.Delly/{sample}/{sample}.delly.cnv.vcf.gz",
    log:
        out = snakedir+"/logs/D02.Delly/{sample}.o",
        err = snakedir+"/logs/D02.Delly/{sample}.e",
    threads:  2
    resources:
        mem  = '8g',
        extra = ' --gres=lscratch:20 ',
    shell:
        """
        module load {config[modules][delly]} {config[modules][bcftools]} {config[modules][samtools]}
        delly call \
          -g {config[references][fasta]} \
          -x {config[references][delly]}/human.hg38.excl.tsv \
          -o {output.sv} \
          {input.bam} > {log.out} 2> {log.err}
        delly cnv \
          -g {config[references][fasta]} \
          -m {config[references][delly]}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz \
          -l {output.sv} \
          -o {output.cnv} \
          {input.bam} >> {log.out} 2>> {log.err}
        bcftools view {output.sv}|bgzip > {output.svgz} 2>>{log.err}
        bcftools view {output.cnv}|bgzip > {output.cvgz} 2>>{log.err}
        tabix -p vcf {output.svgz}
        tabix -p vcf {output.cvgz}
        """
        
rule Canvas:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
        vgz="03.Germline.Strelka/{sample}/results/variants/variants.vcf.gz",
    output:
        vgz="05.CNV.Canvas/SingleSample/{sample}/CNV.vcf.gz",
    params:
        outdir=directory("05.CNV.Canvas/SingleSample/"),
    log:
        out = snakedir+"/logs/E01.Canvas/{sample}.o",
        err = snakedir+"/logs/E01.Canvas/{sample}.e",
    threads:  32
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:50 ',
    shell:
        """
        module load {config[modules][canvas]} {config[modules][samtools]}
        cd {params.outdir}
        Canvas.dll SmallPedigree-WGS \
          -b {config[workdir]}/{input.bam} \
          -r /fdb/Canvas/GRCh38/kmer.fa \
          -g /fdb/Canvas/GRCh38/Sequence/WholeGenomeFasta \
          -f /fdb/Canvas/GRCh38/filter13.bed \
          --sample-b-allele-vcf {config[workdir]}/{input.vgz} \
          -o {wildcards.sample} >> {log.out} 2>> {log.err}
        ls {wildcards.sample} >> {log.out} 2>> {log.err}
        tabix -p vcf {wildcards.sample}/cnv.vcf.gz >> {log.out} 2>> {log.err}
        rm -rf {wildcards.sample}/TempCNV >> {log.out} 2>> {log.err}
        """
        
rule Windowdepth:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        gz="02.Alignment/WindowDepth5K/{sample}/{sample}.bin.5K.depth.gz",
    log:
        out = snakedir+"/logs/E02.WindowDepth/{sample}.o",
        err = snakedir+"/logs/E02.WindowDepth/{sample}.e",
    threads:  8
    resources:
        mem  = '8g',
        extra = ' --gres=lscratch:50 ',
    shell:
        """
        module load {config[modules][sambamba]}
        sambamba depth \
          window \
          -w 5000 \
          -F "mapping_quality>10" \
          -t {threads} \
          -a {input.bam}|gzip -c - > {output.gz} 2>> {log.err}
        """

rule Scramble:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
        vgz="03.Germline.Strelka/{sample}/results/variants/variants.vcf.gz",
    output:
        dir=directory("06.MEI.Scramble/{sample}/"),
        cluster="06.MEI.Scramble/{sample}/{sample}.cluster.txt",
        mei    ="06.MEI.Scramble/{sample}/{sample}.scramble_MEIs.txt",
    params:
        mei = "06.MEI.Scramble/{sample}/{sample}.scramble"
    log:
        out = snakedir+"/logs/F01.Scramble/{sample}.o",
        err = snakedir+"/logs/F01.Scramble/{sample}.e",
    threads:  16
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:50 ',
    shell:
        """
        module load {config[modules][scramble]} {config[modules][singularity]}
        singularity exec -e \
            -B /gpfs,/gs9,/data,/home \
            {config[simg][scramble]} \
            cluster_identifier \
            {config[workdir]}/{input.bam} \
            > {output.cluster} 2>> {log.err}
        echo "Step 2" > {log.out} 2>>{log.err}
        singularity exec -e \
            -B /gpfs,/gs9,/data,/home \
            {config[simg][scramble]} \
        Rscript \
            --vanilla /app/cluster_analysis/bin/SCRAMble.R \
            --out-name {config[workdir]}/{params.mei} \
            --cluster-file {config[workdir]}/{output.cluster} \
            --install-dir /app/cluster_analysis/bin \
            --mei-refs /app/cluster_analysis/resources/MEI_consensus_seqs.fa \
            --ref {config[references][scramblefa]} \
            --eval-dels \
            --eval-meis \
            --no-vcf >> {log.out} 2>> {log.err}
        """


