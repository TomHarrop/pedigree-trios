#!/usr/bin/env python3

import pandas

'''
cnv_map:

tail -n +2 data/trios_gt.csv \
| cut -d',' -f 1,5 --output-delimiter=' '

'''

###########
# GLOBALS #
###########

honeybee_genotype_pipeline = (
    'shub://TomHarrop/'
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.12')
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'
r = 'shub://TomHarrop/r-containers:r_4.0.0'


sample_info = 'data/trios_gt.csv'
ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'


########
# MAIN #
########

# read sample data
sample_df = pandas.read_csv(sample_info,
                            index_col='sample')

all_samples = sorted(set(sample_df.index))


rule target:
    input:
        'output/030_trios/mendelian.txt'

# check mendelian patterns (re run after filtering)
rule test_trios:
    input:
        vcf = 'output/010_genotypes/calls.vcf.gz',
        trios = 'output/030_trios/all_trios.txt',
        rules = 'output/030_trios/rules.txt'
    output:
        'output/030_trios/mendelian.txt'
    log:
        'output/logs/test_trios.log'
    singularity:
        samtools
    shell:
        'bcftools +mendelian '
        '{input.vcf} '
        '-T {input.trios} '
        '-R {input.rules} '
        '-c '
        '> {output} '
        '2> {log}'


rule generate_trios:
    input:
        sample_info = sample_info,
        fai = f'{ref}.fai'
    output:
        trios = 'output/030_trios/all_trios.txt',
        rules = 'output/030_trios/rules.txt'
    log:
        'output/logs/generate_trios.log'
    singularity:
        r
    script:
        'src/generate_trios.R'


# genotype
checkpoint genotype:
    input:
        csv = sample_info,
        ref = ref,
        cnv_map = 'output/tmp/cnv_map.txt',
        # reads = expand('output/000_tmp/reads/{sample}_R{r}.fq.gz',
        #                sample=all_samples,
        #                r=[1, 2])
    output:
        cutoffs = 'output/010_genotypes/040_stats/ldepth.mean_cutoffs.csv',
        vcf = 'output/010_genotypes/calls.vcf.gz',
        bam = 'output/010_genotypes/merged.bam',
        ref = 'output/010_genotypes/015_ref/ref.fasta',
        fai = 'output/010_genotypes/015_ref/ref.fasta.fai'
    params:
        wd = 'output/010_genotypes',
    log:
        'output/logs/genotype.log'
    threads:
        workflow.cores
    container:
        honeybee_genotype_pipeline
    shell:
        'honeybee_genotype_pipeline '
        '--ref {input.ref} '
        '--samples_csv {input.csv} '
        '--outdir {params.wd} '
        '--cnv_map {input.cnv_map} '
        '--threads {threads} '
        '--csd '
        '--restart_times 1 '
        '&>> {log}'

rule cnv_map:
    input:
        sample_info
    output:
        'output/tmp/cnv_map.txt'
    container:
        honeybee_genotype_pipeline
    shell:
        'tail -n +2 {input} '
        '| cut -d\',\' -f 1,5 --output-delimiter=\' \' '
        '> {output}'


rule index_fa:
    input:
        '{path}/{file}.{ext}'
    output:
        '{path}/{file}.{ext}.fai'
    wildcard_constraints:
        ext = 'fasta|fa|fna'
    singularity:
        samtools
    shell:
        'samtools faidx {input}'
