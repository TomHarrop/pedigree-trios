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
        'output/030_trios/mendelian.pdf'

# check mendelian patterns (re run after filtering)
# filter for biallelec, SNPs only, check maf etc.
rule plot_trios:
    input:
        trios = 'output/030_trios/mendelian.txt',
    output:
        plot = 'output/030_trios/mendelian.pdf'
    log:
        'output/logs/plot_trios.log'
    singularity:
        r
    script:
        'src/plot_trios.R'

rule test_trios:
    input:
        vcf = 'output/020_filtered/filtered.vcf.gz',
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


# filter
rule filter_vcf:
    input:
        vcf = 'output/010_genotypes/calls.vcf.gz',
    output:
        temp('output/020_filtered/filtered.vcf')
    params:
        min_maf = 0.05,
        f_missing = 0.2,
        min_qual = 30,
        max_dp = 12
    log:
        'output/logs/filter_vcf.log'
    singularity:
        samtools
    shell:
        'bcftools view '
        '-v snps '          # SNPs only
        '-m2 -M2 '          # biallelic sites only
        '--min-af {params.min_maf}:nonmajor '
        '--exclude '
        '"F_MISSING>{params.f_missing} '
        '|| FMT/DP>{params.max_dp} '
        '|| QUAL<{params.min_qual}" '
        '{input.vcf} '
        '> {output} '
        '2> {log}'


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

# generics
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


rule generic_index_vcf:
    input:
        Path('{folder}', '{file}.vcf')
    wildcard_constraints:
        folder = 'output/(?!010).*'
    output:
        gz = Path('{folder}', '{file}.vcf.gz'),
        tbi = Path('{folder}', '{file}.vcf.gz.tbi')
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output.gz} '
        '; '
        'tabix -p vcf {output.gz}'
