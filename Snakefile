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


sample_info = 'data/trios_gt.csv'
ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'


########
# MAIN #
########

# read sample data
sample_df = pandas.read_csv(sample_info,
                            index_col='sample')

all_samples = sorted(set(sample_df.index))

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

