#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/vcf2maf-nf
========================================================================================
 lifebit-ai/vcf2maf-nf Single process nextflow component for converting VCF to MAF
 #### Homepage / Documentation
 https://github.com/lifebit-ai/vcf2maf-nf
----------------------------------------------------------------------------------------
*/

Channel.fromPath(params.somatic_vcf, type: 'file')
       .set { somatic_vcf_channel }

Channel.fromPath(params.fasta, type: 'file')
       .set { fasta_channel }

vcf2maf_channel = somatic_vcf_channel.combine(fasta_channel)

process Vcf2maf {
    tag "$vcf"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set file(vcf), file(fasta) from vcf2maf_channel

    output:
    file maf into vcf_variant_eval

    script:
    """
    perl /opt/vcf2maf/vcf2maf.pl \
    --input-vcf $vcf \
    --output-maf maf  \
    --tumor-id H46126 \
    --normal-id H06530 \
    --ref-fasta /vepdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --ncbi-build ${params.genome} \
    --filter-vcf /vepdata/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
    --vep-path /opt/variant_effect_predictor_${params.vep_cache_version}/ensembl-tools-release-89/scripts/variant_effect_predictor \
    --vep-data /vepdata/ \
    --vep-forks ${params.vep_forks} \
    --buffer-size ${params.buffer_size} \
    --species ${params.species} \
    --cache-version ${params.vep_cache_version}
    """
}