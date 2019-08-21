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

// check if user set the params tumour_id & normal_id, if not get it from the VCF filename
if (params.tumour_id == false && params.normal_id == false) {
    tumour_regex= "(.+)_vs"
    normal_regex= "vs_(.+).vcf*"
    tumour_id = (params.somatic_vcf =~ tumour_regex)[0][1]
    normal_id = (params.somatic_vcf =~ normal_regex)[0][1]
} else {
    tumour_id = params.tumour_id
    normal_id = params.normal_id
}
ids = Channel.value(["${tumour_id}","${normal_id}"])

process Vcf2maf {
    tag "$vcf"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(tumour_id), val(normal_id) from ids
    each file(vcf) from somatic_vcf_channel
    each file(fasta) from fasta_channel

    output:
    file("${tumour_id}_vs_${normal_id}.maf") into vcf_variant_eval

    script:
    """
    perl /opt/vcf2maf/vcf2maf.pl \
    --input-vcf $vcf \
    --output-maf ${tumour_id}_vs_${normal_id}.maf  \
    --tumor-id $tumour_id \
    --normal-id $normal_id \
    --ref-fasta /vepdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --ncbi-build ${params.genome} \
    --filter-vcf /vepdata/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
    --vep-path /opt/variant_effect_predictor_${params.vep_cache_version}/ensembl-tools-release-${params.vep_cache_version}/scripts/variant_effect_predictor \
    --vep-data /vepdata/ \
    --vep-forks ${params.vep_forks} \
    --buffer-size ${params.buffer_size} \
    --species ${params.species} \
    --cache-version ${params.vep_cache_version}
    """
}