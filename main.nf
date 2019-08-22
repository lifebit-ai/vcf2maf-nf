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

// setup optional params (already in the container for GRCh37)
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "FASTA file not found: ${params.fasta}" }
           .into { fasta }
}
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "FASTA index file not found: ${params.fai}" }
           .into { fai }
}
params.filter_vcf = params.genome ? params.genomes[ params.genome ].filter_vcf ?: false : false
if (params.filter_vcf) {
    Channel.fromPath(params.filter_vcf)
           .ifEmpty { exit 1, "Filter VCF file not found: ${params.filter_vcf}" }
           .into { filter_vcf }
}
params.filter_vcf = params.genome ? params.genomes[ params.genome ].filter_vcf ?: false : false
if (params.filter_vcf) {
    Channel.fromPath(params.filter_vcf)
           .ifEmpty { exit 1, "Filter VCF file not found: ${params.filter_vcf}" }
           .into { filter_vcf }
}
params.vep_path = params.genome ? params.genomes[ params.genome ].vep_path ?: false : false
if (params.vep_path) {
    Channel.fromPath(params.vep_path)
           .ifEmpty { exit 1, "VEP path folder not found: ${params.vep_path}" }
           .into { vep_path }
}

// check if user set the params tumour_id & normal_id, if not get it from the VCF filename
if (params.tumour_id == false && params.normal_id == false) {
    vcf = params.somatic_vcf.split("/").last()
    tumour_regex= "(.+)_vs"
    normal_regex= "vs_(.+).vcf*"
    try {
        tumour_id = (vcf =~ tumour_regex)[0][1]
        normal_id = (vcf =~ normal_regex)[0][1]
    } catch(Exception ex) {
        exit 1, "Please specify either --tumour_id & --normal_id OR rename your VCF to `tumourID_vs_normalID.vcf`"
    }
} else {
    tumour_id = params.tumour_id
    normal_id = params.normal_id
}
ids = Channel.value(["${tumour_id}","${normal_id}"])

if (params.genome == "GRCh37") {
    process Vcf2maf_GRCh37 {
    tag "$vcf"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(tumour_id), val(normal_id) from ids
    each file(vcf) from somatic_vcf_channel

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
} else {
    process Vcf2maf {
    tag "$vcf"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(tumour_id), val(normal_id) from ids
    each file(vcf) from somatic_vcf_channel
    each file(fasta) from fasta
    each file(fai) from fai
    each file(filter_vcf) from filter_vcf
    each file(vep_path) from vep_path

    output:
    file("${tumour_id}_vs_${normal_id}.maf") into vcf_variant_eval

    script:
    """
    mv $fasta /vepdata/ && mv $fai /vepdata/ && mv $filter_vcf /vepdata
    perl /opt/vcf2maf/vcf2maf.pl \
    --input-vcf $vcf \
    --output-maf ${tumour_id}_vs_${normal_id}.maf  \
    --tumor-id $tumour_id \
    --normal-id $normal_id \
    --ref-fasta /vepdata/${fasta} \
    --ncbi-build ${params.genome} \
    --filter-vcf /vepdata/${filter_vcf} \
    --vep-path $vep_path \
    --vep-data /vepdata/ \
    --vep-forks ${params.vep_forks} \
    --buffer-size ${params.buffer_size} \
    --species ${params.species} \
    --cache-version ${params.vep_cache_version}
    """
    }
}