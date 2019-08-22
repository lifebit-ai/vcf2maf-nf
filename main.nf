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
fasta = file(params.fasta)
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
fai = file(params.fai)
params.filter_vcf = params.genome ? params.genomes[ params.genome ].filter_vcf ?: false : false
filter_vcf = file(params.filter_vcf)
params.vep_path = params.genome ? params.genomes[ params.genome ].vep_path ?: false : false
vep_path = file(params.vep_path)

// check if user set the params tumour_id & normal_id, if not get it from the VCF filename
if (params.tumour_id == false && params.normal_id == false) {
    vcf = params.somatic_vcf.split("/").last()
    tumour_regex= "(.+)_vs"
    normal_regex= "vs_(.+).vcf*"
    tumour_id = (vcf =~ tumour_regex)[0][1]
    normal_id = (vcf =~ normal_regex)[0][1]
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
    each file(fasta_opt) from fasta
    each file(fai_opt) from fai
    each file(filter_vcf_opt) from filter_vcf
    each file(vep_path_opt) from vep_path

    output:
    file("${tumour_id}_vs_${normal_id}.maf") into vcf_variant_eval

    script:
    def fasta = fasta.name != 'NO_FILE' ? "$fasta_opt" : 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
    def fai = fai.name != 'NO_FILE1' ? "$fai_opt" : ''
    def filter_vcf = filter_vcf.name != 'NO_FILE2' ? "$filter_vcf_opt" : 'ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz'
    def vep_path = vep_path.name != 'NO_FILE3' ? "$vep_path_opt" : "/opt/variant_effect_predictor_${params.vep_cache_version}/ensembl-tools-release-${params.vep_cache_version}/scripts/variant_effect_predictor"
    def command = params.genome != 'GRCh37' ? "mv $fasta /vepdata/ && mv $fai /vepdata/ && mv $filter_vcf /vepdata" : ''
    """
    $command
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