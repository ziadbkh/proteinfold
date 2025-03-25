/*
 * Run Boltz
 */
process RUN_BOLTZ {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/nf-core/proteinfold_boltz:dev"
    
    input:
    tuple val(meta), path(fasta)
    path ('boltz1_conf.ckpt')
    path ('ccd.pkl')
    
    output:
    path ("boltz_results_*/processed/msa/*.npz"), emit: msa
    path ("boltz_results_*/processed/structures/*.npz"), emit: structures
    path ("boltz_results_*/predictions/*/confidence*.json"), emit: confidence
    path ("*"), emit: plddt
    path ("boltz_results_*/predictions/*/*.pdb"), emit: pdb
    
    script:
    """
    boltz predict --output_format pdb "./${fasta.name}" --cache ./
    """
}
