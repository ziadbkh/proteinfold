/*
 * Run Boltz
 */
process RUN_BOLTZ {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/nf-core/proteinfold_boltz:dev"
    
    input:
    tuple val(meta), path(fasta)
    path (files)
    path ('boltz1_conf.ckpt')
    path ('ccd.pkl')
    
    output:
    tuple val(meta), path ("boltz_results_*/processed/msa/*.npz"), emit: msa
    tuple val(meta), path ("boltz_results_*/processed/structures/*.npz"), emit: structures
    tuple val(meta), path ("boltz_results_*/predictions/*/confidence*.json"), emit: confidence
    tuple val(meta), path ("*"), emit: plddt
    tuple val(meta), path ("boltz_results_*/predictions/*/*.pdb"), emit: pdb

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''

    """
    boltz predict --output_format pdb ${args} "${fasta}" --cache ./
    """
    stub:
    """
    mkdir -p boltz_results_S1/processed/msa/
    mkdir -p boltz_results_S1/processed/structures/
    mkdir -p boltz_results_S1/predictions/S1/
    
    touch boltz_results_S1/processed/msa/S1.npz
    touch boltz_results_S1/processed/structures/S1.npz
    touch boltz_results_S1/predictions/S1/confidence_S1.json
    touch boltz_results_S1/predictions/S1/S1.pdb
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
