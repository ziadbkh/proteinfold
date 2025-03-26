/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { CREATE_SAMPLESHEET_YAML } from '../modules/local/create_samplesheet'
include { MMSEQS_COLABFOLDSEARCH } from '../modules/local/mmseqs_colabfoldsearch'
include { PREPARE_INTERACTIONS } from '../modules/local/prepare_interactions'
//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_proteinfold_pipeline'

//
// MODULE: Boltz
//
include { RUN_BOLTZ } from '../modules/local/run_boltz'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INTERACTIONS {
    
    take:
    ch_samplesheet_1  // channel: samplesheet read from --input
    samplesheet_2  // channel: samplesheet read from --input
    ch_versions     // channel: [ path(versions.yml) ]
    ch_boltz_ccd    // channel: [ path(boltz_ccd) ]
    ch_boltz_model  // channel: [ path(model) ]
    ch_colabfold_params // channel: [ path(colabfold_params) ]
    ch_colabfold_db // channel: [ path(colabfold_db) ]
    ch_uniref30     // channel: [ path(uniref30) ]

    main:
    ch_multiqc_files = Channel.empty()
    
    Channel.fromPath(samplesheet_2)
    .splitCsv(header: true)
    .map{row -> [["id": row.id], file(row.fasta, checkIfExists: true)]}
    .set{ch_samplesheet_2}
    
    ch_samplesheet_1.map{it[0]}.combine(ch_samplesheet_2.map{it[0]}).map{["id": [it[0]["id"], it[1]["id"]].min() + "-" + [it[0]["id"], it[1]["id"]].max()]}.unique().set{ch_unique_pairs}

    MMSEQS_COLABFOLDSEARCH (
                ch_samplesheet_1.mix(ch_samplesheet_2).unique(),
                ch_colabfold_params,
                ch_colabfold_db,
                ch_uniref30
            )
    ch_versions = ch_versions.mix(MMSEQS_COLABFOLDSEARCH.out.versions)
    ch_left_list = ch_samplesheet_1.join(MMSEQS_COLABFOLDSEARCH.out.a3m)
    ch_right_list = ch_samplesheet_2.join(MMSEQS_COLABFOLDSEARCH.out.a3m)

    ch_left_list
    .combine(ch_right_list)
    .map{[
            ["id": [it[0]["id"], it[3]["id"]].min() + "-" + [it[0]["id"], it[3]["id"]].max()], 
            [it[1].name, it[4].name], 
            [it[2].name, it[5].name], 
            [it[1], it[2], it[4], it[5]].unique{it.toUriString()}
    ]}.set{ch_pairs}

    //ch_pairs.view()
    //ch_pairs.main.view()
    ch_pairs
    .join(ch_unique_pairs)
    .set{ch_interaction_input}
    
    PREPARE_INTERACTIONS(
        ch_interaction_input
    )
    ch_interaction_input.view()

    PREPARE_INTERACTIONS.out.versions.view()
    
    PREPARE_INTERACTIONS.out.fasta.view()
    
    PREPARE_INTERACTIONS.out.fasta
    .join(ch_interaction_input)
    .set{ch_boltz_input}
    
/*
//ch_interaction_input.view()
    // RUN_BOLTZ 
    //PREPARE_INTERACTIONS.out.fasta
    //.join(ch_interaction_input)
    //.map { it[0]["model"] = "interactions"; it }
    //.set{ch_boltz_input}
    
    //ch_boltz_input.view()
    RUN_BOLTZ(
        ch_boltz_input.map{[it[0], it[1]]},
        ch_boltz_input.map{[it[0], it[2]]},
        ch_boltz_model,
        ch_boltz_ccd
    )

    
    emit:
    versions   = ch_versions
    msa        = RUN_BOLTZ.out.msa
    structures = RUN_BOLTZ.out.structures
    confidence = RUN_BOLTZ.out.confidence
    plddt      = RUN_BOLTZ.out.plddt
    pdb        = ch_pdb
    */
    emit:
    versions   = ch_versions
} 
