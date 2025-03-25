process PREPARE_INTERACTIONS {
    tag   "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), val(fasta), val(a3m), path(files)
    output:
    tuple val(meta), path ("*.fasta"), emit: fasta
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    #!/usr/bin/env python3
    import os, sys
    import string
    #from itertools import product
    
    #single_letters = list(string.ascii_uppercase)
    #two_letter_combinations = [''.join(p) for p in product(string.ascii_uppercase, repeat=2)]
    #all_combinations = single_letters + two_letter_combinations
    all_combinations = list(string.ascii_uppercase) + list(string.ascii_lowercase) + [str(x) for x in range(0, 10)]
    fasta_files = ["${fasta.join('", "')}"]
    a3m_files = ["${a3m.join('", "')}"]
    if len(fasta_files) != len(a3m_files):
        raise ValueError("FASTA and A3M file lists must be of the same length and order.")

    output_file = "${meta.id}_boltz_interaction_input.fasta"
    counter = 0
    with open(output_file, "w") as outfile:
        for fasta, a3m in zip(fasta_files, a3m_files):
            with open(fasta, "r") as f:
                lines = f.readlines()

            if not lines:
                continue  # Skip empty FASTA files

            #header = lines[0].strip()
            header = f">{all_combinations[counter]}"
            counter += 1
            body = lines[1:]
            if header[-1] == "|":
                new_header = f"{header}protein|{os.path.basename(a3m)}\\n"
            else:
                new_header = f"{header}|protein|{os.path.basename(a3m)}\\n"

            outfile.write(new_header)
            outfile.writelines(body)

    with open ("versions.yml", "w") as version_file:
	    version_file.write("\\"${task.process}\\":\\n    python: {}\\n".format(sys.version.split()[0].strip()))
    """

    stub:
    """
    touch "${meta.id}_boltz_interaction_input.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        generate_comparison_report.py: \$(python3 --version)
    END_VERSIONS
    """
}
