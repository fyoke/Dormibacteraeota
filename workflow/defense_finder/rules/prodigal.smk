rule prodigal:
    input:
        "../../input_folder/genomes/{sample}.fasta"
    output:
        genes = "results/prodigal/{sample}/{genetic_code}/{sample}_gc{genetic_code}.gff",
        aa = "results/prodigal/{sample}/{genetic_code}/{sample}_gc{genetic_code}.faa"
    log:
        "results/log/prodigal/{sample}_gc{genetic_code}.log"
    conda:
        "../envs/prodigal.yaml"
    threads: 8
    shell:
        """
        prodigal -i {input} -o {output.genes} -g {wildcards.genetic_code} -a {output.aa} -f gff 2> {log}
        """