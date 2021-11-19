rule defensefinder:
    input:
        rules.prodigal.output.aa
    output:
        "results/defense_finder/{sample}/gc{genetic_code}/defense_finder_systems.tsv",
    log:
        "results/log/defense_finder/{sample}_gc{genetic_code}.log"
    conda:
        "../envs/defense.yaml"
    params:
        tmpdir = "results/.tmp/defensefinder/{sample}/gc{genetic_code}"
    threads: 8
    shell:
        """
        outdir=$(dirname {output})
        mkdir -p {params.tmpdir}
        defense-finder run -o {params.tmpdir} -w {threads} {input} 2> {log}
        mv $(find {params.tmpdir} -type f) $outdir
        rm -rf {params.tmpdir}
        """