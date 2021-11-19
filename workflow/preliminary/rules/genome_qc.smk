rule checkm_genome:
	input:
		"input_folder/genomes/{id}.fna"	
	output:
		tmpout = temp("results/.tmp/checkm/{id}/{id}.fasta"),
		out1 = "results/checkm/Bacteria_{id}/storage/bin_stats_ext.tsv"
	conda:
		"../envs/checkm.yaml"
	params:
		tmpdir = "tmp/{id}",
		ext = "fasta",
		rank = "domain",
		taxon = "Bacteria",
	threads: 2
	log:
		"results/log/checkm/{id}.log"
	message:
		"Running checkm on sample {wildcards.id}"
	shell:
		"""
		cp {input} {output.tmpout}
		
		tmpdir=$(dirname {output.tmpout})
		outdir=$(dirname $(dirname {output.out1}))

		checkm taxonomy_wf -t {threads} -x {params.ext} \
		 {params.rank} \
		 {params.taxon} \
		 $tmpdir \
		 $outdir &> {log}
		"""

rule compile_checkm_summary:
	input:
		metadata = config["metadata"],
		checkms = expand(rules.checkm_genome.output.out1, id = ids)
	output:
		"results/summary_checkm.tsv"
	shell:
		"""
		python script/extract_checkm.py {input.metadata} > {output}
		"""
