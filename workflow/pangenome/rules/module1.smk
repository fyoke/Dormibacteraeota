rule XXX:
	input:
		""
	output:
		""
	log:
		""
	conda:
		"../envs/tool1.yaml"
	threads: 1
	params:
		a = "",
		b = ""
	shell:
		"""
		command -i {input} \
			-o {output} \
			-p1 {params.a} \
			-p2 {params.b} \
			--threads {threads} 2> {log}
		"""

rule YYY:
	input:
		rules.XXX.output
	output:
		""
	shell:
		"""

		"""
