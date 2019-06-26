metadata = open(config['metadata'])

# first dict hold SAMPLE (key) FASTQ files (values)
SAMPLE_FASTQ = {}
# second dict holds SAMPLE (key) and PairedEnd/SingleEnd status (value)
SAMPLE_END = {}
for line in metadata:
	line = line.split()
	if line[0] == 'file_name':
		continue
	fastq = line[0]
	sample = line[1]
	if sample in SAMPLE_FASTQ:
		old_vals = SAMPLE_FASTQ[sample]
		new_vals = old_vals + [fastq]
		SAMPLE_FASTQ[sample] = new_vals
	else:
		SAMPLE_FASTQ[sample] = [fastq]

# also make sure values (fastq) in SAMPLE_FASTQ are name sorted
for key in SAMPLE_FASTQ.keys():
	# if R2 in fastq, then it is paired end
	if [x for x in SAMPLE_FASTQ[key] if re.search('R2', x)]:
		SAMPLE_END[key] = 'PE'
		SAMPLE_FASTQ[key].sort()
	else:
		SAMPLE_END[key] = 'SE'
		SAMPLE_FASTQ[key].sort()

rule all:
	input:
		aligned = expand('bams/{sample}.bam', sample = SAMPLE_FASTQ.keys()),
		multiqc = 'fastqc/multiqc/multiqc_report.html'


rule align:
	input:
		lambda wildcards: expand('fastq/{sample}/{fastq_files}', \
				sample = wildcards.sample, \
				fastq_files = SAMPLE_FASTQ[wildcards.sample])
	conda:
		'envs/mapping.yaml'
	output:
		'bams/{sample}.bam'
	threads: 
		12
	shell:	
		"""
		ls {input} | grep R1_001 | xargs zcat | gzip > {output}R1
		if [ `ls {input} | grep R2_001 | wc -l` -eq 0 ]; then
			bwa mem -t {threads} -B 4 -O 6 -E 1 -M {config[fasta]} \
				{output}R1 | \
				samtools view -1 - > {output}
			rm {output}R1
		else
			ls {input} | grep R2_001 | xargs zcat | gzip > {output}R2
			bwa mem -t {threads} -B 4 -O 6 -E 1 -M {config[fasta]} \
				{output}R1 \
				{output}R2 | \
				samtools view -1 - > {output}
			rm {output}R1
			rm {output}R2
		fi
		""" 
	
rule fastqc:
	input:
		expand('bams/{sample}.bam', sample = SAMPLE_FASTQ.keys())
	conda:
		'envs/fastqc.yaml'
	output:
		directory('fastqc/{sample}')
	shell:
		"""
		mkdir -p {output}
			fastqc -t {threads} -o {output} {input}
		"""

rule multiqc:
	input:
		expand('fastqc/{sample}', sample = SAMPLE_FASTQ.keys())
	conda:
		'envs/fastqc.yaml'
	output:
		'fastqc/multiqc/multiqc_report.html'
	shell:
		"""
		multiqc fastqc/ -o fastqc/multiqc
		"""
