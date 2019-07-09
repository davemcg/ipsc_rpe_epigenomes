metadata = open(config['metadata'])

# first dict hold SAMPLE (key) FASTQ files (values)
SAMPLE_FASTQ = {}
# second dict holds SAMPLE (key) and PairedEnd/SingleEnd status (value)
SAMPLE_END = {}
# which samples are paired ended
PE = []
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
	if 'R2' in fastq:
		PE.append(line[1])	

# also make sure values (fastq) in SAMPLE_FASTQ are name sorted
for key in SAMPLE_FASTQ.keys():
	# if R2 in fastq, then it is paired end
	if [x for x in SAMPLE_FASTQ[key] if re.search('R2', x)]:
		SAMPLE_END[key] = 'PE'
		SAMPLE_FASTQ[key].sort()
	else:
		SAMPLE_END[key] = 'SE'
		SAMPLE_FASTQ[key].sort()

# build experiment <-> control dict
exp_dict = {}
experiment_ctrl = open(config['input_exp_key'])
for line in experiment_ctrl:
	line = line.split()
	if line[2] == 'input':
		continue
	exp_dict[line[0]] = line[2]

# build antibody <-> replicate num dict
ab_replicate = {}
for exp in list(exp_dict.keys()):
	line = exp.split('_')
	if line[0] not in ab_replicate:
		ab_replicate[line[0]] = list(line[1])
	else:
		existing = ab_replicate[line[0]]
		existing.append(line[1])
		ab_replicate[line[0]] = existing

wildcard_constraints:
	sample = '|'.join(list(SAMPLE_FASTQ.keys()))

rule all:
	input:
		'/data/mcgaugheyd/datashare/ipsc_rpe_epigenomes/hg19/trackDb.txt',
		expand('bigWig/{sample}.bw', sample = SAMPLE_FASTQ.keys()),
		expand('crams/{sample}.cram', sample = SAMPLE_FASTQ.keys()),
		'fastqc/multiqc/multiqc_report.html',
		expand('macs2/{sample}_peaks.xls', sample = list(exp_dict.keys())),
		expand('macs2/{type}_common_peaks.bb', type = ab_replicate.keys())

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

rule bam_to_cram:
	input:
		'bams/{sample}.bam'
	conda:
		'envs/mapping.yaml'
	output:
		cram = 'crams/{sample}.cram',
		crai = 'crams/{sample}.crai'
	threads:
		8
	shell:
		"""
		samtools sort -O bam -l 0 --threads {threads} -T /lscratch/$SLURM_JOB_ID {input} | \
			samtools view -T {config[fasta_with_index]} --threads {threads} -C -o {output.cram} -
		samtools index {output.cram} {output.crai}
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

rule filter_bam:
	input:
		'bams/{sample}.bam'
	output:
		metrics = 'metrics/{sample}.picard.metrics',
		bam = 'bams/{sample}.q5.rmdup.bam'
	conda:
		'envs/mapping.yaml'
	shell:
		"""
		samtools view -O bam -F 4 -q 5 -b {input} | \
		samtools sort - -O bam -o {output.bam}TEMP
		samtools index {output.bam}TEMP
		java -Xmx4g -XX:ParallelGCThreads={threads} -jar picard.jar \
			MarkDuplicates \
			INPUT={output.bam}TEMP \
			OUTPUT={output.bam} \
			REMOVE_DUPLICATES=true \
			METRICS_FILE={output.metrics}
		samtools index {output.bam}
		rm {output.bam}TEMP*
		"""

# to 15e6 reads
rule downsample:
	input:
		'bams/{sample}.q5.rmdup.bam'
	output:
		'bams/{sample}.ds.q5.rmdup.bam'
	conda:
		'envs/mapping.yaml'
	shell:		
		"""
		frac=$( samtools idxstats {input} | cut -f3 | \
			awk 'BEGIN {{total=0}} {{total += $1}} \
				END {{frac=15000000/total;\
			if (frac > 1) {{print "1.0"}} else {{print frac}}}}' )
		samtools view -s $frac -b {input} > {output}
		samtools index {output}
		"""

rule macs2_peaks:
	input:
		'bams/{sample}.ds.q5.rmdup.bam'
	output:
		'macs2/{sample}_peaks.narrowPeak'
	run:
		ctrl = exp_dict[wildcards.sample]
		ctrl_bam = 'bams/' + ctrl + '.ds.q5.rmdup.bam'
		if wildcards.sample in PE:
			type = 'BAMPE'
		else:
			type = 'BAM'
		shell("module load macs/2.1.2; macs2 callpeak -f " + type + " -g \"hs\" -t {input} \
				-c " + ctrl_bam + " -q 0.01 -n {wildcards.sample} \
				--outdir macs2 --scale-to small")

rule union_peaks:
	input:
		lambda wildcards: expand('macs2/{{type}}_{replicate}_peaks.narrowPeak', replicate = list(ab_replicate[wildcards.type])) 
	output:
		'macs2/{type}_union_peaks.narrowPeak'
	shell:
		"""
		cat {input} | sort -k1,1 -k2,2n > {output}
		"""
	
rule merge_peaks:
	input:
		'macs2/{type}_union_peaks.narrowPeak'
	output:
		'macs2/{type}_merge_peaks.narrowPeak'
	conda:
		'envs/bedtools.yaml'
	shell:
		"""
		bedtools merge -d 50 -i {input} -c 5,7,10 -o mean > {output}
		"""
		
rule common_peaks:
	input:
		union = 'macs2/{type}_union_peaks.narrowPeak',
		merge = 'macs2/{type}_merge_peaks.narrowPeak',
		blacklist = 'ENCFF001TDO.bed.gz'
	output:
		'macs2/{type}_common_peaks.narrowPeak'
	run:
		shell("module load bedtools; \
				bedtools intersect -a {input.merge} -b {input.union} -c -f 0.4 | \
				bedtools intersect -v -a - -b {input.blacklist} | \
				awk '$7>1 {{print $0}}' > {output}T")
		tsv = open(output[0] + 'T')
		out = open(output[0], 'w')
		if wildcards.type == 'H3K27ac':
			color = '100,0,0'
		elif wildcards.type == 'H3K27me3':
			color = '200,0,0'
		elif wildcards.type == 'MITF':
			color = '0,100,0'
		elif wildcards.type == 'PAX6':
			color = '0,0,100'
		elif wildcards.type == 'SOX10':
			color = '100,100,100'
		else:
			color = '180,180,180'
		for line in tsv:
			line = line.split()
			line[3] = str(min(999, round(float(line[3])))) # round for bigBed, no more than 999
			if int(line[3]) < 20:
				continue
				# don't keep peaks with average score under q 0.01
				# -10 * log10(qvalue)
			new_line = '\t'.join(line) + '\t1\t.\t' + line[1] + '\t' + line[2] + '\t' + color + '\n'	
			out.write(new_line)
		tsv.close()
		out.close()
		shell("rm {output}T")
		
rule blacklist:
	output:
		'ENCFF001TDO.bed.gz'
	shell:
		'wget https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz'

rule bam_to_bigWig:
	input:
		'bams/{sample}.ds.q5.rmdup.bam'
	output:
		bedgraph = temp('bigWig/{sample}.bG'),
		bw = 'bigWig/{sample}.bw'
	threads:
		16
	conda:
		'envs/bedtools.yaml'
	shell:
		"""
		module load ucsc
		bamCoverage --bam {input} -o {output.bedgraph} \
			--numberOfProcessors {threads} \
			--binSize 10 \
			--normalizeUsing RPGC \
			--effectiveGenomeSize 2864785220 \
			--outFileFormat bedgraph
		sort -k1,1 -k2,2n {output.bedgraph} > {output.bedgraph}TEMP 
		bedGraphToBigWig {output.bedgraph}TEMP /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {output.bw}
		rm {output.bedgraph}TEMP
		"""

localrules: bed_to_bigBed
rule bed_to_bigBed:
	input:
		'macs2/{type}_common_peaks.narrowPeak'
	output:
		'macs2/{type}_common_peaks.bb'
	shell:
		"""
		module load ucsc
		cut -f1,2,3,4 {input} | sort -k1,1 -k2,2n > {input}TEMP
		bedToBigBed {input}TEMP /data/mcgaugheyd/genomes/hg19/hg19.chrom.sizes {output}
		rm {input}TEMP
		"""
localrules: trackDb_maker
rule trackDb_maker:
	input:
		bigBeds = expand('macs2/{type}_common_peaks.bb', type = ab_replicate.keys()),
		bigWigs = expand('bigWig/{sample}.bw', sample = SAMPLE_FASTQ.keys())
	output:
		trackDb = '/data/mcgaugheyd/datashare/ipsc_rpe_epigenomes/hg19/trackDb.txt',
		hub = '/data/mcgaugheyd/datashare/ipsc_rpe_epigenomes/hub.txt',
		genome = '/data/mcgaugheyd/datashare/ipsc_rpe_epigenomes/genomes.txt'
	run:
		types = list(ab_replicate.keys())
		samples = list(SAMPLE_FASTQ.keys())
		container_names =  ','.join(ab_replicate.keys())
		container_names = container_names + ',INP,IgG'
		big_wigs = ','.join([x + '.bw' for x in samples])
		parents = ','.join([x.split('_')[0] for x in samples])
		big_wig_colors = []
		for sample in [x.split('_')[0] for x in samples]:
			if sample == 'H3K27ac':
				big_wig_colors.append('100,0,0')
			elif sample == 'H3K27me3':
				big_wig_colors.append('200,0,0')
			elif sample == 'MITF':
				big_wig_colors.append('0,100,0')
			elif sample == 'PAX6':
				big_wig_colors.append('0,0,100')
			elif sample == 'SOX10':
				big_wig_colors.append('100,100,100')
			else:
				big_wig_colors.append('180,180,180')
		big_wig_colors = '\;'.join(big_wig_colors)
		bigBeds = ','.join([x + '_common_peaks.bb' for x in types])
		bigBed_colors = []
		for sample in [x.split('_')[0] for x in types]:
			if sample == 'H3K27ac':
				bigBed_colors.append('100,0,0')
			elif sample == 'H3K27me3':
				bigBed_colors.append('200,0,0')
			elif sample == 'MITF':
				bigBed_colors.append('0,100,0')
			elif sample == 'PAX6':
				bigBed_colors.append('0,0,100')
			elif sample == 'SOX10':
				bigBed_colors.append('100,100,100')
			else:
				bigBed_colors.append('180,180,180')
		bigBed_colors = '\;'.join(bigBed_colors)
		command = "module load python/3.6; \
			python3 /home/mcgaugheyd/git/ipsc_rpe_epigenomes/src/trackdb.py " + \
			container_names + ' ' + \
			big_wigs + ' ' + \
			parents + ' ' + \
			big_wig_colors + ' ' + \
			bigBeds + ' ' + \
			bigBed_colors + " > {output.trackDb}"
		shell(command)
		command = ''
		for file in input.bigBeds:
			command += "cp " + file + ' /data/mcgaugheyd/datashare/ipsc_rpe_epigenomes/hg19/; '
		for file in input.bigWigs:
			command += "cp " + file + ' /data/mcgaugheyd/datashare/ipsc_rpe_epigenomes/hg19/; '
		shell(command)
		# make hub.txt
		shell("printf \"hub ipsc_rpe_epigenomes\n\
shortLabel hufEpigenomes\n\
longLabel Hufnagel iPSC RPE Epigenomes\n\
genomesFile genomes.txt\n\
email mcgaugheyd@mail.nih.gov\n\
descriptionUrl ucscHub.html\n\
itemRgb on \" > {output.hub}")
		shell("printf \"genome hg19\ntrackDb hg19/trackDb.txt \" > {output.genome}")
