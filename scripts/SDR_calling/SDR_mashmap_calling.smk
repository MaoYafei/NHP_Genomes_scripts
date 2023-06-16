import os
import numpy as np
import pandas as pd
import json
import random
import tempfile

SAMPLES=[line.strip() for line in open("samples.fofn")]

localrules:all

rule all:
	input:
		expand("{sample}.merge.bed", sample = SAMPLES),

rule mashmap_qry:
	input:
		REF="/net/eichler/vol26/projects/chm13_t2t/nobackups/assemblies/chm13.draft_v1.0.fasta",
		QRY="{sample}",
	output:
		"{sample}_qry.out.bed",
	params:
		sge_opts="-l mfree=10G -l h_rt=5:24:00:00 -pe serial 12",
	shell:"""
		module load minimap2 bedtools samtools/1.10 gcc/8.2.0 
		mashmap -r {input.REF} -q {input.QRY} -s 10000 -f one-to-one --pi 85 -t 12 -o {input.QRY}_qry.out
		cut -f 6,8,9 -d ' ' {input.QRY}_qry.out|tr " " "\\t" |bedtools subtract -b - -a chm13.draft_v1.0.fasta.bed \
			|awk '$3-$2>10000'|bedtools merge -i - >{input.QRY}_qry.out.bed
		awk '{{print $3-$2+1}}' {output}|~/bin/sum_awk.sh 1 >{output}.length
"""

rule mashmap_ref:
	input:
		REF="/net/eichler/vol26/projects/chm13_t2t/nobackups/assemblies/chm13.draft_v1.0.fasta",
		QRY="{sample}",
	output:
		"{sample}_ref.out.bed"
	params:
		sge_opts="-l mfree=10G -l h_rt=5:24:00:00 -pe serial 12",
	shell:"""
		module load minimap2 bedtools samtools/1.10 gcc/8.2.0
		mashmap -r {input.QRY} -q {input.REF} -s 10000 -f one-to-one --pi 85 -t 12 -o {input.QRY}_ref.out
		cut -f 1,3,4 -d ' ' {input.QRY}_ref.out|tr " " "\\t" |bedtools subtract -b - -a chm13.draft_v1.0.fasta.bed \
			|awk '$3-$2>10000'|bedtools merge -i - >{input.QRY}_ref.out.bed
		awk '{{print $3-$2+1}}' {output}|~/bin/sum_awk.sh 1 >{output}.length
"""

rule merge:
	input:
		"{sample}_qry.out.bed",
		"{sample}_ref.out.bed",
	output:
		"{sample}.merge.bed",
	params:
		sge_opts="-l mfree=2G -l h_rt=5:24:00:00 -pe serial 2",
	shell:"""
		module load minimap2 bedtools samtools/1.10
		cat {input}|sort -k1,1 -k2,2n|bedtools merge -i - >{output}
		awk '{{print $3-$2+1}}' {output}|~/bin/sum_awk.sh 1 >{output}.length
		bedtools intersect -wo -a {output} -b ../CHM13_reference/matrix/chm13.draft_v1.0.exon.bed \
			|grep protein_coding|cut -f 1-3,7,8|sort -k 1,1 -k2,2n -k8gr|uniq |sort -k 1,1 -k2,2n >{output}.anno
"""
