import os
import numpy as np
import pandas as pd
import json
import random
import tempfile

SAMPLES=[line.strip() for line in open("samples.fofn")]
SAMPLES_REF=[line.strip() for line in open("samples_ref.fofn")]

localrules:all

rule all:
	input:
		expand("{sample}/{sample_ref}_{sample}_h1.png", sample = SAMPLES, sample_ref = SAMPLES_REF),
		expand("{sample}/{sample_ref}_{sample}_h2.png", sample = SAMPLES, sample_ref = SAMPLES_REF),

rule minimap2_qry_h1:
	input:
		REF="REF_ID/{sample_ref}",
		QRY="{sample}",
	output:
		"{sample}/{sample_ref}_{sample}_h1.png",
	params:
		sge_opts="-l mfree=10G -l h_rt=5:24:00:00 -pe serial 4",
	shell:"""
		module load minimap2 bedtools samtools/1.10 gcc/8.2.0 R/3.6.1
		minimap2 -K 8G -t 4 -c --secondary=no ../../{input.QRY}/h1.fa  {input.REF}|awk '$12>30 && $11>2000'|cut -f 1-12 \
			>{input.QRY}/`echo {input.REF}|cut -f 2 -d '/'`_h1.paf
		cd `echo {input.QRY}`
		Rscript ../pafCoordsDotPlotyly.R -p 5 -i `echo {input.REF}|cut -f 2 -d '/'`_h1.paf -o \
			`echo {input.REF}|cut -f 2 -d '/'`_{input.QRY}_h1 -q 1 -m 1
		rm `echo {input.REF}|cut -f 2 -d '/'`_{input.QRY}_h1.html
		cd ..
#		ps2pdf {input.QRY}/`echo {input.REF}|cut -f 2 -d '/'`_h1.out {output}
"""

rule minimap2_qry_h2:
	input:
		REF="REF_ID/{sample_ref}",
		QRY="{sample}",
	output:
		"{sample}/{sample_ref}_{sample}_h2.png",
	params:
		sge_opts="-l mfree=10G -l h_rt=5:24:00:00 -pe serial 4",
	shell:"""
		module load minimap2 bedtools samtools/1.10 gcc/8.2.0 R/3.6.1

		minimap2 -K 8G -t 4 -c --secondary=no ../../{input.QRY}/h2.fa  {input.REF}|awk '$12>30 && $11>2000'|cut -f 1-12 \
			>{input.QRY}/`echo {input.REF}|cut -f 2 -d '/'`_h2.paf
		cd `echo {input.QRY}`
		Rscript ../pafCoordsDotPlotyly.R -p 5 -i `echo {input.REF}|cut -f 2 -d '/'`_h2.paf -o \
                        `echo {input.REF}|cut -f 2 -d '/'`_{input.QRY}_h2 -q 1 -m 1
		rm `echo {input.REF}|cut -f 2 -d '/'`_{input.QRY}_h2.html
		cd ..
#		ps2pdf {input.QRY}/`echo {input.REF}|cut -f 2 -d '/'`_h2.ps {output}
"""
