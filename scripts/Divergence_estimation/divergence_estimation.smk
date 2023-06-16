SAMPLES=[line.strip() for line in open("samples.fofn")]

rule all:
        input:
                expand("{sample}/aln.sam", sample = SAMPLES),

rule BAC_mapping:
        input:
              	"{sample}",
        params:
               	sge_opts="-l mfree=5G -l h_rt=5:24:00:00 -pe serial 12",
        output:
                "{sample}/aln.sam"
        shell:"""
module load minimap2/2.16
module load samtools/1.9
minimap2 -s 10000 -t 12 --eqx -ax asm20 --secondary=no -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
 {input}/ref.fa {input}/qry.fa > {output} 
python samIdentity.py --header --bed {output} >{output}.bed
"""
