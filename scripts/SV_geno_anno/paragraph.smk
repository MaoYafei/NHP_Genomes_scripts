SAMPLES=[line.strip() for line in open("samples.fofn")]

rule all:
        input:
                expand("{sample}_output.txt",sample = SAMPLES),

rule IDX:
        input:
                "{sample}"
        params:
                sge_opts="-l mfree=3G -l h_rt=1:20:00:00 -pe serial 4",
        output:
                "{sample}.json"
        shell:"""
/net/eichler/vol26/projects/primate_sv/nobackups/Tools/paragraph/bin/idxdepth \
-b {input} -r /net/eichler/vol26/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa --threads 4 \
-o {output} --altcontig 1
"""


rule Mainest:
        input:
              	"{sample}.json"
        params:
               	sge_opts="-l mfree=2G -l h_rt=20:00:00 -pe serial 1",
        output:
               	"{sample}_ind.txt"
        shell:"""
echo -E "id\tpath\tidxdepth" > {input}.header.txt

echo {input}.header.txt|cut -f 2 -d '/'|cut -f 1 -d '.' > {input}.1.txt
sed -i 's/\//_/g' {input}.1.txt

echo /net/eichler/vol26/projects/primate_sv/nobackups/NHP_Pop_bam/`echo {input}|sed 's/.json//g'` > {input}.2.txt

echo /net/eichler/vol26/projects/primate_sv/nobackups/NHP_Pop_bam/{input} >{input}.3.txt

paste {input}.1.txt {input}.2.txt {input}.3.txt >{input}.4.txt
cat {input}.header.txt {input}.4.txt>{output}
rm {input}.1.txt {input}.2.txt {input}.3.txt {input}.4.txt {input}.header.txt
"""

rule Paragraph:
        input:
                TXT="{sample}_ind.txt",
                VCF="/net/eichler/vol26/projects/primate_sv/nobackups/sv/vcf/merged/asmsv_svpop_panPri.vcf.gz"
        params:
               	sge_opts="-l mfree=10G -l h_rt=25:20:00:00 -pe serial 12",
        output:
               	"{sample}_output.txt"
        shell:"""
python3 /net/eichler/vol26/projects/primate_sv/nobackups/Tools/paragraph/bin/multigrmpy.py -M 700 \
-i {input.VCF} -m {input.TXT} -o `echo {input.TXT}|cut -f 1 -d '.'`_GT_paragraph -t 12 -r /net/eichler/vol26/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa 
echo 'done'> {output}
"""
