#!/usr/bin/env  bash
#
# The remainder are options passed to the script
#$ -S /bin/bash -V
#$ -P eichlerlab
#$ -l mfree=3G
#$ -l h_rt=25:06:00:00
#$ -pe serial 12
#$ -cwd
#$ -q eichler-short.q


module load bwakit/0.7.15 miniconda/4.8.3
module load samtools/1.9
module load htslib/1.9-20
/net/eichler/vol27/projects/assemblies/nobackups/macaque/NDD_analysis/per_gene_anlaysis/software/ensembl-vep-release-99/vep -i ../asmsv_svpop_panPri.vcf.gz -o asmsv_svpop_panPri.vcf.anno --cache --pick --pick_order biotype,ccds,rank,canonical  --force_overwrite --fork 12 --dir_cache /net/eichler/vol26/home/smurali/.vep --vcf 
