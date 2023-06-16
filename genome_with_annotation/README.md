Due to size limitation of GitHub, all data are deposited in https://eichlerlab.gs.washington.edu/help/yafmao/NHP_Genomes_scripts/genome_with_annotation/{sample}_scaffolding
1. genome_with_annotation
This folder contains 16 genome assemblies and each assembly are named as a subfolder. each subfolder contains 3 files and 2 sufolders (detailed below).
	#ALL FILES CAN BE UPLOAD TO IGV/JBBROSWER FOR VISUALIZATION
	*_scaffolding_ragtag.scaffold.fasta: This FASTA file contains chromosome level NHP HiFi genome assembly (https://eichlerlab.gs.washington.edu/help/yafmao/NHP_Genomes_scripts/genome_with_annotation/{sample}_scaffolding/{sample}_scaffolding_ragtag.scaffold.fasta)
	*_scaffolding_ragtag.liftoff.gff: This gff file contains gene model prediction based on the "*_scaffolding_ragtag.scaffold.fasta" assembly (https://eichlerlab.gs.washington.edu/help/yafmao/NHP_Genomes_scripts/genome_with_annotation/{sample}_scaffolding/{sample}_scaffolding_ragtag.liftoff.gff) 
	*_h1_scaffolding_ragtag.scaffold.agp: This APG file contains the apg path for scaffolding the assembly
	rpmask_bed_tracks: This folder contains all repeatmasker annotation tracks
	trf_bed_tracks: This folder contains all tandem repeat finder annotation tracks
2. syntenic_comparison_with_T2T-CHM13
This folder contains 16 bed files as input for SAFFILE website (https://mrvollger.github.io/SafFire/#dataset=USER&ref=USER_REF&query=USER_QUERY).
Each NHP genome assembly and T2T-CHM13 can be visulaized in SAFFILE website.

3. lineage_specific_SV_with_annotation
This folder contains 2 files related to the lineage specific SVs and their annotation.
	#TWO FILES CAN BE VISULAIZED IN IGV/JBBROSWER
	*all_genotypes.vcf.gz: This file refers to the genotyping raw result from paragraph
	*asmsv_svpop_panPri.vcf.anno.bed: This file contains SV annoation for all genotyped SVs
4. scripts
