for i in `cat samples.fofn`; do mkdir $i ; cp ../../NHP_HiFi_RagTag_scaffolding/${i}/*.fasta ${i}/; cp ../../NHP_HiFi_RagTag_scaffolding/${i}/*.gff ${i}/; cp ../../NHP_HiFi_RagTag_scaffolding//${i}/ragtag.scaffold.agp ${i}/${i}_ragtag.scaffold.agp; cp -r ../../NHP_HiFi_RagTag_scaffolding//${i}/*scaffold.fasta_rp/rmskClass/ ${i}/rpmask_bed_tracks; mkdir ${i}/trf_bed_tracks; cp ../../NHP_HiFi_RagTag_scaffolding//${i}/*trf.filter.* trf_bed_tracks; done
 for i in `cat samples.fofn`; do mkdir -p ${i}/trf_bed_tracks && cp ../../NHP_HiFi_RagTag_scaffolding//${i}/TRF.bb ${i}/trf_bed_tracks; done
