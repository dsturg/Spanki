ECHO "[***************] Starting female spankyjunc runs"
spankyjunc -m all -o female_rep1 -i testdata/female_r1.bam -g testdata/annotation/genemodels.gtf -f testdata/fasta/myref.fa
spankyjunc -m all -o female_rep2 -i testdata/female_r2.bam -g testdata/annotation/genemodels.gtf -f testdata/fasta/myref.fa
ECHO "[***************] Starting male spankyjunc runs"
spankyjunc -m all -o male_rep1 -i testdata/male_r1.bam -g testdata/annotation/genemodels.gtf -f testdata/fasta/myref.fa
spankyjunc -m all -o male_rep2 -i testdata/male_r2.bam -g testdata/annotation/genemodels.gtf -f testdata/fasta/myref.fa

ECHO "[***************] Merging the replicates together"
merge_jtabs_all female_rep1/juncs.all,female_rep2/juncs.all > female_repsmerged.juncs
merge_jtabs_all male_rep1/juncs.all,male_rep2/juncs.all > male_repsmerged.juncs

ECHO "[***************] Parsing asta output"
spankysplice -o female_repsmerged_asta -j female_repsmerged.juncs -c testdata/female_cuff/isoforms.fpkm_tracking -f testdata/fasta/myref.fa -g testdata/annotation/genemodels.gtf -a testdata/annotation/genemodels_splices.out 
spankysplice -o male_repsmerged_asta -j male_repsmerged.juncs -c testdata/male_cuff/isoforms.fpkm_tracking -f testdata/fasta/myref.fa -g testdata/annotation/genemodels.gtf -a testdata/annotation/genemodels_splices.out 

ECHO "[***************] Doing pairwise comparison"
splicecomp -a female_repsmerged_asta/asta.out -b male_repsmerged_asta/asta.out -o F_vs_M_splicecomp
junccomp -a female_repsmerged.juncs -b male_repsmerged.juncs -o F_vs_M_junccomp
