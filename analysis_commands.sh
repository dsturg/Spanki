# Splicing Analysis Toolkit (Spanki)
# Example code for an analysis of splicing
# Uses testdata availble from:
# http://www.cbcb.umd.edu/software/spanki/testdata.tar.gz

# Instructions:
# Make a working directory, eg "spankitest"
# go the the spankitest directory, and extract the testdata there
# copy the "analysis_commands.sh" file from the Spanki directory to here
# Run the commands:  ./analysis_commands.sh  or manually copy and paste lines of code

echo "[***************] Starting female spankijunc runs"
spankijunc -m all -o female_rep1 -i testdata/female_r1.bam -g testdata/annotation/genemodels.gtf -f testdata/fasta/myref.fa
spankijunc -m all -o female_rep2 -i testdata/female_r2.bam -g testdata/annotation/genemodels.gtf -f testdata/fasta/myref.fa
echo "[***************] Starting male spankijunc runs"
spankijunc -m all -o male_rep1 -i testdata/male_r1.bam -g testdata/annotation/genemodels.gtf -f testdata/fasta/myref.fa
spankijunc -m all -o male_rep2 -i testdata/male_r2.bam -g testdata/annotation/genemodels.gtf -f testdata/fasta/myref.fa

# Note:  If you want all your tables to be symmetrical (have the same number of rows and the same junctions),
# then create a list of junctions (either curated, or a merge of all junctions detected), and run make_curated_jtab
# on each sample. Example:
	
echo "[***************] Making a table of all junctions detected in any sample"
cat female_rep1/juncs.all female_rep2/juncs.all male_rep1/juncs.all male_rep2/juncs.all > all_tables.txt
head -1 female_rep1/juncs.all | awk '{print $1"\t"$15"\t"$8"\t"$5"\t"$16"\t"$9"\t"$10"\t"$3}' > all_juncs.txt 
awk '{if ($1 != "juncid") {print $1"\t"$15"\t"$8"\t"$5"\t"$16"\t"$9"\t"$10"\t"$3}}' all_tables.txt | sort -u >> all_juncs.txt

echo "[***************] Getting additional data from the female runs"
make_curated_jtab -o female_rep1_curated -i  testdata/female_r1.bam -jlist all_juncs.txt -jtab female_rep1/juncs.all
make_curated_jtab -o female_rep2_curated -i  testdata/female_r2.bam -jlist all_juncs.txt -jtab female_rep2/juncs.all
echo "[***************] Getting additional data from the male runs"
make_curated_jtab -o male_rep1_curated -i  testdata/male_r1.bam -jlist all_juncs.txt -jtab male_rep1/juncs.all
make_curated_jtab -o male_rep2_curated -i  testdata/male_r2.bam -jlist all_juncs.txt -jtab male_rep2/juncs.all

echo "[***************] Merging the replicates together"
merge_jtabs female_rep1_curated/juncs.all,female_rep2_curated/juncs.all > female_repsmerged.juncs
merge_jtabs male_rep1_curated/juncs.all,male_rep2_curated/juncs.all > male_repsmerged.juncs

# Alternatively, just merge the tables without the extra step of making the tables symmetrical:
#merge_jtabs female_rep1/juncs.all,female_rep2/juncs.all > female_repsmerged.juncs
#merge_jtabs male_rep1/juncs.all,male_rep2/juncs.all > male_repsmerged.juncs

echo "[***************] Parsing events output"
spankisplice -o female_repsmerged_events -j female_repsmerged.juncs -c testdata/female_cuff/isoforms.fpkm_tracking -f testdata/fasta/myref.fa -g testdata/annotation/genemodels.gtf -a testdata/annotation/genemodels_splices.out 
spankisplice -o male_repsmerged_events -j male_repsmerged.juncs -c testdata/male_cuff/isoforms.fpkm_tracking -f testdata/fasta/myref.fa -g testdata/annotation/genemodels.gtf -a testdata/annotation/genemodels_splices.out 

echo "[***************] Doing pairwise comparison"
splicecomp -a female_repsmerged_events/events.out -b male_repsmerged_events/events.out -o F_vs_M_splicecomp
junccomp -a female_repsmerged.juncs -b male_repsmerged.juncs -o F_vs_M_junccomp

########################
# Other functions      #
########################

echo "[***************] Annotating a junction set (From a junction table) "
annotate_junctions -o female_rep1_annotated_junctab -jtab female_rep1/juncs.all -f testdata/fasta/myref.fa -g testdata/annotation/genemodels.gtf
# Annotating a junction set (From a junction list)
annotate_junctions -o female_rep1_annotated_junclist -jlist female_rep1/juncs.list -f testdata/fasta/myref.fa -g testdata/annotation/genemodels.gtf
# Annotating annotated junctions
annotate_junctions -o annotated_reference_juncs -f testdata/fasta/myref.fa -g testdata/annotation/genemodels.gtf



