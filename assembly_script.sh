#! /bin/bash

# data directory
dir=/home/users/omt1027/data


# activate genomics environment
source activate genomics 


# calculate average coverage
echo Calculating average coverage
for gz in $(ls $dir/gzs); do # iterates through each gz file
  num_reads=$(zgrep -c "^@" $dir/gzs/$gz) # counts the number of reads 
  echo "Average coverage for " $gz ":"  $(( num_reads * 250 * 2 / 7000000)) # multiplies num_reads by the number of basepairs in each then divides by the genome length
  #: # colon is so that bash doesn't give errors for an empty for loop 
done  



# run fastqc to generate quality scores for the raw reads
echo Running fastqc
mkdir $dir/fastqc_out # fastqc requires the output directory to already be made
fastqc $dir/gzs/* -o $dir/fastqc_out
# to view the results, open sftp, navigate to fastqc_out, then get *.html



# run trimmomatic through the trim/scriptV2.sh script 
echo running trimmomatic
/usr/local/bin/trim_scriptV2.sh $dir/gzs/04*
/usr/local/bin/trim_scriptV2.sh $dir/gzs/05*
mkdir $dir/trimmed_reads # create output directory
mv ./trimmed-reads/* $dir/trimmed_reads #move reads to output directory


### run fastq on trimmed reads
echo running fastq on the trimmed reads
mkdir $dir/trimmed_fastqc
fastqc $dir/trimmed_reads/* -o $dir/trimmed_fastqc



# run SPAdes on trimmed reads to assemble them
echo running SPAdes 
mkdir $dir/spades_output_04
mkdir $dir/spades_output_05
t=$dir/trimmed_reads #trimmed directory

# this is so that SPAdes will take the file 
for file in $t/* ; do 
  mv $file $file.fastq.gz
  #:
done

# actually running SPAdes
spades.py -1 $t/04*R1* -2 $t/04*R2* -s $t/unpaired-04*R1* -s $t/unpaired-04*R2* -o $dir/spades_output_04
spades.py -1 $t/05*R1* -2 $t/05*R2* -s $t/unpaired-05*R1* -s $t/unpaired-05*R2* -o $dir/spades_output_05



# run quast to look at the contiguity of the assembly
echo running QUAST
mkdir $dir/quast_out_04
mkdir $dir/quast_out_05
quast.py $dir/spades_output_04/contigs.fasta -o $dir/quast_out_04
quast.py $dir/spades_output_05/contigs.fasta -o $dir/quast_out_05 


# run busco to assess the completeness of the assembly
echo running BUSCO
mkdir $dir/busco_out_04
mkdir $dir/busco_out_05
busco -i $dir/spades_output_04/contigs.fasta -m genome -o busco_out_04 --out_path $dir -l bacteria -f
busco -i $dir/spades_output_05/contigs.fasta -m genome -o busco_out_05 --out_path $dir -l bacteria -f


# run prokka for genome annotation
echo running PROKKA
prokka $dir/spades_output_04/contigs.fasta --outdir $dir/prokka_out_04 --cpus 24 --mincontiglen 200
prokka $dir/spades_output_05/contigs.fasta --outdir $dir/prokka_out_05 --cpus 24 --mincontiglen 200



# BLAST 
echo running BLAST
mkdir $dir/contigs_db_04
mkdir $dir/contigs_db_05 
# make the blastdb
makeblastdb -in $dir/spades_output_04/contigs.fasta -dbtype nucl -out $dir/contigs_db_04/contigs_db
makeblastdb -in $dir/spades_output_05/contigs.fasta -dbtype nucl -out $dir/contigs_db_05/contigs_db
# run blast
blast-ncbi-nt.sh $dir/spades_output_04/contigs.fasta
mv ./contigs.fasta.vs.nt.cul5.maxt10.1e5.megablast.out $dir/megablast_04.out
blast-ncbi-nt.sh $dir/spades_output_05/contigs.fasta
mv ./contigs.fasta.vs.nt.cul5.maxt10.1e5.megablast.out $dir/megablast_05.out



# read mapping to de novo assembly
echo read mapping
#index reference genome
bwa index $dir/spades_output_04/contigs.fasta 
bwa index $dir/spades_output_05/contigs.fasta
#map reads and construct a SAM file
bwa mem -t 24 $dir/spades_output_04/contigs.fasta $dir/trimmed_reads/04*R1* $dir/trimmed_reads/04*R2* > $dir/raw_mapped_04.sam
bwa mem -t 24 $dir/spades_output_05/contigs.fasta $dir/trimmed_reads/05*R1* $dir/trimmed_reads/05*R2* > $dir/raw_mapped_05.sam


#construct a coverage table
echo construct a coverage table
#remove sequencing reads that did not map to the assembly and convert the SAM to a BAM 
samtools view -@ 24 -Sb $dir/raw_mapped_04.sam | samtools sort -@ 24 -o $dir/sorted_mapped_04.bam 
samtools view -@ 24 -Sb $dir/raw_mapped_05.sam | samtools sort -@ 24 -o $dir/sorted_mapped_05.bam
# examine how many reads mapped with samtools
samtools flagstat $dir/sorted_mapped_04.bam
samtools flagstat $dir/sorted_mapped_05.bam
# index the new .bam files
samtools index $dir/sorted_mapped_04.bam
samtools index $dir/sorted_mapped_05.bam 
# create coverage tables 
bedtools genomecov -ibam $dir/sorted_mapped_04.bam > $dir/coverage_04.out
bedtools genomecov -ibam $dir/sorted_mapped_05.bam > $dir/coverage_05.out


# copy the python script and give myself permission to run it
cp /tmp/genome_back/gen_input_table.py ./
chmod u+rwx gen_input_table.py
# calculate coverage per contig
./gen_input_table.py --isbedfiles $dir/spades_output_04/contigs.fasta $dir/coverage_04.out > $dir/coverage_table_04.tsv
./gen_input_table.py --isbedfiles $dir/spades_output_05/contigs.fasta $dir/coverage_05.out > $dir/coverage_table_05.tsv


# run blob tools
echo running blob tools
blobtools create -i $dir/spades_output_04/contigs.fasta -b $dir/sorted_mapped_04.bam -t $dir/megablast_04.out -o $dir/blob_out_04
blobtools create -i $dir/spades_output_05/contigs.fasta -b $dir/sorted_mapped_05.bam -t $dir/megablast_05.out -o $dir/blob_out_05

blobtools view -i $dir/*04*.json -r all -o $dir/blob_taxonomy
blobtools view -i $dir/*05*.json -r all -o $dir/blob_taxonomy

blobtools plot -i $dir/blob_out_04.blobDB.json -r species -o $dir/blob_plot_species_04
blobtools plot -i $dir/blob_out_05.blobDB.json -r species -o $dir/blob_plot_species_05

echo running blob plot
mv $dir/blob_plot*04*bam* $dir/blobplot_04.bam0.png
mv $dir/blob_plot*04*cov* $dir/blobplot_04.read_cov.png
mv $dir/*04*stats* $dir/blobplot_04.stats.txt
mv $dir/blob_plot*05*bam* $dir/blobplot_05.bam0.png
mv $dir/blob_plot*05*cov* $dir/blobplot_05.read_cov.png
mv $dir/*05*stats* $dir/blobplot_05.stats.txt


#filtering by length > 500 and coverage > 20
echo filter genome assembly
contigs_to_keep_04=$(grep -v '#' $dir/blob_tax*04* | awk -F'\t' '$2 < 500' | awk -F'\t' '$5 > 20' | cut -f1)
for contig in $contigs_to_keep_04; do
	#:
	grep -A1 "$contig" $dir/spades_output_04/contigs.fasta >> $dir/filtered_04.fasta
done

contigs_to_keep_05=$(grep -v '#' $dir/blob_tax*05* | awk -F'\t' '$2 < 500' | awk -F'\t' '$5 > 20' | cut -f1)
for contig in $contigs_to_keep_05; do
	#:
	grep -A1 "$contig" $dir/spades_output_05/contigs.fasta >> $dir/filtered_05.fasta
done


# calculate average coverage (before filter)
echo calculating average coverage
# 215.88
cat $dir/blob_tax*04* | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}'
# 173.588
cat $dir/blob_tax*05* | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}'



















