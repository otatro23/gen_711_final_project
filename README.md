### Gen 711 Final Project
Author: Olivia Tatro

# Assembly and Identififcation of Mystery Bacteria Genomes 
Goals: 
1. Create a shell script that assembles both bacterial genomes using loops
2. Identify the species
3. Create visualizations of the genome assembly 

# Background
- The data are in 4 .fastq files, containing the forward and reverse raw 250bp Illumina reads for two different mystery bacteria.
- The goal is to use these files to do a de novo genome assembly, annotate the genome, and identify the bacteria species. 

# Methods
Input data: the initial data is from samples of unknown bacteria species. They are 250bp long paired-end reads from the Illumina sequencer in fastq file format. 

**FastQC** 
This program takes the forward and reverse reads for a sample and summarizes the base-call qualities. It produced plots that show the base-call qualities, which were in an HTML report that included GC content, the total number of sequences, and the total number of bases.

**Trimmomatic**
Trimmomatic took the forward and reverse reads and trimmed the poor quality parts of the reads off to create better quality reads. Typically, the ends of raw reads have lower base-call quality. Trimmomatic output trimmed reads in the .fastq file format. I called this program through the trim_scriptV2.sh wrapper script on the UNH ron teaching server. 

**FastQC**
At this point in the pipeline, I ran FastQC again on the trimmed reads from Trimmomatic. 

**SPAdes**
This program assembled the raw reads into contigs. So, it took trimmed reads (forward, reverse, paired, and unpaired) in the format of .fastq, and the output was labelled contigs in .fasta file format. SPAdes uses De Bruijin Graphs where each node is a kmer, and the kmers are different sizes. 

**QUAST**
QUAST took the contigs.fasta file (output from SPAdes) and generated reports on the contiguity of the genome in various formats. These reports include the number of contigs, largest contig, Contig N50, GC content, and total length in base pairs. 

**BUSCO**
BUSCO assessed the completeness of the genome assembly. It does this by comparing how many core genes are present in the assembly to the amount of core genes in broad taxonomic groups. The input to this program is the contigs.fasta file from SPAdes, and it also requires specifying that this project has a genome assembly and it is a bacterium. The output is summary sheets and tables of the results. 

**PROKKA**
This program annotated the assembled genome. It takes the contigs.fasta file as input and allows you to specify the minimum contig length. It outputs several different files, one of which is a genome feature file (.gff), and ffn and fna files for nucleotides and amino acid sequences. Basicaly, this program predicts the genes in the assembly and annotates them by function based on a reference database. 

**BLAST**
This is the Basic Local Alignment and Search Tool. This is a powerful program that allows you to "blast" (or search for) a specific sequence in a reference sequence. For This project, BLAST helped identify the species. It takes the contigs.fasta file from SPAdes and outputs a megablast file that contains the contigs and information about them such as which species they are from. I used a pre-written script called blob_blast.sh on UNH's ron teaching server to call BLAST. 

**BWA**
Burrow's Wheeler Aligner is a program used to align reads to reference sequences. In this case, I mapped the reads onto the de novo assembly. BWA can be used for many different purposes, sucha s identifying SNPs between a sample and reference. The input to this program is the contigs.fasta file from SPAdes, and the forward and reverse reads, and it outputs a .bam file. 

**Samtools**
Samtools is a program that allows conversion between sam, bam and cram files, allows you to sort and index files, and calculates statistics and quality checks. I used samtools to remove reads that did not match to the assembly, convert the .bam to a .sam, and sort and index the .bam. 

**Bedtools**
Bedtools can also be used for many different things, but I used it to make coverage tables from the sorted and indexed .bam. I used these tables to create the coverage histograms in Figure 2. 

**Blobtools**
Blobtools is a program that plots the contigs as circles, with the size being proportional to their length, the x axis being GC proportion, and the y axis being log-transformed coverage. It requires the contigs from the assembly, the bam file, and megablast file. The results of this are shown in Figure 1. Blobtools also lets you specify which taxonomic grouping to categorize the contigs by. I ran it with both species and genus arguments. 

Lastly, I filtered blob_taxonomy files to exclude contigs that had less than 500 base pairs or coverage less than 20. 

I used ggplot() in R to create the histograms in Figure 2 from the bedtools coverage table. 


# Results
- Sample 04
  - Species identification: _Eschericia coli_
  - Number of contigs: 533
  - Total length of contigs: 5,278,143 bp
  - Average coverage: 216
  - N50: 393,671
  - GC: 51.3%
- Sample 05
  - Species identification: _Eschericia coli_
  - Number of contigs: 1711
  - Total length of contigs: 5,874,153 bp
  - Average coverage: 174
  - N50: 392,646
  - GC: 50.5%
 
![Sample 04 blobplot](images/blob_plot_species_04.bam0.png)
![Sample 05 blobplot](images/blob_plot_species_05.bam0.png)
Figure 1: Blobplots of the genome assemblies for sample 04 (top) and sample 05 (bottom). Both samples are identified as _Eschericia coli_, but have contamination from other organisms, mostly other bacteria species and humans. These figures plot coverage against GC proportion so that contigs from contamination are separate from the sample contigs. The size of the dots represent the length of the contig, so both plots show relatively larger contigs with high coverage. I used blobtools to create these plots. 

![Sample 04 coverage](images/cov_04.png)
![Sample 05 coverage](images/cov_05.png)
Figure 2: Histograms showing the number of contigs with each amount of coverage in each of the 2 bacterial genome assemblies: sample 04 (top) and sample 05 (bottom). Contigs with less than 20x coverage were filtered out before the construction of these plots. I used ggplot() in R to create these plots. 


# References
- Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403–410. https://doi.org/10.1016/S0022-2836(05)80360-2
- Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021
- Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170
- Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. Bioinformatics (Oxford, England), 29(8), 1072–1075. https://doi.org/10.1093/bioinformatics/btt086
- Laetsch, D. R., & Blaxter, M. L. (2017). BlobTools: Interrogation of genome assemblies. F1000Research, 6, 1287. https://doi.org/10.12688/f1000research.12232.1
- Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics (Oxford, England), 25(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324
- Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics (Oxford, England), 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
- Quinlan, A. R., & Hall, I. M. (2010). BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics (Oxford, England), 26(6), 841–842. https://doi.org/10.1093/bioinformatics/btq033
- Seemann, T. (2014). Prokka: Rapid prokaryotic genome annotation. Bioinformatics (Oxford, England), 30(14), 2068–2069. https://doi.org/10.1093/bioinformatics/btu153
- Simão, F. A., Waterhouse, R. M., Ioannidis, P., Kriventseva, E. V., & Zdobnov, E. M. (2015). BUSCO: Assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31(19), 3210–3212. https://doi.org/10.1093/bioinformatics/btv351
