### Gen 711 Final Project
*I will be working on this project alone*

# Assembly and Identififcation of Mystery Bacteria Genomes 
Goals: 
1. Create a shell script that assembles both bacterial genomes using loops
2. Identify the species
3. Create visualizations of the genomes and phylogeny 

# Background
- The data are in 4 .fastq files, containing the forward and reverse raw 250bp Illumina reads for two different mystery bacteria.
- The goal is to use these files to do a de novo genome assembly, annotate the genome, and identify the bacteria species. 

# Methods
Input data: the initial data is from samples of unknown bacteria species. They are 250bp long paired-end reads from the Illumina sequencer in fastq file format. 

**FastQC** 
This program takes the forward and reverse reads for a sample and summarizes the base-call qualities. It produced plots that show the base-call qualities, which were in an HTML report that included GC content, the total number of sequences, and the total number of bases.

**Trimmomatic**
Trimmomatic took the forward and reverse reads and trimmed the poor quality parts of the reads off to create better quality reads. Typically, the ends of raw reads have lower base-call quality. Trimmomatic output trimmed reads in the .fastq file format. I called this program through the trim_scriptV2.sh wrapper script on the UNH ron teaching server. 

**SPAdes**
This program assembled the raw reads into contigs. So, it took trimmed reads (forward, reverse, paired, and unpaired) in the format of .fastq, and the output was labelled contigs in .fasta file format. SPAdes uses De Bruijin Graphs where each node is a kmer, and the kmers are different sizes. 



# Results
- Sample 04 and sample 05 were both identified as _Eschericia coli_.
- The average coverage for sample 04 was 216 and for sample 05 was 174. 

![Sample 04 blobplot](images/blob_plot_species_04.bam0.png)
![Sample 05 blobplot](images/blob_plot_species_05.bam0.png)
Figure 1: Blobplots of the genome assemblies for sample 04 (top) and sample 05 (bottom). Both samples are identified as _Eschericia coli_, but have contamination from other organisms, mostly other bacteria species and humans. These figures plot coverage against GC proportion so that contigs from contamination are separate from the sample contigs. The size of the dots represent the length of the contig, so both plots show relatively larger contigs with high coverage. I used blobtools to create these plots. 

![Sample 04 coverage](images/cov_04.png)
![Sample 05 coverage](images/cov_05.png)
Figure 2: Histograms showing the number of contigs with each amount of coverage in each of the 2 bacterial genome assemblies: sample 04 (top) and sample 05 (bottom). Contigs with less than 20x coverage were filtered out before the construction of these plots. I used ggplot() in R to create these plots. 
