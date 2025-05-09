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

# Results
- Species Identification: Both species were _Eschericia coli_.

![Sample 04 blobplot](images/blob_plot_species_04.bam0.png)
![Sample 05 blobplot](images/blob_plot_species_05.bam0.png)
Figure 1: Blobplots of the genome assemblies for sample 04 (top) and sample 05 (bottom). Both samples are identified as _Eschericia coli_, but have contamination from other organisms, mostly bacteria and _Homo sapiens_. 

![Sample 04 coverage](images/cov_04.png)
![Sample 05 coverage](images/cov_05.png)
Figure 2: Histograms showing the coverage of each contig in each of the 2 bacterial genome assemblies. Contigs with less than 20x coverage were filtered out. 
