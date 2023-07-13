## tetR RNAseq:
RNA-Seq analysis of *Bacteroides fragilis* NCTC 9343 samples for tetR project. Two strains, WT and ΔtetR were sequenced to look at gene expression differences between these two treatment groups. Samples were sent to SeqCoast for sequencing, 12M reads, 2x150bp sequencing on an Illumina NextSeq2000.  

Read demultiplexing, trimming, and analytics were performed by SeqCoast using DRAGEN v3.10.12, and fastq metrics were reported along with the raw sequencing files. 

### Software Needed:  
[fastp](https://github.com/OpenGene/fastp) - For trimming the raw reads  
[bwa](https://github.com/lh3/bwa) - For alignment to the *B. fragilis* genome  
[HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/count.html) -  For gene counts per sample  
[Snakemake 6](https://snakemake.readthedocs.io/en/v6.0.0/getting_started/installation.html) - For automating each of the steps involved in pre-processing,and maintaining the process in one place  
[Python 3](https://www.python.org/) - To make conda environment specific to these processes  
[R](https://www.r-project.org/about.html) - For Statistical analysis and graphing  
[tidyverse](https://www.tidyverse.org/) - For data analysis in R  
[ggplot2](https://ggplot2.tidyverse.org/) - For making figures

### Sample Outline:  
The initial naming convention used by SeqCoast is OrderNumber_SeqCoastTubeID_IlluminaSampleSheetID_Read1orRead2. This can be mapped back to the sample sheet provided in the reported data. The first thing we did was rename the samples though, so that they made more sense to us and minimized the length of the name. The strain information for our samples is as follows:  

Sample 1: WT-1  
Sample 2: ΔtetR-1  
Sample 3: WT-2  
Sample 4: ΔtetR-2  
Sample 5: WT-3  
Sample 6: ΔtetR-3  

We thus renamed our samples with the following convention: StrainReplicate_R1orR2.fastq.gz. For example, Sample 1 Read 1 would be WT1_R1.fastq.gz, and Sample 2 Read 1 would be tetR1_R1.fastq.gz. 