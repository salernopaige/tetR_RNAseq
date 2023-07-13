## tetR RNAseq:
RNA-Seq analysis of *Bacteroides fragilis* NCTC 9343 samples for tetR project. Two strains, WT and ΔtetR were sequenced to look at gene expression differences between these two treatment groups. Samples were sent to SeqCoast for sequencing, 12M reads, 2x150bp sequencing. 

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
Sample 1: WT-1  
Sample 2: ΔtetR-1  
Sample 3: WT-2  
Sample 4: ΔtetR-2  
Sample 5: WT-3  
Sample 6: ΔtetR-3