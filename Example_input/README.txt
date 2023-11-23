# Details about the provided input data: 

Organism: Mycobacterium tuberculosis H37Rv
sRNA-seq data: Strand Specific
SRA ID: SRR7058126 [1]
Sample: Rich medium

#Procedure to create per base coverage file:

a. BAM file converted into a sorted BAM file.
Command: samtools sort <bam file> -o <output sorted bam file>
b. Strand-specific sRNA-seq data requires splitting the BAM file into forward and reverse files. Non-strand-specific RNA-seq data uses a single BAM file.
c. BAM files need to be converted into BED files using the BEDtools command.
Command: bedtools genomcov -ibam <bam file> -d <output bed file>
d. The per base (.bed) file consists of three columns: Chromosome, Position, and Read Count. This file format serves as the input for POSIF. In the case of strand-specific data, two separate input files are required: Forward.bed and Reverse.bed. However, for non-strand-specific data, a single .bed file is used as the input for POSIF.

Reference:
1. Gerrick, E. R et al. (2018) Small RNA profiling in Mycobacterium tuberculosis identifies MrsI as necessary for an anticipatory iron sparing response.  Proc Natl Acad Sci U S A, 115(25), 6464–6469. https://doi.org/10.1073/pnas.1718003115


