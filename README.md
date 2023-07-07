# TUM Bioinformatics Java scripts

Scripts from my undergraduate degree in bioinformatics.

The `algorithmic-bioinformatics` folder contains snippets created in the scope of coursework on **algorithms** that 
played a central role in the development of the field of **dynamic programming** and thus bioinformatics. These include 
the Boyer-Moore and Knuth-Morris-Pratt string/pattern searching algorithms as well as code for finding the maximal 
scoring subsequence from a finite sequence of real numbers.

The `src` directories contain more complex Java code for tackling various challenges related to **genome-oriented 
bioinformatics** of RNA-seq data.
- _Exon skipping_: Implement a program to extract exon-skipping events (a form of alternative splicing) between 
  coding DNA sequences for a .gtf file. Create an overview as well as two cumulative plots to analyse the found 
  distributions.
- _Read simulation_: Implement a paired-end read simulator, allowing for a variety of settings. An exemplary simulation 
  is to be generated including a visual analysis of the result distribution.
- _BAM feature counting_: (A BAM file (.bam) is the binary version of a SAM file.  A SAM file (.sam) is a 
  tab-delimited text file that contains sequence alignment data.) Implement a program extracting information on all 
  mapped read pairs, genic and intergenic, 
  sensitive to experiment strandness, annotating several read pair features as defined. Include plots analysing the 
  RPKM (reads per kilobase of transcript per million reads mapped) distribution.
- _Gene enrichment_: Implement an algorithm to analyse Gene Ontology (GO) DAG overlap properties and to perform 
  enrichment analysis on simulated data, including a variety of statistical tests.

### One second --
-- If any current Bioinformatics students should find this repository, I would kindly advise you to not risk
plagiarism - most that do get caught. Besides, the extent to which my code is doing what your specific challenge 
requires might be limited. Good luck with GoBi & Algo I! 🤞


