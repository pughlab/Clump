# README #

### Overview ###
CluMP (Clustering of Mate Pairs) is a tool for inferring the presence of structural rearrangements from targeted DNA sequencing data. A two-step algorithm is used to detect these events. The tool first isolates read pairs aligned to different chromosomes or separated by large genomic distances (indicated by large insert sizes) that may indicate a genome rearrangement. Read pairs supporting a similar rearrangement are grouped into clusters and the minimum and maximum genome coordinates are used to define two candidate break point regions, one for each end of a mate pair. The tool next searches within each of these candidate regions for single reads with high base quality (20 default), high mapping quality (30), and a large fraction (33%) of soft-clipped (i.e. unaligned) bases that may indicate possible genomic break-points. A combination of large-insert read-pairs and break-point-spanning reads is taken as indicative of a rearrangement. In addition to the number of supporting events, each candidate breakpoint region is annotated with genome coordinates, a gene annotation (currently derived from the UCSC Gene track) and the number of soft-clipped reads in the region. Note that the tool does not check that the soft-clipped reads support a common breakpoint nor does it compare the soft-clipped region to the paired candidate breakpoint region. However, all reads are written to a new bam file to support manual review of the underlying data. In summary, CluMP extracts large-insert and soft-clipped reads indicative of a rearrangement in targeted DNA sequence data and presents a summary of these counts for further filtering to arrive at a high-confidence list of rearrangements for manual review and confirmation.

### Execution ###

CluMP is written in R and uses the rbamtools package (http://cran.r-project.org/web/packages/rbamtools/rbamtools.pdf) for reading and writing bam files and the “intervals” package (http://cran.r-project.org/web/packages/intervals/intervals.pdf) for analysis of overlapping genome coordinates. CluMP is implemented as a command-line tool that accepts three arguments: config file, input bam file, and output directory. Here is an example call to CluMP:

~/R-2.15.1/bin/Rscript CluMP.R CluMP.cfg bamfile.bam outdir

### Who do I talk to? ###

* Trevor Pugh, trevor.pugh@utoronto.ca