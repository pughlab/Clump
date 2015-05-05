# README #

## Overview ##
CluMP (Clustering of Mate Pairs) is a tool for inferring the presence of structural rearrangements from targeted DNA sequencing data. A two-step algorithm is used to detect these events. The tool first isolates read pairs aligned to different chromosomes or separated by large genomic distances (indicated by large insert sizes) that may indicate a genome rearrangement. Read pairs supporting a similar rearrangement are grouped into clusters and the minimum and maximum genome coordinates are used to define two candidate break point regions, one for each end of a mate pair. The tool next searches within each of these candidate regions for single reads with high base quality (20 default), high mapping quality (30), and a large fraction (33%) of soft-clipped (i.e. unaligned) bases that may indicate possible genomic break-points. A combination of large-insert read-pairs and break-point-spanning reads is taken as indicative of a rearrangement. In addition to the number of supporting events, each candidate breakpoint region is annotated with genome coordinates, a gene annotation (currently derived from the UCSC Gene track) and the number of soft-clipped reads in the region. Note that the tool does not check that the soft-clipped reads support a common breakpoint nor does it compare the soft-clipped region to the paired candidate breakpoint region. However, all reads are written to a new bam file to support manual review of the underlying data. In summary, CluMP extracts large-insert and soft-clipped reads indicative of a rearrangement in targeted DNA sequence data and presents a summary of these counts for further filtering to arrive at a high-confidence list of rearrangements for manual review and confirmation.

## Execution ##

CluMP is written in R and uses the rbamtools package (http://cran.r-project.org/web/packages/rbamtools/rbamtools.pdf) for reading and writing bam files and the “intervals” package (http://cran.r-project.org/web/packages/intervals/intervals.pdf) for analysis of overlapping genome coordinates. CluMP is implemented as a command-line tool that accepts three arguments: config file, input bam file, and output directory. Here is an example call to CluMP from the command line:


```
#!bash

~/R-2.15.1/bin/Rscript CluMP.R CluMP.cfg bamfile.bam outdir
```

## Configuration file ##
interval_list.file     <- "/xchip/cga_home/tpugh/AKT3_MAGI3_cell_lines/AKT3_MAGI3.interval_list"
ucsc.table.file        <- "/xchip/cga_home/tpugh/AKT3_MAGI3_cell_lines/refgene.txt"
remove.duplicates      <- TRUE #Set to TRUE to ignore duplicate-marked reads
min.interval.length    <- 100
min.mapping.quality    <- 30
min.base.quality       <- 20 #Used to isolate soft-clipped reads, discards reads with min.base.quality.frac bases below min.base.quality
min.soft.clipped.bases <- 33 #Number of	softclipped base pairs in read required	to keep	read as	evidence of translocation
min.insert.size        <- 5000
max.cluster.distance   <- 300
max.mate.partners      <- 5
min.cluster.size       <- 2

## Interval list format ##
A simple 5-column, tab-separated text file with no header.

|Chromosome|Start|End|Strand|Description|
|1|243649535|244008584|-|AKT3|
|1|113931475|114230545|+|MAGI3|

## Output file descriptions ##
Output files and directories are prefixed with a sample name assigned by the Profile bioinformatics pipeline. Often these names can include an identifier such as BL-13-A12345 as well as run information such as L_00012340_HC_0000056_78. In the table below, this prefix is denoted by “NAME”.

|Filename|Purpose|Description|
|NAME.CluMP_run_info.tsv|Record-keeping|A list of parameters and program versions used to generate CluMP output as well as the run date and time.|
|NAME.rearrangement_counts.tsv|Primary output|A table of all candidate rearrangements, the genes affected, the genome coordinates of the two breakpoint regions, and the number of large-insert and soft-clipped reads supporting each lesion.|
|NAME.rearrangement_regions.bed|Manual review|A BED file containing the genome coordinates of all candidate rearrangements. This file can be read by the IGV “Region Navigator” to help navigate between breakpoint regions when reviewing the primary or CluMP-generated bam files.|
|NAME.rearrangements.bam|Manual review|A bam file containing only reads supporting the candidate rearrangements identified by CluMP. |
|NAME.rearrangements.bai|Manual review|Index for bam file, required by IGV and other analysis software.|

## Pilot data ##

From the Center for Advanced Molecular Diagnostics, Brigham and Women's Hospital (Neal Lindeman)

To evaluate the ability of the CluMP algorithm to detect known rearrangements, an early version of the software was run on 14 cases with translocations detected using fluorescent in situ hybridization. Using analysis of large-insert and soft-clipped reads, we uncovered evidence supporting the known event in 10 of 14 cases (blank “Read_pair” entries denote 4 additional samples that were not analyzed):

![clump_pilot_data.png](https://bitbucket.org/repo/LjrKnq/images/4010912555-clump_pilot_data.png) 

Notably, all four cases in which a translocation could not be found had IGH rearrangements. This is likely due to incomplete tiling of baits across the IGH region. Furthermore, the full IGH region is not annotated by the UCSC gene track currently used for gene annotation. Therefore the IGH breakpoints are named "intergenic" using the default download of this table. This can be solved by adding additional lines to the ucsc.refgene.txt file that list the genome coordinate for this non-genic region. Alternately, an external program could be used to map the genome coordinates of each breakpoint region to a specific gene or annotation.

## Data review ##

In addition to the table of rearrangements, CluMP outputs a bam file containing only the reads supporting each candidate alteration. Viewing this bam alongside the primary bam file may provide some context to the otherwise isolated reads extracted by CluMP. A second useful visualization method is to view both breakpoint regions simultaneously using IGV’s “View mate region in split screen” option. The following examples show single breakpoint regions:

### Known BCL2-IGH translocation ###
![clump_bcl2_IGH_example.png](https://bitbucket.org/repo/LjrKnq/images/70566353-clump_bcl2_IGH_example.png)
The above screenshot shows reads from the original bam file and reads extracted by CluMP. In the CluMP track, reads colored orange are the large-insert reads with mates mapped to chromosome 14 containing IGH. The multi-colored reads contain soft-clipped bases. In this case, the soft-clipped bases begin at the exact same genome coordinate, suggest this is the exact translocation breakpoint in this sample. Analysis of the large-insert mate-pairs alignments identified two breakpoint regions in this case, indicative of an event more complex than a simple, balanced translocation. In this case, paired-reads mapped to two locations in IGH separated by over 30,000 bp of sequence (see below). As there is no sequence captured and sequenced from this region, there is no way to tell whether these are two different translocations in the same sample or a single, complex event resulting from further rearrangement of the IGH locus. Further investigation of the mapping location of the soft-clipped reads may help to resolve the exact IGH breakpoints in this case.

![clump_bcl2_IGH_example2.png](https://bitbucket.org/repo/LjrKnq/images/4009201026-clump_bcl2_IGH_example2.png)
![clump_bcl2_IGH_example3.png](https://bitbucket.org/repo/LjrKnq/images/610153287-clump_bcl2_IGH_example3.png)

### EML4-ALK ###
![clump_EML4-ALK_example.png](https://bitbucket.org/repo/LjrKnq/images/413535777-clump_EML4-ALK_example.png)
The EML4-ALK fusion arises from a large inversion on chromosome 2. Therefore, there are no interchromosomal mate-pairs in this example. Rather, the large-insert reads are colored grey as both mates map to the same chromosome. Most informative in this case are the soft-clipped reads that indicate a potential short stretch of retained sequence between the two soft-clipped start points. This is very likely a short stretch of sequence homology between the two inversion breakpoints, a common mediator of structural rearrangement in cancer. At the top of the read track, there is one read with two stretches of soft-clipped bases that are inconsistent with the pattern of soft-clipping in the other seven soft-clipped reads. This is likely a mapping artifact and it does not support or refute the presence of an EML4-ALK rearrangement.

### EWS-FLI ###
![clump_EWS-FLI_example.png](https://bitbucket.org/repo/LjrKnq/images/1367458212-clump_EWS-FLI_example.png)
In this example, the large-insert reads are found only on one side of the stack of soft-clipped reads. The soft-clipped reads were found as CluMP searches within a defined window around the cluster of large-insert alignments.

### PML-RARA ###
![clump_PML-RARA_example.png](https://bitbucket.org/repo/LjrKnq/images/2972833986-clump_PML-RARA_example.png)
This appears to be a balanced translocation with excellent read support from large-insert pairs and soft-clipped reads. One of the large-insert reads appears to have a single base crossing the breakpoint and may be counted twice in the CluMP output.

### BCR-ABL replicates ###
![clump_BCR-ABL_example.png](https://bitbucket.org/repo/LjrKnq/images/587266467-clump_BCR-ABL_example.png)
Megan Hanna from the Center for Cancer Genome Discovery provided primary bam files from three replicate runs of a single cancer sample. The above screenshot contains three split-screen views of the ABL1 and BCR breakpoint regions detected in this these three replicates. In this case, there appear to be multiple soft-clipping start sites suggesting there may be complex breakpoint structure or multiple clones in the cancer specimen. Across the replicates, the breakpoint regions and soft-clipping boundaries are consistent, suggesting these are reflective of the underlying tumor DNA sequence and not technical artifacts.

## Common false positive calls ##
### ZNF708 ###
![clump_false_positives_ZNF708.png](https://bitbucket.org/repo/LjrKnq/images/1325591711-clump_false_positives_ZNF708.png)
This repetitive region is homologous to several other locations in the genome. As a result, reads mapped to this location have mates aligned to many different chromosomes as indicated by the myriad read colors in the above screenshot (each color represents the chromosomal to which the read-mate was aligned). There are also several imperfect alignments, as evidenced the many reads with soft-clipped bases. This region may be a good candidate for a blacklist of frequent false positive calls.

### LINC000486 ###
![clump_false_positives_LINC000486.png](https://bitbucket.org/repo/LjrKnq/images/4081660520-clump_false_positives_LINC000486.png)
Many polyG reads are aligned to this region due to this long stretch of G-rich sequence (orange in the reference sequence track). Since the sequence identity is good, these reads have high mapping quality and are therefore used by CluMP. This region may be a good candidate for a blacklist of frequent false positive calls.

## Who do I talk to? ##

* Trevor Pugh, trevor.pugh@utoronto.ca