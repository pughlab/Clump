#Mate-pair finder
#By Trevor Pugh, June 2013

library("intervals")
library("rbamtools") #requires R-2.15
library("Rsamtools") 
library("GenomicAlignments")
library(logging)
svn.revision <- "$Id$"

###########
#Arguments#
###########

arguments <- commandArgs(trailingOnly = TRUE)
config.file <- arguments[1]
bam.file <- arguments[2]
out.dir  <- arguments[3]

#Test data sets for development, override command-line arguments
#config.file <- "/xchip/cga_home/tpugh/AKT3_MAGI3_cell_lines/CluMP.cfg"
#bam.file    <- "/seq/picard_aggregation/C1513/SUM159_DNA/v1/SUM159_DNA.bam"
#bam.file    <- "/seq/picard_aggregation/C1513/T47D_DNA/v1/T47D_DNA.bam"
#out.dir     <- "/xchip/cga_home/tpugh/AKT3_MAGI3_cell_lines"

if(!file.exists(config.file))  {
    stop(paste("Config file not found:", config.file))
} else {
    source(config.file)
}

# start logging to STDOUT
basicConfig(level='FINEST')
addHandler(writeToConsole)
setLevel(30, getHandler('basic.stdout'))


###########
#Functions#
###########

load_interval_names <- function(interval_list.file) {
    #input:  interval_list.file = string containing a file path
    #output: interval_list = matrix with rownames listing genome coordinates, colnames listing data types
   interval_list <- read.table(interval_list.file, header=FALSE, comment.char = "@", stringsAsFactors=FALSE)
   colnames(interval_list) <- c("chr", "start", "end", "strand", "interval_name")
   rownames(interval_list) <- paste(interval_list$chr, ":", interval_list$start, "-", interval_list$end, sep = "")
   return(interval_list)
}

lookup_reference_number <- function(refs, interval_list) {
    #Rather than string matching, rbamtools uses numbered references, 0-indexed in the
    #order that they appear in the bam file. Therefore, we need to convert our chromosome
    #numbers to reference numbers.
    #input=refs derived from bam file using rbamtools getRefData.
    #output=reference numbers determined through simple lookup.
    refNumbers <- refs[match(interval_list[,"chr"], refs$SN), "ID"]
    return(refNumbers)
}

get_mate_coord <- function(reads.bamRange) {
    rewind(reader)
    coords <- c()
    align <- getNextAlign(reads.bamRange)
    while(!is.null(align)) {
        mate.refID      <- mateRefID(align)
        mate.position   <- matePosition(align)
        mate.comparison <- (mate.refID * 10^9) + mate.position
        coord <- c(mate.refID, mate.position, mate.comparison)
        coords <- append(coords, coord)
        align <- getNextAlign(reads.bamRange)
    }
    coords.mat <- matrix(data=coords, ncol=3, byrow=TRUE, dimnames=list(c(), c("mate.refid", "mate.position", "mate.genome.location")))
    return(coords.mat)
}

get_mate_read <- function(align, reader) {
    read.mate <- c()
    mate.refID      <- mateRefID(align)
    mate.position   <- matePosition(align)
    mate.range <- bamRange(reader, c(mate.refID, mate.position, mate.position+1))
    mate.align <- getNextAlign(mate.range)
    while(!is.null(mate.align)) {
        if(name(mate.align) == name(align)) {
            read.mate <- mate.align
            break
        }
        mate.align <- getNextAlign(mate.range)
    }        
    return(read.mate)
}

get_candidate_translocation_reads <- function(reader, interval_list, min.insert.size, min.interval.length, remove.duplicates, min.mapping.quality) {
    require("rbamtools")
    #Isolate read pairs with intrachromosomal mates or large insert size
    #Input:
    #   reader           = bam object
    #   interval_list    = matrix of interval coordinates and names derived from interval_list file
    #   min.insert.size  = single integer used as threshold to identify read pairs with large inserts
    #Output:
    #   reads.candidates = bamRange object containing read pairs
    
    #Initialize tracking variable, then loop through interval_list entries
    reads.candidates <- bamRange()
    for(i in 1:dim(interval_list)[1]) {
        #Isolate coordinates, skip coordinates below min.interval.length
        coords <- as.integer(interval_list[i, c("ref", "start", "end")])
        if((coords[3] - coords[2]) < min.interval.length) {
            write(paste(sep=", ",
                        paste(sep=" ", "Skipping interval", i, "of", dim(interval_list)[1]),
                        paste(coords[3] - coords[2], "bp", sep=" "),
                        rownames(interval_list)[i]),
                  stdout())
            next
        } else {
	    loginfo("Processing interval %d of %d", i, dim(interval_list)[1])
#            write(paste(sep=", ",
#                        paste(sep=" ", "Processing interval", i, "of", dim(interval_list)[1]),
#                        paste(coords[3] - coords[2], "bp", sep=" "),
#                        rownames(interval_list)[i]),
#                  stdout())
        }        

        #Reset bam file and initialize tracking variables
        rewind(reader)
        reads.interval <- bamRange(reader, coords)
        align <- getNextAlign(reads.interval)

        #Read alignments one-by-one to identify read pairs that may support translocations
        while(!is.null(align)) {
            #Require reads to have adequate mapping quality and to be a non-duplicate (if requested)
            if( (mapQuality(align) >= min.mapping.quality) &
                !((remove.duplicates=="TRUE") & (pcrORopt_duplicate(align)=="TRUE"))) { #Read only non-duplicate reads
                    #align.mate <- get_mate_read(align, reader)
                    mate.refID      <- mateRefID(align)
                    mate.position   <- matePosition(align)

                    if( (mate.refID != refID(align)) | (insertSize(align) >= min.insert.size) ) {
                        push_back(reads.candidates, align)
                    }
            }
            align <- getNextAlign(reads.interval)
        }
    }
    return(reads.candidates)
}

remove_singleton_outliers <- function(reads.candidates, max.cluster.distance, min.cluster.size) {
    require("rbamtools")
    #Covert bamRange object to data.frame (df)
    df.reads <- as.data.frame(reads.candidates)
    df.reads <- cbind(df.reads, get_mate_coord(reads.candidates))

    #Cluster read starts and store reads that pass filters
    df.reads.filtered <- c()
    for(ref in unique(df.reads$refid)) { #Query df.reads on same chromosome
        df.reads.ref <- df.reads[df.reads$refid == ref,]
        #Count df.reads near one another and remove if cluster is too small
        df.reads.within.cluster.distance <- sapply(df.reads.ref$position, function(x) sum(abs((df.reads.ref$position - x)) <= max.cluster.distance)-1)
        mates.within.cluster.distance <- sapply(df.reads.ref$mate.genome.location, function(x) sum(abs((df.reads.ref$mate.genome.location - x)) <= max.cluster.distance)-1)
        clustered.pairs <- (df.reads.within.cluster.distance >= min.cluster.size-1) & (mates.within.cluster.distance >= min.cluster.size-1)
        df.reads.filtered <- rbind(df.reads.filtered, df.reads.ref[clustered.pairs,])
    }

    #Store reads as bamRange object
    reads.filtered <- bamRange()
    rewind(reader)
    align <- getNextAlign(reads.candidates)
    #Read alignments one-by-one to identify read pairs that may support translocations
    for(i in 0:(rbamtools::size(reads.candidates)-1)) {
        if(i %in% rownames(df.reads.filtered)) {
            push_back(reads.filtered, align)
        }
        align <- getNextAlign(reads.candidates)        
    }
    return(reads.filtered)
}

get_softclipped_reads <- function(candidate.interval_list, reader, min.soft.clipped.bases, min.base.quality, remove.duplicates=TRUE) {
    reads.softclipped <- bamRange()
    for(i in 1:dim(candidate.interval_list)[1]) {
        loginfo("outer for loop iteration number %d", i)
        coords <- as.integer(candidate.interval_list[i, c("ref", "start", "end")])
        rewind(reader)
        reads.interval <- bamRange(reader, coords, complex=TRUE) #extract reads with complex cigar strings
        align <- getNextAlign(reads.interval)
    
        #Read alignments one-by-one to identify read pairs that may support translocations
	count=0
        while(!is.null(align)) {
            if(!((remove.duplicates=="TRUE") & (pcrORopt_duplicate(align)=="TRUE"))) { #Read only non-duplicate reads
                #Extract high quality soft-clipped reads nearby candidate translocations
                cigar.data         <- cigarData(align)
                base.qualities     <- as.numeric(charToRaw(alignQual(align)))[1:sum(cigar.data$Length)]
                softclip <- match("S", cigar.data$Type)
                if(!is.na(softclip)) {
                    if((cigar.data$Length[softclip] >= min.soft.clipped.bases) & (min(base.qualities, na.rm=TRUE) >= min.base.quality)) {
                            push_back(reads.softclipped, align)
                    }
                }
            }
            align <- getNextAlign(reads.interval)
	    count = count +1
        }
	loginfo("while loop had %d iterations", count)
    }
    return(reads.softclipped)
}

add_genes_to_reads <- function(reads, interval_list, ucsc.table) {
    reads.filtered.genes <- c()
    for(i in 1:dim(reads)[1]) {
        read <- reads[i,]
        interval <- interval_list[(interval_list$chr %in% read["chr"]) &
                                  (as.numeric(read["position"]) >= interval_list$start) &
                                  (as.numeric(read["position"]) <= interval_list$end),"interval_name"]
        mate.interval <- interval_list[(interval_list$chr %in% read["mate.chr"]) &
                                       (as.numeric(read["mate.position"]) >= interval_list$start) &
                                       (as.numeric(read["mate.position"])  <= interval_list$end),"interval_name"]
        #Lookup mate.interval in UCSC table if not found on interval list
        if(length(interval) == 0) {
            interval <- paste(unique(ucsc.table[(ucsc.table$chrom %in% paste("chr", read["chr"], sep="")) &
                               (as.numeric(read["position"]) >= ucsc.table$txStart) &
                               (as.numeric(read["position"]) <= ucsc.table$txEnd),"name2"]), collapse=",")
           loginfo("found %s for %s %d %d", interval, paste("chr", read["chr"], sep=""), as.numeric(read["position"]), as.numeric(read["position"]))
            if((length(interval) == 0) | (interval == "")) {
                    interval <- "intergenic"
            }
        }

        if(length(mate.interval) == 0) {
            mate.interval <- paste(unique(ucsc.table[(ucsc.table$chrom %in% paste("chr", read["mate.chr"], sep="")) &
                               (as.numeric(read["mate.position"]) >= ucsc.table$txStart) &
                               (as.numeric(read["mate.position"]) <= ucsc.table$txEnd),"name2"]), collapse=",")

           loginfo("found %s for %s %d %d", interval, paste("chr", read["mate.chr"], sep=""), as.numeric(read["mate.position"]), as.numeric(read["mate.position"]))
            if((length(mate.interval) == 0) | (mate.interval == "")) {
                    mate.interval <- "intergenic"
            }
        }
        chrs <- c(read["chr"], read["mate.chr"])
        chrs.order <- order(factor(chrs, levels=c(1:22, "X", "Y")))
        ints <- sapply(c(interval, mate.interval), function(x) sapply(strsplit(x, "_"), "[", 1))
        rearrangement.partners <- paste(paste(chrs[chrs.order], collapse=";"), paste(ints[chrs.order], collapse=";"), sep="/")
        read <- cbind(read, interval, mate.interval, rearrangement.partners, stringsAsFactors=FALSE)
        reads.filtered.genes <- rbind(reads.filtered.genes, read)
    }
    return(reads.filtered.genes)
}

count_reads <- function(reads, ref, start, end, mate.ref=FALSE, mate.start=FALSE, mate.end=FALSE) {
    df.reads <- as.data.frame(reads)
    #Count reads in range
    df.reads.in.range <- (df.reads$refid == ref) & (df.reads$position >= start) & (df.reads$position <= end)

    #Check mate pairs, if coordinates provided
    if(!(FALSE %in% c(mate.ref, mate.start, mate.end))) {
        #Annotate mate pairs
        df.reads <- cbind(df.reads, get_mate_coord(reads))        
        
        #Check that both reads are in range
        df.reads.in.range <- df.reads.in.range & ((df.reads$mate.refid == mate.ref) & (df.reads$mate.position >= mate.start) & (df.reads$mate.position <= mate.end))
    }
    count <- sum(df.reads.in.range)
    return(count)
}

lookup_interval_name <- function(candidate, ucsc.table, prefix="") {
    #Lookup mate.interval in UCSC table
#    interval <- paste(unique(ucsc.table[(ucsc.table$chrom %in% paste("chr", candidate[paste(prefix,"chr",sep="")], sep="")) &
    interval <- paste(unique(ucsc.table[(ucsc.table$chrom %in% candidate[paste(prefix,"chr",sep="")]) &
                           (as.numeric(candidate[paste(prefix,"start",sep="")]) >= ucsc.table$txStart) &
                           (as.numeric(candidate[paste(prefix,"end",sep="")]) <= ucsc.table$txEnd),"name2"]), collapse=",")

    #List "intergenic" if no gene found
    if((length(interval) == 0) | (interval == "")) {
        interval <- "intergenic"
    }

    return(interval)
}

summarize_rearrangements <- function(ucsc.table, candidate.intervals, reads.large_insert, reads.softclipped) {
    summary.all <- c()
    for(i in 1:dim(candidate.intervals)[1]) {
        candidate <- candidate.intervals[i,]

        #Lookup mate.interval in UCSC table
        interval      <- lookup_interval_name(candidate, ucsc.table)
        mate.interval <- lookup_interval_name(candidate, ucsc.table, prefix="mate.")

        #Prepare chromosome and gene/interval strings in order
        chrs <- c(candidate["chr"], candidate["mate.chr"])
        chrs.order <- order(factor(chrs, levels=c(1:22, "X", "Y")))
        ints <- sapply(c(interval, mate.interval), function(x) sapply(strsplit(x, "_"), "[", 1)) #Splits out gene name from interval string description e.g. GENE_Exon_1.
        rearrangement.partners <- paste(paste(chrs[chrs.order], collapse=";"), paste(ints[chrs.order], collapse=";"), sep="/")

        #Count read pairs in each interval
        count.large_insert <- count_reads(reads.large_insert, candidate["ref"], candidate["start"], candidate["end"], candidate["mate.ref"], candidate["mate.start"], candidate["mate.end"])
        count.softclipped.BP1  <- count_reads(reads.softclipped, candidate["ref"], candidate["start"], candidate["end"])
        count.softclipped.BP2  <- count_reads(reads.softclipped, candidate["mate.ref"], candidate["mate.start"], candidate["mate.end"])
        count.total <- count.large_insert + count.softclipped.BP1 + count.softclipped.BP2
        summary.line <- c(rearrangement.partners,
                          interval,
                          candidate["interval_name"],
                          mate.interval,
                          candidate["mate.interval_name"],
                          count.total,
                          count.large_insert,
                          count.softclipped.BP1,
                          count.softclipped.BP2)
        names(summary.line) <- c("Rearrangement_partners",
                                 "Breakpoint_1_Gene",
                                 "Breakpoint_1_Read_range",
                                 "Breakpoint_2_Gene",
                                 "Breakpoint_2_Read_range",
                                 "Total_support",
                                 "Large_insert_read_pairs",
                                 "Breakpoint_1_Soft-clipped_reads",
                                 "Breakpoint_2_Soft-clipped_reads")

        #Add to master tracking matrix
        if(is.null(summary.all)) {
            summary.all <- summary.line
        } else {
            summary.all <- rbind(summary.all, summary.line)
        }
    }
    rownames(summary.all) <- c()

    #Re-order rearrangements from highest to lowest confidence
    summary.ordered <- summary.all[order(decreasing=TRUE,
                                      as.numeric(summary.all[,"Total_support"]),
                                      as.numeric(summary.all[,"Large_insert_read_pairs"]),
                                      as.numeric(summary.all[,"Breakpoint_1_Soft-clipped_reads"]),
                                      as.numeric(summary.all[,"Breakpoint_2_Soft-clipped_reads"]),
                                      as.character(summary.all[,"Rearrangement_partners"])
                                    ),]

    return(summary.ordered)
}

convert_refids_to_chr <- function(reads, refs) {
    reads$chr      <- refs$SN[match(reads$refid, refs$ID)]
    reads$mate.chr <- refs$SN[match(reads$mate.refid, refs$ID)]
    reads.trimmed <- reads[,c("chr", "position", "mate.chr", "mate.position")]
    return(reads.trimmed)
}

format_coords_as_interval_list <- function(candidate.intervals) {
    colnames(candidate.intervals) <- c("ref", "start", "end")
    candidate.intervals <- cbind(candidate.intervals, strand=rep("+", dim(candidate.intervals)[1]))
    candidate.intervals <- cbind(candidate.intervals, chr=refs$SN[match(candidate.intervals[,"ref"], refs$ID)])
    candidate.intervals <- cbind(candidate.intervals, interval_name=
                                     paste(candidate.intervals[,"chr"], ":",
                                           candidate.intervals[,"start"], "-",
                                           candidate.intervals[,"end"], sep=""))
    rownames(candidate.intervals) <- candidate.intervals[,"interval_name"]
    
    #Reorder columns to be consistent with other interval_list files
    candidate.interval_list <- candidate.intervals[,c("chr", "start","end", "strand", "interval_name", "ref"), drop=FALSE]

    return(candidate.interval_list)
}

make_intervals_from_reads <- function(reads.large_insert, max.cluster.distance) {
    require("intervals")
    candidate.coords <- c()
    df.reads <- as.data.frame(reads.large_insert)
    df.reads$end <- df.reads$position + nchar(df.reads$seq)
    for(ref in unique(df.reads$refid)) { #Query df.reads on same chromosome
        df.reads.ref <- df.reads[df.reads$refid == ref,]
        df.reads.intervals <- matrix(c(df.reads.ref$position - max.cluster.distance,
                                       df.reads.ref$end + max.cluster.distance), ncol=2)
        df.reads.unions <- unlist(interval_union(Intervals(df.reads.intervals)))
        df.reads.unions <- matrix(data=df.reads.unions, ncol=2)
        ref.candidate.coords <- matrix(data=c(rep(ref, dim(df.reads.unions)[1]),
                                                 df.reads.unions[,1],
                                                 df.reads.unions[,2]),
                                           ncol=3)
        #Merge unions into master tracking matrix
        if(is.null(candidate.coords)) {
            candidate.coords <- ref.candidate.coords
        } else {
            candidate.coords <- rbind(candidate.coords, ref.candidate.coords)
        }
    }

    #Format intervals as interval_list
    candidate.interval_list <- format_coords_as_interval_list(candidate.coords)

    return(candidate.interval_list)
}

map_unmapped_reads <- function(reads.softclipped, reader, candidate.intervals, remove.duplicates) {
    #TODO: This function is under construction
    #Read alignments one-by-one to identify read pairs that may support translocations
    
    #Extract all unmapped reads
    #/usr/local/samtools-0.1.18-default/samtools view -b -h -f 4 31468-3_L_000286_BC23_HC-000026_49.dedup.cleaned.bam > ~/31468-3_L_000286_BC23_HC-000026_49.dedup.cleaned.unmapped.bam
    reads.unmapped <- bamRange()
    rewind(reader)
    align <- getNextAlign(reader)
    while(!is.null(align)) {
        if(!((remove.duplicates=="TRUE") & (pcrORopt_duplicate(align)=="TRUE"))) { #Read only non-duplicate reads
            if(unmapped(align)=="TRUE") {
                push_back(reads.unmapped, align)
            }
        }
        align <- getNextAlign(reader)
    }

    #Filter unmapped reads by base quality

    #Align unmapped reads to breakpoint reference sequences
    candidate.intervals.with.count.of.new.maps <- c()
    for(i in seq(1:dim(candidate.intervals)[1])) {
        reads.softclipped.interval  <- isolate_reads_in_interval(reads.softclipped)
        breakpoint.refs             <- make_fasta_from_overlapping_reads(reads.softclipped, interval_list)
        reads.unmapped.aligned.bp   <- bowtie_align_reads_to_fasta(reads.softclipped.interval, breakpopint.refs)
        reads.unmapped.aligned.hg19 <- bowtie_align_reads_to_fasta(reads.unmapped.aligned.bp, hg19.fasta)
        candidate.intervals.with.count.of.new.maps <- append(candidate.intervals[i,], rbamtools::size(reads.unmapped.aligned.hg19))
    }    
    return(candidate.intervals.with.remapped.counts)
}

make_bam_from_reads <- function(out.dir, out.bam, bam.file, ...) {
    #Trailing arguments can be any number of bamRange objects

    #Need to change working directory as bamWriter doesn't support full paths
    #wd.orig <- getwd()
    #setwd(out.dir)  

    #write new bam file, using header from original bam file
    out.bam.unsorted <- paste(out.dir,sub(".bam$", ".unsorted.bam", out.bam), sep="/")
    orig.bam <- BamFile(bam.file)
    reads.list <- list(...)
    # want contains the names of the reads that we want to keep
    want <- c()
    for(reads in reads.list) {
    	# write.table(as.data.frame(reads), file="reads2.txt", sep="\t", append=T)
	want <- c(want,as.data.frame(reads)$name)
    }
    #loginfo("writing want to text file")
    #write.table(want, file=paste(out.dir,"want.txt",sep="/"), sep="\t", append=F)

    # methods to filter the original bam file to get the desired reads in a new bam
    filter_factory <- function(want) {
    		   list(KeepQname = function(x) x$qname %in% want)
		   }
    filter <- FilterRules(filter_factory(want))
    loginfo("trying filterBam")
    dest <- filterBam(orig.bam, out.bam.unsorted, filter=filter,
                  param=ScanBamParam(what="qname"))
    readGAlignments(dest)

    # sort and index the resulting file
    loginfo("trying sort")
    # oddly, sortBam adds '.bam' to the filename, while the other tools don't 
    sortBam(out.bam.unsorted,paste(out.dir,sub(".bam$", "", basename(out.bam)),sep="/"))
    loginfo("trying index")
    indexBam(paste(out.dir,out.bam,sep="/"))
    loginfo("finished making bam")

    #Delete unsorted bam file
    #file.remove(out.bam.unsorted)
    #setwd(wd.orig)
}

convert_bamRange_to_annotated_data_frame <- function(reads.bamRange) {
    df.reads <- as.data.frame(reads.bamRange)
    df.reads <- cbind(df.reads, get_mate_coord(reads.bamRange))

    #Convert refIDs to chromosomes
    df.reads$chr <- refs$SN[match(df.reads$refid, refs$ID)]
    df.reads$mate.chr <- refs$SN[match(df.reads$mate.refid, refs$ID)]

    #Add translocation partner information
    df.reads.annotated <- add_genes_to_reads(df.reads, interval_list, ucsc.table)

    return(df.reads.annotated)
}

get_mate_pair_reads <- function(reads, reader) {
    reads.mates <- bamRange()
    loginfo("rewind")
    rewind(reader)
    loginfo("align")
    align <- getNextAlign(reads)

    #Read alignments one-by-one to identify read pair
    loginfo("start while")
    count = 0
    while(!is.null(align)) {
        align.mate <- get_mate_read(align, reader)
        if(class(align.mate)[1] == "bamAlign") {    
            push_back(reads.mates, align.mate)
        }
        align <- getNextAlign(reads)
    count = count + 1
    }
    loginfo("while loop had %d iterations", count)
    return(reads.mates)
}

join_bamRanges <- function(...) {
    joined.reads <- bamRange()
    reads.list <- list(...)
    for(reads in reads.list) {
        align <- getNextAlign(reads)
        while(!is.null(align)) {
            push_back(joined.reads, align)
            align <- getNextAlign(reads)
        }
    }
    return(joined.reads)
}

add_mate_intervals <- function(candidate.intervals.discovery, reads.large_insert.mates, max.cluster.distance, max.mate.partners) {
    candidate.mate.coords <- c() #master tracking variable
    
    #Convert reads to data.frame for parsing, add mate-pair information    
    df.reads          <- as.data.frame(reads.large_insert.mates)
    df.reads$end      <- df.reads$position + nchar(df.reads$seq)
    df.reads          <- cbind(df.reads, get_mate_coord(reads.large_insert.mates))
    df.reads$mate.end <- df.reads$mate.position + nchar(df.reads$seq)

    loginfo("start for loop")
    for(i in 1:dim(candidate.intervals.discovery)[1]) {
    	loginfo("outer for loop iteration number %d", i)
        interval <- candidate.intervals.discovery[i,]
        df.reads.interval <- df.reads[(df.reads$mate.refid == interval["ref"]) & 
                                      (df.reads$mate.position >= interval["start"]) &
                                      (df.reads$mate.position <= interval["end"]),]

        #Define left- and right-most boundaries
        for(ref in unique(df.reads.interval$refid)) { #Query df.reads on same chromosome
	    loginfo("inner for loop, ref is %s", ref)
            df.reads.ref <- df.reads.interval[df.reads.interval$refid == ref,]
            df.reads.ranges <- matrix(c(df.reads.ref$position - max.cluster.distance,
                                           df.reads.ref$end + max.cluster.distance), ncol=2)
            df.reads.unions <- unlist(interval_union(Intervals(df.reads.ranges)))
            df.reads.unions <- matrix(data=df.reads.unions, ncol=2)
            ref.candidate.coords <- matrix(data=c(rep(ref, dim(df.reads.unions)[1]),
                                                     df.reads.unions[,1],
                                                     df.reads.unions[,2]),
                                               ncol=3)

            #Merge unions into master tracking matrix
            mate.count <- dim(ref.candidate.coords)[1]
            intervals.to.join <-  matrix(data=rep(interval, mate.count), nrow=mate.count, byrow=TRUE, dimnames=list(c(),names(interval)))
            ref.candidates <- cbind(intervals.to.join, format_coords_as_interval_list(ref.candidate.coords))
            colnames(ref.candidates) <- c(colnames(intervals.to.join), paste("mate",colnames(intervals.to.join), sep="."))
            rownames(ref.candidates) <- c()
            if(dim(ref.candidates)[1] <= max.mate.partners) {
                if(is.null(candidate.mate.coords)) {
                    candidate.mate.coords <- ref.candidates
                } else {
                    candidate.mate.coords <- rbind(candidate.mate.coords, ref.candidates)
                }
            }
        }
    }
    return(candidate.mate.coords)
}

make_bed_from_rearrangement_counts <- function(rearrangement.counts) {
    #Parse chr, start, end into separate columns for use by IGV and UCSC genome browser
    coords <- c(rbind(rearrangement.counts[,"Breakpoint_2_Read_range"],
                      rearrangement.counts[,"Breakpoint_1_Read_range"]))
    bed.coord <- matrix(data=unlist(sapply(coords, strsplit, ":|-")), ncol=3, byrow=TRUE)

    #Parse out descriptions and use to annotate coordinates
    descr <- c(rbind(rearrangement.counts[,"Rearrangement_partners"],
                     rearrangement.counts[,"Rearrangement_partners"]))
    bed <- cbind(bed.coord, descr)

    #Prefix chromosome names with "chr"
    bed[,1] <- paste("chr", bed[,1], sep="")
    return(bed)
}

######
#Main#
######


#if (FALSE) {   # start commented out code
#Read input files
loginfo("reading ucsc table and interval list")
ucsc.table <- read.table(ucsc.table.file, header=TRUE, comment.char="")
interval_list <- load_interval_names(interval_list.file)
#interval_list <- interval_list[grep("BCL2", interval_list[,"interval_name"]),]  #For development only. Restrict to BCL2.
#TODO: Add IGH to UCSC file
loginfo("done reading")

#Read bam and index files, extract reads with mates in interval.pair
write(paste("Reading bam file", bam.file, sep=" "), stdout())
idx.file <- sub(".bam$", ".bai", bam.file) #bam index file
reader<-bamReader(bam.file, idx.file)
if(!isOpen(reader)) {
    stop(paste("ERROR: Cannot open bam file", bam.file, sep=" "))
} 
#load.index(reader, idx.file)
#if(!index.initialized(reader)) {
#    warning(paste("ERROR: Cannot initialize index", idx.file, sep=" "))
#    next
#}

#Parse reference names and provided intervals
loginfo("getRefData")
refs <- getRefData(reader)
interval_list <- cbind(interval_list, ref=lookup_reference_number(refs, interval_list))

#Isolate read pairs with intrachromosomal mates or large insert size
loginfo("reads.candidates")
reads.candidates <- get_candidate_translocation_reads(reader, interval_list, min.insert.size, min.interval.length, remove.duplicates, min.mapping.quality)

#Filter candidates that are distant from other reads and add in mate pairs
loginfo("remove singleton outliers")
reads.large_insert <- remove_singleton_outliers(reads.candidates, max.cluster.distance, min.cluster.size)
loginfo("get mate pair reads")
reads.large_insert.mates <- get_mate_pair_reads(reads.large_insert, reader)

#} # end commented out code



loginfo("optional save point")
#Optional save/load point for development
#save.image(paste(out.dir, "temp.rData", sep="/"))
#load(paste(out.dir, "temp.rData", sep="/"))

loginfo("create output dir")
#Create output directory, before writing bam file
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

#Generate intervals for breakpoint regions and find soft-clipped reads
if(dim(as.data.frame(reads.large_insert))[1] > 0) {
    loginfo("make intervals from reads")
    candidate.intervals.discovery <- make_intervals_from_reads(reads.large_insert, max.cluster.distance)
    loginfo("add mate intervals")
    candidate.intervals <- add_mate_intervals(candidate.intervals.discovery, reads.large_insert.mates, max.cluster.distance, max.mate.partners)

    #Get soft-clipped reads near the intrachromsomal or large insert size mates
    loginfo("get softclipped reads")
    reads.softclipped <- get_softclipped_reads(candidate.intervals, reader, min.soft.clipped.bases, min.base.quality, remove.duplicates)

    #Add count of reads supporting each interval
    loginfo("summarize rearrangements")
    rearrangement.counts  <- summarize_rearrangements(ucsc.table, candidate.intervals, reads.large_insert, reads.softclipped)
        loginfo("make bed from rearrangement counts")
	rearrangement.regions <- make_bed_from_rearrangement_counts(rearrangement.counts)

	#Write out single bam file for manual review and close source bam file
	loginfo("write out single bamfile")
	out.bam <- sub(".bam$", ".rearrangements.bam", basename(bam.file))
	loginfo("out.bam is %s", out.bam)
	#write.table(as.data.frame(reads.large_insert),file="reads1.txt",sep="\t")
	make_bam_from_reads(out.dir, out.bam, bam.file, reads.large_insert, reads.large_insert.mates, reads.softclipped)
	bamClose(reader)

} else {
    loginfo("no rearrangements found")
    rearrangement.counts <- c("No rearrangement candidates found.")
	rearrangement.regions <- ""
}

loginfo("write out to files")
#Write out to files
base.out.file <- paste(out.dir, sub(".bam$", "", basename(bam.file)), sep="/")
write.table(rearrangement.counts, paste(base.out.file, ".rearrangement_counts.tsv", sep=""), row.names=FALSE, quote=FALSE, sep="\t")

#Write out list of regions for import to IGV for rapid review
write.table(rearrangement.regions, paste(base.out.file, ".rearrangement_regions.bed", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

#Write log file
run_info_table <- matrix(ncol = 2, byrow=TRUE, data = c(
        "Date",                                         date(),
        "CluMP command",                                paste(commandArgs(), collapse=" "),
        "Subversion revision information",              svn.revision,
        "Config file",                                  config.file,      
        "Bam file",                                     bam.file,
        "Output directory",                             out.dir,              
        "Interval file",                                interval_list.file,
        "UCSC gene table file",                         ucsc.table.file,
        "Remove read duplicates",                       remove.duplicates,
        "Minimum interval length",                      min.interval.length,
        "Minimum read mapping quality",                 min.mapping.quality,
        "Minimum base quality across entire read",      min.base.quality,
        "Minimum number of soft-clipped bases in read", min.soft.clipped.bases,
        "Minimum read pair insert size",                min.insert.size,
        "Maximum bp between reads to form cluster",     max.cluster.distance,
        "Maximum mate regions before discarding",       max.mate.partners,
        "Minimum number of reads in cluster",           min.cluster.size
        ))
run_info_outfile <- paste(base.out.file, ".CluMP_run_info.tsv", sep="")
write.table(run_info_table, run_info_outfile, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

loginfo("run complete")