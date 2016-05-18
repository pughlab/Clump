### Simple filter of CluMP results

# 1. Hard filter of number of Large insert reads
# 2. Filter known false positive partners
# 3. Output snapshot script for each site in the filtered pairs

## User input
CluMP_data_file = "/mnt/work1/users/pughlab/projects/Meningioma_targeted/GeneFusion/output/rearrangement_count_files/ALL_rearrangements.txt"
out_dir = "/mnt/work1/users/pughlab/projects/Meningioma_targeted/GeneFusion/output/"
min_large_insert_pairs = 10

### #### ###
### MAIN ###
### #### ###

clump_results <- read.table(CluMP_data_file,
           sep="\t",header=TRUE)

## Parse sample names

partner_sample <- unlist(strsplit(as.character(clump_results$Rearrangement_partners),":"))

clump_results$Rearrangement_partners <- partner_sample[seq(2,length(partner_sample),by=2)]
clump_results$Sample <- gsub(".processed.rearrangement_counts.tsv","",
                                partner_sample[seq(1,length(partner_sample),by=2)])

## Filter CluMP results

clump_results_filtered <- subset(clump_results,
                                    !(grepl("LINC00486",Rearrangement_partners)) &
                                    !(grepl("ZNF708",Rearrangement_partners)) &
                                    Large_insert_read_pairs>min_large_insert_pairs)


## output IGV_snapshoter commands
snap_Breakpoint_2 <- as.data.frame(paste("/mnt/work1/users/pughlab/src/IGVSnapshot/snap ",
                                         clump_results_filtered$Sample,".processed.rearrangements.bam ",
                                         clump_results_filtered$Breakpoint_2_Read_range,sep=""))
      

snap_Breakpoint_1 <- paste("/mnt/work1/users/pughlab/src/IGVSnapshot/snap ",
                           clump_results_filtered$Sample,".processed.rearrangements.bam ",
                           clump_results_filtered$Breakpoint_1_Read_range,sep="")

## output filtered breakpoint list
bas <- basename(CluMP_data_file)
write.table(clump_results_filtered,file=paste(out_dir,"/",bas,".filtered",sep=""),
            quote=FALSE,row.names=FALSE,col.names=TRUE)


## output IGVSnapshot scripts (to be run on samwise) for each breakpoint region in the filtered list

write.table(snap_Breakpoint_1,file=paste(out_dir,"/snap_Breakpoints_1.sh",sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(snap_Breakpoint_2,paste(out_dir,"/snap_Breakpoints_2.sh"),
            quote=FALSE,row.names=FALSE,col.names=FALSE)