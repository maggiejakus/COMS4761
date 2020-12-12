#reads in clumping data so we only get SNPs that are (probably not) in LD 
#reads in alkes group annotations. alkes group in genome build GrCh38, while everything else in GrCh37.
#reads in conversion file (converting SNPs from GrCh38 to GrCh37) + converts alkes files
#check sizes/for duplicates + sort by order often
all_combos <- list()

for (i in 1:22){
  print(i)

  alkes <- read.table(paste0("~/Documents/CGNMX/baseline_v1.2/baseline.",i, ".annot.gz"), stringsAsFactors=F, header=T)
  convert_pos <- read.table(paste0("~/Documents/CGNMX/conversions/conv", i, "_38to37.txt"), stringsAsFactors=F, header=T, fill = T)
  
  convert_pos <- convert_pos[convert_pos[,2] != "" & !is.na(convert_pos[,1]) & !is.na(convert_pos[,2]),]
  convert_pos <- convert_pos[!(convert_pos$source_start %in% convert_pos$source_start[duplicated(convert_pos$source_start)]),]
  #convert_pos <- convert_pos[!is.na(convert_pos[,2]),]
  alkes <- alkes[alkes$BP %in% convert_pos[,1],]

  #add GrCh37 BP location to alkes data
  alkes <- alkes[order(alkes$BP),]
  convert_pos <- convert_pos[order(convert_pos[,1]),]
  alkes$hg19_bp <- convert_pos[,2]


  #num should be the bp pos, indexing at 1
  bp_pos = 2 
  ss <- read.table(paste0("keep.ss.", i), stringsAsFactors=F, header=F)
  ss <- ss[ss[,bp_pos] %in% alkes$hg19_bp,]
  ss <- ss[!(ss[,bp_pos] %in% ss[,bp_pos][duplicated(ss[,bp_pos])]),]
  #all the 3s were 2 before

  sub_alkes <- alkes[alkes$hg19_bp %in% ss[,bp_pos],]
  sub_alkes <- sub_alkes[order(sub_alkes$hg19_bp),]
  ss <- ss[order(ss[,bp_pos]),]

  #combine clumped SNPs with their annotations, all in GrCh37/hg19
  all_combos[[i]] <- cbind(ss, sub_alkes)
}

combo <- do.call("rbind", all_combos)
write.table(combo, "converted_annot.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ' ')
