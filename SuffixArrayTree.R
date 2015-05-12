# ./art_454 chr1.fa chr1_reads 2
#./art_illumina -i chr1.fa -l 50 -f 5 -o chr1_reads_illumina
require(BSgenome.Hsapiens.UCSC.hg19)
seqChr = unmasked(Hsapiens$chr1)
#sa <- sort(DNAStringSet(Views(seqChr1, start=seq_len(nchar(seqChr1)), length(seqChr1))))

seqLength <- length(seqChr)
sampleSize <- 100
runs <- round(seqLength/sampleSize)
matchSize <- 100

library("ShortRead") 
library(rbenchmark)
library(doParallel)
registerDoParallel(cores = 4)
library(foreach)

# reads <- readFastq('chr1_reads.fq')
reads <- readFastq('~/Documents/CMSC702/Project/chr1_reads.fq')
rs <- sread(reads)

# for (i in 1:runs) {
for (i in 1:1) {
  if (i != runs) {
    sa <- DNAStringSet(Views(seqChr, start=seq(i*sampleSize, (i+1)*sampleSize-1), length(seqChr)))
    
  } else {
    sa <- DNAStringSet(Views(seqChr, start=seq(i*sampleSize, length(seqChr)), length(seqChr)))
    
  }
  
  print(benchmark(replications = 1, order = "elapsed",
                  non_parallel = {
                    for (j in 1:matchSize) {
                      foreach(ref = sa) %do% pmatch(rs[[j]], substr(ref, 1, length(rs[[j]]))) 
                    }
                  }
                  , parallel = {
                    test <- foreach(ref = sa, .combine = cbind) %dopar%  for (j in 1:matchSize) {pmatch(rs[[j]], substr(ref, 1, length(rs[[j]])))}
                    # pairwiseAlignment(rs[[j]], substr(ref, 1, length(rs[[j]])))
                  }
  ))
}

for (i = 1; i <= runs; i++) {
  if (i != runs) {
    sa <- sort(DNAStringSet(Views(seqChr, start=seq(i*sampleSize, (i+1)*sampleSize), length(seqChr))))
    
  } else {
    sa <- sort(DNAStringSet(Views(seqChr, start=seq(i*sampleSize, length(seqChr)), length(seqChr))))
    
  }
  # do binary search
  benchmark(replications = 1, order = "elapsed",
            non_parallel = {
              test1 <- for (j = 1 ; j <= length(rs); j++) {
                for (k = 1; k <= length(sa); k++) {
                  referencePrefix <- substr(sa[[k]], 1, length(rs[[j]]))
                  pairwiseAlignment(rs[[j]], referencePrefix)
                }
              }
            },
            
            parallel = {
              test2 <- foreach(i = sa, .combine = cbind) %dopar% pairwiseAlignment(rs[[j]], substr(i, 1, length(rs[[j]])))
            }  
  )
}


# sa <- sort(DNAStringSet(Views(seqChr1, start=seq(1, 10000000), length(seqChr1))))

# Load fq file generated from art program

library("ShortRead")
# reads <- readFastq('chr1_reads.fq')
reads <- readFastq('~/Documents/CMSC702/Project/chr1_reads.fq')
rs <- sread(reads)
# pairwiseAlignment(rs[[1]], rs[[1]])

# library("ape")
# reads <- read.dna("~/Documents/CMSC702/Project/chr1_reads_illumina.aln", format = "clustal")
# sequences <- sread(reads)

#pairWiseAlignment

# library("Rsamtools")

library("PST")
testString <- toString(seqChr1[1:100000])
testString <- gsub("([A-Z])", "-\\1", testString)
testString <- substr(testString, 2, nchar(testString))
seq1 <- seqdef(testString)
tree <- pstree(seq1)
