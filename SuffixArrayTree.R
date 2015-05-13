

# Art code
# ./art_454 chr1.fa chr1_reads 2
#./art_illumina -i chr1.fa -l 50 -f 5 -o chr1_reads_illumina

# Suffix Array based alignment algorithm created by S. Karthik Badam

library(ShortRead) 
library(rbenchmark)
library(doParallel)
registerDoParallel(cores = 4)
library(foreach)
library(gtools)

# Reference genome
require(BSgenome.Hsapiens.UCSC.hg19)
seqChr = unmasked(Hsapiens$chr1)[100000:110000]
seqLength <- length(seqChr)

# Size of the suffix array loaded into memory
sampleSize <- 1000
runs <- round(seqLength/sampleSize)

# Number of reads to be aligned
matchSize <- 100

# reads <- readFastq('~/Documents/CMSC702/Project/chr1_reads.fq')
# Reads from the ART simulator
reads <- readFastq('chr1_reads.fq')
rs <- sort (sread(reads), decreasing=TRUE)

# Looping through the entire suffix array --currently lets just do 1 segment of suffix array
# for (i in 1:runs) {
for (i in 1:1) {
  
  # Load the part of the suffix array for alignment
  if (i != runs) {
    sa <- DNAStringSet(Views(seqChr, start=seq(i*sampleSize, (i+1)*sampleSize-1), length(seqChr)))
    
  } else {
    sa <- DNAStringSet(Views(seqChr, start=seq(i*sampleSize, length(seqChr)), length(seqChr)))
    
  }
  
  # benchmark code -- comparing linear search, parallel linear search, and parallel binary search
  print(benchmark(replications = 1, order = "elapsed",
                  non_parallel = {
                    for (j in 1:matchSize) {
                      foreach(ref = sa) %do% compareStrings(rs[[j]], substr(ref, 1, length(rs[[j]]))) 
                    }
                  }
                  , linear_parallel = {
                    # Matching parallel
                    test <- foreach(ref = sa, .combine = cbind) %dopar%  for (j in 1:matchSize) {
                      compareStrings(rs[[j]], substr(ref, 1, length(rs[[j]])))
                    }
                  }
                  , binary_parallel = {
                   
                    sa <- sort(sa)
                    
                    foreach(search = rs[1:10], .combine = cbind) %dopar%  {  
                      
                      # Binary search
                      range <- c(1, length(sa))
                      index <- 1
                      while (range[2] > range[1]) {
                        index <- as.integer(sum(range)/2)
                        boolArray <- as.integer(search) < as.integer(substr(sa[[index]], 1, length(search)))
                        if (!is.na(sum(boolArray))) {
                          boolArray <- boolArray*1
                          if (median(boolArray) == 1) {
                            range[2] = index - 1
                          } else if (median(boolArray) == 0) {
                            range[1] = index + 1
                          } else {
                            break
                          }
                        } else {
                          range[1] = index - 1
                        }
                      }
                      compareStrings(search, substr(sa[[index]], 1, length(search)))
                    }
                  }
  ))
}

# Suffix Tree created by S. Karthik Badam

# Load fq file generated from art program
library("ShortRead")
library("PST")

reads <- readFastq('chr1_reads.fq')
rs <- sread(reads)

# Creating a suffix tree for a test string
testString <- toString(seqChr[100000:101000])
testString <- gsub("([A-Z])", "-\\1", testString)
testString <- substr(testString, 2, nchar(testString))
seq1 <- seqdef(testString)
tree <- pstree(seq1)
