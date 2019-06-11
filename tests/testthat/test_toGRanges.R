library(regioneR)
context("toGRanges")

#Define some GRanges to use in the tests
starts <- c(1, 15, 24)
ends <- c(10, 20, 30)
chrs <- c("1", "1", "2")
A <- data.frame(chr=chrs, start=starts, end=ends,  x=c(1,2,3), y=c("a", "b", "c"))
gr1 <- toGRanges(A)



#from data.frame
test_that(" toGRanges works for a data.frame", {
  
  expect_is(gr1, "GRanges")
  expect_length(gr1, 3)
  expect_equal(names(mcols(gr1)), c("x", "y"))
  expect_is(gr1$x, "numeric")
  expect_is(gr1$y, "factor")  
  expect_equal(GenomicRanges::start(gr1), starts)
  expect_equal(GenomicRanges::end(gr1), ends)
  expect_equal(as.character(seqnames(gr1)), chrs)
  
  gr1.2 <- toGRanges(data.frame(1, starts, ends))
  expect_equal(as.character(seqnames(gr1.2)), c("1", "1", "1"))

  expect_equal(toGRanges(1, starts, ends), gr1.2)
    
  expect_equal(toGRanges(toDataframe(gr1)), gr1)
  
  
  #check the genome
  expect_length(genome(gr1), 2)
  expect_length(genome(gr1.2), 1)

  gr1.3 <- toGRanges(A, genome = "hg19")
  expect_equal(seqlengths(gr1.3),   seqlengths(getGenome("hg19")))
  expect_equal(as.character(seqnames(gr1.3)), c("chr1", "chr1", "chr2"))
  
})


     



#' A <- data.frame(chr=1, start=c(1, 15, 24), end=c(10, 20, 30),  x=c(1,2,3), y=c("a", "b", "c"))
#' gr1 <- toGRanges(A)
#' 
#' #No need to give the data.frame columns any specific name
#' A <- data.frame(1, c(1, 15, 24), c(10, 20, 30),  x=c(1,2,3), y=c("a", "b", "c"))
#' gr2 <- toGRanges(A)
#' 
#' #We can pass the data without building the data.frame
#' gr3 <- toGRanges("chr9", 34229289, 34982376, x="X")
#' 
#' #And each argument can be a vector (they will be recycled as needed)
#' gr4 <- toGRanges("chr9", c(34229289, 40000000), c(34982376, 50000000), x="X", y=c("a", "b"))
#' 
#' #toGRanges will automatically convert the second and third argument into numerics
#' gr5 <- toGRanges("chr9", "34229289", "34982376") 
#' 
#' #It can be a file from disk
#' bed.file <- system.file("extdata", "my.special.genes.bed", package="regioneR")
#' gr6 <- toGRanges(bed.file)
#' 
#' #Or a URL to a valid file
#' gr7 <- toGRanges("http://molb7621.github.io/workshop/_downloads/lamina.bed")
#' 
#' #It can also parse genomic location strings
#' gr8 <- toGRanges("chr9:34229289-34982376")
#' 
#' #more than one
#' gr9 <- toGRanges(c("chr9:34229289-34982376", "chr10:1000-2000"))
#' 
#' #even with mixed strange and mixed syntaxes
#' gr10 <- toGRanges(c("chr4:3873-92928", "chr4:3873,92928", "chr5:33,444-45,555"))
#' 
#' #if the genome is given it is used to annotate the resulting GRanges
#' gr11 <- toGRanges(c("chr9:34229289-34982376", "chr10:1000-2000"), genome="hg19")
#' 
#' 
#' #and the genome is added to the GRanges even if A is a GRanges
#' gr12 <- toGRanges(gr6, genome="hg19")
#' 
#' #And it will change the chromosome naming of the GRanges to match that of the
#' #genome if it is possible (using GenomeInfoDb::seqlevelsStyle)
#' gr2
#' gr13 <- toGRanges(gr2, genome="hg19")
#' 
#' #in addition, it can convert other objects into GRanges such as the 
#' #result of GenomicRanges::coverage
#' 
#' gr14 <- toGRanges(c("1:1-20", "1:5-25", "1:18-40"))
#' cover <- GenomicRanges::coverage(gr14)
#' gr15 <- toGRanges(cover)

