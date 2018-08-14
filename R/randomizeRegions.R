#' Randomize Regions
#' 
#' @description 
#' Given a set of regions A and a genome, this function returns a new set of regions randomly distributted in the genome. 
#' 
#' @details
#' The new set of regions will be created with the same sizes of the original ones, and optionally placed in the same chromosomes.
#' 
#' In addition, they can be made explicitly non overlapping and a mask can be provided so no regions fall in an undesirable part of the genome.
#' 
#' 
#' @usage 
#' randomizeRegions(A, genome="hg19", mask=NULL, allow.overlaps=TRUE, per.chromosome=FALSE, ...)
#' 
#' @param A The set of regions to randomize. A region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param genome The reference genome to use. A valid genome object. Either a \code{\link{GenomicRanges}} or \code{\link{data.frame}} containing one region per whole chromosome or a character uniquely identifying a genome in \code{\link{BSgenome}} (e.g. "hg19", "mm10",... but not "hg"). Internally it uses \code{\link{getGenomeAndMask}}.
#' @param mask The set of regions specifying where a random region can not be (centromeres, repetitive regions, unmappable regions...). A region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}},\code{\link{data.frame}}, ...). If \code{\link{NULL}} it will try to derive a mask from the genome (currently only works if the genome is a character string). If \code{\link{NA}} it gives, explicitly, an empty mask.
#' @param allow.overlaps A boolean stating whether the random regions can overlap (FALSE) or not (TRUE).
#' @param per.chromosome Boolean. If TRUE, the regions will be created in a per chromosome maner -every region in A will be moved into a random position at the same chromosome where it was originally-.
#' @param ... further arguments to be passed to or from methods.
#' 
#' 
#' @return
#' It returns a \code{\link{GenomicRanges}} object with the regions resulting from the randomization process.
#' 
#' @seealso  \code{\link{toDataframe}}, \code{\link{toGRanges}}, \code{\link{getGenome}}, \code{\link{getMask}}, \code{\link{getGenomeAndMask}}, \code{\link{characterToBSGenome}}, \code{\link{maskFromBSGenome}}, \code{\link{resampleRegions}}, \code{\link{createRandomRegions}}, \code{\link{circularRandomizeRegions}}
#' 
#' @examples
#' A <- data.frame("chr1", c(1, 10, 20, 30), c(12, 13, 28, 40))
#' 
#' mask <- data.frame("chr1", c(20000000, 100000000), c(22000000, 130000000))
#' 
#' genome <- data.frame(c("chr1", "chr2"), c(1, 1), c(180000000, 20000000))
#' 
#' randomizeRegions(A)
#' 
#' randomizeRegions(A, genome=genome, mask=mask, per.chromosome=TRUE, allow.overlaps=FALSE)
#' 
#' @export randomizeRegions
#' 
#' @importFrom IRanges IRanges




randomizeRegions <- function(A, genome="hg19",  mask=NULL, allow.overlaps=TRUE, per.chromosome=FALSE, ...) { 
  
  if(!hasArg(A)) stop("A is missing")
  if(!is.logical(allow.overlaps)) stop("allow.overlaps must be logical")
  if(!is.logical(per.chromosome)) stop("per.chromosome must be logical")
  
  A <- toGRanges(A)
  
  if(any(end(A) < start(A))) {
    stop("There are regions with negative width. Start must always be smaller than end.")
  }
  
  #The randomization of an empty region set, is an empty region set
  if(length(A)==0) {
    return(A)
  }
  
  #Prepare the genome and mask
  gam <- getGenomeAndMask(genome=genome, mask=mask)
  genome <- gam$genome
  mask <- gam$mask
  
  
  #If per.chromosome, split by chromosomes and recall ourselves recursively
  if(per.chromosome == TRUE) { 
    
    #Check if there's any chromosome in A not in the genome
    missing.chrs <- setdiff(levels(seqnames(A)), levels(seqnames(genome)))
    if(length(missing.chrs > 0)) {
      stop(paste("The chromosomes **", paste(missing.chrs, sep=", "), "** of A are not present in the genome. It is not possible to use \"per.chromosome\" with this dataset.", sep=""))
    }
    
    #Split the regions in A per chromosomes and create the random regions separately per each one
    random.regions <- copySeqLevels(to=GRanges(), from=genome)    
    for(chr in seqlevels(A)) {  #TODO: Should we change it to an apply and possibly parallelize?
      chr.A <- A[seqnames(A)==chr]
      if(length(chr.A)>0) {
        chr.genome <- genome[seqnames(genome)==chr]
        chr.mask <- mask[seqnames(mask)==chr]
        new.regions <- randomizeRegions(chr.A, chr.genome, chr.mask, per.chromosome=FALSE, allow.overlaps=allow.overlaps)
        random.regions <- c(random.regions, new.regions)
      }
    }
    return(random.regions)
  }
  
  #Else, If we are here, per.chromosome==FALSE. Place all the regions
  
  #We have three different algorithms, 
  #    1 - a fast but not guaranteed to completely suceed
  #    2 - a linear one that does not account for overlaps between regions
  #    3 - a quadratic, that produces non-overlapping regions
  
  # 1 - FAST ALGORITHM 
  #To place the regions, start with the fastest one, and then, if necessary, finish the work
  # with the slower ones
  
  
  A <- copySeqLevels(to=A, from=genome)
  
  ntimes <- 0
  
  random.regions <- A[0]
  pending <- A
  #TODO: Adjust these parameters
  while(ntimes < 10) {
    orig.pending.len <- length(pending)
    #print(paste0("pending regions: ", orig.pending.len))
    
    rand.reg.list <- randomizeFastWithNoOverlappingControl(pending, genome=genome, mask=mask, ...)
    
    tt <- do.call(rbind, rand.reg.list[!is.na(rand.reg.list)])
    random.regs <- GRanges(seqnames=as.character(tt[,1]), ranges=IRanges::IRanges(start=as.numeric(tt[,2]), end=as.numeric(tt[,3])))
    random.regs <- copySeqLevels(to=random.regs, from=genome)
    pend <- pending[is.na(rand.reg.list)] #if there is any NA in rand.reg.list, these are pending regions for which we have not a valid random place
    
    
    if(allow.overlaps==FALSE) { #In addition, if overlapping regions are not allowed, any overlapping region needs to be replaced
      rr <- removeOverlapping(random.regs)
      random.regs <- rr$rs
      pend <- c(pend, rr$overlapping, ignore.mcols=TRUE)
    }
    
    random.regions <- c(random.regions, random.regs, ignore.mcols=TRUE)
    pending <- pend
   
    if(length(pending)==0) {
      return(random.regions)
    } else { #Update and continue
      if(allow.overlaps==FALSE) {
        mask <- c(mask, copySeqLevels(to=random.regions, from=genome), ignore.mcols=TRUE)
      }
      ntimes <- ntimes + 1
      if(length(pending)>orig.pending.len*0.7) { #if we placed less than a third of the regions, move to the slower algorithms
        ntimes <- 99 #break
      }
    }
  }
  
 
  # --------------------------------------------------------------------------------
  
  # 2 - SLOWER ALGORITHM
  
  #At this point, random.regs contain the already placed regions and pending those not yet placed
  
  #It's time to call the slower algorithms to place the remaining regions.
  
  if(allow.overlaps==TRUE) { #simply pass to the actual function
    random.regions <- c(random.regions, private_randomizeRegions(A=pending, genome=genome, mask=mask, allow.overlaps=TRUE,...), ignore.mcols=TRUE)    

    return(random.regions)
  } else {
    
    #add the already placed regions to the mask (so there's no overlap)
    mask <- c(mask, copySeqLevels(to=random.regions, from=genome), ignore.mcols=TRUE)
    
    #Create a set of POSSIBLY OVERLAPPING regions first and then remove the overlapping ones
    
    ntimes <- 0
    
    #TODO: Adjust these parameters
    while(ntimes < 10 && length(pending)>100) {
      pending.len <- length(pending)
      
      #Create a set of regions allowing overlaps
      ran.regs <- private_randomizeRegions(A=pending, genome=genome, mask=mask, allow.overlaps=TRUE, ...)
      rr <- removeOverlapping(rs=ran.regs)
      random.regions <- c(random.regions, rr$rs, ignore.mcols=TRUE) #add the newly placed regions to the final GRanges
      pending <- rr$overlapping
      ntimes <- ntimes + 1
    }
    if(length(pending)==0) { #if all regions have been placed, return
      return(random.regions)
    } else {
      #Launch the quadratic algorithm
      ran.regs <- private_randomizeRegions(A=pending, genome=genome, mask=mask, allow.overlaps=FALSE, ...)
      return(c(random.regions, ran.regs, ignore.mcols=TRUE))
    }
  }
  warning("We should not reach this point (end of randomizeRegions), this is a bug")
  return(NA)
}


#This is by far the fastest implementation we have of randomize regions. 
#However, it does not contron for overlaps and is not guaranteed to return a complete
#set of regions (some of them might not be randomized)
randomizeFastWithNoOverlappingControl <- function(A, genome,  mask, max.retries=5, ...) { 
  gam <- getGenomeAndMask(genome=genome, mask=mask)
  genome <- toDataframe(genome)
  genome$chr <- as.character(genome$chr)
  
  mask.chr <- as.character(seqnames(mask))
  mask.start <- start(mask)
  mask.end <- end(mask)
  
  #Prepare chromsome lengths
  chrs <- as.numeric(genome$end)
  chrs <- cumsum(chrs)   
  genome.length <- chrs[length(chrs)]
  
  #This function gets region width and creates a random region in the genome. 
  #If the region is valid (does not span between chromosomes and does not overlap the mask), returns the region
  #else, if it's not valid, it return NA
  getValidReg <- function(reg.width) {
    #Get a random position on the whole genome
    reg.start <- round(runif(1)*(genome.length - reg.width))
    
    num.chr.start <- which(chrs >= reg.start)[1]  #We can be sure al least one will be >= than pos
    num.chr.end <- which(chrs >= reg.start+reg.width)[1]  #We can be sure al least one will be >= than pos
    
    if(num.chr.start != num.chr.end) {
      return(NA)
    }
    
    if(num.chr.start>1) reg.start <- reg.start - chrs[num.chr.start-1]
    reg.end <- reg.start + reg.width
    chr <- genome[num.chr.start, 1]
    
    #Does it overlap with any of the masked regions?
    if(any(mask.chr == chr & mask.start < reg.end  & mask.end > reg.start )) {
      return(NA)
    }
    return(list(chr=chr, start=reg.start, end=reg.end))
  }
  
  
  reg.widths <- width(A) - 1
  #Create a set of valid random of regions. For some regions it might succeed, and for others it will be NA
  rand.reg.list <- lapply(reg.widths, getValidReg)
  
  #While there are still NAs (and for max retries), try to recreate them
  retries <- 0
  nas <- which(is.na(rand.reg.list))
  while(length(nas)>0 & retries < max.retries) {
    rand.reg.list[nas] <- lapply(reg.widths[nas], getValidReg)
    retries <- retries + 1
    nas <- which(is.na(rand.reg.list))
  }
  return(rand.reg.list) #WARNING! This list might contain NAs for regions not randomly placed!
}





#This function creates the random regions in the valid regions and it's capable of taking into account
#the overlap with other regions
private_randomizeRegions <- function(A, genome, mask, allow.overlaps=TRUE, max.retries=5, ...) { 
  if(length(A)<1) {
    return(A)
  }
  
  #use subtractRegions to get the masked genome
  valid.regions <- toDataframe(subtractRegions(genome, mask))[,c(1,2,3)]
  valid.regions[,1] <- as.character(valid.regions[,1])
  names(valid.regions) <- c("chr", "start", "end")
  
  
  #Convenience function to map back the position in the vector of valid positions to the original regions
  map.back <- function(pos) {
    num.region <- which(len.valid >= pos & len.valid != 0)[1]  #We can be sure al least one len.valid will be >= than pos
    if(num.region>1) pos <- pos - len.valid[num.region-1]
    #start <- data.frame(chr=valid.regions[num.region, 1], start=valid.regions[num.region, 2]+pos)
    #new.pos <- data.frame(chr=as.character(seqnames(valid.regions)[num.region]), start=start(valid.regions)[num.region]+pos)
    new.pos <- list(num.region=num.region, chr=as.character(valid.regions[num.region, 1]), start=valid.regions[num.region, 2]+pos)
    return(new.pos)
  }
  
  
  #TODO: Should we really retry? or just giveup immediately?
  #Set up a retry system in case the region set does not fit into the genome with a given configuration
  original.valid.regions <- valid.regions
  for(ntry in 1:max.retries) {
    
    #reset the valid.regions
    valid.regions <- original.valid.regions
    
    #TODO: Sort from largest to smallest to improve the probability of finding a place for all regions in highly fragmented genomes
    
    failed <- FALSE
    first <- TRUE
    #Bernat
    new.chr <- new.start <- new.end <- c()
    for(i in 1:length(A)) {
      
      len <- width(A)[i] - 1
      #CHECK: a -1 is not necessary when using the dataframe fr valid.regions
      len.valid <- valid.regions[,3]-valid.regions[,2] - len -1 #The region cannot start in the last len positions of the valid region. It would not fit
      #len.valid <- width(valid.regions) - len - 1 #The region cannot start in the last len positions of the valid region. It would not fit
      len.valid[len.valid<0] <- 0
      if(length(len.valid)>1) {
        for(j in 2:length(len.valid)) {
          len.valid[j] <- len.valid[j] + len.valid[j-1]
          #print(paste(j,len.valid[j]))
        }
      }
      
      if(max(len.valid)==0) { #If theres no space left to put this region, stop this try
        failed = TRUE
        break
      }
      
      rand.num <- round(runif(1)*len.valid[length(len.valid)])
      
      new.pos <- map.back(rand.num)
      
      #New without dataframe
      new.chr <- c(new.chr, new.pos[["chr"]])
      new.start <- c(new.start, new.pos[["start"]])
      new.end <- c(new.end, new.pos[["start"]] + len)
      
      #Finally, if non overlapping regions are required, remove the last region from the valid regions set
      #Note: the non overlapping processing changes the order of the regions, but that does not have any impact on the randomness
      if(allow.overlaps==FALSE) {
        #Substitute the current region (old, num.region) by two new ones, substracting the new random region just created
        
        old.end <- valid.regions[new.pos[["num.region"]],]
        old.end[1,"start"] <- new.pos[["start"]]+len+1
        
        valid.regions[new.pos[["num.region"]], "end"] <- new.pos[["start"]]-1
        
        valid.regions <- rbind(valid.regions, old.end)  
        
        #If the newly created region was at the exact start position or end position of a valid region, this will create a valid region with size -1. Remove them
        #TODO: we could check this condition before instead of creating and removing.
        negative.lengths <- (valid.regions[,"end"]-valid.regions[,"start"])<0
        if(any(negative.lengths)) {
          valid.regions <- valid.regions[-negative.lengths,]
        }
        
      }
      
    }
    
    if(!failed) {
      random.regions <- toGRanges(data.frame(chr=new.chr, start=new.start, end=new.end))
      random.regions <- copySeqLevels(to=random.regions, from=genome)
      return(random.regions)
    } else {
      warning("It was not possible to create the random region set because there was no space available. Retrying.")
    }
  }
  #If are here, all the retries have failed to produce a valid region set. Raise an error and stop.
  stop(paste("It was not possible to create a valid random region set after ", max.retries, " retries. This might be due to the regions in A covering most of the available genome. Allowing overlapping random regions could solve this problem.", sep=""))
}



#helper functions
copySeqLevels <- function(to, from) {
  seqlevels(to)<-seqlevels(from)
  levels(seqnames(to))<-seqlevels(from)
  return(to)
}

#Detect the overlaps between different regions in rs and remove the from the rs
removeOverlapping <- function(rs) {
  ff <- findOverlaps(rs, rs)
  overlapping.regs <- unique(queryHits(ff)[((queryHits(ff) - subjectHits(ff)) != 0)]) #which region overlaps any region that is not itself?
  
  #and move them from random.regs to pending
  if(length(overlapping.regs)>0) {
    pending <- rs[overlapping.regs]
    rs <- rs[-overlapping.regs]  
  } else {
    pending <- rs[0]
  }
  return(list(rs=rs, overlapping=pending))
}