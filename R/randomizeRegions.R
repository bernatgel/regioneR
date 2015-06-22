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
#' randomizeRegions(A, genome="hg19", mask=NULL, non.overlapping=TRUE, per.chromosome=FALSE, ...)
#' 
#' @param A The set of regions to randomize. A region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}}, \code{\link{data.frame}}, etc...)
#' @param genome The reference genome to use. A valid genome object. Either a \code{\link{GenomicRanges}} or \code{\link{data.frame}} containing one region per whole chromosome or a character uniquely identifying a genome in \code{\link{BSgenome}} (e.g. "hg19", "mm10",... but not "hg"). Internally it uses \code{\link{getGenomeAndMask}}.
#' @param mask The set of regions specifying where a random region can not be (centromeres, repetitive regions, unmappable regions...). A region set in any of the accepted formats by \code{\link{toGRanges}} (\code{\link{GenomicRanges}},\code{\link{data.frame}}, ...). If \code{\link{NULL}} it will try to derive a mask from the genome (currently only works if the genome is a character string). If \code{\link{NA}} it gives, explicitly, an empty mask.
#' @param non.overlapping A boolean stating whether the random regions can overlap (FALSE) or not (TRUE).
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
#' randomizeRegions(A, genome=genome, mask=mask, per.chromosome=TRUE, non.overlapping=TRUE)
#' 
#' @export randomizeRegions
#' 



randomizeRegions <- function(A, genome="hg19",  mask=NULL, non.overlapping=TRUE, per.chromosome=FALSE, ...) { 

  if(!hasArg(A)) stop("A is missing")
  if(!is.logical(non.overlapping)) stop("non.overlapping must be logical")
  if(!is.logical(per.chromosome)) stop("per.chromosome must be logical")
  
  
  A <- toGRanges(A)
  
  #The randomization of an empty region set, is an empty region set
  if(length(A)<1) {
    return(A)
  }
  
  gam <- getGenomeAndMask(genome=genome, mask=mask)
  genome <- gam[["genome"]]
  mask <- gam[["mask"]]
  
  #use subtractRegions to get the masked genome
  valid.regions <- toDataframe(subtractRegions(genome, mask))[,c(1,2,3)]
  valid.regions[,1] <- as.character(valid.regions[,1])
  names(valid.regions) <- c("chr", "start", "end")
  
  

  if(per.chromosome) { 
    missing.chrs <- setdiff(levels(seqnames(A)), valid.regions[,1])
    if(length(missing.chrs > 0)) {
      stop(paste("The chromosomes **", paste(missing.chrs, sep=", "), "** of A are not present in the genome or are completely masked. It is not possible to use \"per.chromosome\" with this dataset.", sep=""))
    }
    #Split the regions in A and the valid regions per chromosomes and create the random regiosn separately per each chromosome
    random.regions <- toGRanges(data.frame(chr=character(), start=numeric(), end=numeric(), stringsAsFactors=FALSE))
    seqlevels(random.regions)<-seqlevels(A)
    levels(seqnames(random.regions))<-seqlevels(A)
       
    
    for(chr in seqlevels(A)) {
      chr.A <- A[seqnames(A)==chr]
      chr.valid <- valid.regions[valid.regions[,1]==chr,]
      new.regions<-hybrid_randomizeRegions(chr.A, chr.valid, non.overlapping=non.overlapping)
      seqlevels(new.regions)<-seqlevels(A)
      levels(seqnames(new.regions))<-seqlevels(A)
      random.regions <- c(random.regions, new.regions)
    }
  } else {
    random.regions <- hybrid_randomizeRegions(A, valid.regions=valid.regions, non.overlapping=non.overlapping)
  }
  
  return(random.regions)
    
}


#Creating non-overlapping regions has quadratic cost, while allowing overlaps has linear cost. 
#To speed up the non.overlapping option, we will use a hybrid approach. Create the the regions allowing overlaps,
#detect the overlapping regions and then use the quadratic algorithm to recreate the overlapping ones.

hybrid_randomizeRegions <- function(A, valid.regions, non.overlapping=FALSE, max.retries=5) {
  
  if(length(A)<1) {
    return(A)
  }
  
  if(non.overlapping==FALSE) { #simply pass to the actual function
    rr <- private_randomizeRegions(A=A, valid.regions=valid.regions, non.overlapping=non.overlapping, max.retries=max.retries)
  } else {
    #Create a set of regions allowing overlaps
      rr <- private_randomizeRegions(A=A, valid.regions=valid.regions, non.overlapping=FALSE, max.retries=max.retries)
            
      
    #Detect the overlaps between different regions
      ov <- overlapRegions(A=rr, B=rr)
      dups <- ov$type != "equal"   #TODO: take into account the situation where, by chance, two regions are exactly the same
      
      if(length(which(dups))>0) {
        pending <- toGRanges(ov[dups, c(1,2,3)]) #the regions taht overlapped something are still pending
            
        toRemove <- which(paste0(seqnames(rr), "#",start(rr),"#", end(rr)) %in%  paste0(seqnames(pending), "#",start(pending),"#", end(pending)))
      
        rr <- rr[-toRemove,] #the regions that overlapped nothing, are kept in the random regions set
      
        #Add the already placed regions into the mask
          valid.regions <- toDataframe(subtractRegions(valid.regions, rr))
        #and place the remaining regions using the quadratic algorithm
          #if the number of regions to place is high (>1000), call hybrid_randomizeRegions recursively, else, call the quadratic algorithm
            #TODO: What if the regions do not fit into the genome? We should place a guard against taht situation limiting the number of recursive calls
          if(length(pending)>1000) {
            rr2 <- hybrid_randomizeRegions(A=pending, valid.regions=valid.regions, non.overlapping=TRUE, max.retries=max.retries)
          } else {
            rr2 <- private_randomizeRegions(A=pending, valid.regions=valid.regions, non.overlapping=TRUE, max.retries=max.retries)
          }
          suppressWarnings(rr <- append(rr, rr2)) #supress warnings, since it will warn of potential different reference genome if we have just a few regions
      }
  }
  
  return(rr)
  
}



#This is the function that actually creates the random regions in the valid regions  
#Idea: Generar una "estructura" amb les posicions POSSIBLES (depen de la regi? a tirar). 
#Triar aleatoriament i remapar enrera.
#TODO: Vigilar amb els possible -1? O no cal pq estem parlant de longituds?

private_randomizeRegions <- function(A, valid.regions, non.overlapping=FALSE, max.retries=5) {
  
  if(length(A)<1) {
    return(A)
  }
  
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
     
        len <- width(A)[i]
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
        if(non.overlapping) {
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
        random.regions <- data.frame(chr=new.chr, start=new.start, end=new.end)
        return(toGRanges(random.regions))
      } else {
        warning("It was not possible to create the random region set because there was no space available. Retrying.")
      }
  }
  #If are here, all the retries have failed to produce a valid region set. Raise an error and stop.
  stop(paste("It was not possible to create a valid random region set after ", max.retries, " retries. This might be due to the regions in A covering most of the available genome. Allowing overlapping random regions could solve this problem.", sep=""))
}

