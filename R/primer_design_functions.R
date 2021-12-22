#' Primer design
#'
#' This function allows you to design qPCR/ddPCR primer and probe sets that fit to a large set of related nucleic acid sequences. 
#' Tm is calculated using the nearest neighbor method.
#' @param sequences Input sequences (DNAStringSet object) that primers should be designed for
#' @param length_max Maximum length of primers to be considered. Defaults to 30 nt.
#' @param length_min Minimum length of primers to be considered. Defaults to 16 nt.
#' @param ROI Numeric vector of format c(start, end) that specifies region of interest within target sequence. Defaults to NA (whole sequence is examined)
#' @param min_dist Minimum distance between primers and probes. Defaults to 2 nt.
#' @param max_dist Maximum distance between primers and probes. Defaults to 50 nt. WARNING: strongly increases computation time
#' @param gc_thres Minimum GC content (in %) that any sequence must have to be considered.
#' @param amplicon_length_opt Optimal amplicon length. Defaults to 90 nt.
#' @param mismatch_tolerance How many ambiguous bases will be tolerated in any primer/probe sequence. Increase only for highly heterogeneous data sets.
#' @param Tm_in Target melting temperature for primers. Defaults to NA, then any primers with Tm >50 degree C and Tm<65 degree C will be considered.
#' @param Tm_delta_max Maximum Tm difference between forward and reverse primer that will be tolerated.
#' @param Tm_delta_probe Minimum difference between primer Tm and probe Tm in any given primer-probe set. Defaults to 5 degree C.
#' @param n_sets How many primer-probe sets to return. Defaults to 200.
#' @param Na Sodium concentration in mM to be assumed for Tm calculation. Defaults to 300 mM. If set to a non-numeric value, Tm will be calculated based on GC content.
#' @param primer_conc Primer concentration in microM to be assumed for Tm calculation. Defaults to 0.5 microM. If set to a non-numeric value, Tm will be calculated based on GC content.
#' @param mismatch_threshold Threshold of which fraction of the input sequences can disagree at one specific base for it still to be considered a clear consensus. Defaults to 0.05.
#' @param max_repetitions Maximum number of repetitions of the same (di-) nucleotide in any primer/probe sequence. Defaults to 3.
#' @param max_complementarity Maxmium alignment score with self or other primers of the same set. Defaults to 4.
#' @param seqs_aligned Whether input sequences are aligned. Setting to TRUE if sequences are already aligned saves computation time.
#' @param position_probe If given (as nt from start of the complete sequence), only probes covering this position in the sequence will be considered
#' @return list of two elements:
#' @return  all_seqs: data.frame with unique ID, start position within sequence (pos), sequence information (seq), melting Temperature (Tm), GC content (GC), number of ambivalent bases (mismatches), self complementarity score (self_complementary) and self 3prime complementarity score (self_3prime_complementarity) 
#' @keywords primer
#' @details This function searches for primer-probe sets within the input sequences that fit without any ambiguous base pairings to as many of the sequences as possible. It first aligns the input sequences and determines the consensus sequence, then it produces a list of all possible primers that could be found in the target sequence. This list is then filtered to eliminate any oligonucleotides that would not be suitable as primers(too many ambiguous bases, wrong melting temperature, too low GC content ...). Then, for each primer combinations with a probe and another primer are determined and examined in accordance with the input parameters (difference in melting temperature, position, ...). The function then gives a list of potential primer-probe sets with a score that gives an approximate measurement of the performance of the primer-probe set in qPCR/ddPCR.
#' 

#TO DO: investigate whether porting loops to C++ could increase speed!
#TO DO: if not, get rid of loops!

# function to get primer candidates
design_primers<-function(sequences, length_max=30, length_min=16, ROI=NA, min_dist=2, max_dist=50, gc_thres=30, amplicon_length_opt=90, mismatch_tolerance=0, Tm_in=NA, Tm_delta_max=1.5, Tm_delta_probe=5, n_sets=200, Na=300, primer_conc=0.5, mismatch_threshold=0.05, max_repetitions=3, max_complementarity=4, check_spec=NA, position_probe=NA, seqs_aligned=FALSE){
  suppressPackageStartupMessages(library(DECIPHER))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(TmCalculator))
  suppressPackageStartupMessages(library(purrr))
  
  # set Tm values, if no input is given determine upper and lower threshold according to general guidelines
  if(is.na(Tm_in)){
    Tm_upper<-65
    Tm_lower<-50
  }else{
    Tm_upper<-Tm_in+Tm_delta_max
    Tm_lower<-Tm_in-Tm_delta_max
  }
  #adjust position of probe to be in relation to start of ROI
  if(!is.na(position_probe)){
    if(!is.na(ROI)){
      position_probe<-position_probe-ROI[[1]]
    }
  }
  
  if(length(sequences)>1){
    if(seqs_aligned==FALSE){
      #align sequences
      alignment<-sequences %>% RemoveGaps() %>% AlignSeqs(., verbose=FALSE) # align sequences
    }
    
    # get consensus sequence
    consensus<-ConsensusSequence(alignment, threshold=mismatch_threshold)
  }else{
    consensus<-sequences
  }
  if(length(which(is.na(ROI)))==0){
    if(ROI[[1]]<1){
      ROI[[1]]<-1
    }
    if(ROI[[2]]>width(consensus)-length_max){
      ROI[[2]]<-width(consensus)-length_max
    }
  }
  
  # if probe position is set, set ROI around its position
  if(!is.na(position_probe)){
    
    ROI<-c(position_probe-150, position_probe+150)
    position_probe<-151
    
  }
  
  if(!length(ROI)==2){
    ROI<-c(1, width(consensus)-length_max)
  }
  
  # crop sequence to ROI
  consensus<-subseq(consensus, ROI[[1]], ROI[[2]]+length_max)
  # go through all positions and find primers of different lengths, make df with Tm, GC, mismatches
  # find all possible primer and probe sequences in the given ROI with Tm, GC content etc
  n=(width(consensus))*(length_max-length_min+1)
  df<-data.frame(ID=numeric(length = n), start=numeric(length = n), end=numeric(length=n), seq=character(length=n), Tm=numeric(length=n), GC=numeric(length=n), mismatches=numeric(length=n), self_complementary=numeric(length=n), self_3prime_complementary=numeric(length=n))
  ind=0
  
  # IDEA: for-loop only to create df with positions and IDs; then use apply/maps to calculate Tm etc.?
  for(i in 1:(width(consensus))){
    for(j in length_min:length_max){
      ind<-ind+1
      df$start[[ind]]<-as.numeric(i)
      df$end[[ind]]<-as.numeric(i+j-1)
      df$ID[[ind]]<-ind
    }
  }
  
  # get rid of primer sequences extending over the ROI
  df <- df %>% filter(end<width(consensus))
  df <- df %>% filter(start>0)
  
  subset<-function(vec, names){
    temp<-vec[names]
    return(temp[which(!is.na(temp))])
  }
  
  generate_primer<-function(vec, seq){
    vec<-as.numeric(vec)
    sub<-subseq(seq, start=vec[[2]], end=vec[[3]]) %>% as.character()
    #get GC content
    vec[[6]]<-GC(sub)
    #get sequence
    vec[[4]]<-sub
    
    #find number of mismatches
    vec[[7]]<-sub %>% s2c() %>% table() %>% subset(., c("M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N", "-", "+", ".")) %>% sum()
    # estimate self-complementarity
    return(unlist(vec))
  }
  
  # restore proper dataframe structure
  df_new<-apply(df, 1, generate_primer, seq=consensus) %>% t() 
  df_new<- df_new %>% data.frame()
  colnames(df_new)<-colnames(df)
  rownames(df_new)<-df_new$ID
  
  df<-df_new
  df$start<-as.numeric(df$start)
  df$end<-as.numeric(df$end)
  df$ID<-as.numeric(df$ID)
  df$GC<-as.numeric(df$GC)
  df$mismatches<-as.numeric(df$mismatches)
  df$Tm<-as.numeric(df$Tm)
  df$self_complementary<-as.numeric(df$self_complementary)
  df$self_3prime_complementary<-as.numeric(df$self_3prime_complementary)
  
  # filter out primers with too extreme GC content or too many mismatches to save computation time
  df <- df %>% filter(GC>gc_thres)  %>% filter(mismatches<mismatch_tolerance+1)
  
  #calculate all melting temperatures
  if(!Tm_method=="TibMolBiol"){
    if(mismatch_tolerance==0 && length(which(is.na(c(Na,primer_conc))))==0){
      # NN only works with non-ambiguous bases and with salt concentrations given
      df$Tm<-df$seq %>% lapply(., FUN=Tm_NN, Mg=1, Na=Na,dnac1=primer_conc, dnac2=0.1, nn_table="DNA_NN3") 
    }else{
      # otherwise use Wallace method
      df$Tm<-df$seq %>% lapply(., FUN=Tm_Wallace, ambiguous=TRUE)
    }
  }else{
    dnac1=0.25
    dnac2=0.01
    Na=0.1
    Mg=0.5
    K=0.1
    Tris=0.1
    df$Tm<-df$seq %>% lapply(., FUN=Tm_NN, dnac1 = dnac1, dnac2=dnac2, Na=Na, Mg=Mg, K=K, Tris=Tris, nn_table="DNA_NN1")
  }
  
  
  # filter out primers with Tm outside of range to be considered
  df <- df %>% filter(Tm>Tm_lower) %>% filter(Tm<Tm_upper+Tm_delta_probe)
  
  # estimate self-complementarity
  df$self_complementary<-lapply(df$seq, check_self_complementarity)
  df$self_3prime_complementary<-lapply(df$seq, check_3prime_selfcomplementarity)
  
  # filter out those with more than (max_repetitions) repetitions of same nucleotide or dinucleotide
  failed<-unique(c(grep(paste("(.)\\1{", as.character(max_repetitions), ",}", sep=""), df$seq), grep(paste("(..)\\1{", as.character(max_repetitions), ",}", sep=""), df$seq)))
  df<-df[which(!df$ID %in% failed),]
  
  # filter primers and probes according to input parameters
  df_primer<- df %>% filter(Tm>Tm_lower) %>% filter(Tm<Tm_upper) %>% filter(self_complementary<max_complementarity+1) %>% filter(self_3prime_complementary<max_complementarity+1)
  
  df_probe<- df %>% filter(Tm>Tm_lower+Tm_delta_probe) %>% filter(Tm<Tm_upper+Tm_delta_probe) %>% filter(GC>gc_thres)  %>% filter(mismatches<mismatch_tolerance+1)  %>% filter(self_complementary<max_complementarity+1)  %>% filter(self_3prime_complementary<max_complementarity+1)
  
  df_sets<-data.frame(ID_fwd=numeric(), ID_probe=numeric(), ID_rev=numeric(), deltaTm=numeric())
  # find potential matching primer and probe 
  for(i in 1:dim(df_primer)[[1]]){
    if(!is.na(position_probe)){
      if(df_primer$end[[i]]<position_probe-length_max*2){
        next
      }
      if(df_primer$start[[i]]>position_probe){
        next
      }
    }
    
    #find potential probe positions
    range_probe=c(df_primer$start[[i]]+nchar(df_primer$seq[[i]])+min_dist+1, df_primer$end[[i]]+nchar(df_primer$seq[[i]])+max_dist+1)
    
    if(!is.na(position_probe)){
      if(!(position_probe > range_probe[[1]] && position_probe < range_probe[[2]])){
        next
      }
    }
    
    probes<-df_probe %>% filter(start>range_probe[[1]]) %>% filter(end<range_probe[[2]]) %>% filter(Tm>df_primer$Tm[[i]]+Tm_delta_probe) %>% filter(Tm<df_primer$Tm[[i]]+Tm_delta_probe+5)
    # filter out probes that start with a G (quenches fluorophore!)
    probes<-probes[grep("^G\\w+", unlist(probes$seq), invert=TRUE),]
    #if(!is.na(position_probe)){
    #  probes_temp <- probes %>% filter(!pos>position_probe) %>% filter(!(pos+nchar(seq))< position_probe)
    #}
    if(dim(probes)[[1]]>0){
      for(j in 1:dim(probes)[[1]]){
        #find potential reverse primer positions
        range_rev<-c(probes$start[[j]]+nchar(probes$seq[[j]])+min_dist+1, probes$end[[j]]+nchar(probes$seq[[j]])+max_dist+1)
        reverse<-df_primer %>% filter(start>range_rev[[1]]) %>% filter(end<range_rev[[2]]) %>% filter(Tm>df_primer$Tm[[i]]-Tm_delta_max) %>% filter(Tm<df_primer$Tm[[i]]+Tm_delta_max)
        
        if(dim(reverse)[[1]]>0){
          df_sets<-rbind(df_sets, data.frame(ID_fwd=rep(df_primer$ID[[i]], times=dim(reverse)[[1]]), ID_probe=rep(probes$ID[[j]], times=dim(reverse)[[1]]), ID_rev= reverse$ID))
        }
      }
    }
  }
  # calculate deviances from optimum values and score + order primer-probe sets according toadherence to optimum values
  if(dim(df_sets)[[1]]==0){
    stop("Error: No primer-probe combinations were found for your search parameters. Please repeat the search with less stringent requirements.")
  }
  
  df_sets$mismatches<-NA
  df_sets$deltaTm<-NA
  df_sets$delta_GC<-NA
  df_sets$amplicon_size<-NA
  df_sets$gc_3prime<-NA
  df_sets$Tm_amplicon<-NA
  df_sets$cross_complementarity<-NA
  df_sets$score<-NA
  df_sets$fwd<-NA
  df_sets$probe<-NA
  df_sets$rev<-NA
  df_sets$amplicon_seq<-NA
  df_sets$comp_mod<-NA
  df_sets$GC_ratio<-NA
  
  temp<-list(df[as.character(unlist(df_sets[,1])),], df[as.character(unlist(df_sets[,2])),], df[as.character(unlist(df_sets[,3])),])
  
  # get difference in Tm
  df_sets$deltaTm<-abs(unlist(temp[[1]]$Tm) - unlist(temp[[3]]$Tm))
  
  # get difference from optimal GC
  df_sets$delta_GC<-rowSums(abs(cbind(temp[[1]]$GC-50, temp[[2]]$GC-50, temp[[3]]$GC-50)))
  
  # get amplicon size
  df_sets$amplicon_size<-temp[[3]]$end-temp[[1]]$start
  
  # get number of mismatches
  df_sets$mismatches<-rowSums(cbind(temp[[1]]$mismatches, temp[[2]]$mismatches, temp[[3]]$mismatches))
  
  # check for GC clamp
  gc_end<-numeric(length=dim(temp[[1]])[[1]])
  
  gc_fwd<-c(grep("G$", temp[[1]]$seq), grep("C$", temp[[1]]$seq))
  gc_rev<-c(grep("^G", temp[[3]]$seq), grep("^C", temp[[3]]$seq))
  gc_end[gc_fwd]<-gc_end[gc_fwd]+1
  gc_end[gc_rev]<-gc_end[gc_rev]+1
  df_sets$gc_3prime<-gc_end
  
  df_sets$fwd<-temp[[1]]$seq
  df_sets$probe<-temp[[2]]$seq
  df_sets$rev<-temp[[3]]$seq %>% DNAStringSet() %>% reverseComplement() %>% as.character()
  
  #calculate score
  df_sets$score<-1/(df_sets$deltaTm^2+0.001)-(df_sets$delta_GC^2)/2-df_sets$mismatches^4*300-(abs(df_sets$amplicon_size-amplicon_length_opt)^2*3)+df_sets$gc_3prime*400#-df_sets$comp_mod
  
  #check position of probe in relation to target sequence
  if(!is.na(position_probe)){
    df_sets$position_probe_delta<-abs(temp[[2]]$pos+nchar(temp[[2]]$seq)/2-position_probe)
    df_sets$score<-df_sets$score-df_sets$position_probe_delta*50
  }
  
  if(n_sets>dim(df_sets)[[1]]){
    n_sets<-dim(df_sets)[[1]]
  }
  
  df_sets<-df_sets[order(df_sets$score, decreasing=TRUE),]
  df_sets<-df_sets[1:n_sets,]
  
  for(i in 1:dim(df_sets)[[1]]){
    sub<-df[as.character(unlist(df_sets[i,c(1:3)])),]
    
    # get Tmand sequence of amplicon
    amplicon<-substr(consensus, sub$start[[1]], sub$end[[3]])
    df_sets$amplicon_seq[[i]]<-amplicon
    if(mismatch_tolerance==0 && length(which(is.na(c(Na, primer_conc))))==0){
      # NN only works with non-ambiguous bases and with salt concentrations given
      df_sets$Tm_amplicon[[i]]<-Tm_NN(amplicon, Mg=1, Na=Na,dnac1=primer_conc, dnac2=0.1, nn_table="DNA_NN3")
    }else{
      # otherwise use Wallace method
      df_sets$Tm_amplicon[[i]]<-Tm_Wallace(amplicon, ambiguous=TRUE)
    }
    
    # get maximum cross-complementarity between any of the primers and/or probes:
    df_sets$cross_complementarity[[i]]<-check_cross_complementarity(sub$seq[[1]], sub$seq[[2]], sub$seq[[3]])
    comp_mod<-0
    if(df_sets$cross_complementarity[[i]]>max_complementarity){
      comp_mod<-3000
    }
    df_sets$comp_mod[[i]]<-comp_mod
    
    t<- df_sets$probe[i] %>% s2c() %>%  table()
    if(t["G"]<t["C"]){
      df_sets$GC_ratio[[i]]<-1000
    }else{
      df_sets$GC_ratio[[i]]<-0
    }
    
  }
  df_sets$comp_mod<-as.numeric(unlist(df_sets$comp_mod))
  df_sets$score<-df_sets$score-df_sets$comp_mod
  df_sets$score<-df_sets$score-df_sets$GC_ratio
  df_sets<-df_sets[order(df_sets$score, decreasing=TRUE),] 
  
  if(class(check_spec)=="DNAStringSet"){
    df_sets <- check_specificity(df_sets, check_spec)
  }
  
  out<-list(all_seqs=df, best_candidates=df_sets[1:n_sets,])
  
  if(length(sequences)>1){
    out<-validate_primers(out, sequences)
  }
  
  return(out)
}

check_3prime_selfcomplementarity<-function(p){
  # make reverse complement
  p2<-reverseComplement(p) %>% as.character() %>% s2c()
  p<-as.character(p) %>% s2c()
  
  #loop through sequence and find whether 5 or more nucleotides are complementary
  score=0
  i=1
  while(i<length(p)){
    sub1=p[(length(p)-i+1):length(p)]
    sub2=p2[1:i]
    temp=length(which(sub1==sub2))-length(which(!sub1==sub2))
    i=i+1
    if(temp>score){
      score<-temp
    }
  }
  return(score)
}

check_cross_complementarity<-function(p1, p2, p3){
  p1 <- p1 %>% DNAStringSet()
  p2 <- p2 %>% DNAStringSet()
  p3 <- p3 %>% DNAStringSet() %>% reverseComplement()
  
  c1 <- append(p1, reverseComplement(p2)) %>% AlignSeqs(. , verbose = FALSE, gapExtension = -20) %>% ConsensusSequence(., includeTerminalGaps = TRUE) %>% as.character()
  c2 <- append(p1, reverseComplement(p3)) %>% AlignSeqs(. , verbose = FALSE, gapExtension = -20) %>% ConsensusSequence(., includeTerminalGaps = TRUE) %>% as.character()
  c3 <- append(p2, reverseComplement(p3)) %>% AlignSeqs(. , verbose = FALSE, gapExtension = -20) %>% ConsensusSequence(., includeTerminalGaps = TRUE) %>% as.character()
  
  l<-list(c1, c2, c3) %>% lapply(., function(x){return(gsub("^\\++", "", x))}) %>% lapply(., function(x){return(gsub("\\++$", "", x))}) %>% lapply(., s2c)
  scores<-unlist(lapply(l, function(x){return(length(which(x %in% c("A", "C", "G", "T"))))})) - unlist(lapply(l, function(x){return(length(which(!x %in% c("A", "C", "G", "T"))))}))
  
  return(max(scores))
}

check_self_complementarity<-function(p1){
  p1 <- p1 %>% DNAStringSet()
  
  c1 <- append(p1, reverseComplement(p1)) %>% AlignSeqs(. , verbose = FALSE, gapExtension = -20) %>% ConsensusSequence(., includeTerminalGaps = TRUE) %>% as.character()

  l<-c1 %>% gsub("^\\++", "", .) %>% gsub("\\++$", "", .) %>% s2c
  score<-length(which(l %in% c("A", "C", "G", "T"))) - length(which(!l %in% c("A", "C", "G", "T")))
  
  return(score)
}

#' Primer design
#'
#' This function searches for matches with other sequences than the target sequence.
#' @param primers A data.frame as produced by design_primers()
#' @param reference Input sequences (DNAStringSet object) that primers should be compared to.
#' @param cutoff Threshold of the Hamming distance below which sequences are treated as part of the same cluster.
#' @return data.frame as from design_primers() with additional information on matching sequences.
#' @keywords primer
#' @details This function searches for sequences in the data set that match with the primer-probe set and stores the sequence name and minimum alignment score of any mismatching sequence in the primers data.frame. To this end it first splits the input data into multiple clusters and compares the primer sequences to each of the cluster consensus sequences.

check_specificity<-function(primers, reference, cutoff=0.005){ # save=FALSE,
    if(length(reference)>10){
      clusters<- DistanceMatrix(reference, verbose = FALSE) %>% IdClusters(type="clusters", cutoff=cutoff, method = "UPGMA", showPlot = TRUE, verbose = FALSE)
      cluster_cons<-DNAStringSet()
      nclust<-max(clusters$cluster)
      for(i in 1:nclust){
        temp<-reference[filter(clusters, cluster==i) %>% rownames()]
        if(length(temp)>1){
          cons<-temp %>% AlignSeqs(., verbose=FALSE) %>% ConsensusSequence() 
          names(cons)[[1]]<-names(temp)[[1]]
        }else{
          cons<-temp
        }
        cluster_cons<-append(cluster_cons, cons)
      }
      seqs<-cluster_cons 
    }else{
      seqs<-reference
    }
    
    # get amplicon sequences from output of previous function
    #amplicon_seqs<-primers$best_candidates %>% pull(amplicon_seq) %>% unlist() %>% DNAStringSet()
    fwd<-primers$best_candidates %>% pull(fwd) %>% DNAStringSet()
    probe<-primers$best_candidates %>% pull(probe) %>% DNAStringSet()
    rev<-primers$best_candidates %>% pull(rev) %>% DNAStringSet() %>% reverseComplement()
    names_seq<-seqs %>% names()
     
    # align each amplicon with the reference, find reference sample with lowest distance to it
    for(i in 1:dim(primers$best_candidates)[[1]]){
      # for every cluster within the data
      al_fwd<-pairwiseAlignment(seqs, fwd[i], type="local")
      al_probe<-pairwiseAlignment(seqs, probe[i], type="local")
      al_rev<-pairwiseAlignment(seqs, rev[i], type="local")
      
      sc_fwd<-al_fwd %>% score()
      name_fwd<-names_seq[which(sc_fwd==max(sc_fwd))[[1]]]
      primers$best_candidates$fwd_pot_match[[i]]<-name_fwd
      primers$best_candidates$fwd_al_score[[i]]<-max(sc_fwd)
      
      sc_probe<-al_probe %>% score()
      name_probe<-names_seq[which(sc_probe==max(sc_probe))[[1]]]
      primers$best_candidates$probe_pot_match[[i]]<-name_probe
      primers$best_candidates$probe_al_score[[i]]<-max(sc_probe)
      
      sc_rev<-al_rev %>% score()
      name_rev<-names_seq[which(sc_rev==max(sc_rev))[[1]]]
      primers$best_candidates$rev_pot_match[[i]]<-name_rev
      primers$best_candidates$rev_al_score[[i]]<-max(sc_rev)
    }

  return(primers)
}

#' Primer design
#'
#' This function searches for mismatches of primer-probe sets with a set of sequences.
#' @param primers A data.frame as produced by design_primers()
#' @param seqs Input sequences (DNAStringSet object) that primers should be compared to.
#' @return data.frame as from design_primers() with additional information on mismatching sequences.
#' @keywords primer
#' @details This function searches for sequences in the data set that have mismatches with any of the primer-probe sets and stores the sequence name and minimum alignment score of any mismatching sequence in the primers data.frame. 

validate_primers <- function(primers, seqs){
  fwd<-primers$best_candidates %>% pull(fwd) %>% DNAStringSet()
  probe<-primers$best_candidates %>% pull(probe) %>% DNAStringSet() 
  rev<-primers$best_candidates %>% pull(rev) %>% DNAStringSet() %>% reverseComplement()
  names_seq<-seqs %>% names()
  
  primers$best_candidates$fwd_pot_mismatch<-character(length(dim(primers$best_candidates)[[1]]))
  primers$best_candidates$fwd_mismatch_score<-numeric(length(dim(primers$best_candidates)[[1]]))
  primers$best_candidates$rev_pot_mismatch<-character(length(dim(primers$best_candidates)[[1]]))
  primers$best_candidates$rev_mismatch_score<-numeric(length(dim(primers$best_candidates)[[1]]))
  primers$best_candidates$probe_pot_mismatch<-character(length(dim(primers$best_candidates)[[1]]))
  primers$best_candidates$probe_mismatch_score<-numeric(length(dim(primers$best_candidates)[[1]]))
  
  for( i in 1:dim(primers$best_candidates)[[1]]){
    al_fwd<-pairwiseAlignment(seqs, fwd[i], type="local")
    al_probe<-pairwiseAlignment(seqs, probe[i], type="local")
    al_rev<-pairwiseAlignment(seqs, rev[i], type="local")
    
    sc_fwd<-al_fwd %>% score()
    if(!min(sc_fwd) == max(sc_fwd)){
      lim<-max(sc_fwd)
      names_fwd<-names_seq[which(sc_fwd<lim)]
      primers$best_candidates$fwd_pot_mismatch[[i]]<-paste(names_fwd, collapse="_AND_")
      primers$best_candidates$fwd_mismatch_score[[i]]<-min(sc_fwd)
    }
        
    sc_probe<-al_probe %>% score()
    if(!min(sc_probe) == max(sc_probe)){
      lim<-max(sc_probe)
      names_probe<-names_seq[which(sc_probe<lim)]
      primers$best_candidates$probe_pot_mismatch[[i]]<-paste(names_probe, collapse="_AND_")
      primers$best_candidates$probe_mismatch_score[[i]]<-min(sc_probe)
    }
    
    sc_rev<-al_rev %>% score()
    if(!min(sc_rev) == max(sc_rev)){
      lim<-max(sc_rev)
      names_rev<-names_seq[which(sc_rev<lim)]
      primers$best_candidates$rev_pot_mismatch[[i]]<-paste(names_rev, collapse="_AND_")
      primers$best_candidates$rev_mismatch_score[[i]]<-min(sc_rev)
    }
  }
  
  seq_names_fwd <-primers$best_candidates$fwd_pot_mismatch %>% strsplit(., split="_AND_") #%>% unlist() %>% unique()
  seq_names_rev <-primers$best_candidates$rev_pot_mismatch %>% strsplit(., split="_AND_") #%>% unlist() %>% unique()
  seq_names_probe <-primers$best_candidates$probe_pot_mismatch %>% strsplit(., split="_AND_") #%>% unlist() %>% unique()
  
  nseq <- c()
  for(i in 1:length(seq_names_fwd)){
    nseq[i] <- c(seq_names_fwd[[i]], seq_names_rev[[i]], seq_names_probe[[i]]) %>% unlist() %>% unique() %>% length()
  }
  primers$best_candidates$n_mismatched_seqs<-nseq
  primers$best_candidates$score<-primers$best_candidates$score-(((primers$best_candidates$n_mismatched_seqs/length(seqs))^2)*10000)
  primers$best_candidates<-primers$best_candidates[order(primers$best_candidates$n_mismatched_seqs),]
  
  return(primers)
}

#' Primer design
#'
#' This function displays the sequences of a primer-probe set together with any sequences in the data set that it has any mismatches with.
#' @param primers A line of a data.frame as produced by design_primers()
#' @param sequences Input sequences (DNAStringSet object) that primers were validated against (call to validate_primers() before using this function!).
#' @param ref_seq Additional sequences that should always be displayed.
#' @return NULL. Opens browser window to inspect sequences.
#' @keywords primer
#' @details This function aligns the primer-probe set with sequences that previously mismatches were detected with. It then displays the alignment for manual inspection. Additional sequences that should always be displayed, regardless of whether there are mismatches with the primers, can be passed as ref_seq.

display_mismatched_sequences<-function(primers, sequences, ref_seq=NA){
  seq_names<-c(primers$fwd_pot_mismatch, primers$rev_pot_mismatch, primers$probe_pot_mismatch) %>% strsplit(., split="_AND_") %>% unlist() %>% unique()
  
  prim_seq<-DNAStringSet(c(primers$fwd, primers$probe, primers$rev))
  names(prim_seq)<-c("forward", "probe", "reverse")
  prim_seq["reverse"] <- prim_seq["reverse"] %>% reverseComplement()
  seqs<-sequences[seq_names]
  
  seqs<-append(seqs, prim_seq)
  
  if(!is.na(ref_seq)){
    seqs<-append(seqs, ref_seq)
  }
  
  al<- seqs %>% RemoveGaps() %>% AlignSeqs(., verbose=FALSE)
  BrowseSeqs(al)
  
}

#' Primer design
#'
#' This function displays the sequences of a primer-probe set together with any sequences in the data set that it has significant sequence similarity to.
#' @param primers A line of a data.frame as produced by design_primers()
#' @param sequences Input sequences (DNAStringSet object) that primers were checked for unspecific binding (call to check_specificity() before using this function!).
#' @param ref_seq Additional sequences that should always be displayed.
#' @return NULL. Opens browser window to inspect sequences.
#' @keywords primer
#' @details This function aligns the primer-probe set with sequences that previously matches were detected with. It then displays the alignment for manual inspection. Additional sequences that should always be displayed, regardless of whether there are mismatches with the primers, can be passed as ref_seq.


display_matched_sequences<-function(primers, sequences, ref_seq=NA){
  seq_names<-unique(c(primers$fwd_pot_match, primers$rev_pot_match, primers$probe_pot_match))
  prim_seq<-DNAStringSet(c(primers$fwd, primers$probe, primers$rev))
  names(prim_seq)<-c("forward", "probe", "reverse")
  prim_seq["reverse"] <- prim_seq["reverse"] %>% reverseComplement()
  seqs<-sequences[unlist(seq_names)]
  
  seqs<-append(seqs, prim_seq)
  
  if(!is.na(ref_seq)){
    seqs<-append(seqs, ref_seq)
  }
  
  al<- seqs %>% RemoveGaps() %>% AlignSeqs(., verbose=FALSE)
  BrowseSeqs(al)
}

#' Primer design for large and heterogeneous sets of sequences
#'
#' This function allows you to design qPCR/ddPCR primer and probe sets that fit to a large and very heterogeneous set of nucleic acid sequences. Multiple sets are designed for distinct clusters of sequences.
#' @param sequences Input sequences (DNAStringSet object) that primers should be designed for
#' @param cutoff Maximum Hamming distance between two sequences to still be considered part of the same cluster. 
#' @param min_cluster_size Minimum number of samples in a cluster to warrant designing a specialised primer-probe set for
#' @param ROI as in design_primers()
#' @param position_probe as in design_primers()
#' @param ... additional parameters passed to design_primers()
#' @return list of data.frames as produced by design_primers()
#' @keywords primer
#' @details This function first calls the design_primer() function to design a primer-probe set that fits for as many of the input sequences as possible. Then, the remaining sequences are split into clusters according to their sequence similarity. If any of the clusters reach a size of min_cluster_size, an additional primer-probe set is designed for it.


automated_primer_design<-function(sequences, cutoff=0.005, min_cluster_size=5, ROI=NA, position_probe=NA, ...){
  print("Designing primers for whole dataset")
  
  cluster_seqs<-list()
  cluster_primers<-list()
  
  m<-length(sequences)
  
  primers<-design_primers(sequences, ROI=ROI, position_probe=position_probe, ...)
  print("Identifying clusters in sequences with mismatches for more specific primer design")
  val<-primers#validate_primers(primers, sequences)
  seq_names<-c(val$best_candidates$fwd_pot_mismatch[[1]], val$best_candidates$rev_pot_mismatch[[1]], val$best_candidates$probe_pot_mismatch[[1]]) %>% strsplit(., split="_AND_") %>% unlist() %>% unique()
  matched_seq<-sequences[which(!names(sequences) %in% seq_names)]
  sequences <- sequences[seq_names] %>% RemoveGaps() %>% AlignSeqs(., verbose=FALSE)
  
  if(!is.na(ROI)){
    subs<-subseq(sequences, start=ROI[1], end=ROI[2])
  }
  
  if(!is.na(position_probe)){
    coords<-c(position_probe-100, position_probe+100)
    coords[which(coords<1)]<-1
    coords[which(coords>width(sequences))]<-width(sequences)
    subs<-subseq(sequences, start=coords[1], end=coords[2])
  }
  clusters<- DistanceMatrix(subs, verbose = FALSE) %>% IdClusters(type="clusters", cutoff=cutoff, method = "UPGMA", showPlot = TRUE, verbose = FALSE)
  nclust<-max(clusters$cluster)

  for(i in 1:nclust){
    c("Processing cluster", as.character(i), "of", as.character(nclust)) %>% paste(., collapse=" ") %>% print()
    temp<-sequences[filter(clusters, cluster==i) %>% rownames()]
    if(length(temp)>(min_cluster_size-1)){
      temp_primers<-design_primers(temp, ROI=ROI, position_probe=position_probe, seqs_aligned=TRUE, ...)
      cluster_primers[[i]]<-temp_primers
      
    }else{
      c("Cluster", as.character(i), "was ignored due to size < min_cluser_size") %>% paste(., collapse=" ") %>% warning()
    }
    cluster_seqs[[i]]<-temp
  }
  cluster_primers[[length(cluster_primers)+1]]<-val
  cluster_seqs[[length(cluster_seqs)+1]]<-matched_seq
  vec<-c()
  for(i in 1:length(cluster_seqs)){vec<-c( vec, length(cluster_seqs[[i]]))}
  k<-sum(vec)
  cluster_seqs<-cluster_seqs[which(!is.null(cluster_seqs))]
  cluster_primers<-cluster_primers[which(!is.null(cluster_primers))]
  warning(m-k, " out of ", m, " samples have been ignored due to having too many sequence differences")
  return(list(cluster_primers, cluster_seqs))
}

eval_primer_combination<-function(primers, seqs){
  val<-validate_primers(primers, seqs)
  
  eval<-function(primer){
    seq_names<-primer %>% strsplit(., split="_AND_") %>% unlist() %>% table()
    seq_names<-seq_names[which(seq_names==2)] %>% names()
    return(seq_names)
  }
  
  ind<-which(colnames(primers) %in% c("fwd_pot_mismatch", "rev_pot_mismatch", "probe_pot_mismatch"))
  sequences<-apply(primers[,ind], 2, eval) %>% unlist() %>% unique()

  prim_seq<-DNAStringSet(c(primers$fwd, primers$probe, primers$rev))
  names(prim_seq)<-c(rep("forward", dim(primers)[[1]]), rep("probe", dim(primers)[[1]]), rep("reverse", dim(primers)[[1]]))
  prim_seq[grep("reverse", names(prim_seq))] <- prim_seq["reverse"] %>% reverseComplement()
  seqs<-seqs[sequences]
  
  seqs<-append(seqs, prim_seq)
  
  if(!is.na(ref_seq)){
    seqs<-append(seqs, ref_seq)
  }
  
  al<- seqs %>% RemoveGaps() %>% AlignSeqs(., verbose=FALSE)
  BrowseSeqs(al)
  
  #, primers$rev_pot_mismatch, primers$probe_pot_mismatch)
}


