#Note: If using an AIRR formatted file or output file from pipeline v1.11 or later, change every instance of "Chain_Type" to "locus", "Molecule_Count" to "duplicate_count" and ALL locus / chain calls.
#TCR locus / chain calls must be modified as follows: "TCR_Alpha" -> "TRA", "TCR_Beta" -> "TRB", "TCR_Delta" -> "TRD", "TCR_Gamma" -> "TRG"
#BCR locus / chain calls must be modified as follows: "BCR_Heavy" -> "IGH", "BCR_Kappa" -> "IGK", "BCR_Lambda" -> "IGL"

demo <- read.table(file = "~/gitrepos/rhapsody-analysis-notebook/data/vdj/copied-FASTQs-expected-cells-30K-v2-rerun-MC_VDJ_Dominant_Contigs_AIRR.tsv", sep = ",", header = T)
table(demo$Chain_Type)
dat <- demo
new_df <- as.data.frame(matrix(nrow = nrow(dat), ncol = ncol(dat) ))
colnames(new_df) <- colnames(dat)

#Enter proof ratio: Enter X where X:1 is the ratio of one cell type (AB or GD) to overcome the call of another in conflicted locus calls
proof_ratio <- 2

`%notin%` <- Negate(`%in%`)
newline_counter <- 1
for(i in 1:length(rownames(dat))){
  len <- length(which(dat[i,1] == dat[,1]))
  lines <- dat[which(dat[i,1] == dat[,1]),]
  
##Deal with case of T and B cell chains present in same cell:
  if(("TCR_Alpha" %in% lines$Chain_Type | "TCR_Beta" %in% lines$Chain_Type | "TCR_Delta" %in% lines$Chain_Type | "TCR_Gamma" %in% lines$Chain_Type) & 
     ("BCR_Heavy" %in% lines$Chain_Type | "BCR_Lambda" %in% lines$Chain_Type | "BCR_Kappa" %in% lines$Chain_Type)){
    
    TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
    TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
    TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
    TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
    
    BCR_Lambda_rec <- lines[which(lines$Chain_Type == "BCR_Lambda"),]
    BCR_Kappa_rec <- lines[which(lines$Chain_Type == "BCR_Kappa"),]
    BCR_Heavy_rec <- lines[which(lines$Chain_Type == "BCR_Heavy"),]
    
    TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
    if(identical(TCR_Alpha_count, integer(0))){TCR_Alpha_count <- 0}
    TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
    if(identical(TCR_Beta_count, integer(0))){TCR_Beta_count <- 0}
    TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
    if(identical(TCR_Delta_count, integer(0))){TCR_Delta_count <- 0}
    TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
    if(identical(TCR_Gamma_count, integer(0))){TCR_Gamma_count <- 0}
    
    BCR_Heavy_count <- BCR_Heavy_rec$Molecule_Count
    if(identical(BCR_Heavy_count, integer(0))){BCR_Heavy_count <- 0}
    BCR_Lambda_count <- BCR_Lambda_rec$Molecule_Count
    if(identical(BCR_Lambda_count, integer(0))){BCR_Lambda_count <- 0}
    BCR_Kappa_count <- BCR_Kappa_rec$Molecule_Count
    if(identical(BCR_Kappa_count, integer(0))){BCR_Kappa_count <- 0}
    
    ################################################################################################################################################
    #TCR Analysis if T counts outweigh B counts
    if((TCR_Alpha_count + TCR_Beta_count + TCR_Delta_count + TCR_Gamma_count) > 
       ((BCR_Heavy_count + BCR_Lambda_count + BCR_Kappa_count)*proof_ratio)){
      lines <- lines[-which(lines$Chain_Type == "BCR_Heavy"),]
      lines <- lines[-which(lines$Chain_Type == "BCR_Lambda"),]
      lines <- lines[-which(lines$Chain_Type == "BCR_Kappa"),]
      len_counter <- length(which(c(BCR_Heavy_count, BCR_Lambda_count, BCR_Kappa_count) == 0))
      
      len <- len - 3 + len_counter
      #ABGD
      if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Beta" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
        
        TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
        TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
        TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
        TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
        TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
        TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
        TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
        TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
        
        if((TCR_Alpha_count + TCR_Beta_count) > ((TCR_Gamma_count + TCR_Delta_count)*proof_ratio)){
          lines <- lines[-(which(lines$Chain_Type == "TCR_Gamma")),]
          lines <- lines[-(which(lines$Chain_Type == "TCR_Delta")),]
          len <- len - 2
        }
        else if ((TCR_Gamma_count + TCR_Delta_count) > ((TCR_Alpha_count + TCR_Beta_count)*proof_ratio)){
          lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
          lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
          len <- len - 2
        }
        else{
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
        }
      } 
      #ABG
      if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Beta" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type){
        TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
        TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
        TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
        
        TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
        TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
        TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
        
        if((TCR_Alpha_count + TCR_Beta_count) > ((TCR_Gamma_count)*proof_ratio)){
          lines <- lines[-(which(lines$Chain_Type == "TCR_Gamma")),]
          len <- len - 1
        }
        else if ((TCR_Gamma_count) > ((TCR_Alpha_count + TCR_Beta_count)*proof_ratio)){
          lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
          lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
          len <- len - 2
        }
        else{
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
        }
      }
      #ABD
      if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Beta" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
        TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
        TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
        TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
        
        TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
        TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
        TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
        
        if((TCR_Alpha_count + TCR_Beta_count) > ((TCR_Delta_count)*proof_ratio)){
          lines <- lines[-(which(lines$Chain_Type == "TCR_Delta")),]
          len <- len - 1
        }
        else if ((TCR_Delta_count) > ((TCR_Alpha_count + TCR_Beta_count)*proof_ratio)){
          lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
          lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
          len <- len - 2
        }
        else{
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
        }
      }
      #ADG
      if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
        TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
        TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
        TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
        
        TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
        TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
        TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
        
        if((TCR_Gamma_count + TCR_Delta_count) > ((TCR_Alpha_count)*proof_ratio)){
          lines <- lines[-(which(lines$Chain_Type == "TCR_Alpha")),]
          len <- len - 1
        }
        else if ((TCR_Alpha_count) > ((TCR_Gamma_count + TCR_Delta_count)*proof_ratio)){
          lines <- lines[-which(lines$Chain_Type == "TCR_Delta"),]
          lines <- lines[-which(lines$Chain_Type == "TCR_Gamma"),]
          len <- len - 2
        }
        else{
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
        }
      }
      #BDG
      if("TCR_Beta" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
        TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
        TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
        TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
        
        TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
        TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
        TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
        
        if((TCR_Gamma_count + TCR_Delta_count) > ((TCR_Beta_count)*proof_ratio)){
          lines <- lines[-(which(lines$Chain_Type == "TCR_Beta")),]
          len <- len - 1
        }
        else if ((TCR_Beta_count) > ((TCR_Gamma_count + TCR_Delta_count)*proof_ratio)){
          lines <- lines[-which(lines$Chain_Type == "TCR_Gamma"),]
          lines <- lines[-which(lines$Chain_Type == "TCR_Delta"),]
          len <- len - 1
        }
        else{
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
        }
      }
      #AD
      if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
        TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
        TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
        
        TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
        TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
        
        if((TCR_Alpha_count) > ((TCR_Delta_count)*proof_ratio)){
          lines <- lines[-(which(lines$Chain_Type == "TCR_Delta")),]
          len <- len - 1
        }
        else if ((TCR_Delta_count) > ((TCR_Alpha_count)*proof_ratio)){
          lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
          len <- len - 1
        }
        else{
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
        }
      }
      #AG
      if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type){
        TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
        TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
        
        TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
        TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
        
        if((TCR_Alpha_count) > ((TCR_Gamma_count)*proof_ratio)){
          lines <- lines[-(which(lines$Chain_Type == "TCR_Gamma")),]
          len <- len - 1
        }
        else if ((TCR_Gamma_count) > ((TCR_Alpha_count)*proof_ratio)){
          lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
          len <- len - 1
        }
        else{
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
        }
      }
      #BG
      if("TCR_Beta" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type){
        TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
        TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
        
        TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
        TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
        
        if((TCR_Beta_count) > ((TCR_Gamma_count)*proof_ratio)){
          lines <- lines[-(which(lines$Chain_Type == "TCR_Gamma")),]
          len <- len - 1
        }
        else if ((TCR_Gamma_count) > ((TCR_Beta_count)*proof_ratio)){
          lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
          len <- len - 1
        }
        else{
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
        }
      }
      #BD
      if("TCR_Beta" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
        TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
        TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
        
        TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
        TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
        
        if((TCR_Beta_count) > ((TCR_Delta_count)*proof_ratio)){
          lines <- lines[-(which(lines$Chain_Type == "TCR_Delta")),]
          len <- len - 1
        }
        else if ((TCR_Delta_count) > ((TCR_Beta_count)*proof_ratio)){
          lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
          len <- len - 1
        }
        else{
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
        }
      }
    }
    ###########################################################################################################################################
    #Perform BCR Analysis if B counts outweigh T counts
    else if((BCR_Heavy_count + BCR_Lambda_count + BCR_Kappa_count) > 
       ((TCR_Alpha_count + TCR_Beta_count + TCR_Delta_count + TCR_Gamma_count)*proof_ratio)){
      lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
      lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
      lines <- lines[-which(lines$Chain_Type == "TCR_Delta"),]
      lines <- lines[-which(lines$Chain_Type == "TCR_Gamma"),]
      
      len_counter <- length(which(c(TCR_Alpha_count, TCR_Beta_count, TCR_Delta_count, TCR_Gamma_count) == 0))
      len <- len - 4 + len_counter
      #INSERT B CELL ANALYSIS BELOW
      if("BCR_Lambda" %in% lines$Chain_Type & "BCR_Kappa" %in% lines$Chain_Type){
        
        BCR_Lambda_rec <- lines[which(lines$Chain_Type == "BCR_Lambda"),]
        BCR_Kappa_rec <- lines[which(lines$Chain_Type == "BCR_Kappa"),]
        BCR_Lambda_count <- BCR_Lambda_rec$Molecule_Count
        BCR_Kappa_count <- BCR_Kappa_rec$Molecule_Count
        
        if(BCR_Lambda_count > (BCR_Kappa_count*proof_ratio)){
          lines <- lines[-(which(lines$Chain_Type == "BCR_Kappa")),]
          len <- len - 1
        }
        else if (BCR_Kappa_count > (BCR_Lambda_count*proof_ratio)){
          lines <- lines[-which(lines$Chain_Type == "BCR_Lambda"),]
          len <- len - 1
        }
        else{
          lines$Chain_Type[which(lines$Chain_Type == "BCR_Kappa")] <- "BCR_Kappa_UNDETERMINED"
          lines$Chain_Type[which(lines$Chain_Type == "BCR_Lambda")] <- "BCR_Lambda_UNDETERMINED"
        }
      }
    
    }
    #If no clear winner between BCR and TCR, change everything to undetermined
    else{
      lines$Chain_Type[which(lines$Chain_Type == "BCR_Kappa")] <- "BCR_Kappa_UNDETERMINED"
      lines$Chain_Type[which(lines$Chain_Type == "BCR_Lambda")] <- "BCR_Lambda_UNDETERMINED"
      lines$Chain_Type[which(lines$Chain_Type == "BCR_Heavy")] <- "BCR_Heavy_UNDETERMINED"
      lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
      lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
      lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
      lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
      
    }
  }
  
  #Deal with only TCR chains present in same cell
  else if(("TCR_Alpha" %in% lines$Chain_Type | "TCR_Beta" %in% lines$Chain_Type | "TCR_Delta" %in% lines$Chain_Type | "TCR_Gamma" %in% lines$Chain_Type)){
    
    #ABGD
    if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Beta" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
      
      TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
      TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
      TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
      TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
      TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
      TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
      TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
      TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
      
      if((TCR_Alpha_count + TCR_Beta_count) > ((TCR_Gamma_count + TCR_Delta_count)*proof_ratio)){
        lines <- lines[-(which(lines$Chain_Type == "TCR_Gamma")),]
        lines <- lines[-(which(lines$Chain_Type == "TCR_Delta")),]
        len <- len - 2
      }
      else if ((TCR_Gamma_count + TCR_Delta_count) > ((TCR_Alpha_count + TCR_Beta_count)*proof_ratio)){
        lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
        lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
        len <- len - 2
      }
      else{
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
      }
    } 
    #ABG
    if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Beta" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type){
      TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
      TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
      TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
      
      TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
      TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
      TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
      
      if((TCR_Alpha_count + TCR_Beta_count) > ((TCR_Gamma_count)*proof_ratio)){
        lines <- lines[-(which(lines$Chain_Type == "TCR_Gamma")),]
        len <- len - 1
      }
      else if ((TCR_Gamma_count) > ((TCR_Alpha_count + TCR_Beta_count)*proof_ratio)){
        lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
        lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
        len <- len - 2
      }
      else{
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
      }
    }
    #ABD
    if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Beta" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
      TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
      TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
      TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
      
      TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
      TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
      TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
      
      if((TCR_Alpha_count + TCR_Beta_count) > ((TCR_Delta_count)*proof_ratio)){
        lines <- lines[-(which(lines$Chain_Type == "TCR_Delta")),]
        len <- len - 1
      }
      else if ((TCR_Delta_count) > ((TCR_Alpha_count + TCR_Beta_count)*proof_ratio)){
        lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
        lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
        len <- len - 2
      }
      else{
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
      }
    }
    #ADG
    if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
      TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
      TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
      TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
      
      TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
      TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
      TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
      
      if((TCR_Gamma_count + TCR_Delta_count) > ((TCR_Alpha_count)*proof_ratio)){
        lines <- lines[-(which(lines$Chain_Type == "TCR_Alpha")),]
        len <- len - 1
      }
      else if ((TCR_Alpha_count) > ((TCR_Gamma_count + TCR_Delta_count)*proof_ratio)){
        lines <- lines[-which(lines$Chain_Type == "TCR_Delta"),]
        lines <- lines[-which(lines$Chain_Type == "TCR_Gamma"),]
        len <- len - 2
      }
      else{
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
      }
    }
    #BDG
    if("TCR_Beta" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
      TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
      TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
      TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
      
      TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
      TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
      TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
      
      if((TCR_Gamma_count + TCR_Delta_count) > ((TCR_Beta_count)*proof_ratio)){
        lines <- lines[-(which(lines$Chain_Type == "TCR_Beta")),]
        len <- len - 1
      }
      else if ((TCR_Beta_count) > ((TCR_Gamma_count + TCR_Delta_count)*proof_ratio)){
        lines <- lines[-which(lines$Chain_Type == "TCR_Gamma"),]
        lines <- lines[-which(lines$Chain_Type == "TCR_Delta"),]
        len <- len - 1
      }
      else{
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
      }
    }
    #AD
    if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
      TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
      TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
      
      TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
      TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
      
      if((TCR_Alpha_count) > ((TCR_Delta_count)*proof_ratio)){
        lines <- lines[-(which(lines$Chain_Type == "TCR_Delta")),]
        len <- len - 1
      }
      else if ((TCR_Delta_count) > ((TCR_Alpha_count)*proof_ratio)){
        lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
        len <- len - 1
      }
      else{
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
      }
    }
    #AG
    if("TCR_Alpha" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type){
      TCR_Alpha_rec <- lines[which(lines$Chain_Type == "TCR_Alpha"),]
      TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
      
      TCR_Alpha_count <- TCR_Alpha_rec$Molecule_Count
      TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
      
      if((TCR_Alpha_count) > ((TCR_Gamma_count)*proof_ratio)){
        lines <- lines[-(which(lines$Chain_Type == "TCR_Gamma")),]
        len <- len - 1
      }
      else if ((TCR_Gamma_count) > ((TCR_Alpha_count)*proof_ratio)){
        lines <- lines[-which(lines$Chain_Type == "TCR_Alpha"),]
        len <- len - 1
      }
      else{
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Alpha")] <- "TCR_Alpha_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
      }
    }
    #BG
    if("TCR_Beta" %in% lines$Chain_Type & "TCR_Gamma" %in% lines$Chain_Type){
      TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
      TCR_Gamma_rec <- lines[which(lines$Chain_Type == "TCR_Gamma"),]
      
      TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
      TCR_Gamma_count <- TCR_Gamma_rec$Molecule_Count
      
      if((TCR_Beta_count) > ((TCR_Gamma_count)*proof_ratio)){
        lines <- lines[-(which(lines$Chain_Type == "TCR_Gamma")),]
        len <- len - 1
      }
      else if ((TCR_Gamma_count) > ((TCR_Beta_count)*proof_ratio)){
        lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
        len <- len - 1
      }
      else{
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Gamma")] <- "TCR_Gamma_UNDETERMINED"
      }
    }
    #BD
    if("TCR_Beta" %in% lines$Chain_Type & "TCR_Delta" %in% lines$Chain_Type){
      TCR_Beta_rec <- lines[which(lines$Chain_Type == "TCR_Beta"),]
      TCR_Delta_rec <- lines[which(lines$Chain_Type == "TCR_Delta"),]
      
      TCR_Beta_count <- TCR_Beta_rec$Molecule_Count
      TCR_Delta_count <- TCR_Delta_rec$Molecule_Count
      
      if((TCR_Beta_count) > ((TCR_Delta_count)*proof_ratio)){
        lines <- lines[-(which(lines$Chain_Type == "TCR_Delta")),]
        len <- len - 1
      }
      else if ((TCR_Delta_count) > ((TCR_Beta_count)*proof_ratio)){
        lines <- lines[-which(lines$Chain_Type == "TCR_Beta"),]
        len <- len - 1
      }
      else{
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Beta")] <- "TCR_Beta_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "TCR_Delta")] <- "TCR_Delta_UNDETERMINED"
      }
      complete_line <- apply(lines[1:len,] , 2, paste , collapse = ";" )
      complete_line[1] <- dat[i,1]
      if(complete_line[1] %notin% new_df[,1]){
        new_df[newline_counter,] <- complete_line
        newline_counter <- newline_counter + 1
      }
    }
  }
 
  #Deal with only BCR chains present in same cell
  else if(("BCR_Heavy" %in% lines$Chain_Type | "BCR_Lambda" %in% lines$Chain_Type | "BCR_Kappa" %in% lines$Chain_Type)){
    if("BCR_Lambda" %in% lines$Chain_Type & "BCR_Kappa" %in% lines$Chain_Type){
      
      BCR_Lambda_rec <- lines[which(lines$Chain_Type == "BCR_Lambda"),]
      BCR_Kappa_rec <- lines[which(lines$Chain_Type == "BCR_Kappa"),]
      BCR_Lambda_count <- BCR_Lambda_rec$Molecule_Count
      BCR_Kappa_count <- BCR_Kappa_rec$Molecule_Count
      
      if(BCR_Lambda_count > (BCR_Kappa_count*proof_ratio)){
        lines <- lines[-(which(lines$Chain_Type == "BCR_Kappa")),]
        len <- len - 1
      }
      else if (BCR_Kappa_count > (BCR_Lambda_count*proof_ratio)){
        lines <- lines[-which(lines$Chain_Type == "BCR_Lambda"),]
        len <- len - 1
      }
      else{
        lines$Chain_Type[which(lines$Chain_Type == "BCR_Kappa")] <- "BCR_Kappa_UNDETERMINED"
        lines$Chain_Type[which(lines$Chain_Type == "BCR_Lambda")] <- "BCR_Lambda_UNDETERMINED"
      }
    }
  }
  
  complete_line <- apply(lines[1:len,] , 2, paste , collapse = ";" )
  complete_line[1] <- dat[i,1]
  if(complete_line[1] %notin% new_df[,1]){
    new_df[newline_counter,] <- complete_line
    newline_counter <- newline_counter + 1
  }
}
    
demo <- new_df[-c(grep("NA", new_df$Chain_Type)),]


vals <- as.data.frame(table(demo$Chain_Type))
vals$Percent <- vals$Freq/sum(vals$Freq)*100
vals <- vals[order(vals$Percent, decreasing = T),]
par(mar = c(12,5,2,1))
barplot(vals$Percent, names.arg = vals$Var1, las = 2, main = "BCR/TCR Chain Calls per Cell Demo -- Cleaned", ylab = "Percent of Total Cells", col = "seagreen", cex.names = 0.75)
