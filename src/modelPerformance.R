# Title        : Simulation study
# Objective    : Compute the bias and coefficient of variation of the estimates,
#                and the estimates of prevalence
# Created by   : Christian Tsoungui Obama
# Created on   : 05.05.22
# Last modified: 05.01.2023

# Importing libraries
library(wordspace)

# Functions
psi <- function (inp){
  inp / (1 - exp(-inp))
}

#################################
# Function hapl(arch) takes as input a vector of the number of alleles per locu and outputs
#  a nm x 2 matrix of all possible haplotypes in geadic representation
#################################
hapl <- function(arch){
  H <- array(0,c(prod(arch),2))
  H[,1] <- rep(0:(arch[1]-1), each=arch[2])
  H[,2] <- rep(0:(arch[2]-1), arch[1])
  H
}

gen_func <- function(x, lambd){
  (exp(x*lambd)-1)/(exp(lambd) - 1)
}

setUh <- function(haplo, cardUh, arch){
  uh  <- t(array(rep(haplo, cardUh), dim = c(2, cardUh)))
  idx <- which(arch==max(arch))
  if(length(idx) > 1){
    idx <- idx[1]
  }
  v1 <- 0:(arch[idx]-1)
  v  <- v1[-which(v1==uh[1,idx])]
  uh[2:(length(v)+1), idx] <- v
  idx  <- which(1:2 !=idx)
  idx2 <- length(v)+2
  v2   <- 0:(arch[idx]-1)
  v    <- v2[-which(v2==uh[1,idx])]
  uh[idx2:cardUh, idx] <- v
  uh
}

bias     <- function(sim_Param, name){#(estim_Param, sim_Param, name){
  # This function implements bias of mean MOI and haplotype frequencies as defined in the manuscript
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from Microatellite data"

  moibiasloc <- vector(mode = "list", length = n_Sim_Loci)
  freqbiasloc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    estim_Param        <- readRDS(paste0(path, "dataset/full_modelEstimates_", l, name, ".rds"))
    freqbiasSamp <- vector(mode = "list", length = n_Sampl)
    moibiasSamp  <- vector(mode = "list", length = n_Sampl)

    # Number of haplotypes
    num_Hapl   <- n_Hapl[l]
    num_Hapl_PlusOne <- num_Hapl + 1

    for (k in 1:n_Sampl){  # For each sample size
      freqbias_lamb <- vector(mode = "list", length = n_Freq_Distr)
      moibias_lamb  <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        freq_bias <- vector(mode = "list", length = n_Freq_Distr)
        moi_bias  <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          # Access each of the 10000 estimates
          tmp1 <- estim_Param[[k]][[j]][, , i] #[[l]][[k]][[j]][, , i]

          # True frequencies
          tmp2 <- t(sim_Param[[1]][[l]])
          tmp3 <- tmp2[, i]

          # Relative bias of haplotype freq. in percent
          freqbias <- rowMeans((tmp1[2:num_Hapl_PlusOne,]/tmp3 - 1), na.rm = TRUE) * 100

          # Remove the estimates of Lambda that yield an "Infinite mean MOI"
          mean_moi_estim <- psi(tmp1[1, ])
          
          # Relative bias of MOI in percent
          bias_moi <- (mean_moi_estim/true_Mean_MOI[j] - 1)
          bias_moi <- mean(bias_moi[!is.infinite(bias_moi)], na.rm = TRUE) * 100

          # Save frequency and MOI bias in lists
          freq_bias[[i]] <- freqbias
          moi_bias[[i]]  <- bias_moi
        }
        freqbias_lamb[[j]] <- freq_bias
        moibias_lamb[[j]]  <- moi_bias
      }
      freqbiasSamp[[k]] <- freqbias_lamb
      moibiasSamp[[k]]  <- moibias_lamb
    }
    freqbiasloc[[l]] <- freqbiasSamp
    moibiasloc[[l]]  <- moibiasSamp
  }

  # Saving the estimates
  saveRDS(freqbiasloc, file = paste0(path, "dataset/freqbias",name ,".rds"))
  saveRDS(moibiasloc, file = paste0(path, "dataset/moibias",name ,".rds"))
}

coefvar  <- function(sim_Param, name){#(estim_Param, sim_Param, name){
  # This function implements coefficient of variation of mean MOI as defined in the manuscript
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from Microatellite data"

  moicvloc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    estim_Param        <- readRDS(paste0(path, "dataset/full_modelEstimates_", l, name, ".rds"))
    moicvSamp  <- vector(mode = "list", length = n_Sampl)

    # Number of haplotypes
    num_Hapl   <- n_Hapl[l]
    num_Hapl_PlusOne <- num_Hapl + 1

    for (k in 1:n_Sampl){  # For each sample size
      moicv_lamb  <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        moi_cv  <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          # Access each of the 10000 estimates
          tmp1 <- estim_Param[[k]][[j]][, , i]#estim_Param[[l]][[k]][[j]][, , i]

          # Remove the estimates of Lambda that yield an "Infinite mean MOI"
          mean_moi_estim <- psi(tmp1[1, ])
          mean_moi_estim <- mean_moi_estim[!is.infinite(mean_moi_estim)]
          
          # Dimensionless coefficient of variation of MOI
          moicv <- sd(mean_moi_estim, na.rm = TRUE) / true_Mean_MOI[j]

          # Save moi coefficient of variation in list
          moi_cv[[i]]  <- moicv
        }
        moicv_lamb[[j]]  <- moi_cv
      }
      moicvSamp[[k]]  <- moicv_lamb
    }
    moicvloc[[l]]  <- moicvSamp
  }

  # Saving the estimates
  saveRDS(moicvloc, file = paste0(path, "dataset/moicv", name, ".rds"))
}

coefvar_LD  <- function(estimLD, true_LD, ldtype, jj, name){
  # This function implements coefficient of variation of Linkage disequilibrium adhoc measures
  #  as defined in the manuscript "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from Microatellite data"

  moicvloc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    moicvSamp  <- vector(mode = "list", length = n_Sampl)

    for (k in 1:n_Sampl){  # For each sample size
      moicv_lamb  <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        moi_cv  <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          # Access each of the 10000 estimates
          if (ldtype %in% c('D', 'Q', 'r')){
            tmp1   <- estimLD[[l]][[k]][[j]][[i]] 
          }else{
           tmp1   <- estimLD[[l]][[k]][[j]][, , i] 
          }
          
          # Dimensionless coefficient of variation of MOI
          moicv <- sd(tmp1, na.rm = TRUE)

          # Save moi coefficient of variation in list
          moi_cv[[i]]  <- moicv
        }
        moicv_lamb[[j]]  <- moi_cv
      }
      moicvSamp[[k]]  <- moicv_lamb
    }
    moicvloc[[l]]  <- moicvSamp
  }

  # Saving the estimates
  saveRDS(moicvloc, file = paste0(path, "dataset/LDcv_", ldtype, "_", jj, name, ".rds"))
}

true_LD <- function(haplfreq, gen, name){
  # This function implements the true linkage disequilibrium measures (D', r^2) as defined in the manuscript of Tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from STR data"

  D_loc <- vector(mode = "list", length = n_Sim_Loci)
  r_loc <- vector(mode = "list", length = n_Sim_Loci)
  Q_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (i in 1:n_Sim_Loci){
    D_loc[[i]] <- vector(mode = "list", length = n_Freq_Distr)
    r_loc[[i]] <- vector(mode = "list", length = n_Freq_Distr)
    Q_loc[[i]] <- vector(mode = "list", length = n_Freq_Distr)
    arch <- gen[i,]
  
    for (j in 1:n_Freq_Distr){
      freq <- haplfreq[[i]][j,]
      pmat  <- matrix(freq, arch, byrow = TRUE)
      Afreq <- rowSums(pmat)
      Bfreq <- colSums(pmat)

      pij   <- Afreq%*%t(Bfreq)        # pipj
      picjc <- (1-Afreq)%*%t(1-Bfreq)  # (1-pi)(1-pj)
      picj  <- (1-Afreq)%*%t(Bfreq)    # (1-pi)pj
      pijc  <- Afreq%*%t(1-Bfreq)      # pi(1-pj)

      D     <- round(pmat - pij, 3)    # Dij
      pick  <- D > 0

      D_pos       <- array(0, arch)
      D_pos[pick] <- D[pick]

      D_neg        <- array(0, arch)
      D_neg[!pick] <- D[!pick]

      Dmax_pos  <- D
      Dmax_neg  <- D

      Dmax_pos  <- do.call(pmin, list(picj, pijc))  # if Dij > 0
      Dmax_neg  <- do.call(pmin, list(pij, picjc))  # if Dij < 0

      pij_pos       <- array(0, arch)
      pij_pos[pick] <- pij[pick]

      pij_neg        <- array(0, arch)
      pij_neg[!pick] <- pij[!pick]
    
      D_loc[[i]][[j]]   <- sum(pij_pos*D_pos/Dmax_pos) + sum(pij_neg*abs(D_neg)/Dmax_neg)       # D'
      r_loc[[i]][[j]]  <- sum(D^2)/((1-sum(Afreq**2))*(1-sum(Bfreq**2))) # D* # tmp_sum/min(arch-1)     # r^2
      Q_loc[[i]][[j]]  <- sum(D^2/pij)/prod(arch-1)    # Q*
    }
  }
  # Save LD estimates
  saveRDS(D_loc, file = paste0(path, "dataset/true_LD_D", name, ".rds"))
  saveRDS(r_loc, file = paste0(path, "dataset/true_LD_r", name, ".rds"))
  saveRDS(Q_loc, file = paste0(path, "dataset/true_LD_Q", name, ".rds"))
}

estim_LD <- function(gen, name){
  # This function implements the true linkage disequilibrium measures (D', r^2, Q*) as defined in the manuscript of Tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from STR data"

  Dloc <- vector(mode = "list", length = n_Sim_Loci)
  rloc <- vector(mode = "list", length = n_Sim_Loci)
  Qloc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){
    estim    <- readRDS(paste0(path, "dataset/full_modelEstimates_", l, name, ".rds"))
    DSamp    <- vector(mode = "list", length = n_Sampl)
    rSamp   <- vector(mode = "list", length = n_Sampl)
    QSamp   <- vector(mode = "list", length = n_Sampl)
    arch     <- gen[l,]
    nhapl    <- prod(arch)
    nhaplo   <- nhapl+1

    for (k in 1:n_Sampl){  # For each true sample size
      Dlamb <- vector(mode = "list", length = n_Freq_Distr)
      rlamb <- vector(mode = "list", length = n_Freq_Distr)
      Qlamb <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        Dfreq  <- vector(mode = "list", length = n_Freq_Distr)
        rfreq <- vector(mode = "list", length = n_Freq_Distr)
        Qfreq <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          Dp  <- matrix(0, ncol = n_Sampl_Gen)
          r <- matrix(0, ncol = n_Sampl_Gen)
          Q <- matrix(0, ncol = n_Sampl_Gen)

          ## Access each of the estimates frequencies
          estfreq <- estim[[k]][[j]][, , i]
          estfreq <- estfreq[2:nhaplo,]

          for(p in 1:n_Sampl_Gen){
            freq  <- estfreq[,p]
            pmat  <- matrix(freq, arch)
            Afreq <- rowSums(pmat)
            Bfreq <- colSums(pmat)

            pij   <- Afreq%*%t(Bfreq)        # pipj
            picjc <- (1-Afreq)%*%t(1-Bfreq)  # (1-pi)(1-pj)
            picj  <- (1-Afreq)%*%t(Bfreq)    # (1-pi)pj
            pijc  <- Afreq%*%t(1-Bfreq)      # pi(1-pj)

            D     <- round(freq - pij, 3)               # Dij
            pick  <- D > 0

            D_pos       <- array(0, arch)
            D_pos[pick] <- D[pick]

            D_neg        <- array(0, arch)
            D_neg[!pick] <- D[!pick]

            Dmax_pos  <- D
            Dmax_neg  <- D

            Dmax_pos  <- do.call(pmin, list(picj, pijc))  # if Dij > 0
            Dmax_neg  <- do.call(pmin, list(pij, picjc))  # if Dij < 0

            pij_pos       <- array(0, arch)
            pij_pos[pick] <- pij[pick]

            pij_neg        <- array(0, arch)
            pij_neg[!pick] <- pij[!pick]

            Dp[,p]   <- sum(pij_pos*D_pos/Dmax_pos) + sum(pij_neg*abs(D_neg)/Dmax_neg)       # D'
            r[,p]    <- sum(D^2)/((1-sum(Afreq**2))*(1-sum(Bfreq**2))) # D* # tmp_sum/min(arch-1)     # r^2
            Q[,p]    <- sum(D^2/pij)/prod(arch-1)    # Q*
          }
          Dfreq[[i]] <- rowMeans(Dp, na.rm = TRUE) #Dp# 
          rfreq[[i]] <- rowMeans(r, na.rm = TRUE) #r# 
          Qfreq[[i]] <- rowMeans(Q, na.rm = TRUE) #Q# 
        }
        Dlamb[[j]]  <- Dfreq
        rlamb[[j]] <- rfreq
        Qlamb[[j]] <- Qfreq
      }
      DSamp[[k]]  <- Dlamb
      rSamp[[k]] <- rlamb
      QSamp[[k]] <- Qlamb
    }
    Dloc[[l]]  <- DSamp
    rloc[[l]] <- rSamp
    Qloc[[l]] <- QSamp
  }
  # Save LD estimates
  saveRDS(Dloc,  file = paste0(path, "dataset/estim_LD_D", name, ".rds"))
  saveRDS(rloc, file = paste0(path, "dataset/estim_LD_r", name, ".rds"))
  saveRDS(Qloc, file = paste0(path, "dataset/estim_LD_Q", name, ".rds"))
}

estim_LDtot <- function(gen, name){
  # This function implements the true linkage disequilibrium measures (D', r^2) as defined in the manuscript of Tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from STR data"

  Dloc <- vector(mode = "list", length = n_Sim_Loci)
  rloc <- vector(mode = "list", length = n_Sim_Loci)
  Qloc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){
    estim    <- readRDS(paste0(path, "dataset/full_modelEstimates_", l, name, ".rds"))
    DSamp    <- vector(mode = "list", length = n_Sampl)
    rSamp   <- vector(mode = "list", length = n_Sampl)
    QSamp   <- vector(mode = "list", length = n_Sampl)
    arch     <- gen[l,]
    nhapl    <- prod(arch)
    nhaplo   <- nhapl+1

    for (k in 1:n_Sampl){  # For each true sample size
      Dlamb <- vector(mode = "list", length = n_Freq_Distr)
      rlamb <- vector(mode = "list", length = n_Freq_Distr)
      Qlamb <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        Dfreq  <- vector(mode = "list", length = n_Freq_Distr)
        rfreq <- vector(mode = "list", length = n_Freq_Distr)
        Qfreq <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          Dp  <- matrix(0, ncol = n_Sampl_Gen)
          r <- matrix(0, ncol = n_Sampl_Gen)
          Q <- matrix(0, ncol = n_Sampl_Gen)

          ## Access each of the estimates frequencies
          estfreq <- estim[[k]][[j]][, , i]#estim[[l]][[k]][[j]][, , i]
          estfreq <- estfreq[2:nhaplo,]

          for(p in 1:n_Sampl_Gen){
            freq  <- estfreq[,p]
            pmat  <- matrix(freq, arch, byrow = TRUE)
            Afreq <- rowSums(pmat)
            Bfreq <- colSums(pmat)

            pij   <- Afreq%*%t(Bfreq)        # pipj
            picjc <- (1-Afreq)%*%t(1-Bfreq)  # (1-pi)(1-pj)
            picj  <- (1-Afreq)%*%t(Bfreq)    # (1-pi)pj
            pijc  <- Afreq%*%t(1-Bfreq)      # pi(1-pj)

            D     <- round(freq - pij, 3)               # Dij
            pick  <- D > 0

            D_pos       <- array(0, arch)
            D_pos[pick] <- D[pick]

            D_neg        <- array(0, arch)
            D_neg[!pick] <- D[!pick]

            Dmax_pos  <- D
            Dmax_neg  <- D

            Dmax_pos  <- do.call(pmin, list(picj, pijc))  # if Dij > 0
            Dmax_neg  <- do.call(pmin, list(pij, picjc))  # if Dij < 0

            pij_pos       <- array(0, arch)
            pij_pos[pick] <- pij[pick]

            pij_neg        <- array(0, arch)
            pij_neg[!pick] <- pij[!pick]

            Dp[,p]   <- sum(pij_pos*D_pos/Dmax_pos) + sum(pij_neg*abs(D_neg)/Dmax_neg)       # D'
            r[,p]    <- sum(D^2)/((1-sum(Afreq**2))*(1-sum(Bfreq**2))) # D* # tmp_sum/min(arch-1)     # r^2
            Q[,p]    <- sum(D^2/pij)/prod(arch-1)    # Q*
          }
          Dfreq[[i]]  <- Dp
          rfreq[[i]] <- r
          Qfreq[[i]] <- Q
        }
        Dlamb[[j]]  <- Dfreq
        rlamb[[j]] <- rfreq
        Qlamb[[j]] <- Qfreq
      }
      DSamp[[k]]  <- Dlamb
      rSamp[[k]] <- rlamb
      QSamp[[k]] <- Qlamb
    }
    Dloc[[l]]  <- DSamp
    rloc[[l]] <- rSamp
    Qloc[[l]] <- QSamp
  }

  # Save LD estimates
  saveRDS(Dloc,  file = paste0(path, "dataset/estim_LD_Dtot", name, ".rds"))
  saveRDS(rloc, file = paste0(path, "dataset/estim_LD_rtot", name, ".rds"))
  saveRDS(Qloc, file = paste0(path, "dataset/estim_LD_Qtot", name, ".rds"))
}

main <- function(sim_Param, reshap_Sim_Param, name, gen){

  # Loading estimates using the ad-hoc model and for each of the types of frequencies estimates
  adhoc_estim_Param1 <- list(readRDS(paste0(path, "dataset/full_adhocModelEstimates1_1", name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelEstimates1_2", name, ".rds")))
  adhoc_estim_Param2 <- list(readRDS(paste0(path, "dataset/full_adhocModelEstimates2_1", name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelEstimates2_2", name, ".rds")))
  adhoc_estim_Param3 <- list(readRDS(paste0(path, "dataset/full_adhocModelEstimates3_1", name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelEstimates3_2", name, ".rds")))

  # Loading estimates of Linkage Disequilibrium for each type of frequency estimates

  adhoc_Dptmp1 <- list(readRDS(paste0(path, "dataset/full_adhocModelDpEstimates1_1",  name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelDpEstimates1_2",  name, ".rds")))
  adhoc_Dptmp2 <- list(readRDS(paste0(path, "dataset/full_adhocModelDpEstimates2_1",  name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelDpEstimates2_2",  name, ".rds")))
  
  adhoc_rtmp1  <- list(readRDS(paste0(path, "dataset/full_adhocModelrEstimates1_1",  name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelrEstimates1_2",  name, ".rds")))
  adhoc_rtmp2  <- list(readRDS(paste0(path, "dataset/full_adhocModelrEstimates2_1",  name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelrEstimates2_2",  name, ".rds")))
  
  adhoc_Qtmp1  <- list(readRDS(paste0(path, "dataset/full_adhocModelQEstimates1_1",  name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelQEstimates1_2",  name, ".rds")))
  adhoc_Qtmp2  <- list(readRDS(paste0(path, "dataset/full_adhocModelQEstimates2_1",  name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelQEstimates2_2",  name, ".rds")))
  
  #adhoc_r_ld1  <- list(readRDS(paste0(path, "dataset/full_adhocModelrEstimates1_1",   name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelrEstimates1_2",   name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelrEstimates1_3",   name, ".rds")))
  #adhoc_Q_ld1  <- list(readRDS(paste0(path, "dataset/full_adhocModelQEstimates1_1",   name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelQEstimates1_2",   name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelQEstimates1_3",   name, ".rds")))

  #adhoc_D_ld2  <- list(readRDS(paste0(path, "dataset/full_adhocModelDpEstimates2_1",  name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelDpEstimates2_2",  name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelDpEstimates2_3",  name, ".rds")))
  #adhoc_r_ld2  <- list(readRDS(paste0(path, "dataset/full_adhocModelrEstimates2_1",   name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelrEstimates2_2",   name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelrEstimates2_3",   name, ".rds")))
  #adhoc_Q_ld2  <- list(readRDS(paste0(path, "dataset/full_adhocModelQEstimates2_1",   name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelQEstimates2_2",   name, ".rds")), readRDS(paste0(path, "dataset/full_adhocModelQEstimates2_3",   name, ".rds")))

  # Bias of frequencies and MOI
  bias(sim_Param, name)

  # Coefficient of variation of MOI
  coefvar(sim_Param, name)
  
  # True linkage disequilibrium measure
  true_LD(sim_Param[[1]], gen, name)

  # Estimated linkage disequilibrium measure
  estim_LD(gen, name)

  # Merging adhoc estimates of each LD measure type for different frequencies estimates
  adhoc_estim_Dp <- list(adhoc_Dptmp1, adhoc_Dptmp2, readRDS(paste0(path, "dataset/estim_LD_D", name, ".rds")))
  adhoc_estim_R  <- list(adhoc_rtmp1, adhoc_rtmp2, readRDS(paste0(path, "dataset/estim_LD_r", name, ".rds")))
  adhoc_estim_Q  <- list(adhoc_Qtmp1, adhoc_Qtmp2, readRDS(paste0(path, "dataset/estim_LD_Q", name, ".rds")))
  adhoc_LD_estim <- list(adhoc_estim_Dp, adhoc_estim_R, adhoc_estim_Q)

  # Coefficient of variation of LD
  vec <- c('D', 'r', 'Q')

  for(i in 1:3){
    trueLD <- readRDS(paste0(path, "dataset/true_LD_", vec[i], ".rds")) 
    for(j in 1:3){
      if(j != 3){
        ad <- 'adhoc'
      }else{
        ad <- ''
      }
      estimLD <- adhoc_LD_estim[[i]][[j]]
      coefvar_LD(estimLD, trueLD, paste0(ad,vec[i]), j, name)
    }
  }

  # Estimated ambiguous prevalence
  # estim_amb_prevalence(true_Amb_Prev, name)

  # # Estimated conditional prevalence
  # estim_conditional_prevalence(sim_Param, true_Conditional_Prev, name, gen)

  # # Estimated relative prevalence (adhoc Model)
  # for(idx in 1:3){
  #   estim_relative_prevalence(name, idx)
  # }
  #estim_relative_prevalence(adhoc_estim_Param1, name, 1)
  #estim_relative_prevalence(adhoc_estim_Param2, name, 2)
  #estim_relative_prevalence(adhoc_estim_Param3, name, 3)


}
 
# Relative path
path <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/MultiAllelicBilociModel/"

name <- c('')
print(paste0('Ongoing simulation for: ', name, ' data!'))

# Loading true haplotype parameters for the simulation
sim_Param <- readRDS(paste0(path, "dataset/full_true_Parameters", name, ".rds"))

# Loading extra parameters
extra_Sim_Param  <- readRDS(paste0(path, "dataset/full_extra_Parameters", name, ".rds"))

# Simulation parameters
n_Lbda        <- extra_Sim_Param[[1]]
n_Sim_Loci    <- extra_Sim_Param[[2]]
n_Hapl        <- extra_Sim_Param[[3]]
n_Sampl       <- extra_Sim_Param[[4]]
n_Sampl_Gen   <- extra_Sim_Param[[5]]
n_Freq_Distr  <- extra_Sim_Param[[6]]
genArch       <- extra_Sim_Param[[7]]
true_Mean_MOI <- psi(sim_Param[[2]])

# Reformatting true parameters to compute true prevalence
reshap_Sim_Param <- vector(mode='list', length=n_Sim_Loci)

for (i in 1:n_Sim_Loci){
  numb_Row <- n_Hapl[i]+1
  reshap_Sim_Param[[i]] <- vector(mode='list', length=n_Freq_Distr)

  for (j in 1:n_Freq_Distr){
    reshap_Sim_Param[[i]][[j]]              <- array(0, c(numb_Row, n_Lbda))
    reshap_Sim_Param[[i]][[j]][1,]          <- sim_Param[[2]]
    reshap_Sim_Param[[i]][[j]][2:numb_Row,] <- sim_Param[[1]][[i]][j,]
  }
}

main(sim_Param, reshap_Sim_Param, name, genArch)
