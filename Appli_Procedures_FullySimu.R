#Functions
appli_procedures_ind_simu <- function(N_iteration, pvalues_list, covariates_list, alpha, m1, subset_m1_1, subset_m1_2, subset_m1_3,limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3){
  # Allows to apply the procedures on a independent simulation dataset (with 500 iterations) and to compute the FDR and the powers (global and subgroup)

  all_FDP <- c()
  all_Power <- c()
  all_SubPower <- c()

  for(i in c(1:N_iteration)){
    pvalues <- pvalues_list[[i]]
    covariates <- covariates_list[[i]]

    # wBHa procedure
    res_wBHa <- wBHa(pvalues, covariates, alpha, B_size, K)
    rejects_wBHa <- which(res_wBHa$adjusted_pvalues<=alpha)
    res_wBHa <- power_FDR_calc_ind_simu(rejects_wBHa, m1, subset_m1_1, subset_m1_2, subset_m1_3, 
                                        limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

    # BH procedure
    rejects_BH <- which(p.adjust(pvalues, method="BH")<=alpha)
    res_BH <- power_FDR_calc_ind_simu(rejects_BH, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

    # IHW procedure
    res_IHW <- ihw(pvalues~covariates, alpha=alpha)
    reject_IHW <- which(adj_pvalues(res_IHW)<alpha)
    res_IHW <- power_FDR_calc_ind_simu(reject_IHW, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

    # wBH without a_opt
    rejects_wBH <- which(p.adjust(pvalues/((length(pvalues)/sum(1/covariates))*(1/covariates)), method="BH")<=alpha)
    res_wBH <- power_FDR_calc_ind_simu(rejects_wBH, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

    # Qvalue procedure
    res_qvalue <- qvalue(pvalues)
    reject_qvalue <- which(res_qvalue$qvalues<alpha)
    res_qvalue <- power_FDR_calc_ind_simu(reject_qvalue, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

    # Swfdr procedure
    res_swfdr <- lm_qvalue(pvalues, covariates)
    reject_swfdr <- which(res_swfdr$qvalue<alpha)
    res_swfdr <- power_FDR_calc_ind_simu(reject_swfdr, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

    # FDRreg procedure
    pvalues2 <- pvalues
    pvalues2[pvalues2==1] <- (1-10^-7)
    zscores <- qnorm(pvalues2)
    res_fdrreg <- FDRreg(zscores, as.matrix(covariates))
    reject_fdrreg <- which(res_fdrreg$FDR<alpha)
    res_fdrreg <- power_FDR_calc_ind_simu(reject_fdrreg, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

    # CAMT procedure
    res_camt <- camt.fdr(pvals=pvalues, pi0.var=covariates)
    reject_camt <- which(c(res_camt$fdr<alpha))
    res_camt <- power_FDR_calc_ind_simu(reject_camt, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

    # Results
    all_FDP <- cbind(all_FDP, rbind(res_BH[[1]], res_wBH[[1]], res_wBHa[[1]], res_IHW[[1]], res_qvalue[[1]], res_swfdr[[1]], res_fdrreg[[1]], res_camt[[1]]))
    all_Power <- cbind(all_Power, rbind(res_BH[[2]],res_wBH[[2]],res_wBHa[[2]],res_IHW[[2]], res_qvalue[[2]],res_swfdr[[2]],res_fdrreg[[2]],res_camt[[2]]))
    all_SubPower <- rbind(all_SubPower,rbind(res_BH[[3]],res_wBH[[3]],res_wBHa[[3]],res_IHW[[3]],res_qvalue[[3]],res_swfdr[[3]],res_fdrreg[[3]],res_camt[[3]]))
  }

  # N/A or 0 management
  for(i in c(1:8)){
    all_FDP[i, which(is.na(all_FDP[i,]))] <- 0
    all_Power[i, which(is.na(all_Power[i,]))] <- 0
  }

  # Final results editing
  res_FDR_power <- data.frame(row.names=c("FDR", "Power"),
                              rbind( apply(all_FDP, 1, mean), apply(all_Power, 1, mean) ))
  colnames(res_FDR_power) <- c("BH","wBH","wBHa","IHW","qvalue","swfdr","fdrreg","CAMT")
  res_Subpower <- data.frame(row.names=c("BH_subpower", "wBH_subpower", "wBHa_subpower",
                                         "IHW_subpower","Qvalue_subpower","Swfdr_subpower","FDRreg_subpower", "CAMT_subpower"),
                             rbind(apply(all_SubPower[seq(1,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(2,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(3,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(4,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(5,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(6,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(7,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(8,(N_iteration*8),8),], 2, mean)
                             ))
  res <- list(res_FDR_power, res_Subpower)
  return(res)
}
power_FDR_calc_ind_simu <- function(rejects_vect, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3){
  # Allows to compute the FDP and Powers in independent simulation data
  R <- length(rejects_vect)
  V <- sum(rejects_vect>m1)
  FDP <- c(V/R)
  Power <- c((R-V)/m1)

  rejects_subset_m1_1 <- sum(rejects_vect<=limit_subset_m1_1)
  rejects_subset_m1_2 <- sum(rejects_vect<=limit_subset_m1_2)-rejects_subset_m1_1
  rejects_subset_m1_3 <- sum(rejects_vect<=limit_subset_m1_3)-rejects_subset_m1_1-rejects_subset_m1_2
  Power_subgroup <- c(rejects_subset_m1_1/subset_m1_1, rejects_subset_m1_2/subset_m1_2, rejects_subset_m1_3/subset_m1_3)

  return(list(FDP, Power, Power_subgroup))

}

appli_procedures_corr_simu <- function(N_iteration, pvalues_list, covariates_list, SNPs_Drawn_Corr_list,
                                       corr_draw_all, draw_list, alpha, m1,
                                       limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3,
                                       subset_m1_1,subset_m1_2,subset_m1_3,corr_threshold){
  # Allows to apply the procedures on a correlated simulated dataset (with 500 iterations) and to compute the FDR and the powers (global and subgroup)
  all_FDP <- c()
  all_Power <- c()
  all_SubPower <- c()
  for(i in c(1:N_iteration)){
    pvalues <- pvalues_list[[i]]
    covariates <- covariates_list[[i]]
    SNPs_Drawn_Correlated <- SNPs_Drawn_Corr_list[[i]]
    correlation_all <- corr_draw_all
    SNPs_Dranwn_gp1 <- draw_list[[i]][c(1:limit_subset_m1_1)]
    SNPs_Dranwn_gp2 <- draw_list[[i]][c((limit_subset_m1_1+1):limit_subset_m1_2)]
    SNPs_Dranwn_gp3 <- draw_list[[i]][c((limit_subset_m1_2+1):limit_subset_m1_3)]
    correlation_draw<-subset(correlation_all,select = c(SNPs_Dranwn_gp1,SNPs_Dranwn_gp2,SNPs_Dranwn_gp3))
    
    # wBHa procedure
    res_wBHa <- wBHa(pvalues, covariates, alpha)
    rejects_wBHa <- which(res_wBHa$adjusted_pvalues<alpha)
    res_wBHa <- power_FDR_calc_corr_simu(rejects_wBHa, SNPs_Drawn_Correlated, correlation_draw, correlation_all,
                                         m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, 
                                         subset_m1_1,subset_m1_2,subset_m1_3,
                                         SNPs_Dranwn_gp1,SNPs_Dranwn_gp2,SNPs_Dranwn_gp3,corr_threshold)
    
    # BH procedure
    rejects_BH <- which(p.adjust(pvalues, method="BH")<=alpha)
    res_BH <- power_FDR_calc_corr_simu(rejects_BH, SNPs_Drawn_Correlated, correlation_draw, correlation_all,
                                       m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, 
                                       subset_m1_1,subset_m1_2,subset_m1_3,
                                       SNPs_Dranwn_gp1,SNPs_Dranwn_gp2,SNPs_Dranwn_gp3,corr_threshold)
    
    # IHW procedure
    res_IHW <- ihw(pvalues~covariates, alpha=alpha)
    reject_IHW <- which(adj_pvalues(res_IHW)<alpha)
    res_IHW <- power_FDR_calc_corr_simu(reject_IHW, SNPs_Drawn_Correlated, correlation_draw, correlation_all,
                                        m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, 
                                        subset_m1_1,subset_m1_2,subset_m1_3,
                                        SNPs_Dranwn_gp1,SNPs_Dranwn_gp2,SNPs_Dranwn_gp3,corr_threshold)
    
    # wBH procedure
    rejects_wBH <- which(p.adjust(pvalues/((length(pvalues)/sum(1/covariates))*(1/covariates)), method = "BH")<=alpha)
    res_wBH <- power_FDR_calc_corr_simu(rejects_wBH, SNPs_Drawn_Correlated, correlation_draw, correlation_all,
                                        m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, 
                                        subset_m1_1,subset_m1_2,subset_m1_3,
                                        SNPs_Dranwn_gp1,SNPs_Dranwn_gp2,SNPs_Dranwn_gp3,corr_threshold)
    
    # Qvalue procedure
    res_qvalue <- qvalue(pvalues)
    reject_qvalue <- which(res_qvalue$qvalues<alpha)
    res_qvalue <- power_FDR_calc_corr_simu(reject_qvalue, SNPs_Drawn_Correlated, correlation_draw, correlation_all,
                                           m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, 
                                           subset_m1_1,subset_m1_2,subset_m1_3,
                                           SNPs_Dranwn_gp1,SNPs_Dranwn_gp2,SNPs_Dranwn_gp3,corr_threshold)
    
    # Swfdr procedure
    res_swfdr <- lm_qvalue(pvalues, covariates)
    reject_swfdr <- which(res_swfdr$qvalue<alpha)
    res_swfdr <- power_FDR_calc_corr_simu(reject_swfdr, SNPs_Drawn_Correlated, correlation_draw, correlation_all,
                                          m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, 
                                          subset_m1_1,subset_m1_2,subset_m1_3,
                                          SNPs_Dranwn_gp1,SNPs_Dranwn_gp2,SNPs_Dranwn_gp3,corr_threshold)
    
    # FDRreg procedure
    pvalues2 <- pvalues
    pvalues2[pvalues2==1] <- (1-10^-7)
    zscores <- qnorm(pvalues2)
    res_fdrreg <- FDRreg(zscores, as.matrix(covariates))
    reject_fdrreg <- which(res_fdrreg$FDR<alpha)
    res_fdrreg <- power_FDR_calc_corr_simu(reject_fdrreg, SNPs_Drawn_Correlated, correlation_draw, correlation_all,
                                           m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, 
                                           subset_m1_1,subset_m1_2,subset_m1_3,
                                           SNPs_Dranwn_gp1,SNPs_Dranwn_gp2,SNPs_Dranwn_gp3,corr_threshold)
    
    # CAMT procedure
    res_camt <- camt.fdr(pvals=pvalues, pi0.var=covariates)
    reject_camt <- which(c(res_camt$fdr<alpha))
    res_camt <- power_FDR_calc_corr_simu(reject_camt, SNPs_Drawn_Correlated, correlation_draw, correlation_all,
                                         m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, 
                                         subset_m1_1,subset_m1_2,subset_m1_3,
                                         SNPs_Dranwn_gp1,SNPs_Dranwn_gp2,SNPs_Dranwn_gp3,corr_threshold)
    
    # Results
    all_FDP <- cbind(all_FDP, rbind(res_BH[[1]], res_wBH[[1]], res_wBHa[[1]], res_IHW[[1]], res_qvalue[[1]], res_swfdr[[1]], res_fdrreg[[1]], res_camt[[1]]))
    all_Power <- cbind(all_Power, rbind(res_BH[[2]],res_wBH[[2]],res_wBHa[[2]],res_IHW[[2]],res_qvalue[[2]],res_swfdr[[2]],res_fdrreg[[2]],res_camt[[2]]))
    all_SubPower <- rbind(all_SubPower,rbind(res_BH[[3]],res_wBH[[3]],res_wBHa[[3]],res_IHW[[3]],res_qvalue[[3]],res_swfdr[[3]],res_fdrreg[[3]],res_camt[[3]]))
  }
  
  # N/A or 0 management
  for(i in c(1:8)){
    all_FDP[i, which(is.na(all_FDP[i,]))] <- 0
    all_Power[i, which(is.na(all_Power[i,]))] <- 0
  }
  for (j in c(1:c(N_iteration*8))) {
    all_SubPower[j, which(is.na(all_SubPower[j,]))] <- 0
  }
  
  # Final results editing
  res_FDR_power <- data.frame(row.names=c("FDR", "Power"),
                              rbind( apply(all_FDP, 1, mean), apply(all_Power, 1, mean) ))
  colnames(res_FDR_power) <- c("BH","wBH","wBHa","IHW","qvalue","swfdr","fdrreg","CAMT")
  res_Subpower <- data.frame(row.names=c("BH_subpower", "wBH_subpower", "wBHa_subpower",
                                         "IHW_subpower","Qvalue_subpower","Swfdr_subpower",
                                         "FDRreg_subpower", "CAMT_subpower"),
                             rbind(apply(all_SubPower[seq(1,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(2,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(3,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(4,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(5,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(6,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(7,(N_iteration*8),8),], 2, mean),
                                   apply(all_SubPower[seq(8,(N_iteration*8),8),], 2, mean)
                             ))
  res <- list(res_FDR_power, res_Subpower)
  return(res)
}
power_FDR_calc_corr_simu <- function(rejects_vect, list_SNPs_Drawn_Correlated, Matrix_correlation_draw, Matrix_correlation_all,
                                     m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, 
                                     subset_m1_1, subset_m1_2, subset_m1_3, 
                                     SNPs_Dranwn_gp1, SNPs_Dranwn_gp2, SNPs_Dranwn_gp3, corr_threshold){
  # Allows to compute the FDP and Powers in correlation simulation data
  if(length(rejects_vect)!=0){
    if(length(which(rejects_vect%in%unlist(list_SNPs_Drawn_Correlated)))!=0){
      H0 <- rejects_vect[-which(rejects_vect%in%unlist(list_SNPs_Drawn_Correlated))]
      S <- sum(unlist(lapply(list_SNPs_Drawn_Correlated, Corr_drawn_rej, rej_list=rejects_vect))>0)
      H0_sbgp<-lapply(H0,Corr_of_SNPs,Matrix_corr=Matrix_correlation_all,threshold=corr_threshold)
      H0_sbgp[unlist(lapply(H0_sbgp, rm_subgp_H0,H0_subgpe_list=H0_sbgp))]<-NULL
      V=length(H0_sbgp)
      R=V+S
    }else{
      H0 <- rejects_vect
      S <- sum(unlist(lapply(list_SNPs_Drawn_Correlated, Corr_drawn_rej, rej_list=rejects_vect))>0)
      H0_sbgp<-lapply(H0,Corr_of_SNPs,Matrix_corr=Matrix_correlation_all,threshold=corr_threshold)
      H0_sbgp[unlist(lapply(H0_sbgp, rm_subgp_H0,H0_subgpe_list=H0_sbgp))]<-NULL
      V=length(H0_sbgp)
      R=V+S
    }
    
    FDP=V/R
    Power <- c((S)/m1)
    
    Power_gp1 <- Corr_drawn_rej(rejects_vect,SNPs_Dranwn_gp1)/subset_m1_1
    Power_gp2 <- Corr_drawn_rej(rejects_vect,SNPs_Dranwn_gp2)/subset_m1_2
    Power_gp3 <- Corr_drawn_rej(rejects_vect,SNPs_Dranwn_gp3)/subset_m1_3
    Power_subgroup<-c(Power_gp1, Power_gp2, Power_gp3)
    Power_subgroup[is.na(Power_subgroup)]=0
    
  }else{
    FDP <- NA
    Power <- NA
    Power_subgroup <- c(NA,NA,NA)
    
  }
  
  
  return(list(FDP, Power, Power_subgroup))
}



Corr_drawn_rej<-function(corr_drawn_list, rej_list){
  # Allows to compute how many elements of corr_drawn_list are also in rej_list
  S <- sum((corr_drawn_list%in%rej_list))
  return(S)
}

Corr_of_SNPs<-function(SNPs_draw, Matrix_corr, threshold){
  # Allows to obtain elements in the SNPs_draw columns of Matrix_corr having a value higher than the threshold
  res<-which(Matrix_corr[,SNPs_draw]>=threshold)
  return(res)
}

rm_subgp_H0 <- function(x,H0_subgpe_list){
  
  if(length(x)==1){
    if(x%in%unlist(H0_subgpe_list[lapply(H0_subgpe_list,length)!=1])){
      return(T)
    }else{
      return(F)
    }
  }else{
    return(F)
  }
}


#Library
library(IHW)
library(qvalue)
library(swfdr)
library(FDRreg)
library(CAMT)
library(wBHa)

#Parameters
size_n <- 2000 #nb of individuals
nb_iteration <- 2
r2 <- 0.2 #coefficient of determination
correlation_threshold <- 0.8 #coefficient correlation threshold
alpha <- 0.05 #significance level
path_in <- c("~/")
path_out <- c("~/")


for (Case in c("independent","correlation")) {
  if(Case==c("independent")){ #case without correlations between SNPs
    vrho <- 0 #rho value
    blsize <- 0 #bloc size
    for (scenario in c("reference","inverse","constant")) {
      if(scenario=="reference"){ #corresponding to scenario 1
        binary_beta=c("B_1.8_1.5_1.3_1.rda")
        quantitative_beta=c("B3210.rda")
      }
      if(scenario=="inverse"){ #corresponding to scenario 2
        binary_beta=c("B_1.3_1.5_1.8_1.rda")
        quantitative_beta=c("B1230.rda")
      }
      if(scenario=="constant"){ #corresponding to scenario 3
        binary_beta=c("B_1.5_1.5_1.5_1.rda")
        quantitative_beta=c("B2220.rda")
      }
      for (design in c("quantitative","case_control")) {
        for (size_m in c(8000,14000,20000)) {
          #size_m is the total number of null hypotheses tested
          K=100
          B_size <- size_m/100
          for(m1 in c(5,10,15,20,25,50,100,150)){
            #m1 is the number of causal variants
            subset_m1_1 <- m1%/%3 ; limit_subset_m1_1<-(1+subset_m1_1-1) # Nb of hypotheses in subgroup 1 and cumulated sum
            subset_m1_2 <- m1%/%3 ; limit_subset_m1_2<-(limit_subset_m1_1+subset_m1_2) # Nb of hypotheses in subgroup 2 and cumulated sum
            subset_m1_3 <- (m1%%3)+(m1%/%3) ; limit_subset_m1_3<-(limit_subset_m1_2+subset_m1_3) # Nb of hypotheses in subgroup 3 and cumulated sum
            if (design=="quantitative") {
              linear_file_in <- paste("LineR","n", size_n,"m", size_m,"m1", m1,"rho", vrho,"tb", blsize,"r2", r2, quantitative_beta, sep="_")
              if(file.exists(paste(path_in, linear_file_in, sep=""))){
                load(paste(path_in, linear_file_in, sep=""))

                results <- appli_procedures_ind_simu(nb_iteration, reg_pvalues, reg_covariates, alpha, m1, subset_m1_1, subset_m1_2,
                                                     subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)
                res_FDR_power <- results[[1]]
                res_Subpower <- results[[2]]

                # Final results analysis
                file_out <- paste("Res", linear_file_in, sep="_")
                save(res_FDR_power, res_Subpower, file = paste(path_out, file_out, sep=""))
              }
            }
            if(design=="case_control"){
              logit_file_in <- paste("Logit","n", size_n,"m", size_m,"m1", m1,"rho", vrho,"tb", blsize, binary_beta, sep="_")
              if(file.exists(paste(path_in, logit_file_in, sep=""))){
                load(paste(path_in, logit_file_in, sep=""))

                results <- appli_procedures_ind_simu(nb_iteration, reg_pvalues, reg_covariates, alpha, m1, subset_m1_1, subset_m1_2,
                                                     subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)
                res_FDR_power <- results[[1]]
                res_Subpower <- results[[2]]

                # Final results analysis
                file_out <- paste("Res", logit_file_in, sep="_")
                save(res_FDR_power, res_Subpower, file = paste(path_out, file_out, sep=""))
              }
            }
          }
        }
      }
    }
  }
  if(Case==c("correlation")){ #case with correlations between SNPs
    size_m <- 8000 #total number of null hypotheses tested
    blsize <- 10 #bloc size
    K=100 #nb of boostrap sample
    B_size <- size_m/100 #size of each bootstrap sample
    for (scenario in c("reference","inverse","constant")) {
      if(scenario=="reference"){ #corresponding to scenario 1
        binary_beta=c("B_1.8_1.5_1.3_1.rda")
        quantitative_beta=c("B3210.rda")
      }
      if(scenario=="inverse"){ #corresponding to scenario 2
        binary_beta=c("B_1.3_1.5_1.8_1.rda")
        quantitative_beta=c("B1230.rda")
      }
      if(scenario=="constant"){ #corresponding to scenario 3
        binary_beta=c("B_1.5_1.5_1.5_1.rda")
        quantitative_beta=c("B2220.rda")
      }
      for (design in c("case_control","quantitative")) {
        for(m1 in c(5,10,15,20,25,50,100,150)){
          #m1 is the number of causal variants
          subset_m1_1 <- m1%/%3 ; limit_subset_m1_1 <- (1+subset_m1_1-1) # Nb of hypotheses in subgroup 1 and cumulated sum
          subset_m1_2 <- m1%/%3 ; limit_subset_m1_2 <- (limit_subset_m1_1+subset_m1_2) # Nb of hypotheses in subgroup 2 and cumulated sum
          subset_m1_3 <- (m1%%3)+(m1%/%3) ; limit_subset_m1_3 <- (limit_subset_m1_2+subset_m1_3) # Nb of hypotheses in subgroup 3 and cumulated sum
          for (design in c("case_control","quantitative")) {
            for (vrho in c(0.10,0.20,0.35,0.5,0.75)) {
              #vrho is the correlation value between variants
              if (design=="quantitative") {
                linear_file_in <- paste("LineR","n", size_n,"m", size_m,"m1", m1,"rho", vrho,"tb", blsize, "r2", r2, quantitative_beta, sep="_")
                if(file.exists(paste(path_in, linear_file_in, sep=""))){
                  load(paste(path_in, linear_file_in, sep=""))
                  
                  results <- appli_procedures_corr_simu(nb_iteration, reg_pvalues, reg_covariates, reg_SNPs_Drawn_Correlated, 
                                                        reg_correlation_all, reg_draw, alpha, m1, 
                                                        limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3,
                                                        subset_m1_1,subset_m1_2,subset_m1_3,correlation_threshold)
                  res_FDR_power <- results[[1]]
                  res_Subpower <- results[[2]]
                  
                  # Final results analysis
                  file_out <- paste("Res", linear_file_in, sep="_")
                  save(res_FDR_power, res_Subpower, file = paste(path_out, file_out, sep=""))
                }
              }
              if(design=="case_control"){
                logit_file_in <- paste("Logit","n", size_n,"m", size_m,"m1", m1,"rho", vrho,"tb", blsize, binary_beta, sep="_")
                if(file.exists(paste(path_in, logit_file_in, sep=""))){
                  load(paste(path_in, logit_file_in, sep=""))
                  
                  results <- appli_procedures_corr_simu(nb_iteration, reg_pvalues, reg_covariates, reg_SNPs_Drawn_Correlated, 
                                                        reg_correlation_all, reg_draw, alpha, m1, 
                                                        limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3,
                                                        subset_m1_1,subset_m1_2,subset_m1_3,correlation_threshold)
                  res_FDR_power <- results[[1]]
                  res_Subpower <- results[[2]]
                  
                  # Final results analysis
                  file_out <- paste("Res", logit_file_in, sep="_")
                  save(res_FDR_power, res_Subpower, file=paste(path_out, file_out, sep=""))
                }
              }
            }
          }
        }
      }
    }
  }
}
