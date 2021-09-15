#Functions
appli_procedures_ind_simu <- function(N_iteration, pvalues_list, covariates_list, alpha, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3){
  # Allows to apply the procedures on a independent simulation dataset (with 500 iterations) and to compute the FDR and the powers (global and subgroup)

  all_FDP <- c()
  all_Power <- c()
  all_SubPower <- c()

  for(i in c(1:N_iteration)){
    pvalues <- pvalues_list[[i]]
    covariates <- covariates_list[[i]]

    # wBHa procedure
    res_wBHa <- wBHa(pvalues, covariates, alpha)
    rejects_wBHa <- which(res_wBHa$adjusted_pvalues<alpha)
    res_wBHa <- power_FDR_calc_ind_simu(rejects_wBHa, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

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
    pvalues[pvalues==1] <- (1-10^-7)
    zscores <- qnorm(pvalues)
    res_fdrreg <- FDRreg(zscores, as.matrix(covariates))
    reject_fdrreg <- which(res_fdrreg$FDR<alpha)
    res_fdrreg <- power_FDR_calc_ind_simu(reject_fdrreg, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

    # CAMT procedure
    res_camt <- camt.fdr(pvals=pvalues, pi0.var=covariates)
    reject_camt <- which(c(res_camt$fdr<alpha))
    res_camt <- power_FDR_calc_ind_simu(reject_camt, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)

    # Results
    all_FDP <- cbind(all_FDP, rbind(res_BH[[1]], res_wBH[[1]], res_wBHa[[1]], res_IHW[[1]], res_qvalue[[1]], res_swfdr[[1]], res_fdrreg[[1]], res_camt[[1]]))
    all_Power <- cbind(all_Power, rbind(res_BH[[2]], res_wBH[[2]], res_wBHa[[2]], res_IHW[[2]],  res_qvalue[[2]], res_swfdr[[2]], res_fdrreg[[2]], res_camt[[2]]))
    all_SubPower <- rbind(all_SubPower, rbind(res_BH[[3]], res_wBH[[3]], res_wBHa[[3]], res_IHW[[3]],  res_qvalue[[3]], res_swfdr[[3]], res_fdrreg[[3]], res_camt[[3]]))
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

appli_procedures_corr_simu <- function(N_iteration, pvalues_list, covariates_list, sNPs_Drawn_Corr_list, corr_draw_list, alpha, m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m){
  # Allows to apply the procedures on a correlated simulated dataset (with 500 iterations) and to compute the FDR and the powers (global and subgroup)
  all_FDP <- c()
  all_Power <- c()
  all_SubPower <- c()
  for(i in c(1:N_iteration)){
    pvalues <- pvalues_list[[i]]
    covariates <- covariates_list[[i]]
    SNPs_Drawn_Correlated <- SNPs_Drawn_Corr_list[[i]]
    correlation_draw <- corr_draw_list[[i]]

    # wBHa procedure
    res_wBHa <- wBHa(pvalues, covariates, alpha)
    rejects_wBHa <- which(res_wBHa$adjusted_pvalues<alpha)
    res_wBHa <- power_FDR_calc_corr_simu(rejects_wBHa, sNPs_Drawn_Correlated, correlation_draw,
                                         m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m)

    # BH procedure
    rejects_BH <- which(p.adjust(pvalues, method="BH")<=alpha)
    res_BH <- power_FDR_calc_corr_simu(rejects_BH, sNPs_Drawn_Correlated, correlation_draw,
                                       m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m)

    # IHW procedure
    res_IHW <- ihw(pvalues~covariates, alpha=alpha)
    reject_IHW <- which(adj_pvalues(res_IHW)<alpha)
    res_IHW <- power_FDR_calc_corr_simu(reject_IHW, sNPs_Drawn_Correlated, correlation_draw,
                                        m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m)

    # wBH procedure
    rejects_wBH <- which(p.adjust(pvalues/((length(pvalues)/sum(1/covariates))*(1/covariates)), method = "BH")<=alpha)
    res_wBH <- power_FDR_calc_corr_simu(rejects_wBH, sNPs_Drawn_Correlated, correlation_draw,
                                        m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m)

    # Qvalue procedure
    res_qvalue <- qvalue(pvalues)
    reject_qvalue <- which(res_qvalue$qvalues<alpha)
    res_qvalue <- power_FDR_calc_corr_simu(reject_qvalue, SNPs_Drawn_Correlated, correlation_draw,
                                           m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m)

    # Swfdr procedure
    res_swfdr <- lm_qvalue(pvalues, covariates)
    reject_swfdr <- which(res_swfdr$qvalue<alpha)
    res_swfdr <- power_FDR_calc_corr_simu(reject_swfdr, sNPs_Drawn_Correlated, correlation_draw,
                                          m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m)

    # FDRreg procedure
    pvalues[pvalues==1] <- (1-10^-7)
    zscores <- qnorm(pvalues)
    res_fdrreg <- FDRreg(zscores, as.matrix(covariates))
    reject_fdrreg <- which(res_fdrreg$FDR<alpha)
    res_fdrreg <- power_FDR_calc_corr_simu(reject_fdrreg, SNPs_Drawn_Correlated, correlation_draw,
                                           m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m)

    # CAMT procedure
    res_camt <- camt.fdr(pvals=pvalues, pi0.var=covariates)
    reject_camt <- which(c(res_camt$fdr<alpha))
    res_camt <- power_FDR_calc_corr_simu(reject_camt, SNPs_Drawn_Correlated, correlation_draw,
                                         m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m)

    # Results
    all_FDP <- cbind(all_FDP, rbind(res_BH[[1]], res_wBH[[1]], res_wBHa[[1]], res_IHW[[1]], res_qvalue[[1]], res_swfdr[[1]], res_fdrreg[[1]], res_camt[[1]]))
    all_Power <- cbind(all_Power, rbind(res_BH[[2]], res_wBH[[2]], res_wBHa[[2]], res_IHW[[2]],  res_qvalue[[2]], res_swfdr[[2]], res_fdrreg[[2]], res_camt[[2]]))
    all_SubPower <- rbind(all_SubPower, rbind(res_BH[[3]], res_wBH[[3]], res_wBHa[[3]], res_IHW[[3]],  res_qvalue[[3]], res_swfdr[[3]], res_fdrreg[[3]], res_camt[[3]]))
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
power_FDR_calc_corr_simu <- function(rejects_vect, list_SNPs_Drawn_Correlated, Matrix_correlation_draw, m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, m){
  # Allows to compute the FDP and Powers in correlation simulation data
  R <- length(rejects_vect)
  V <- length(rejects_vect[-which(rejects_vect%in%unlist(list_SNPs_Drawn_Correlated))])
  FDP <- c(V/R)
  Power <- c((R-V)/m1)

  # Identification of the True Positives (TP)
  snps_gp1 <- c()
  for (j in c(1:limit_subset_m1_1)) {
    snps_gp1 <- c(snps_gp1, which(Matrix_correlation_draw[,j]>0.8))
  }
  snps_gp1 <- unique(snps_gp1)
  snps_gp2 <- c()
  for (k in c((limit_subset_m1_1+1):limit_subset_m1_2)) {
    snps_gp2 <- c(snps_gp2, which(Matrix_correlation_draw[,k]>0.8))
  }
  snps_gp2 <- unique(snps_gp2)
  snps_gp3 <- c()
  for (l in c((limit_subset_m1_2+1):limit_subset_m1_3)) {
    snps_gp3 <- c(snps_gp3, which(Matrix_correlation_draw[, l]>0.8))
  }
  snps_gp3 <- unique(snps_gp3)
  reject_matrix <-matrix(FALSE, nrow=m, ncol=1)
  reject_matrix[rejects_vect,1]=TRUE

  Power_gp1 <- prop.table(table(as.vector(reject_matrix[snps_gp1,])))["TRUE"]
  Power_gp2 <- prop.table(table(as.vector(reject_matrix[snps_gp2,])))["TRUE"]
  Power_gp3 <- prop.table(table(as.vector(reject_matrix[snps_gp3,])))["TRUE"]
  Power_subgroup <- c(Power_gp1, Power_gp2, Power_gp3)

  return(list(FDP, Power, Power_subgroup))
}


#Library
library(IHW)
library(qvalue)
library(swfdr)
library(FDRreg)
library(CAMT)
library(wBHa)

#Parameters
size_n <- 2000
nb_iteration <- 1
r2 <- 0.2
alpha <- 0.05
path_in <- c("/home/lobry/Documents/Nouveau2021/Test/")
path_out <- c("/home/lobry/Documents/Nouveau2021/Test/")


for (Case in c("independent","correlation")) {
  if(Case==c("independent")){
    vrho <- 0
    blsize <- 0
    for (scenario in c("reference","inverse","constant")) {
      for (design in c("quantitative","case_control")) {
        for (size_m in c(8000,14000,20000)) {
          for(m1 in c(25,50,100,150)){
            subset_m1_1 <- m1%/%3 ; limit_subset_m1_1<-(1+subset_m1_1-1) # Nb of hypotheses in subgroup 1 and cumulated sum
            subset_m1_2 <- m1%/%3 ; limit_subset_m1_2<-(limit_subset_m1_1+subset_m1_2) # Nb of hypotheses in subgroup 2 and cumulated sum
            subset_m1_3 <- (m1%%3)+(m1%/%3) ; limit_subset_m1_3<-(limit_subset_m1_2+subset_m1_3) # Nb of hypotheses in subgroup 3 and cumulated sum
            if (design=="quantitative") {
              if(scenario=="reference"){quantitative_beta <- c("B3210.RData")}
              if(scenario=="inverse"){quantitative_beta <- c("B1230.RData")}
              if(scenario=="constant"){quantitative_beta <- c("B2220.RData")}
              linear_file_in <- paste("LineR","n", size_n,"m", size_m,"m1", m1,"rho", vrho,"tb", blsize,"r2", r2, quantitative_beta, sep="_")
              if(file.exists(paste(path_in, linear_file_in, sep=""))){
                load(paste(path_in, linear_file_in, sep=""))

                results <- appli_procedures_ind_simu(nb_iteration, reg_pvalues, reg_covariates, alpha, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)
                res_FDR_power <- results[[1]]
                res_Subpower <- results[[2]]

                # Final results analysis
                file_out <- paste("Res", linear_file_in, sep="_")
                save(res_FDR_power, res_Subpower, file = paste(path_out, file_out, sep=""))
              }
            }
            if(design=="case_control"){
              if(scenario=="reference"){binary_beta <- c("B_1.8_1.5_1.3_1.RData")}
              if(scenario=="inverse"){binary_beta <- c("B_1.3_1.5_1.8_1.RData")}
              if(scenario=="constant"){binary_beta <- c("B_1.5_1.5_1.5_1.RData")}
              logit_file_in <- paste("Logit","n", size_n,"m", size_m,"m1", m1,"rho", vrho,"tb", blsize, binary_beta, sep="_")
              if(file.exists(paste(path_in, logit_file_in, sep=""))){
                load(paste(path_in, logit_file_in, sep=""))

                results <- appli_procedures_ind_simu(nb_iteration, reg_pvalues, reg_covariates, alpha, m1, subset_m1_1, subset_m1_2, subset_m1_3, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3)
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
  if(Case==c("correlation")){
    size_m <- 8000
    m1 <- 50
    subset_m1_1 <- m1%/%3 ; limit_subset_m1_1 <- (1+subset_m1_1-1) # Nb of hypotheses in subgroup 1 and cumulated sum
    subset_m1_2 <- m1%/%3 ; limit_subset_m1_2 <- (limit_subset_m1_1+subset_m1_2) # Nb of hypotheses in subgroup 2 and cumulated sum
    subset_m1_3 <- (m1%%3)+(m1%/%3) ; limit_subset_m1_3 <- (limit_subset_m1_2+subset_m1_3) # Nb of hypotheses in subgroup 3 and cumulated sum
    blsize <- 10
    for (design in c("case_control","quantitative")) {
      for (vrho in c(0.10,0.20,0.35,0.5,0.75)) {
        if (design=="quantitative") {
          linear_file_in <- paste("LineR","n", size_n,"m", size_m,"m1", m1,"rho", vrho,"tb", blsize, "r2", r2, "B3210.RData", sep="_")
          if(file.exists(paste(path_in, linear_file_in, sep=""))){
            load(paste(path_in, linear_file_in, sep=""))

            results <- appli_procedures_corr_simu(nb_iteration, reg_pvalues, reg_covariates, reg_SNPs_Drawn_Correlated, reg_correlation_draw, alpha, m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m)
            res_FDR_power <- results[[1]]
            res_Subpower <- results[[2]]

            # Final results analysis
            file_out <- paste("Res", linear_file_in, sep="_")
            save(res_FDR_power, res_Subpower, file = paste(path_out, file_out, sep=""))
          }
        }
        if(design=="case_control"){
          logit_file_in <- paste("Logit","n", size_n,"m", size_m,"m1", m1,"rho", vrho,"tb", blsize, "B_1.8_1.5_1.3_1.RData", sep="_")
          if(file.exists(paste(path_in, logit_file_in, sep=""))){
            load(paste(path_in, logit_file_in, sep=""))

            results <- appli_procedures_corr_simu(nb_iteration, reg_pvalues, reg_covariates, reg_SNPs_Drawn_Correlated, reg_correlation_draw, alpha, m1, limit_subset_m1_1, limit_subset_m1_2, limit_subset_m1_3, size_m)
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
