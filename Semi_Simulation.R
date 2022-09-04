#Packages
library(bnstruct)
library(pbapply)
library(CAMT)
library(IHW)
library(qvalue)
library(swfdr)
library(FDRreg)
library(wBHa)

#Parameters
alpha <- 0.05
r2 <- 0.8 #coefficient of determination
correlation_threshold <- 0.8 #coefficient correlation threshold
nb_iteration <- 500
path_out <- ("~/")

load(file="~/HIV_data/matgeno/matgeno6.RData")
load(file="~/HIV_data/matpheno2.RData")
vpheno2<- matpheno2$CVinc

#Functions
coeff_calc_pval_beta0_beta<- function(phenotype,SNP_geno){
  #Allows to obtain different coefficients (pvalues, bêta 0 and bêta) from a regression model
  res <- as.numeric(try(summary(glm(phenotype~SNP_geno))$coefficient))
  pval <-res[8]
  beta0<-res[1]
  beta1<-res[2]
  out<-c(pval,beta0,beta1)
  return(out)
}
maf_calc <- function(SNP) {
  # Allows to compute the MAF of a SNP
  freq_A <- (sum(SNP==1,na.rm=T)*2+sum(SNP==2,na.rm=T))/(2*length(SNP))
  MAF <- min(freq_A, 1-freq_A)
  return(MAF)
}
Corr_of_SNPs<-function(SNPs_draw, Matrix_corr, threshold){
  # Allows to obtain elements in the SNPs_draw columns of Matrix_corr having a value higher than the threshold
  res<-which(Matrix_corr[,SNPs_draw]>threshold)
  return(res)
}
quantitative_association_test <- function(phenotype,SNP_geno){
  # Allows to compute the pvalues from a linear model
  pval <- summary(lm(phenotype~SNP_geno))$coefficient[8]
  return(pval)
}
Corr_drawn_rej<-function(corr_drawn_list, rej_list){
  # Allows to compute how many elements of corr_drawn_list are also in rej_list
  S <- sum((corr_drawn_list%in%rej_list))
  return(S)
}
power_FDR_calc_corr_simu <- function(reject_vector, SNPs_Drawn_Correlated, Matrix_correlation_all, SNPs_Dranwn_gp1, SNPs_Dranwn_gp2, SNPs_Dranwn_gp3, seuil_correlation){
  if(length(reject_vector)!=0){
    if(length(which(reject_vector%in%unlist(SNPs_Drawn_Correlated)))!=0){
      H0 <- reject_vector[-which(reject_vector%in%unlist(SNPs_Drawn_Correlated))]
      S<- length(rm_duplicates(lapply(SNPs_Drawn_Correlated,function(x, y){S<-(x%in%y);return(x[S])},y=reject_vector)))
      
      if(length(H0)==1){      
        Matrix_correlation_H0=matrix(Matrix_correlation_all[H0,],nrow = length(H0))
      }else{
        Matrix_correlation_H0=Matrix_correlation_all[H0,]
      }
      H0_sbgp<-rm_duplicates(lapply(H0,Corr_of_SNPs,Matrix_corr=Matrix_correlation_H0,threshold=seuil_correlation))
      
      V=length(H0_sbgp)
      R=V+S
    }else{
      H0 <- reject_vector
      S<- length(rm_duplicates(lapply(SNPs_Drawn_Correlated,function(x, y){S<-(x%in%y);return(x[S])},y=reject_vector)))
      if(length(H0)==1){      
        Matrix_correlation_H0=matrix(Matrix_correlation_all[H0,],nrow = length(H0))
      }else{
        Matrix_correlation_H0=Matrix_correlation_all[H0,]
      }
      H0_sbgp<-rm_duplicates(lapply(H0,Corr_of_SNPs,Matrix_corr=Matrix_correlation_H0,threshold=seuil_correlation))
      V=length(H0_sbgp)
      R=V+S
    }
    m1=length(SNPs_Drawn_Correlated)
    FDP=V/R
    FDP[is.na(FDP)]=0
    Power <- c((S)/m1)
    Power[is.na(Power)]=0
    
    Power_gp1 <- power_subgroup_calc_corr_simu(SNPs_Dranwn_gp1,reject_vector)
    Power_gp2 <- power_subgroup_calc_corr_simu(SNPs_Dranwn_gp2,reject_vector)
    Power_gp3 <- power_subgroup_calc_corr_simu(SNPs_Dranwn_gp3,reject_vector) 
    Power_subgroup<-c(Power_gp1, Power_gp2, Power_gp3)
    Power_subgroup[is.na(Power_subgroup)]=0
    
  }else{
    FDP <- NA
    Power <- NA
    Power_subgroup <- c(NA,NA,NA)
    
  }
  
  
  return(list(FDP, Power, Power_subgroup))
  
}
power_subgroup_calc_corr_simu<-function(draw_list,reject_vector){
  #Allows to compute the power of subgroups
  return(sum(unlist(lapply(draw_list, Corr_drawn_rej, rej_list=reject_vector))>0)/length(draw_list))
}
rm_duplicates <- function(list_input){
  # Allows to search and remove repeated vectors in an input list
  index_rm=NA
  while (length(index_rm)>0) {
    duplicate_element=unlist(list_input)[which(duplicated(unlist(list_input))==T)]
    index_rm=which(sapply(
      lapply(list_input ,match, duplicate_element), 
      function(y,x) x %in% y,NA)==F)
    if (length(index_rm)>0) {
      list_temp=list_input
      list_temp[[index_rm[1]]]<-NULL
      list_input=list_temp 
    }
  }
  return(list_input)
}

# Transposition :
matgeno <- t(matgeno)

# Remove SNPs with a number > 10 NA :
propNAcol <- colSums(is.na(matgeno))
matgeno <- subset(matgeno,select = -c(which(propNAcol>10)))

# Imputation
matrix_imp<-knn.impute(matgeno,k=1)

## Calculation of MAF
maf <- apply(matrix_imp, 2, maf_calc)

# Parameters
size_m <- dim(matrix_imp)[2]
size_n <- dim(matrix_imp)[1]

# Cutting of MAF
SNPs_maf_gp1 <- which(maf>0.05&maf<0.15)
SNPs_maf_gp2 <- which(maf>0.15&maf<0.30)
SNPs_maf_gp3 <- which(maf>0.30&maf<0.50)


coeff_imp <- pbapply(matrix_imp, 2, coeff_calc_pval_beta0_beta, phenotype=vpheno2)
rej_wBH_imp <- which(p.adjust(coeff_imp[1,]/((length(coeff_imp[1,])/sum(1/maf))*(1/maf)), method="BH")<=alpha)
snp_associated_imp <- coeff_imp[1,rej_wBH_imp]
coeff_snp_associated_imp <- subset(as.data.frame(coeff_imp),select = names(snp_associated_imp))
quant_imp <- quantile(abs(as.numeric(coeff_snp_associated_imp[3,])))

## Matrix of correlation
mcor<-as.data.frame(abs(cor(matrix_imp,matrix_imp)))

#save.image("~/Need_for_SemiSimu2.RData")
#load("~/Need_for_SemiSimu2.RData")
for (m1 in c(15,20,25,50,100,150)) {
  #m1 is the number of causal variants
  m0 <- size_m-m1
  subset_m1_1 <- m1%/%3 ; limit_subset_m1_1<-(1+subset_m1_1-1) # Nb of hypotheses in subgroup 1 and cumulated sum
  subset_m1_2 <- m1%/%3 ; limit_subset_m1_2<-(limit_subset_m1_1+subset_m1_2) # Nb of hypotheses in subgroup 2 and cumulated sum
  subset_m1_3 <- (m1%%3)+(m1%/%3) ; limit_subset_m1_3<-(limit_subset_m1_2+subset_m1_3) # Nb of hypotheses in subgroup 3 and cumulated sum
  
  # Draw for causal variants
  draw1<-sample(SNPs_maf_gp1,m1%/%3)
  draw2<-sample(SNPs_maf_gp2,m1%/%3)
  draw3<-sample(SNPs_maf_gp3,m1-2*(m1%/%3))
  draw<-c(draw1,draw2,draw3)
  
  # Dermination of bêta0 and bêta
  beta <- c(rep(0,size_m))
  
  beta[draw1]=quant_imp[4]
  beta[draw2]=quant_imp[3]
  beta[draw3]=quant_imp[2]
  
  beta <- as.matrix(unlist(beta))
  
  XB <- matrix_imp %*% beta # Epsilon Matrix (quantitative)
  sigma2 <- ( (r2-1) * sum((XB-mean(XB))^2) ) / ( r2*(2-size_n) )
  
  # How to define the True Positives (TP) ?
  
  ## Define the list of drawn correlated 
  List_SNPs_Drawn_Correlated <- rm_duplicates(lapply(draw,
                                                     Corr_of_SNPs, 
                                                     Matrix_corr=mcor, 
                                                     threshold=correlation_threshold))
  draw1_bis<-List_SNPs_Drawn_Correlated[unlist(lapply(
    lapply(List_SNPs_Drawn_Correlated,match,draw1),function(x){
      x[is.na(x)]<-0;return(sum(x))}))>0]
  draw1_bis=rm_duplicates(draw1_bis)
  draw2_bis<-List_SNPs_Drawn_Correlated[unlist(lapply(
    lapply(List_SNPs_Drawn_Correlated,match,draw2),function(x){
      x[is.na(x)]<-0;return(sum(x))}))>0]
  draw2_bis=rm_duplicates(draw2_bis)
  draw3_bis<-List_SNPs_Drawn_Correlated[unlist(lapply(
    lapply(List_SNPs_Drawn_Correlated,match,draw3),function(x){
      x[is.na(x)]<-0;return(sum(x))}))>0]
  draw3_bis=rm_duplicates(draw3_bis)
  
  
  all_FDP<-c()
  all_Power<-c()
  all_subpower<-c()
  for (simu in c(1:nb_iteration)) {
    epsilon <- as.matrix(rnorm(size_n, mean=0, sd=sigma2)) # Epsilon matrix
    Y <- XB + epsilon # Computation of Y (quantitative)
    pvalues <- apply(matrix_imp, 2, quantitative_association_test, phenotype=Y) # Computation of quantitative Pj = Pvalues
    
    #BH procedure
    rejects_BH <- which(p.adjust(pvalues, method="BH")<=alpha)
    res_BH <- power_FDR_calc_corr_simu(rejects_BH,List_SNPs_Drawn_Correlated,mcor,
                                       draw1_bis,draw2_bis,draw3_bis,correlation_threshold)
    #IHW procedure
    res_IHW <- ihw(pvalues~maf,alpha=alpha)
    rejects_IHW <- which(adj_pvalues(res_IHW)<0.05)
    res_IHW <- power_FDR_calc_corr_simu(rejects_IHW,List_SNPs_Drawn_Correlated,mcor,
                                        draw1_bis,draw2_bis,draw3_bis,correlation_threshold)
    #wBH procedure
    rejects_wBH <- which(p.adjust(pvalues/((length(pvalues)/sum(1/maf))*(1/maf)), method = "BH")<=alpha)
    res_wBH <- power_FDR_calc_corr_simu(rejects_wBH,List_SNPs_Drawn_Correlated,mcor,
                                        draw1_bis,draw2_bis,draw3_bis,correlation_threshold)
    # wBHa procedure
    res_wBHa <- wBHa(pvalues, maf, alpha)
    rejects_wBHa <- which(res_wBHa$adjusted_pvalues<alpha)
    res_wBHa <- power_FDR_calc_corr_simu(rejects_wBHa,List_SNPs_Drawn_Correlated,mcor,
                                         draw1_bis,draw2_bis,draw3_bis,correlation_threshold)
    # Qvalue procedure
    res_qvalue <- qvalue(pvalues)
    rejects_qvalue <- which(res_qvalue$qvalues<alpha)
    res_qvalue <- power_FDR_calc_corr_simu(rejects_qvalue, List_SNPs_Drawn_Correlated,mcor,
                                           draw1_bis,draw2_bis,draw3_bis,correlation_threshold)
    # Swfdr procedure
    res_lm_qvalue <- lm_qvalue(pvalues, maf)
    reject_swfdr <- which(res_lm_qvalue$qvalue<alpha)
    res_swfdr <- power_FDR_calc_corr_simu(reject_swfdr,List_SNPs_Drawn_Correlated,mcor,
                                          draw1_bis,draw2_bis,draw3_bis,correlation_threshold)
    # FDRreg procedure
    pvalues2 <- pvalues
    pvalues2[pvalues2==1]<-(1-10^-7)
    zscores <- qnorm(pvalues2)
    res_fdrreg <- FDRreg(zscores, as.matrix(maf))
    rejects_fdrreg <- which(res_fdrreg$FDR<alpha)
    res_fdrreg <- power_FDR_calc_corr_simu(rejects_fdrreg, List_SNPs_Drawn_Correlated,mcor,
                                           draw1_bis,draw2_bis,draw3_bis,correlation_threshold)
    # CAMT procedure
    res_camt <- camt.fdr(pvals=pvalues,pi0.var=maf)
    reject_camt <- which(c(res_camt$fdr<alpha))
    res_camt <- power_FDR_calc_corr_simu(reject_camt, List_SNPs_Drawn_Correlated,mcor,
                                         draw1_bis,draw2_bis,draw3_bis,correlation_threshold)
    #Results
    # Results
    all_FDP <- cbind(all_FDP, rbind(res_BH[[1]], res_wBH[[1]], res_wBHa[[1]], res_IHW[[1]], res_qvalue[[1]], res_swfdr[[1]], res_fdrreg[[1]], res_camt[[1]]))
    all_Power <- cbind(all_Power, rbind(res_BH[[2]],res_wBH[[2]],res_wBHa[[2]],res_IHW[[2]],res_qvalue[[2]],res_swfdr[[2]],res_fdrreg[[2]],res_camt[[2]]))
    all_subpower <- rbind(all_subpower,rbind(res_BH[[3]],res_wBH[[3]],res_wBHa[[3]],res_IHW[[3]],res_qvalue[[3]],res_swfdr[[3]],res_fdrreg[[3]],res_camt[[3]]))
  }
  # N/A or 0 management
  for(i in c(1:8)){
    all_FDP[i, which(is.na(all_FDP[i,]))] <- 0
    all_Power[i, which(is.na(all_Power[i,]))] <- 0
  }
  for (j in c(1:c(simu*8))) {
    all_subpower[j, which(is.na(all_subpower[j,]))] <- 0
  }
  
  # Final results editing
  res_FDR_power <- data.frame(row.names=c("FDR", "Power"),
                              rbind( apply(all_FDP, 1, mean), apply(all_Power, 1, mean) ))
  colnames(res_FDR_power) <- c("BH","wBH","wBHa","IHW","qvalue","swfdr","fdrreg","CAMT")
  res_Subpower <- data.frame(row.names=c("BH_subpower", "wBH_subpower", "wBHa_subpower",
                                         "IHW_subpower","Qvalue_subpower","Swfdr_subpower",
                                         "FDRreg_subpower", "CAMT_subpower"),
                             rbind(apply(all_subpower[seq(1,(nb_iteration*8),8),], 2, mean),
                                   apply(all_subpower[seq(2,(nb_iteration*8),8),], 2, mean),
                                   apply(all_subpower[seq(3,(nb_iteration*8),8),], 2, mean),
                                   apply(all_subpower[seq(4,(nb_iteration*8),8),], 2, mean),
                                   apply(all_subpower[seq(5,(nb_iteration*8),8),], 2, mean),
                                   apply(all_subpower[seq(6,(nb_iteration*8),8),], 2, mean),
                                   apply(all_subpower[seq(7,(nb_iteration*8),8),], 2, mean),
                                   apply(all_subpower[seq(8,(nb_iteration*8),8),], 2, mean)
                             ))
  save(res_FDR_power, res_Subpower,
       file = paste(path_out,"SemiSimu_HIV_m1_",m1,"_threshold_",correlation_threshold,"_R2_",r2,".RData",sep=""))
}