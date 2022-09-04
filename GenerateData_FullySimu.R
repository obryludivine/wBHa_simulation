#Functions

geno_calc <- function(NbInd,vecMAF){
  # Allows to obtain the genotypes of NbInd individuals governed under a probability vecMaf.
  z1 <- rbinom(NbInd, 1, vecMAF)
  z2 <- rbinom(NbInd, 1, vecMAF)
  X <- z1 + z2
  return(X)
}

maf_calc <- function(SNP) {
  # Allows to compute the MAF of a SNP
  freq_A <- (sum(SNP==0)*2+sum(SNP==1))/(2*length(SNP))
  MAF <- min(freq_A, 1-freq_A)
  return(MAF)
}

quantitative_association_test <- function(phenotype,SNP_geno){
  # Allows to compute the pvalues from a linear model
  pval <- summary(lm(phenotype~SNP_geno))$coefficient[8]
  return(pval)
}

binary_association_test <- function(phenotype,SNP_geno){
  # Allows to compute the pvalues from a logistic model
  pval <- summary(glm(phenotype~SNP_geno,family=binomial(link=logit)))$coefficient[8]
  return(pval)
}

generate_ind_data <- function(nb_iteration=500, Design="quantitative",ScenarioBeta="reference", size_n=2000, size_m=8000, m1=50, R2=r2){
  # Allows to generate pvalues and covariates in the case of independent data (without correlation between SNPs)
  
  reg_covariates<-list()
  reg_pvalues<-list()
  
  for(iteration in c(1:nb_iteration)){
    #####################################
    #          Data generation          #
    #####################################
    # Design: Choice between quantitative or case_control design to simulate pvalues
    # size_n: Nb of individuals
    # size_m: Nb of SNPs
    # m1: Nb of causal SNPs
    m0<-size_m-m1 #  m0: Nb of non causal SNPs
    subset_m1_1<-m1%/%3 # Nb of hypotheses in subgroup 1
    subset_m1_2<-m1%/%3 # Nb of hypotheses in subgroup 2
    subset_m1_3<-(m1%%3)+(m1%/%3) # Nb of hypotheses in subgroup 3
    ## Regarding covariables
    theoretical_maf_m0 <- runif(m0, 0.05, 0.5)
    theoretical_maf_subset_m1_1 <- runif(subset_m1_1,0.05,0.15)
    theoretical_maf_subset_m1_2 <- runif(subset_m1_2,0.15,0.25)
    theoretical_maf_subset_m1_3 <- runif(subset_m1_3,0.3,0.4)
    
    ## Regarding genotypes
    geno_matrix_m0 <- apply(as.matrix(theoretical_maf_m0), 1, geno_calc,NbInd=size_n)
    geno_matrix_subset_m1_1 <- apply(as.matrix(theoretical_maf_subset_m1_1), 1, geno_calc,NbInd=size_n)
    geno_matrix_subset_m1_2 <- apply(as.matrix(theoretical_maf_subset_m1_2), 1, geno_calc,NbInd=size_n)
    geno_matrix_subset_m1_3 <- apply(as.matrix(theoretical_maf_subset_m1_3), 1, geno_calc,NbInd=size_n)
    matrix_X <- cbind(geno_matrix_subset_m1_1,geno_matrix_subset_m1_2,geno_matrix_subset_m1_3,geno_matrix_m0)
    
    covariates <- apply(matrix_X, 2, maf_calc) # Computation of the estimated MAF = covariates
    
    ## Regarding p-values
    if(Design=="quantitative"){ # If quantitative design
      R2<-r2 # Percentage explaining the model coefficient of determination
      ### Beta building
      if(ScenarioBeta=="reference"){
        beta <- as.matrix(c(rep(3,subset_m1_1),rep(2,subset_m1_2),rep(1,subset_m1_3), rep(0,m0)))
      }
      if(ScenarioBeta=="inverse"){
        beta <- as.matrix(c(rep(1,subset_m1_1),rep(2,subset_m1_2),rep(3,subset_m1_3), rep(0,m0)))
      }
      if(ScenarioBeta=="constant"){
        beta <- as.matrix(c(rep(2,subset_m1_1),rep(2,subset_m1_2),rep(2,subset_m1_3), rep(0,m0)))
      }
      
      ### Regarding p-values
      XB <- matrix_X %*% beta # Epsilon Matrix (quantitative)
      sigma2 <- ( (r2-1) * sum((XB-mean(XB))^2) ) / ( r2*(2-size_n) )
      epsilon <- as.matrix(rnorm(size_n, mean=0, sd=sqrt(sigma2) ))
      Y <- XB + epsilon # Computation of Y (quantitative)
      pvalues <- apply(matrix_X, 2, quantitative_association_test,phenotype=Y) # Computation of quantitative Pj = Pvalues
    } else if(Design=="case_control"){ # If case_control design
      ### Beta building
      if(ScenarioBeta=="reference"){
        beta <- as.matrix(c(rep(log(1.8),subset_m1_1),rep(log(1.5),subset_m1_2),rep(log(1.3),subset_m1_3), rep(log(1),m0)))
      }
      if(ScenarioBeta=="inverse"){
        beta <- as.matrix(c(rep(log(1.3),subset_m1_1),rep(log(1.5),subset_m1_2),rep(log(1.8),subset_m1_3), rep(log(1),m0)))
      }
      if(ScenarioBeta=="constant"){
        beta <- as.matrix(c(rep(log(1.5),subset_m1_1),rep(log(1.5),subset_m1_2),rep(log(1.5),subset_m1_3), rep(log(1),m0)))
      }
      
      ### Choice of bêta0
      list_ratio <- c()
      possible_values_beta0<-seq(0,-50,-0.1)
      for (beta0 in possible_values_beta0) {
        eXB <- exp(beta0+matrix_X%*%beta)
        Y <- rbinom(size_n, 1, (eXB/ (1+eXB)) )
        Y <- as.factor(Y)
        list_ratio[(length(list_ratio)+1)]<-c(prop.table(table(Y))[2])
      }
      beta0 <- possible_values_beta0[order(abs(list_ratio-0.5))[1]]
      eXB <- exp(beta0+matrix_X%*%beta)
      Y <- rbinom(size_n, 1, (eXB/ (1+eXB)) ) # Computation of Y:
      Y <- as.factor(Y)
      pvalues <- apply(matrix_X, 2, binary_association_test,phenotype=Y) # Computation of case_control Pj = Pvalues
    }
    
    reg_covariates[[iteration]]<-covariates
    reg_pvalues[[iteration]]<-pvalues
  }
  
  return(list(reg_pvalues, reg_covariates))
}

discretise_geno <- function(SNP){
  # Allows to discretize the values according to the Hardy-Weinberg equation
  MAF<-SNP[1]
  genotype<-SNP[-1]
  qp2<-quantile(genotype, probs = MAF^2)
  qq2<-quantile(genotype, probs = (MAF^2) + (2*MAF*(1-MAF)) )
  genotype[between(genotype,qp2,qq2)]=1
  genotype[genotype>=qq2&genotype!=1]=0
  genotype[genotype<=qp2&genotype!=1&genotype!=0]=2
  return(genotype)
}

generate_corr_data <- function(nb_iteration=500, Design="quantitative", ScenarioBeta="reference", size_n=2000, size_m=8000, m1=50, blsize=10, vrho=0.20, R2=0.2){
  # Allows to generates p-values and covariates in the case of correlated data (with correlation between SNPs)
  
  library(Matrix)
  library(dplyr)
  library(mvtnorm)
  
  #####################################
  #          Data generation          #
  #####################################
  # Design: Choice between quantitative or case_control design to simulate pvalues
  # size_n: Nb of individuals
  # size_m: Nb of SNPs
  # m1: Nb of causal SNPs
  # blsize : bloc size
  
  blsize <- 10 
  nblock <- size_m/blsize
  mylist <- list()
  mysubmat <- matrix(vrho,blsize,blsize)
  diag(mysubmat) <- 1
  for(i in 1:nblock) mylist[[i]] <- mysubmat
  vsigma <- as.matrix(bdiag(mylist))
  
  m0<-size_m-m1 #  m0: Nb of non causal SNPs
  subset_m1_1<-m1%/%3 # Nb of hypotheses in subgroup 1
  subset_m1_2<-m1%/%3 # Nb of hypotheses in subgroup 2
  subset_m1_3<-(m1%%3)+(m1%/%3) # Nb of hypotheses in subgroup 3
  
  ## Regarding covariables
  theoretical_maf_m0 <- runif(m0, 0.05, 0.5)
  theoretical_maf_subset_m1_1 <- runif(subset_m1_1,0.05,0.15)
  theoretical_maf_subset_m1_2 <- runif(subset_m1_2,0.15,0.25)
  theoretical_maf_subset_m1_3 <- runif(subset_m1_3,0.3,0.4)
  total_theoretical_maf <- c(theoretical_maf_subset_m1_1,
                             theoretical_maf_subset_m1_2,
                             theoretical_maf_subset_m1_3,
                             theoretical_maf_m0)
  
  mat1<-rmvnorm(size_n,mean=rep(0,size_m/4),sigma = vsigma[1:2000,1:2000])
  mat2<-rmvnorm(size_n,mean=rep(0,size_m/4),sigma = vsigma[2001:4000,2001:4000])
  mat3<-rmvnorm(size_n,mean=rep(0,size_m/4),sigma = vsigma[4001:6000,4001:6000])
  mat4<-rmvnorm(size_n,mean=rep(0,size_m/4),sigma = vsigma[6001:8000,6001:8000])
  
  matbis1 <- rbind(total_theoretical_maf[1:2000],mat1)
  matbis2 <- rbind(total_theoretical_maf[2001:4000],mat2)
  matbis3 <- rbind(total_theoretical_maf[4001:6000],mat3)
  matbis4 <- rbind(total_theoretical_maf[6001:8000],mat4)
  
  matrix_X <- apply(cbind(matbis1, matbis2, matbis3, matbis4),
                    2, discretise_geno)
  
  covariates <- apply(matrix_X, 2, maf_calc) # Computation of the estimated MAF = covariates
  
  # Cutting the MAF into subgroups
  SNPs_maf_subgp_m1_1 <- which(covariates>0.05&covariates<0.15)
  SNPs_maf_subgp_m1_2 <- which(covariates>0.15&covariates<0.25)
  SNPs_maf_subgp_m1_3 <- which(covariates>0.30&covariates<0.50)
  
  # How to define the True Positives (TP) ?
  ## Define the correlation matrix
  mcor<-as.data.frame(abs(cor(matrix_X,matrix_X)))
  mcor_transit<-matrix(0,ncol=size_m,nrow=size_m)
  mcor_transit[which(upper.tri(mcor_transit,T),arr.ind=TRUE)[,2:1]]<-c(mcor[which(upper.tri(mcor,T),arr.ind=TRUE)[,2:1]])
  correlation_matrix<-mcor_transit
  
  reg_covariates<-list()
  reg_pvalues<-list()
  List_SNPs_Drawn_Correlated<-list()
  reg_draw<-list()
  for(iteration in c(1:nb_iteration)){
    
    ## Beta initialization
    beta <- c(rep(0,size_m))
    
    ## Draw among them
    draw1<-sample(SNPs_maf_subgp_m1_1,m1%/%3)
    draw2<-sample(SNPs_maf_subgp_m1_2,m1%/%3)
    draw3<-sample(SNPs_maf_subgp_m1_3,m1-2*(m1%/%3))
    draw<-c(draw1,draw2,draw3)
    
    ## Regarding p-values
    if(Design=="quantitative"){ # If quantitative design
      r2 <- R2 # Percentage explaining the model (coefficient of determination)
      ### Beta building
      if(ScenarioBeta=="reference"){
        beta[draw1] <- 3; beta[draw2] <- 2; beta[draw3] <- 1 # Beta for each group
      }
      if(ScenarioBeta=="inverse"){
        beta[draw1] <- 1; beta[draw2] <- 2; beta[draw3] <- 3 # Beta for each group
      }
      if(ScenarioBeta=="constant"){
        beta[draw1] <- 2; beta[draw2] <- 2; beta[draw3] <- 2 # Beta for each group
      }
      beta <- as.matrix(unlist(beta)) # Final Beta vector
      ### Regarding p-values
      XB <- matrix_X %*% beta # Epsilon Matrix (quantitative)
      sigma2 <- ( (r2-1) * sum((XB-mean(XB))^2) ) / ( r2*(2-size_n) )
      epsilon <- as.matrix(rnorm(size_n, mean=0, sd=sqrt(sigma2) ))
      Y <- XB + epsilon # Computation of Y (quantitative)
      pvalues <- apply(matrix_X, 2, quantitative_association_test,phenotype=Y) # Computation of quantitative Pj = Pvalues
    } else if(Design=="case_control"){ # If case_control design
      ### Beta building
      if(ScenarioBeta=="reference"){
        beta[draw1] <- log(1.8); beta[draw2] <- log(1.5); beta[draw3] <- log(1.3) # Beta for each group
      }
      if(ScenarioBeta=="inverse"){
        beta[draw1] <- log(1.3); beta[draw2] <- log(1.5); beta[draw3] <- log(1.8) # Beta for each group
      }
      if(ScenarioBeta=="constant"){
        beta[draw1] <- log(1.5); beta[draw2] <- log(1.5); beta[draw3] <- log(1.5) # Beta for each group
      }
      beta <- as.matrix(unlist(beta)) # Final Beta vector
      ### Choice of bêta0
      list_ratio <- c()
      possible_values_beta0<-seq(0,-50,-0.1)
      for (beta0 in possible_values_beta0) {
        eXB <- exp(beta0+matrix_X%*%beta)
        Y <- rbinom(size_n, 1, (eXB/ (1+eXB)) )
        Y <- as.factor(Y)
        list_ratio[(length(list_ratio)+1)]<-c(prop.table(table(Y))[2])
      }
      beta0 <- possible_values_beta0[order(abs(list_ratio-0.5))[1]]
      eXB <- exp(beta0+matrix_X%*%beta)
      Y <- rbinom(size_n, 1, (eXB/ (1+eXB)) ) # Computation of Y:
      Y <- as.factor(Y)
      pvalues <- apply(matrix_X, 2, binary_association_test,phenotype=Y) # Computation of case_control Pj = Pvalues
    }
    
    ## Define the list of drawn correlated SNPs
    List_SNPs_Drawn_Correlated[[iteration]] <- lapply(draw,correlationsofSNPs,
                                                      CorrMatrix=mcor, Threshold=(0.8))
    reg_covariates[[iteration]]<-covariates
    reg_pvalues[[iteration]]<-pvalues
    reg_draw[[iteration]]<-draw
  }
  
  return(list(reg_pvalues, reg_covariates, List_SNPs_Drawn_Correlated,correlation_matrix, reg_draw))
}

correlationsofSNPs<-function(List_Drawn_SNPs, CorrMatrix, Threshold){
  # Allows you to select the correlated SNPs with the drawn SNPs
  Corr_SNPs<-which(CorrMatrix[,List_Drawn_SNPs]>=Threshold)
  return(Corr_SNPs)
}

#Parameters
nb_iteration <- 2
size_n <- 2000
r2 <- 0.2 #coefficient of determination
path_out <- c("~/")

#Main
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
          for(m1 in c(5,10,15,20,25,50,100,150)){
            #m1 is the number of causal variants
            Config <- generate_ind_data(nb_iteration=nb_iteration,Design=design,ScenarioBeta=scenario,size_n=size_n,size_m=size_m,m1=m1,R2=r2)
            reg_pvalues <- Config[[1]]
            reg_covariates <- Config[[2]]
            if (design=="quantitative") {
              file_out <- paste("LineR","n",size_n,"m",size_m,"m1",m1,"rho",vrho,"tb",blsize,"r2",r2,quantitative_beta,sep="_")
            }
            if(design=="case_control"){
              file_out <- paste("Logit","n",size_n,"m",size_m,"m1",m1,"rho",vrho,"tb",blsize,binary_beta,sep="_")
            }
            save(reg_covariates,reg_pvalues, file=paste(path_out,file_out,sep=""))
          }
        }
      }
    }
  }
  if(Case==c("correlation")){ #case with correlations between SNPs
    size_m <- 8000 #total number of null hypotheses tested
    blsize <- 10 #bloc size
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
          for (vrho in c(0.10,0.20,0.35,0.5,0.75)) {
            #vrho is the correlation value between variants
            
            Config <- generate_corr_data(nb_iteration=nb_iteration,Design=design,ScenarioBeta=scenario,
                                         size_n=size_n,size_m=size_m,m1=m1,blsize=blsize,vrho=vrho,R2=r2)
            reg_pvalues <- Config[[1]]
            reg_covariates <- Config[[2]]
            reg_SNPs_Drawn_Correlated <- Config[[3]]
            reg_correlation_all <- Config[[4]]
            reg_draw <- Config[[5]]
            
            if (design=="quantitative") {
              file_out <- paste("LineR","n",size_n,"m",size_m,"m1",m1,"rho",vrho,"tb",blsize, "r2",r2, quantitative_beta ,sep="_")
            }
            if(design=="case_control"){
              file_out <- paste("Logit","n",size_n,"m",size_m,"m1",m1,"rho",vrho,"tb",blsize, binary_beta,sep="_")
            }
            save(reg_covariates,reg_pvalues,reg_SNPs_Drawn_Correlated,reg_correlation_all, reg_draw, file=paste(path_out,file_out,sep=""))
          }
        }
      }
    }
  }
}


