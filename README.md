# Description

This repository contains the differents script to reproduce the simulated data in the article presenting the wBHa procedure. Data were produced to compare and evaluate the performance of the different procedures in a GWAS context. For more details on the simulation study and procedures, please refer to the article.

3 cases were considered:

* **Independent case** corresponds to the case where correlations exist between SNPs.
* **Correlation case** corresponds to the case where simple correlations exist between SNPs.
* **Semi-simulation case** corresponds to the case where correlations exist between SNPs and are based on a real dataset (data available at **HIV_data/**). 

For the first two cases, 3 scenarios were considered: 

* **Scenario 1** corresponds to a situation where rare causal variants have a higher effect than common variants (reference scenario).
* **Scenario 2** corresponds to a situation where common variants have a higher effect than less frequent variants (inverse scenario).
* **Scenario 3** corresponds to a situation where all variants have equal effects (constant scenario).

For these cases (independent and correlation cases), we considered two types of phenotype: quantitative and binary (corresponding to studies for quantitative trait and case-control design).

For the independent case, we considered different values for the total number of null hypotheses tested (size_m), and the number of causal variants (m1) such as:

* size_m $\in$ {8000, 14000, 20000} 
* m1 $\in$ {5, 10, 15, 20, 25, 50, 100, 150}

For the correlation case, we considered different values for the number of causal variants (m1), and the correlation values between SNPs (vrho) such as:

* m1 $\in$ {5, 10, 15, 20, 25, 50, 100, 150}
* vrho $\in$ {0.1, 0.2, 0.35, 0.5, 0.75}

Finally, for the semi-simulation, we considered different values for the number of causal variants (m1) such as :

* m1 $\in$ {15, 20, 25, 50, 100, 150}

For each configuration, we simulated 500 datasets in two steps: The first one where the pvalues and covariates were produced. Then a second one where the procedures studied in the article were applied on these datasets. In this step, the overall powers and FDRs for each procedure were computed and saved in a data frame (res_FDR_power). Moreover, the powers of subgroups for each procedure were computed and saved in a data frame (res_Subpower). The subgroups correspond to the 3 distinct subsets of causal variants in which the MAF were generated from the following distributions:

* Group 1 (Rare SNPs): U[0.05,0.15]
* Group 2 (Medium SNPs): U[0.15,0.25]
* Group 3 (Common SNPs): U[0.30,0.40]

The pvalues and covariates were generated from **GenerateData_FullySimu.R** script and the procedures were applied from **Appli_Procedures_FullySimu.R** script for independent and correlation cases. For semi-simulation case, both stage are available at **Semi_Simulation.R** script.

Note that the results from the application of the procedures for all data are available in the R package **wBHa**.
