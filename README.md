# Scripts and material to reproduce the results in *Non-regular Inference for Dynamic Weighted Ordinary Least Squares: Understanding the Impact of Solid Food Intake in Infancy on Childhood Weight*

A dynamic treatment regime (DTR) is a set of decision rules to be applied across multiple stages of treatments. The decisions are tailored to individuals, by inputing an individual?s observed characteristics and outputting a treatment decision at each stage for that individual. Dynamic weighted ordinary least squares (dWOLS) is a theoretically robust and easily implementable method for estimating an optimal DTR. As many related DTR methods, the dWOLS treatment effects estimators can be non-regular when true treatment effects are zero or very small, which results in invalid Wald-type or standard bootstrap confidence intervals. We propose the use the *m*-out-of-*n* bootstrap with dWOLS for valid inferences of optimal DTR. We demonstrate our methodology in a real application, analyzing the effect of diet in infancy on measures of weight and body size in later childhood.

## Simulations
The objective of our simulation is to evaluate the performance of the *m*-out-of-*n* bootstrap with dWOLS as method of analysis. Here, we provide all R scripts necessary to reproduce the results found in the paper, specifically in Table 2 and Figure 1. 

We consider 9 simulation scenarios, each with 1000 simulated datasets of size *n=300*. The number of bootstrap replicates to *B=1000*, and use *B1 = B2 = 500* replications for the double bootstrap procedure. Raw results can be obtained as following:
* For the standard bootstrap, run the script `nn_allscenarios.R`. At completion, 9 CSV files (named *nn_psi10_scenario1.csv*, *nn_psi10_scenario2.csv*,...) are outputted. Each row shows the results for one of the 1000 simulated datasets under that scenario. The first column `psi10` contains the main effect estimate in the simulated dataset, and the next 1000 columns contain the main effect estimate in each bootstrap resample. The CSV files should have dimension 1000 rows x 1001 columns, plus a header.
* For the *m*-out-of-*n* bootstrap with fixed *alpha = 0.05*, run the script `mn0.05_allscenarios.R`. The output also consists of 9 CSV files (named *mn0.05_psi10_scenario1.csv*, *mn0.05_psi10_scenario2.csv*,...). Each row shows the results for one of the 1000 simulated datasets. The first column `psi10` contains the main effect estimate in the simulated dataset, the second column `phat` contains the estimate of non-regularity in that dataset, the third column `m` shows the resample size for the dataset, and the next 1000 columns contain the main effect estimate for each bootstrap resample. The CSV files should have dimension 1000 rows x 1004 columns, plus a header.
* For the *m*-out-of-*n* bootstrap with fixed *alpha = 0.1*, run the script `mn0.1_allscenarios.R`. The output also consists of 9 CSV files (named *mn0.1_psi10_scenario1.csv*, *mn0.1_psi10_scenario2.csv*,...) with the same format as described above.
* For the *m*-out-of-*n* bootstrap with *alpha* chosen in a data-driven manner, run the scripts `mnad_scenarioi_1.R` and `mnad_scenarioi_2.R` for each of the 9 scenarios, changing the value of *i* and the seed (see below). For each scenario, the output is 2 CSV files (named *mnad_psi10_scenario1_1.csv*, *mnad_psi10_scenario1_2.csv*,...), each showing the results for 500 simulated datasets. The 2 CSV files have a similar format as described above, plus one column `alpha` showing the value of *alpha* selected with the double bootstrap procedure. For each scenario, the two resulting files are concatenated into one CSV files (named `mnad_psi10_scenario1.R`, `mnad_psi10_scenario2.R`,...) using the following lines of code 
```{r eval=FALSE}
for(i in 1:9)
{
   name1 <- paste("mnad_psi10_scenario", paste(i), paste("_1.csv"), sep = "")
   name2 <- paste("mnad_psi10_scenario", paste(i), paste("_2.csv"), sep = "")
   dat1 <- read.csv(name1, header = TRUE)
   dat2 <- read.csv(name2, header = TRUE)
   dat <- rbind(dat1, dat2)
   write.csv(dat, file = paste("mnad_psi10_scenario", paste(i), paste(".csv"), sep = ""), row.names=FALSE)
}
```
All resulting CSV files can be found in the folder `Simulations/Simulation results`.

Scenario 	| Seed 1	| Seed 2
-----------  	| ------------	| ------------
1		| 18513	| 543
2		| 2101	| 19002
3		| 1070	| 45671
4		| 1990	| 345675
5		| 6980	| 21
6		| 29340	| 113
7		| 777		| 11200
8		| 10010	| 64
9		| 1		| 18763

To reproduce the results in Table 2, run the script `simulation_results.R`. The matrix *Table 2* shows the results presented in the corresponding table in the paper. The matrices `Table 3` and `Table 4` respectively show the corresponding coverage rates and the average confidence interval width used to produce Figure 1. To reproduce Figure 1, run the script `figure1.R`.

## Application
In the paper, we demonstrate the proposed methodology using the data from the PROmotion of Breastfeeding Intervention Trial (PROBIT). As the PROBIT dataset can not be shared, we create an example dataset and sample outputs to show how to reproduce the analyses presented in the paper.

