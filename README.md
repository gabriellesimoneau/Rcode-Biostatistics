# Scripts and material to reproduce the results in *Non-regular Inference for Dynamic Weighted Ordinary Least Squares: Understanding the Impact of Solid Food Intake in Infancy on Childhood Weight*

A dynamic treatment regime (DTR) is a set of decision rules to be applied across multiple stages of treatments. The decisions are tailored to individuals, by inputing an individual?s observed characteristics and outputting a treatment decision at each stage for that individual. Dynamic weighted ordinary least squares (dWOLS) is a theoretically robust and easily implementable method for estimating an optimal DTR. As many related DTR methods, the dWOLS treatment effects estimators can be non-regular when true treatment effects are zero or very small, which results in invalid Wald-type or standard bootstrap confidence intervals. We propose the use the *m*-out-of-*n* bootstrap with dWOLS for valid inferences of optimal DTR. We demonstrate our methodology in a real application, analyzing the effect of diet in infancy on measures of weight and body size in later childhood.

## Simulations
The objective of our simulation is to evaluate the performance of the *m*-out-of-*n* bootstrap with dWOLS as method of analysis. Here, we provide all R scripts necessary to reproduce the results found in the paper, specifically in Table 2 and Figure 1. 

We consider 9 simulation scenarios, each with 1000 simulated datasets of size $n=300$. The number of bootstrap replicates to $B=1000$, and use $B_1 = B_2 = 500$ replications for the double bootstrap procedure. Raw results can be obtained as following:
* For the standard bootstrap, run the script \texttt{nn_allscenarios.R}. At completion, 9 CSV files (named \texttt{nn_psi10_scenario1.R}, \texttt{nn_psi10_scenario2.R},...) are outputted. Each row shows the results for one of the 1000 simulated datasets under that scenario. The first column \texttt{psi10} contains the main effect estimate in the simulated dataset, and the next 1000 columns contain the main effect estimate in each bootstrap resample. The CSV files should have dimension 1000 rows $\times$ 1001 columns, plus a header.
* For the *m*-out-of-*n* bootstrap with fixed $\alpha = 0.05$, run the script \texttt{mn0.05_allscenarios.R}. The output also consists of 9 CSV files (named \texttt{mn0.05_psi10_scenario1.R}, \texttt{mn0.05_psi10_scenario2.R},...). Each row shows the results for one of the 1000 simulated datasets. The first column \texttt{psi10} contains the main effect estimate in the simulated dataset, the second column \texttt{phat} contains the estimate of non-regularity in that dataset, the third column \texttt{m} shows the resample size for the dataset, and the next 1000 columns contain the main effect estimate for each bootstrap resample. The CSV files should have dimension 1000 rows $\times$ 1004 columns, plus a header.
* For the *m*-out-of-*n* bootstrap with fixed $\alpha = 0.1$, run the script \texttt{mn0.1_allscenarios.R}. The output also consists of 9 CSV files (named \texttt{mn0.1_psi10_scenario1.R}, \texttt{mn0.1_psi10_scenario2.R},...) with the same format as described above.
* For the *m*-out-of-*n* bootstrap with $\alpha$ chosen in a data-driven manner, run the script \texttt{mnad_scenario_i.R} 9 times, changing the value of \texttt{i} and the seed (see below) for each scenario. For each scenario, the output is a CSV file (named \texttt{mnad_psi10_scenario1.R}, \texttt{mnad_psi10_scenario2.R},...) , with a similar format as described above, plus one column \texttt{alpha} showing the value of $\alpha$ selected with the double bootstrap procedure. 

Scenario 	| Seed
-----------  	| ------------
1		| 18513
2		| 2101
3		| 1070
4		| 1990
5		| 6980
6		| 29340
7		| 777
8		| 10010
9		| 1



## Application
In the paper, we demonstrate the proposed methodology using the data from the PROmotion of Breastfeeding Intervention Trial (PROBIT). As the PROBIT dataset can not be shared, we create an example dataset and sample outputs to show how to reproduce the analyses presented in the paper.

