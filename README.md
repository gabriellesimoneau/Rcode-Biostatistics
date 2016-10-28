# Scripts and material to reproduce the results in *Non-regular Inference for Dynamic Weighted Ordinary Least Squares: Understanding the Impact of Solid Food Intake in Infancy on Childhood Weight*

A dynamic treatment regime (DTR) is a set of decision rules to be applied across multiple stages of treatments. The decisions are tailored to individuals, by inputing an individual?s observed characteristics and outputting a treatment decision at each stage for that individual. Dynamic weighted ordinary least squares (dWOLS) is a theoretically robust and easily implementable method for estimating an optimal DTR. As many related DTR methods, the dWOLS treatment effects estimators can be non-regular when true treatment effects are zero or very small, which results in invalid Wald-type or standard bootstrap confidence intervals. We propose the use the *m*-out-of-*n* bootstrap with dWOLS for valid inferences of optimal DTR. We demonstrate our methodology in a real application, analyzing the effect of diet in infancy on measures of weight and body size in later childhood.

## Simulations
The objective of our simulation is to evaluate the performance of the *m*-out-of-*n* bootstrap with dWOLS as method of analysis. Here, we provide all R scripts necessary to reproduce the results found in the paper, specifically Table 2 and Figure 1. 

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

