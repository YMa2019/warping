A Stochastic Process Framework for Time Warping Functions

In this project, we provide a linear framework for time warping functions. This new framework can be easily utilized to generate time warping functions by using different types of stochastic process. We demonstrate the effectiveness of this new framework by using it as a new penalty in the penalized function registration, and propose an efficient gradient method to minimize the cost function. We illustrate the new penalized method using simulations that properly characterize nonuniform and correlated constraints in the time domain. Furthermore, we utilize the new framework to develop boxplot for warping functions. This boxplot can be used to estimate templates and identify warping outliers. 

This repository contains 6 main scripts, 2 function scripts and  1 R script: 

1. Generation.m
Used to simulate time warping functions without prior information using Algorithm 1.

2. fPCA.m
Used to resample time warping functions using Algorithm 2

3. CLRFemale.m
Used to resample of Berkey growth data for female

4. CLRmale.m
Used to resample of Berkey growth data for male

5. Penalty_Diag_lag.m
Main script Penalized functional alignment (diagonal case)

6. penaltyFA.m
The function alignment function which used in the main script "Penalty_Diag_lag.m'

7. fBoxplot.m
Used to build the proposed functional boxplot for time warping functions

8. Tukey Depth Calculation.Rmd
Used to calculate the Tukey depth

9. plot_warping.m
function used to plot the time warping function

Files in supplement folder contains the other necessary functions used in the main scripts. 

1. centroid.m
calculate the centroid of a convex hull

2. clr_inv.m
The inverse clr transformation

3. inhull.m
Return all the point index in the convex hull

4. laprnd.m
Simulate laplacian distributed vectors


