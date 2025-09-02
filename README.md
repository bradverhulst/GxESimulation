# GxE Simulation

This GitHub page is a companion website to the manuscript:

Verhulst, B. (Under Review) Identifying Gene-Environment Interactions across Genome-wide, Twin, and Polygenic Risk Score (PRS) approaches

In the paper, we conduct a series of simualation studies to highlight the consistent ability to detecting GxE in genome-wide, twin, and polygenic risk score (PRS) analyses. While each method requires different types of data, we demonstrate that by simulating genome-wide twin data, we can conduct all the necessary data structures for to compare results across methods.

All of the simulation scripts that were used in the paper are provided, including scripts to construct the figures.

## Description of the function to generate GxE data for  Twin, GWAS, and PRS Data

As we cannot possibly anticipate every combination of simulated GxE data, we are providing researchers with the R functions we used to simulates the data for this project. In full disclosure, the functions we built are include options that are not discussed in the manuscript, which allows those interested in GxE to play around and simulate data and answer their own questions.

Below is a short tutorial demonstrating how to use the functions. The first thing that we need to do is load the require packages and functions. You can download TwinGxEsim_Fun.R and save it in your working directory. All of the information for conducting the analyses and generating the figures can be found in TwinGxE_demo.R.

_As this project is currently under review, any suggestions for improving the simulations would be greatly appreciated._
