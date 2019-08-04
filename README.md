# Quantile_Regression
#Gibbs samplings implementation of quantile regression, following the algorithm proposed by Kozumi and Kobayashi.

#Using the classic dataset Auto mpg.

The file Quantile_Regression.R contains the code explained above. 

Whereas the file SimulationQuantileReg.R is a simulation for a quantile fixed parameter linear model. This time the mixture of Kozumi and 
Kobayashi is not used, therefore resulting a posterior distribution without a known closed form. Thus, a metropolis-hastings with log adaptive proposal is used in this case.
