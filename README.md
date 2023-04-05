# Introduction

I am Thet H. Nyein, a graduate student at the University of Calgary. I have created this GitHub repository to share the code that 
I have developed as part of my Master's thesis. My research project is centered around the implementation of data subset-based techniques to 
infer spatial individual-level epidemic models.

# Individual Level Models

The general form of the ILM epidemic model is presented below. In this research, we use the SIR version of the general ILM. 
This assumes that at any given point in time $t$, either individual $i$ is susceptible to disease $(S)$, individual $i$ is infected with the disease $(I)$, 
or individual $i$ is recovered from the disease $(R)$.Individuals move through the states, $S \rightarrow I \rightarrow R$.
For each time point, $t=1,2,...,\infty$, a susceptible or infected is removed. 
Individual i is said to be in the set $S(t)$, $I(t)$ respectively.

Let $P(i,t)$ be the probability that a susceptible individual $i$ is infected in the continuous interval $[t,t+1)$ under the ILM framework it is given by 

$$P(i,t)=1-\exp[{-\Omega_{S}(i)\sum_{j \in I(t)}\Omega_{T}(j)\kappa(i,j)}+\epsilon(i,t)],\epsilon(i,t)<0$$

where, $I(t)$ is the set of infected and infectious individuals at time $t$, $\Omega_{S}(i)$ 
is a function of risk factors associated with a susceptible individual i acquiring the disease (susceptibility), 
$\Omega_{T}(j)$ is a function of risk factor associated with an infectious individual j (transmissibility), 
$\kappa(i,j)$ is an infection kernel that contains risk factors associated with pairs of susceptible and infectious individuals. $\epsilon(i,t)$ 
represents some infection process otherwise unexplained by the model.It is often set to be zero, or assumed to be constant over time and individuals. 
The likelihood over discrete time points $t=1,...,T$ is given by:

$$f(D \mid \theta)=\prod_{i=1}^{T-1}f_t(S(t),I(t),R(t)) \mid \theta)$$


$$f(D \mid \theta)=\prod_{i\in I(t+1)\backslash I(t)}P(i,t)\prod_{i\in S(t+1)}[1-P(i,t)] $$

where $i\in I(t+1)\backslash I(t)$ is the set of infectious individuals who are newly infected in time $t+1$, and $i\in S(t+1)$ is the set of susceptible individuals at time $t+1$.

## Spatial ILM Epidemic Model
In this research, the following spatial ILM epidemic model is considered. Let $\Omega_{T}(j)=1$, $\Omega_{S}(i)=\alpha$ and $\epsilon(i,t)=0$. We then let $k(i,j)=d_{ij}^{-\beta}$,  resulting in a power-law spatial ILM,

$$P(i,t)=1-\exp[{-\alpha\sum_{j \in I(t)} d_{ij}^{-\beta}}] \qquad \alpha,\beta>0$$

where $d_{ij}$ is the Euclidean distance between a susceptible individual $i$ and an infectious individual $j$, 
$\alpha$ is the susceptibility parameter and $\beta$ is the spatial parameter.

# Bayesian MCMC Approach
In this research, we utilize a Bayesian Markov Chain Monte Carlo (MCMC) framework, specifically the Metropolis-Hastings Random Walk MCMC, 
to conduct inference for the parameters $\alpha$ and $\beta$. To implement this framework, we utilize the "epimcmc" function from the EpiILM package in R.
To ensure the convergence of the MCMC chain, we used Geweke's diagnostic method. 

# Absolute Bias in ILMs 

To evaluate the performance of each subset method, we generate 20 epidemics for each one. 
We then explore the absolute biases of the posterior means $(E(.))$ and standard deviations $(SD(.))$ of the parameters $\alpha$ and $\beta$. 
The absolute bias is used to measure the absolute difference between the parameter estimate and the true parameter.
In particular, we measure the absolute biases under the full epidemic data and under the subsetting method.

To calculate the absolute biases of the posterior means and standard deviations, we use the following expressions:


$$
abs[bias(\alpha)]= abs[E_{SS}(\alpha)-E(\alpha)]
$$

$$
abs[bias(\beta)]= abs[E_{SS}(\beta)-E(\beta)]
$$

$$
abs[bias(SD_{\alpha} )]= abs[SD_{SS}(\alpha)-SD(\alpha) ]
$$

$$
abs[bias(SD_{\beta})]= abs[SD_{SS}(\beta)-SD(\beta) ]
$$

where, $E_{SS}(.)$ and $SD_{SS}(.)$ are the posterior mean and standard deviation under the subsetting method, $E(.)$ and $SD(.)$ are 
the posterior mean and standard deviation under the full epidemic data.

To measure the performance of each method, we generate 20 epidemics under each subset method and calculate the median and third quantile of absolute biases. 
This analysis provides insight into the behavior of absolute biases of $\alpha$ and $\beta$ and allows us to compare the performance of different subset methods.

# Methods

To conduct inference for the parameters $\alpha$ and $\beta$, we decided to subset data spatially or in terms of individuals' position in the square grid and temporally or 
by using infection time points.


## Spatial Subset Methods 

We implemented the **k-percent center** method and observed the relationship between population size and effective subset percentage. We selected a square containing $k%$ of the original epidemic data from the center of the $\sqrt{n}\times\sqrt{n}$ square where $k=10%,20%,30%,40%,50%,64%$, and computed $E_{SS}(.)$ and $SD_{SS}(.)$ using posterior means and standard deviations of the subset.

To investigate the impact of location of the spatial subset data on absolute biases, we developed a new method referred to as \emph{k-percent ll} method. This method involves obtaining a subset from the lower left corner of the square from 0 to $\sqrt{nk}$ where $k=10\%,20\%,30\%,40\%,50\%,64\%$, and computed $E_{SS}(.)$ and $SD_{SS}(.)$ using posterior means and standard deviations of the subset.

We also expanded two more methods named **k-percent corner** method and **k-percent random** method to incorporate additional subset areas of different sizes. We computed $E_{SS}(.)$ and $SD_{SS}(.)$ using posterior means and standard deviations of each subset.

To explore the effect of including different subset percentages on inference, we developed the **average** method, which selects four subsets each containing approximately $25\%$, $36\%$, $49\%$, and $64\%$ of the population. We computed $E_{SS}(.)$ and $SD_{SS}(.)$ using posterior means and standard deviations of each subset and obtained the average of these estimates.

Finally, we examined whether incorporating considering different subset areas in the square grid would improve the **average** method. Therefore, we established two methods named, **average corner** and **average random** methods. We computed the $E_{SS}(.)$ and $SD_{SS}(.)$ using posterior means and standard deviations of each subset, obtained the average of these estimates, and find the mean of those averages. The only difference between **average corner** and **average random** method is that the subset areas are chosen differently.

## Temporal Subset Methods

For our research, we explore two cases related to the temporal inference of posterior parameter estimates of $\alpha$ and $\beta$. We investigate which infection time frame is the most useful and whether infection time frames that contain few epidemic data improve the estimation process. To delve deeper into the first topic, we generate a method called **temporal 1 subset-m**, which involves subsetting epidemic data from the minimum infection time to the desired infection time. We derive $E_SS(.)$ and $SD_{SS}(.)$ via the posterior estimates of the subset data. This method is computationally efficient, requiring only one Bayesian MCMC framework to obtain $E_SS(.)$ and $SD_{SS}(.)$.

Concerning the second scenario, we devise two methods, referred to as **temporal mean subset-h** and **temporal median subset-h** methods. These methods require temporally subsetting the epidemic data into $h$ parts. We then acquire $E_{SS}(.)$ and $SD_{SS}(.)$ by finding the average of posterior parameter estimates of $\alpha$ and $\beta$ across all $h$ groups for the **temporal mean subset-h** method and the median of those posterior estimates for the **temporal median subset-h** method. We provide equations that describe the temporal subsets of the epidemic data based on the value of $h$.

# Code Implementation

For this research, we use R and ARC (slurm cluster computing system from University of Calgary). Hence, we need an R and slurm scripts to implement the simulations. We used cluster computing since it makes the computational process faster and more efficent. For more information about R code functions such as epdidata and epimcmc, you can refer to the epiILM documentation on this channel.
