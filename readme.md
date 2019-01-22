
OVERVIEW:
=========

This repository contains code to reproduce the methodology established in the paper "Data-intensive multidimensional modeling of forest dynamics" in Ecological Modeling and Software, 2015: [link to preprint](https://www.biorxiv.org/content/early/2015/05/02/005009) / [link to paper](http://dx.doi.org/10.1016/j.envsoft.2015.01.010).

The code performs the transition matrix estimation of the Quebec forest inventory database. The main code is located in `script.R`, which calls helper functions in `db_glue.R` and `estimate_gibbs.R`.


Minimal example on synthetic data
===

Below, I also provide a minimal example demonstrating transition matrix inference from generated sequences using Gibbs sampling:

```R
source('estimate_gibbs.R')
require('MCMCpack') # used here to draw random numbers from the Dirichlet distribution with the 'rdirichlet' function 

test_MCMC = function(n=3,Tmax=4,m=5000,pm=0.1,pn=0.1,Hmax=1000,burnin=100)
{
  # this function generates sequences from a (random) known transition matrix and initial states, and returns the transition matrix estimated with MCMC
  # the transition matrix is used to generates m sequences of n states with times 1..Tmax
  # then, data from sequences are removed, with pm probability of missing data, pn probability of noisy data (i.e. a "wrong" state)
  # finally, gibbs sampling is used to estimate the estimated transition matrix

  # generate a random transition matrix n*n
  transition <- matrix(runif(n^2),ncol=n)
  transition <- diag(1/rowSums(transition)) %*% transition
  
  # generate a random state initial distribution
  state_ini <- rgamma(n,shape=1,scale=2)
  state_ini <- state_ini / sum(state_ini)
  
  # generate complete sequences for this initial distribution and transition matrix
  S <- matrix(NA,nrow=m,ncol=Tmax) # complete sequences
  for (k in 1:m) {
    S[k,1] <- which(cumsum(state_ini) > runif(1))[1]
    for (t in 2:Tmax) {
      S[k,t] <- which(cumsum(transition[S[k,t-1],]) > runif(1))[1]
    }
  }
  
  # remove points from sequences with probability pm, and change points from sequences with probability pn
  for (k in 1:m) {
    for (t in 2:Tmax) {
      r = runif(1)
      if (r < pm) {
        S[k,t] <- NA
      } else if (r < pm + pn) {
        S[k,t] <- sample((1:n)[-S[k,t]],1)
      }
    }
  }
  
  print('original')
  print(transition)
  
  transitions_gibbs = gibbs_sampling(S,n,Hmax=Hmax,burnin=burnin)
  median_gibbs = apply(transitions_gibbs, c(2,3), median)
  print('Gibbs deduced')
  print(median_gibbs)
  
  return(data.frame(real_trans=transition, gibbs_trans=median_gibbs))
}

test_MCMC()
```
