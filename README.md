
---
# TraitTrainR: Undertanding trait evolution with flexible, large-scale simulations of complex models
**NOTE See the file [https://github.com/radamsRHA/TraitTrainR/blob/main/TraitTrainR_CommandManual.pdf](https://github.com/radamsRHA/TraitTrainR/blob/main/TraitTrainR_CommandManual.pdf) for detailed instructions on functions**
**NOTE See the file https://github.com/radamsRHA/TraitTrainR/blob/main/TraitTrainR_TutorialManual.pdf for detailed example Run**

## NOTE: This README serves as a quick start guide. Please see the above tutorial files for more indepth instructions and details of options and implementation
## NOTE: TraitTrainR was written in R 4.4.0 ("Puppy cup") and we recommend that version or later for installing TraitTrainR

## QUICK START CODE: simulate under Brownian Motion and Ornstein–Uhlenbeck models (see next section for more detailed information):
```
library(TraitTrainR); library(phytools); library(geiger) 

MyTree <- read.tree(text = "((A:1, B:1):1, C:2);") 
list.SimulationModelSettings <- list() 
NumberOfReps <- 5  
list.Rmatrix <- list(); for (i in 1:NumberOfReps){list.Rmatrix[[i]] <- matrix(1, nrow = 1, ncol = 1)}

list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "BM", 
                                          vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                                          vector.AncestralState = rep(1, NumberOfReps), 
                                          list.Rmatrix = list.Rmatrix)

list.SimulationModelSettings[[2]] <- list(string.SimulationModel = "OU", 
                            vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                            vector.AncestralState = rnorm(NumberOfReps),
                            vector.Alpha = runif(n = NumberOfReps, min = exp(-500), max = exp(1)),
                            list.Rmatrix = list.Rmatrix)

MySimulationResults <- TraitTrain(handle.Phylogeny = MyTree,
                       list.SimulationModelSettings = list.SimulationModelSettings,
                       logical.PIC = TRUE, logical.PROJECT = TRUE)
```

  
    
  
  
  
## Step 0: Installing R package TraitTrainR from github
The R package TraitTrainR is freely available to download and distribute from github <https://github.com/radamsRHA/TraitTrainR/>. To install and load TraitTrainR, you must first install the R package `devtools`, 


```
install.packages("devtools")
```
Now using devtools we can install `TraitTrainR` from github:

```
library(devtools)
install_github("radamsRHA/TraitTrainR")
library(TraitTrainR) # Load package 
```
`TraitTrainR` also requires the following dependencies to be installed:

```
install.packages('phytools'); library(phytools) #make sure the version of phytools is: 2.1-1
install.packages('geiger'); library(geiger) # make sure the version of geiger is 2.0.11
```

To begin using `TraitTrainR` try using the examples associated with the manual file (XXX)

## Step 1: Load TraitTrainR R package

We can use `library(TraitTrainR)` to load the R package `TraitTrainR` into your R working environment

```
################
# Load depends #
################
library(TraitTrainR)
library(phytools)
library(geiger)
```

## Step 2: Define the target phylogeny used for simulations with TraitTrainR run

We can use `read.tree` function from the R package ape to either read an external 


```{r}
MyTree <- read.tree(text = "((((((((human: 6, chimp: 6): 1, gorilla: 7): 7, orangutan: 14): 11, macaque: 25): 64, mouse: 89): 91, opossum: 180): 20, platypus: 200): 110, chicken: 310);")
plot(MyTree) # view phylogeny before TraitTrainR
```

## Step 3: Define a TraitTrainR experimental run using the `list.SimulationModelSettings`

For this demo, we will set an experimental run using all 11 target models (without stacking) with the following commands.

```{r}
list.SimulationModelSettings <- list() # define an empty model list here
NumberOfReps <- 5  #number of replicates
list.Rmatrix <- list(); for (i in 1:NumberOfReps){list.Rmatrix[[i]] <- matrix(1, nrow = 1, ncol = 1)} # this makes a list of R matrices; one matrix for each replicate of R

######################
# First model is BM  #
######################
list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "BM", 
                                          vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                                          vector.AncestralState = rep(1, NumberOfReps), 
                                          list.Rmatrix = list.Rmatrix)

#######################
# Second model is OU  #
#######################
list.SimulationModelSettings[[2]] <- list(string.SimulationModel = "OU", 
                            vector.Sig2 = rexp(n = NumberOfReps, rate = 1), # exponential distribution
                            vector.AncestralState = rnorm(NumberOfReps), # standard normal
                            vector.Alpha = runif(n = NumberOfReps, min = exp(-500), max = exp(1)), #uniform distribution for alpha, with bounds set by fitContinuous function in geiger
                            list.Rmatrix = list.Rmatrix)

######################
# third model is EB  #
######################
list.SimulationModelSettings[[3]] <- list(string.SimulationModel = "EB", 
                                     vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                                     vector.AncestralState = rnorm(NumberOfReps), # standard normal
                                     vector.A = runif(n = NumberOfReps, min = log(10^-5)/310, max = -0.000001),
                                     list.Rmatrix = list.Rmatrix)

##########################
# fourth model is kappa  #
##########################
list.SimulationModelSettings[[4]] <- list(string.SimulationModel = "kappa", 
                                     vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                                     vector.AncestralState = rnorm(NumberOfReps), # standard normal
                                     vector.kappa = runif(n = NumberOfReps, min = exp(-500), max = 1), 
                                     list.Rmatrix = list.Rmatrix)

##########################
# fifth model is lambda  #
##########################
list.SimulationModelSettings[[5]] <- list(string.SimulationModel = "lambda", 
                                      vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                                      vector.AncestralState = rnorm(NumberOfReps), # standard normal
                                      vector.lambda = runif(n = NumberOfReps, min = exp(-500), max = 1),
                                      list.Rmatrix = list.Rmatrix)

#######################
# six model is delta  #
#######################
list.SimulationModelSettings[[6]] <- list(string.SimulationModel = "delta", 
                                     vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                                     vector.AncestralState = rnorm(NumberOfReps), # standard normal
                                     vector.delta = runif(n = NumberOfReps, min = exp(-500), max = 3), 
                                     list.Rmatrix = list.Rmatrix)

###########################
# seventh model is trend  #
###########################
list.SimulationModelSettings[[7]] <- list(string.SimulationModel = "trend", 
                                     vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                                     vector.AncestralState = rnorm(NumberOfReps), # standard normal
                                     vector.slope = runif(n = NumberOfReps, min = -100, max = 100),
                                     list.Rmatrix = list.Rmatrix)

#########################
# eight model is white  #
#########################
list.SimulationModelSettings[[8]] <- list(string.SimulationModel = "white", 
                                     vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                                     vector.AncestralState = rnorm(NumberOfReps), # standard normal
                                     list.Rmatrix = list.Rmatrix)

##########################
# nineth model is depth  #
##########################
list.SimulationModelSettings[[9]] <- list(string.SimulationModel = "depth", 
                                     vector.AncestralState = rnorm(NumberOfReps), # standard normal
                                     list.Rmatrix = list.Rmatrix, 
                                     vector.depth = runif(n = NumberOfReps, min = 0, max = 10))

##########################
# tenth model is lrate  #
##########################
list.SimulationModelSettings[[10]] <- list(string.SimulationModel = "lrate", 
                                      vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                                      vector.AncestralState = rnorm(NumberOfReps), # standard normal
                                      list.Rmatrix = list.Rmatrix, 
                                      matrix.lrate_node = matrix(replicate(n = NumberOfReps, runif(n = 1, min = 0, max = 1)), nrow = T),
                                      matrix.lrate_rate = matrix(replicate(n = NumberOfReps , runif(n = 1, min = 0, max = 100)), nrow = T))

############################
# eleventh model is nrate  #
############################
list.SimulationModelSettings[[11]] <- list(string.SimulationModel = "nrate", 
                                      vector.Sig2 = rexp(n = NumberOfReps, rate = 1), 
                                      vector.AncestralState = rnorm(NumberOfReps), 
                                      list.Rmatrix = list.Rmatrix, 
                                      matrix.nrate_time = matrix(replicate(n = NumberOfReps, runif(n = 1, min = 0, max = 1)), nrow = T), 
                                      matrix.nrate_rate = matrix(replicate(n = NumberOfReps , runif(n = 1, min = 0, max = 100)), nrow = T))

```

### Step 4: Time to Rock & Roll with TraitTrainR!

After setting up the `list.SimulationModelSettings`, we are ready to simulate!
We will use the primary function of the package, which (as you've probably guessed), is called `TraitTrain`.
Check out the following command to generate an object called `MySimulationResults` that contains all of the simulation results of `TraitTrain`.


```{r}
MySimulationResults <- TraitTrain(handle.Phylogeny = MyTree,
                       list.SimulationModelSettings = list.SimulationModelSettings,
                       logical.PIC = TRUE, logical.PROJECT = TRUE)
```

