#' Simulate.TraitTrain.nrate_STACKED_01: Function to simulate trait data under the STACK (+AncShift) nrate model
#'
#' This function returns a Matrix of simulated trait data. Species are rows, traits are columns in matrix.TraitData 
#' @param handle.Phylogeny Phylogeny used to simulate training data under nrate model
#' @param numeric.NumberOfSpecies Numeric counting the number of species
#' @param numeric.Sig2 Numeric Value of the sigsg (evolutionary rate) parmater of the nrate model
#' @param numeric.AncestralState Numerical value of the z0 (ancestral state) parameter of the nrate model
#' @param vector.nrate_time Vector that includes all times (can be multiple) for the nrate model
#' @param vector.nrate_rate Vector that includes all rates (can be multiple) for the nrate model
#' @param vector.STACK_AncShiftNode Vector that includes all nodes for a STACKED AncShift model (on top of nrate)
#' @param vector.STACK_AncShiftValue Vector that includes all shift values for each node in vector.STACK_AncShiftNode for a STACKED AncShift model (on top of nrate)
#' @param matrix.R Matrix specifying the among trait covariance for p traits. Can also be a single value (1) for a single trait
#' @return matrix.TraitData Matrix of simulated trait data. Species are rows, traits are columns
#' @export
#' @examples
#' #####################
#' # load dependencies #
#' #####################
#' library(geiger); library(phytools); library(TraitTrainR)
#' 
#' #########################
#' # get example phylogeny #
#' #########################
#' handle.PrimatePhylogeny <- read.tree(text = "((((((((human: 6, chimp:6): 1, gorilla: 7): 7, orangutan: 14): 11, macaque: 25): 64, mouse: 89): 91, opossum: 180): 20, platypus: 200): 110, chicken: 310);")
#' 
#' #############################
#' # Simulation Model Settings #
#' #############################
#' list.SimulationModelSettings <- list() # define an empty model list
#' 
#' ##########################
#' # SET SIMULATION MODEL 1 #
#' ##########################
#' numeric.NumberTrainingReps <- 10 # same number of replicates for all models in list.SimulationModelSettings
#' matrix.nrate_time <- matrix(replicate(n = numeric.NumberTrainingReps, runif(n = 1, min = 0, max = 1)), nrow = T) # single time slice
#' matrix.nrate_rate <- matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T) # single time slice + rate
#' matrix.STACK_AncShiftNode <- replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F))
#' matrix.STACK_AncShiftValue <- replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000))
#' list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3)} # three traits. Different rates for different traits can be specified here. 
#' #list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- 1} # single trait
#' 
#' #############################
#' # append to simulation list #
#' #############################
#' list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "nrate", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_AncShiftNode = matrix.STACK_AncShiftNode, 
#'                                           matrix.STACK_AncShiftValue = matrix.STACK_AncShiftValue, 
#'                                           matrix.nrate_rate = matrix.nrate_rate, 
#'                                           matrix.nrate_time = matrix.nrate_time, 
#'                                           list.Rmatrix = list.Rmatrix)
#' 
#' #################################################################
#' # simulate under nrate model with STACKED 01 (nrate + AncShift) #
#' #################################################################
#' Simulate.TraitTrain.nrate_STACKED_01(numeric.NumberOfSpecies = length(handle.PrimatePhylogeny$tip.label), 
#'                                      handle.Phylogeny = handle.PrimatePhylogeny, 
#'                                      numeric.Sig2 = list.SimulationModelSettings[[1]]$vector.Sig2[1], 
#'                                      numeric.AncestralState = list.SimulationModelSettings[[1]]$vector.AncestralState[1], 
#'                                      vector.nrate_time = list.SimulationModelSettings[[1]]$matrix.nrate_time[,1], 
#'                                      vector.nrate_rate = list.SimulationModelSettings[[1]]$matrix.nrate_rate[,1], 
#'                                      vector.STACK_AncShiftNode = list.SimulationModelSettings[[1]]$matrix.STACK_AncShiftNode[,1], 
#'                                      vector.STACK_AncShiftValue = list.SimulationModelSettings[[1]]$matrix.STACK_AncShiftValue[,1], 
#'                                      matrix.R = list.SimulationModelSettings[[1]]$list.Rmatrix[[1]])

########################################
# Simulate.TraitTrain.nrate_STACKED_01 #
#######################################
Simulate.TraitTrain.nrate_STACKED_01 <- function(numeric.NumberOfSpecies, handle.Phylogeny, numeric.Sig2, numeric.AncestralState, vector.nrate_time, vector.nrate_rate, vector.STACK_AncShiftNode, vector.STACK_AncShiftValue, matrix.R){
  
  ######################
  # transform via nrate #
  #######################
  handle.Phylogeny_TRANSFORM_nrate <- phytools::rescale(handle.Phylogeny, "nrate")(sigsq = numeric.Sig2, time = vector.nrate_time, rate = vector.nrate_rate)
  
  ##################
  # simulate trait #
  ##################
  matrix.TraitData <- as.matrix(geiger::sim.char(phy = handle.Phylogeny_TRANSFORM_nrate, par = matrix.R, nsim = 1, model = "BM", root = numeric.AncestralState)[,,1])
  
  ####################
  # add shift values #
  ####################
  for (i in 1:length(vector.STACK_AncShiftNode)){
    numeric.ShiftedNode <- vector.STACK_AncShiftNode[i]
    numeric.ShiftedValue <- vector.STACK_AncShiftValue[i]
    vector.Descendents <- getDescendants(tree = handle.Phylogeny, numeric.ShiftedNode)
    vector.Descendents <- vector.Descendents[vector.Descendents <= numeric.NumberOfSpecies]
    matrix.TraitData[vector.Descendents,] <- matrix.TraitData[vector.Descendents,] + numeric.ShiftedValue
  }
  return(matrix.TraitData)
}
