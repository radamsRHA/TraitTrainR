#' Simulate.TraitTrain.depth_STACKED_10: Function to simulate trait data under the STACK (+lrates) depth model
#'
#' This function returns a Matrix of simulated trait data. Species are rows, traits are columns in matrix.TraitData 
#' @param handle.Phylogeny Phylogeny used to simulate training data under depth model
#' @param numeric.AncestralState Numerical value of the z0 (ancestral state) parameter of the depth model
#' @param vector.STACK_lrate_nodes Vector that includes all nodes for a STACKED lrates model (on top of depth)
#' @param vector.STACK_lrate_rates Vector that includes all rate values for each node in vector.STACK_lrate_nodes for a STACKED lrates model (on top of depth)
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
#' matrix.STACK_lrate_nodes <- matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T)
#' matrix.STACK_lrate_rates <- matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T)
#' list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3)} # three traits. Different rates for different traits can be specified here. 
#' #list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- 1} # single trait
#' 
#' #############################
#' # append to simulation list #
#' #############################
#' list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "depth", 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix.STACK_lrate_nodes, matrix.STACK_lrate_rates = matrix.STACK_lrate_rates, 
#'                                           vector.Depth = runif(n = numeric.NumberTrainingReps, min = 0, max = 10),
#'                                           list.Rmatrix = list.Rmatrix)
#' 
#' ############################################################
#' # simulate under depth model with STACKED 10 (BM + lrates) #
#' ############################################################
#' Simulate.TraitTrain.depth_STACKED_10(handle.Phylogeny = handle.PrimatePhylogeny, 
#'                                      numeric.AncestralState = list.SimulationModelSettings[[1]]$vector.AncestralState[1], 
#'                                      numeric.depth = list.SimulationModelSettings[[1]]$vector.Depth[1],
#'                                      vector.STACK_lrate_nodes = list.SimulationModelSettings[[1]]$matrix.STACK_lrate_nodes[,1], 
#'                                      vector.STACK_lrate_rates = list.SimulationModelSettings[[1]]$matrix.STACK_lrate_rates[,1], 
#'                                      matrix.R = list.SimulationModelSettings[[1]]$list.Rmatrix[[1]])

#######################################
# Simulate.TraitTrain.depth_STACKED_10 #
########################################
Simulate.TraitTrain.depth_STACKED_10 <- function(handle.Phylogeny, numeric.depth, numeric.AncestralState, vector.STACK_lrate_nodes, vector.STACK_lrate_rates, matrix.R){
  
  #######################
  # transform via depth #
  #######################
  handle.Phylogeny_TRANSFORM_depth <- rescale(handle.Phylogeny, "depth")(depth = numeric.depth)
  
  ##############################
  # transform via STACK lrates #
  ##############################
  handle.Phylogeny_TRANSFORM_depth_lrates <- rescale(handle.Phylogeny_TRANSFORM_depth, "lrate")(node = vector.STACK_lrate_nodes, rate = vector.STACK_lrate_rates)
  
  ##################
  # simulate trait #
  ##################
  matrix.TraitData <- as.matrix(sim.char(phy = handle.Phylogeny_TRANSFORM_depth_lrates, par = matrix.R, nsim = 1, model = "BM", root = numeric.AncestralState)[,,1])
  
  return(matrix.TraitData)
}
