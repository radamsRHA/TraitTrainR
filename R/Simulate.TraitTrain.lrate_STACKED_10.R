#' Simulate.TraitTrain.lrate_STACKED_11: Function to simulate trait data under the STACK (+lrates) lrate model
#'
#' This function returns a Matrix of simulated trait data. Species are rows, traits are columns in matrix.TraitData 
#' @param handle.Phylogeny Phylogeny used to simulate training data under lrate model
#' @param numeric.NumberOfSpecies Numeric counting the number of species
#' @param numeric.Sig2 Numeric Value of the sigsg (evolutionary rate) parmater of the lrate model
#' @param numeric.AncestralState Numerical value of the z0 (ancestral state) parameter of the lrate model
#' @param vector.lrate_node Vector that includes all nodes (can be multiple) for the lrate model
#' @param vector.lrate_rate Vector that includes all rates (can be multiple) for the lrate model
#' @param vector.STACK_lrate_nodes Vector that includes all nodes for a STACKED lrates model (on top of lrate)
#' @param vector.STACK_lrate_rates Vector that includes all rate values for each node in vector.STACK_lrate_nodes for a STACKED lrates model (on top of lrate)
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
#' matrix.lrate_node <- matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T)
#' matrix.lrate_rate <- matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T)
#' matrix.STACK_lrate_nodes <- matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T)
#' matrix.STACK_lrate_rates <- matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T)
#' list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3)} # three traits. Different rates for different traits can be specified here. 
#' #list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- 1} # single trait
#' 
#' #############################
#' # append to simulation list #
#' #############################
#' list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "lrate", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix.STACK_lrate_nodes, 
#'                                           matrix.STACK_lrate_rates = matrix.STACK_lrate_rates, 
#'                                           matrix.lrate_node = matrix.lrate_node, 
#'                                           matrix.lrate_rate = matrix.lrate_rate, 
#'                                           list.Rmatrix = list.Rmatrix)
#' 
#' ###############################################################
#' # simulate under lrate model with STACKED 10 (lrate + lrates) #
#' ###############################################################
#' Simulate.TraitTrain.lrate_STACKED_10(handle.Phylogeny = handle.PrimatePhylogeny, 
#'                                      numeric.Sig2 = list.SimulationModelSettings[[1]]$vector.Sig2[1], 
#'                                      numeric.AncestralState = list.SimulationModelSettings[[1]]$vector.AncestralState[1], 
#'                                      vector.lrate_node = list.SimulationModelSettings[[1]]$matrix.lrate_node[,1], 
#'                                      vector.lrate_rate = list.SimulationModelSettings[[1]]$matrix.lrate_rate[,1], 
#'                                      vector.STACK_lrate_nodes = list.SimulationModelSettings[[1]]$matrix.STACK_lrate_nodes[,1], 
#'                                      vector.STACK_lrate_rates = list.SimulationModelSettings[[1]]$matrix.STACK_lrate_rates[,1], 
#'                                      matrix.R = list.SimulationModelSettings[[1]]$list.Rmatrix[[1]])

########################################
# Simulate.TraitTrain.lrate_STACKED_10 #
########################################
Simulate.TraitTrain.lrate_STACKED_10 <- function(handle.Phylogeny, numeric.Sig2, numeric.AncestralState, vector.lrate_node, vector.lrate_rate, vector.STACK_lrate_nodes, vector.STACK_lrate_rates, matrix.R){
  
  #######################
  # transform via lrate #
  #######################
  handle.Phylogeny_TRANSFORM_lrate <- phytools::rescale(handle.Phylogeny, "lrate")(sigsq = numeric.Sig2, node = vector.lrate_node, rate = vector.lrate_rate)
  
  ##############################
  # transform via STACK lrates #
  ##############################
  handle.Phylogeny_TRANSFORM_lrate_lrates <- phytools::rescale(handle.Phylogeny_TRANSFORM_lrate, "lrate")(node = vector.STACK_lrate_nodes, rate = vector.STACK_lrate_rates)
  
  ##################
  # simulate trait #
  ##################
  matrix.TraitData <- as.matrix(geiger::sim.char(phy = handle.Phylogeny_TRANSFORM_lrate_lrates, par = matrix.R, nsim = 1, model = "BM", root = numeric.AncestralState)[,,1])
  
  return(matrix.TraitData)
}
