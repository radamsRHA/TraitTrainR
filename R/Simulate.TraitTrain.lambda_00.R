#' Simulate.TraitTrain.lambda_00: Function to simulate trait data under the lambda model
#'
#' This function returns a Matrix of simulated trait data. Species are rows, traits are columns in matrix.TraitData 
#' @param handle.Phylogeny Phylogeny used to simulate training data under lambda model
#' @param numeric.NumberOfSpecies Numeric counting the number of species
#' @param numeric.Sig2 Numeric Value of the sigsg (evolutionary rate) parmater of the lambda model
#' @param numeric.AncestralState Numerical value of the z0 (ancestral state) parameter of the lambda model
#' @param numeric.lambda Numerical value of the lambda parameter of the lambda model
#' @param matrix.R Matrix specifying the among trait covariance for p traits. Can also be a single value (1) for a single trait
#' @return matrix.TraitData Matrix of simulated trait data. Species are rows, traits are columns
#' @export
#' @examples
#' 
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
#' list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3)} # three traits. Different rates for different traits can be specified here. 
#' #list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- 1} # single trait
#' 
#' #############################
#' # append to simulation list #
#' #############################
#' list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "lambda", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           vector.lambda = runif(n = numeric.NumberTrainingReps, min = exp(-500), max = 1), 
#'                                           list.Rmatrix = list.Rmatrix)
#' 
#' ###############################
#' # simulate under lambda model #
#' ###############################
#' Simulate.TraitTrain.lambda_00(handle.Phylogeny = handle.PrimatePhylogeny, 
#'                               numeric.Sig2 = list.SimulationModelSettings[[1]]$vector.Sig2[1], 
#'                               numeric.AncestralState = list.SimulationModelSettings[[1]]$vector.AncestralState[1], 
#'                               numeric.lambda = list.SimulationModelSettings[[1]]$vector.lambda[1],
#'                               matrix.R = list.SimulationModelSettings[[1]]$list.Rmatrix[[1]])

#################################
# Simulate.TraitTrain.lambda_00 #
#################################
Simulate.TraitTrain.lambda_00 <- function(handle.Phylogeny, numeric.Sig2, numeric.lambda, numeric.AncestralState, matrix.R){
  
  ########################
  # transform via lambda #
  ########################
  handle.Phylogeny_TRANSFORM_lambda <- phytools::rescale(handle.Phylogeny, "lambda")(sigsq = numeric.Sig2, lambda = numeric.lambda)
  
  ##################
  # simulate trait #
  ##################
  matrix.TraitData <- as.matrix(geiger::sim.char(phy = handle.Phylogeny_TRANSFORM_lambda, par = matrix.R, nsim = 1, model = "BM", root = numeric.AncestralState)[,,1])
  
  return(matrix.TraitData)
}
