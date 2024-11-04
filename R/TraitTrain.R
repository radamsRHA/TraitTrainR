#' TraitTrain: Function to conduct large-scale, flexible simulations of trait data given a particular tree. 
#'
#' This function returns a list of data.frames that include: RESULTS_TRAIT (raw trait values across simulations); RESULTS_PIC (compute PICS across simulations); RESULTS_PIC_DEPTH (computed PICS across simulation using target tree scaled to depth); RESULTS_PROJECT (projected trait using input tree); RESULTS_PROJECT_DEPTH (projected traits using input tree scaled to unit depth = 1)
#' @param handle.Phylogeny Phylogeny used to simulate traits with TrainTrainR
#' @param list.SimulationModelSettings List containing the model values and parameters to be used during simulation. Must contain (at least) the following for each entry in the list: string.SimulationModel (string defines the model; must be one of c("BM", "OU", "EB", "nrate", "lrate", "trend", "lambda", "kappa", "delta", "white", "depth")); vector.Sig2 (value of the evolutionary rate parameter for each replicate); vector.AncestralState (value of the ancestral state for each replicate)
#' @param logical.PIC True/False logical indicating whether to compute PICs (or not).
#' @param logical.PROJECT True/False logical indicating whether to compute phylogenetic projections (or not)
#' @param numeric.MeasurementError Numerical value representing the variance of the sample error (assumed Normally distributed with mean = 0 , sd = sqrt(numeric.MeasurementError))
#' @return LIST A list of data.frames that include: RESULTS_TRAIT (raw trait values across simulations); RESULTS_PIC (compute PICS across simulations); RESULTS_PIC_DEPTH (computed PICS across simulation using target tree scaled to depth); RESULTS_PROJECT (projected trait using input tree); RESULTS_PROJECT_DEPTH (projected traits using input tree scaled to unit depth = 1)
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
#' numeric.NumberTrainingReps <- 2 # same number of replicates for all models in list.SimulationModelSettings
#' list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3)} # three traits. Different rates for different traits can be specified here. 
#' #list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- matrix(1, nrow = 1, ncol = 1)} # three traits. Different rates for different traits can be specified here. 
#' 
#' ######################
#' # First model is BM  #
#' ######################
#' list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "BM", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                           matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                           matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                           matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                           list.Rmatrix = list.Rmatrix)
#' 
#' #######################
#' # Second model is OU  #
#' #######################
#' list.SimulationModelSettings[[2]] <- list(string.SimulationModel = "OU", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                           matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                           matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                           matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                           vector.Alpha = runif(n = numeric.NumberTrainingReps, min = exp(-500), max = exp(1)),
#'                                           list.Rmatrix = list.Rmatrix)
#' 
#' 
#' ######################
#' # third model is EB  #
#' ######################
#' list.SimulationModelSettings[[3]] <- list(string.SimulationModel = "EB", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                           matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                           matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                           matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                           vector.A = runif(n = numeric.NumberTrainingReps, min = log(10^-5)/310, max = -0.000001),
#'                                           list.Rmatrix = list.Rmatrix)
#' 
#' ##########################
#' # fourth model is kappa  #
#' ##########################
#' list.SimulationModelSettings[[4]] <- list(string.SimulationModel = "kappa", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                           matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                           matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                           matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                           vector.kappa = runif(n = numeric.NumberTrainingReps, min = exp(-500), max = 1), 
#'                                           list.Rmatrix = list.Rmatrix)
#' ##########################
#' # fifth model is lambda  #
#' ##########################
#' list.SimulationModelSettings[[5]] <- list(string.SimulationModel = "lambda", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                           matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                           matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                           matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                           vector.lambda = runif(n = numeric.NumberTrainingReps, min = exp(-500), max = 1), 
#'                                           list.Rmatrix = list.Rmatrix)
#' 
#' #######################
#' # six model is delta  #
#' #######################
#' list.SimulationModelSettings[[6]] <- list(string.SimulationModel = "delta", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                           matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                           matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                           matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                           vector.delta = runif(n = numeric.NumberTrainingReps, min = exp(-500), max = 3), 
#'                                           list.Rmatrix = list.Rmatrix)
#' 
#' 
#' ###########################
#' # seventh model is trend  #
#' ###########################
#' list.SimulationModelSettings[[7]] <- list(string.SimulationModel = "trend", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                           matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                           matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                           matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                           vector.slope = rexp(n = numeric.NumberTrainingReps, rate = 1),
#'                                           list.Rmatrix = list.Rmatrix)
#' #########################
#' # eight model is white  #
#' #########################
#' list.SimulationModelSettings[[8]] <- list(string.SimulationModel = "white", 
#'                                           vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                           matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                           matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                           matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                           list.Rmatrix = list.Rmatrix)
#' 
#' 
#' ##########################
#' # nineth model is depth  #
#' ##########################
#' list.SimulationModelSettings[[9]] <- list(string.SimulationModel = "depth", 
#'                                           vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                           matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                           matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                           matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                           matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                           list.Rmatrix = list.Rmatrix, 
#'                                           vector.depth = runif(n = numeric.NumberTrainingReps, min = 0, max = 10))
#' 
#' ##########################
#' # tenth model is lrate  #
#' ##########################
#' list.SimulationModelSettings[[10]] <- list(string.SimulationModel = "lrate", 
#'                                            vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                            vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                            matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                            matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                            matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                            matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                            list.Rmatrix = list.Rmatrix, 
#'                                            matrix.lrate_node = matrix(replicate(n = numeric.NumberTrainingReps, runif(n = 1, min = 0, max = 1)), nrow = T),
#'                                            matrix.lrate_rate = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T))
#' 
#' ############################
#' # eleventh model is nrate  #
#' ############################
#' list.SimulationModelSettings[[11]] <- list(string.SimulationModel = "nrate", 
#'                                            vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
#'                                            vector.AncestralState = runif(n = numeric.NumberTrainingReps, min = 0, max = 100), 
#'                                            matrix.STACK_lrate_nodes = matrix(replicate(n = numeric.NumberTrainingReps, sample(size = 1, x = 1:16, replace = F)), nrow = T), 
#'                                            matrix.STACK_lrate_rates = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T), 
#'                                            matrix.STACK_AncShiftNode = replicate(n = numeric.NumberTrainingReps, sample(size = 2, x = 1:16, replace = F)), 
#'                                            matrix.STACK_AncShiftValue = replicate(n = numeric.NumberTrainingReps , rnorm(n = 2, mean = 10000)), 
#'                                            list.Rmatrix = list.Rmatrix, 
#'                                            matrix.nrate_time = matrix(replicate(n = numeric.NumberTrainingReps, runif(n = 1, min = 0, max = 1)), nrow = T),
#'                                            matrix.nrate_rate = matrix(replicate(n = numeric.NumberTrainingReps , runif(n = 1, min = 0, max = 100)), nrow = T))
#' 
#' ####################
#' # SIMULATE TRAITS! #
#' ####################
#' handle.RESULTS_TEST <- TraitTrain(handle.Phylogeny = handle.PrimatePhylogeny, 
#'                                         list.SimulationModelSettings = list.SimulationModelSettings, 
#'                                         logical.PIC = T, logical.PROJECT = T)
#' 


##############
# TraitTrain #
##############
TraitTrain <- function(handle.Phylogeny, list.SimulationModelSettings, numeric.MeasurementError = 0, logical.PIC = F, logical.PROJECT = F){
  
  ##########################
  # Check rescale function #
  ##########################
#  if(exists('rescale', where='package:geiger', mode='function')){rescale <- geiger::rescale; print("Found rescale function in the installed package geiger")}
#  if(exists('rescale', where='package:phytools', mode='function')){rescale <- phytools::rescale; print("Found rescale function in the installed package phytools")}
  handle.Phylogeny_DEPTH <- phytools::rescale(handle.Phylogeny, model = "depth")(depth = 1)
  numeric.NumberOfTraits <- ncol(list.SimulationModelSettings[[1]]$list.Rmatrix[[1]])
  numeric.NumberOfSpecies <- length(handle.Phylogeny$tip.label); print(paste0("Found a single phylogeny of ", numeric.NumberOfSpecies, " species that will be used for simulation with TraitTrainR for models with ", numeric.NumberOfTraits, " traits"))
  vector.SpeciesNames <- handle.Phylogeny$tip.label
  matrix.RESULTS_PROJECTION <- matrix.RESULTS_PROJECTION_DEPTH <- NULL
  
  ####################
  # summarize models #
  ####################
  numeric.NumberOfModels <- length(list.SimulationModelSettings); print(paste0("Found a total of ", numeric.NumberOfModels, " models that will be used for simulation with TraitTrainR."))
  numeric.COUNTER <- 0
  numeric.NumberOfTotalReplicates <- 0; for (model in 1:numeric.NumberOfModels){numeric.NumberOfTotalReplicates = numeric.NumberOfTotalReplicates + length(list.SimulationModelSettings[[model]]$vector.AncestralState)}
  
  ###############
  # SET RESULTS #
  ###############
  matrix.RESULTS_TRAIT <- matrix(nrow = numeric.NumberOfTotalReplicates, ncol = (numeric.NumberOfTraits * numeric.NumberOfSpecies) + 1); colnames(matrix.RESULTS_TRAIT) <- c(paste0(rep(paste0("Trait", 1:numeric.NumberOfTraits), each = numeric.NumberOfSpecies), "_", vector.SpeciesNames), "SimulationModelNumber")
  matrix.RESULTS_PIC <- matrix(nrow = numeric.NumberOfTotalReplicates, ncol = (numeric.NumberOfTraits * (numeric.NumberOfSpecies-1)) + 1); colnames(matrix.RESULTS_PIC) <- c(paste0(rep(paste0("Trait", 1:numeric.NumberOfTraits), each = (numeric.NumberOfSpecies-1)), "_PIC", 1:(numeric.NumberOfSpecies-1)), "SimulationModelNumber")
  matrix.RESULTS_PIC_DEPTH <- matrix.RESULTS_PIC
  
  #######################
  # loop through models #
  #######################
  for (i in 1:numeric.NumberOfModels){
    
    #########
    # CLEAN #
    #########
    list.SimulationModel <- string.PrimaryModel <- vector.Sig2 <- vector.AncestralState <- CHECK.STACK_10 <- CHECK.STACK_01 <- CHECK.STACK_11 <- CHECK.STACK_00 <- NULL
    
    #############
    # GET MODEL #
    #############
    list.SimulationModel <- list.SimulationModelSettings[[i]]
    string.PrimaryModel <- list.SimulationModel$string.SimulationModel
    vector.Sig2 <- list.SimulationModel$vector.Sig2
    vector.AncestralState <- list.SimulationModel$vector.AncestralState
    numeric.NumberModelReplicates <- length(vector.AncestralState); print(paste0("Simulating ", numeric.NumberModelReplicates, " replicates for the ", i, " model, which belongs to primary model class: ", string.PrimaryModel))
    if (numeric.NumberModelReplicates != length(vector.Sig2) & string.PrimaryModel != "depth"){stop("The number of rate parameters and the number of ancestral states for the " , i, " model do not mach")}
    
    ##########################
    # CHECK IF STACKED MODEL #
    ##########################
    CHECK.STACK_00 <- is.null(list.SimulationModel$matrix.STACK_AncShiftNode) & is.null(list.SimulationModel$matrix.STACK_lrate_nodes)
    CHECK.STACK_10 <- !is.null(list.SimulationModel$matrix.STACK_lrate_nodes) & is.null(list.SimulationModel$matrix.STACK_AncShiftNode) 
    CHECK.STACK_01 <- is.null(list.SimulationModel$matrix.STACK_lrate_nodes) & !is.null(list.SimulationModel$matrix.STACK_AncShiftNode) 
    CHECK.STACK_11 <- !is.null(list.SimulationModel$matrix.STACK_AncShiftNode) & !is.null(list.SimulationModel$matrix.STACK_lrate_nodes)
    
    ############
    # MODEL BM #
    ############
    if (string.PrimaryModel == "BM"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.BM_00(handle.Phylogeny = handle.Phylogeny, 
                                                            numeric.Sig2 = vector.Sig2[j], 
                                                            numeric.AncestralState = vector.AncestralState[j], 
                                                            matrix.R = list.SimulationModel$list.Rmatrix[[j]])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.BM_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                    numeric.Sig2 = vector.Sig2[j], 
                                                                    numeric.AncestralState = vector.AncestralState[j], 
                                                                    matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                    numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                    vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                    vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.BM_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                    numeric.Sig2 = vector.Sig2[j], 
                                                                    numeric.AncestralState = vector.AncestralState[j], 
                                                                    matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                    vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                    vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.BM_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                    numeric.Sig2 = vector.Sig2[j], 
                                                                    numeric.AncestralState = vector.AncestralState[j], 
                                                                    matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                    numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                    vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                    vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                    vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                    vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
    
    ############
    # MODEL OU #
    ############
    if (string.PrimaryModel == "OU"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.OU_00(handle.Phylogeny = handle.Phylogeny, 
                                                            numeric.Sig2 = vector.Sig2[j], 
                                                            numeric.AncestralState = vector.AncestralState[j], 
                                                            matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                            numeric.Alpha = list.SimulationModel$vector.Alpha[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.OU_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                    numeric.Sig2 = vector.Sig2[j], 
                                                                    numeric.AncestralState = vector.AncestralState[j], 
                                                                    matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                    numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                    vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                    vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                    numeric.Alpha = list.SimulationModel$vector.Alpha[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.OU_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                    numeric.Sig2 = vector.Sig2[j], 
                                                                    numeric.AncestralState = vector.AncestralState[j], 
                                                                    matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                    vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                    vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                    numeric.Alpha = list.SimulationModel$vector.Alpha[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.OU_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                    numeric.Sig2 = vector.Sig2[j], 
                                                                    numeric.AncestralState = vector.AncestralState[j], 
                                                                    matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                    numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                    vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                    vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                    vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                    vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                    numeric.Alpha = list.SimulationModel$vector.Alpha[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
    
    ############
    # MODEL EB #
    ############
    if (string.PrimaryModel == "EB"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.EB_00(handle.Phylogeny = handle.Phylogeny, 
                                                            numeric.Sig2 = vector.Sig2[j], 
                                                            numeric.AncestralState = vector.AncestralState[j], 
                                                            matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                            numeric.A = list.SimulationModel$vector.A[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.EB_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                    numeric.Sig2 = vector.Sig2[j], 
                                                                    numeric.AncestralState = vector.AncestralState[j], 
                                                                    matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                    numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                    vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                    vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                    numeric.A = list.SimulationModel$vector.A[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.EB_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                    numeric.Sig2 = vector.Sig2[j], 
                                                                    numeric.AncestralState = vector.AncestralState[j], 
                                                                    matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                    vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                    vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                    numeric.A = list.SimulationModel$vector.A[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.EB_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                    numeric.Sig2 = vector.Sig2[j], 
                                                                    numeric.AncestralState = vector.AncestralState[j], 
                                                                    matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                    numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                    vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                    vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                    vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                    vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                    numeric.A = list.SimulationModel$vector.A[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
    
    ###############
    # MODEL kappa #
    ###############
    if (string.PrimaryModel == "kappa"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.kappa_00(handle.Phylogeny = handle.Phylogeny, 
                                                               numeric.Sig2 = vector.Sig2[j], 
                                                               numeric.AncestralState = vector.AncestralState[j], 
                                                               matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                               numeric.kappa = list.SimulationModel$vector.kappa[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.kappa_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                       numeric.kappa = list.SimulationModel$vector.kappa[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.kappa_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       numeric.kappa = list.SimulationModel$vector.kappa[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.kappa_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       numeric.kappa = list.SimulationModel$vector.kappa[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
    
    ###############
    # MODEL lambda #
    ###############
    if (string.PrimaryModel == "lambda"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.lambda_00(handle.Phylogeny = handle.Phylogeny, 
                                                                numeric.Sig2 = vector.Sig2[j], 
                                                                numeric.AncestralState = vector.AncestralState[j], 
                                                                matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                numeric.lambda = list.SimulationModel$vector.lambda[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.lambda_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                        numeric.Sig2 = vector.Sig2[j], 
                                                                        numeric.AncestralState = vector.AncestralState[j], 
                                                                        matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                        numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                        vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                        vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                        numeric.lambda = list.SimulationModel$vector.lambda[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.lambda_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                        numeric.Sig2 = vector.Sig2[j], 
                                                                        numeric.AncestralState = vector.AncestralState[j], 
                                                                        matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                        vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                        vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                        numeric.lambda = list.SimulationModel$vector.lambda[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.lambda_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                        numeric.Sig2 = vector.Sig2[j], 
                                                                        numeric.AncestralState = vector.AncestralState[j], 
                                                                        matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                        numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                        vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                        vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                        vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                        vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                        numeric.lambda = list.SimulationModel$vector.lambda[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
    
    ###############
    # MODEL delta #
    ###############
    if (string.PrimaryModel == "delta"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.delta_00(handle.Phylogeny = handle.Phylogeny, 
                                                               numeric.Sig2 = vector.Sig2[j], 
                                                               numeric.AncestralState = vector.AncestralState[j], 
                                                               matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                               numeric.delta = list.SimulationModel$vector.delta[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.delta_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                       numeric.delta = list.SimulationModel$vector.delta[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.delta_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       numeric.delta = list.SimulationModel$vector.delta[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.delta_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       numeric.delta = list.SimulationModel$vector.delta[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
    
    ###############
    # MODEL trend #
    ###############
    if (string.PrimaryModel == "trend"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.trend_00(handle.Phylogeny = handle.Phylogeny, 
                                                               numeric.Sig2 = vector.Sig2[j], 
                                                               numeric.AncestralState = vector.AncestralState[j], 
                                                               matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                               numeric.slope = list.SimulationModel$vector.slope[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.trend_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                       numeric.slope = list.SimulationModel$vector.slope[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.trend_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       numeric.slope = list.SimulationModel$vector.slope[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.trend_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       numeric.slope = list.SimulationModel$vector.slope[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
    
    ###############
    # MODEL white #
    ###############
    if (string.PrimaryModel == "white"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.white_00(handle.Phylogeny = handle.Phylogeny, 
                                                               numeric.Sig2 = vector.Sig2[j], 
                                                               numeric.AncestralState = vector.AncestralState[j], 
                                                               matrix.R = list.SimulationModel$list.Rmatrix[[j]])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.white_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.white_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.white_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
    
    ###############
    # MODEL depth #
    ###############
    if (string.PrimaryModel == "depth"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.depth_00(handle.Phylogeny = handle.Phylogeny, 
                                                               numeric.AncestralState = vector.AncestralState[j], 
                                                               matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                               numeric.depth = list.SimulationModel$vector.depth[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.depth_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j],
                                                                       numeric.depth = list.SimulationModel$vector.depth[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.depth_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j],
                                                                       numeric.depth = list.SimulationModel$vector.depth[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.depth_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       numeric.depth = list.SimulationModel$vector.depth[j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
    
    ###############
    # MODEL lrate #
    ###############
    if (string.PrimaryModel == "lrate"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.lrate_00(handle.Phylogeny = handle.Phylogeny, 
                                                               numeric.Sig2 = vector.Sig2[j], 
                                                               numeric.AncestralState = vector.AncestralState[j], 
                                                               matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                               vector.lrate_node = list.SimulationModel$matrix.lrate_node[,j], vector.lrate_rate = list.SimulationModel$matrix.lrate_rate[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.lrate_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j],
                                                                       vector.lrate_node = list.SimulationModel$matrix.lrate_node[,j], vector.lrate_rate = list.SimulationModel$matrix.lrate_rate[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.lrate_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       vector.lrate_node = list.SimulationModel$matrix.lrate_node[,j], vector.lrate_rate = list.SimulationModel$matrix.lrate_rate[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.lrate_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       vector.lrate_node = list.SimulationModel$matrix.lrate_node[,j], vector.lrate_rate = list.SimulationModel$matrix.lrate_rate[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
    
    ###############
    # MODEL nrate #
    ###############
    if (string.PrimaryModel == "nrate"){
      
      ######
      # 00 #
      ######
      if (CHECK.STACK_00 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.nrate_00(handle.Phylogeny = handle.Phylogeny, 
                                                               numeric.Sig2 = vector.Sig2[j], 
                                                               numeric.AncestralState = vector.AncestralState[j], 
                                                               matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                               vector.nrate_time = list.SimulationModel$matrix.nrate_time[,j], vector.nrate_rate = list.SimulationModel$matrix.nrate_rate[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 01 #
      ######
      if (CHECK.STACK_01 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.nrate_STACKED_01(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j],
                                                                       vector.nrate_time = list.SimulationModel$matrix.nrate_time[,j], vector.nrate_rate = list.SimulationModel$matrix.nrate_rate[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 10 #
      ######
      if (CHECK.STACK_10 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.nrate_STACKED_10(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       vector.nrate_time = list.SimulationModel$matrix.nrate_time[,j], vector.nrate_rate = list.SimulationModel$matrix.nrate_rate[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
      
      ######
      # 11 #
      ######
      if (CHECK.STACK_11 == T){
        for (j in 1:numeric.NumberModelReplicates){
          numeric.COUNTER <- numeric.COUNTER + 1
          
          #############
          # RAW TRAIT #
          #############
          matrix.TraitData_RAW <- Simulate.TraitTrain.nrate_STACKED_11(handle.Phylogeny = handle.Phylogeny, 
                                                                       numeric.Sig2 = vector.Sig2[j], 
                                                                       numeric.AncestralState = vector.AncestralState[j], 
                                                                       matrix.R = list.SimulationModel$list.Rmatrix[[j]], 
                                                                       numeric.NumberOfSpecies = numeric.NumberOfSpecies, 
                                                                       vector.STACK_AncShiftNode = list.SimulationModel$matrix.STACK_AncShiftNode[,j], 
                                                                       vector.STACK_AncShiftValue = list.SimulationModel$matrix.STACK_AncShiftValue[,j], 
                                                                       vector.STACK_lrate_nodes = list.SimulationModel$matrix.STACK_lrate_nodes[,j],
                                                                       vector.STACK_lrate_rates = list.SimulationModel$matrix.STACK_lrate_rates[,j], 
                                                                       vector.nrate_time = list.SimulationModel$matrix.nrate_time[,j], vector.nrate_rate = list.SimulationModel$matrix.nrate_rate[,j])
          matrix.TraitData_RAW <- matrix.TraitData_RAW + rnorm(n = length(matrix.TraitData_RAW), mean = 0, sd = sqrt(numeric.MeasurementError))
          matrix.RESULTS_TRAIT[numeric.COUNTER,] <- c(matrix.TraitData_RAW, i)
          
          ##################
          # CONVERT TO PIC #
          ##################
          if (logical.PIC == T){matrix.PIC_DEPTH <- matrix.PIC <- matrix(nrow = (numeric.NumberOfSpecies-1), ncol = numeric.NumberOfTraits); for (p in 1:numeric.NumberOfTraits){matrix.PIC[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny); matrix.PIC_DEPTH[,p] <- ape::pic(x = matrix.TraitData_RAW[,p], phy = handle.Phylogeny_DEPTH)}; matrix.RESULTS_PIC[numeric.COUNTER,] <- c(matrix.PIC, i); matrix.RESULTS_PIC_DEPTH[numeric.COUNTER,] <- c(matrix.PIC_DEPTH, i)}}
      }
    }
  }
  
  #######################
  # COMPUTE PROJECTIONS #
  #######################
  if (logical.PROJECT == T){
    
    ###################
    # COLLECT RESULTS #
    ###################
    list.RESULTS_PROJECTION <- list()
    list.RESULTS_PROJECTION_DEPTH <- list()
    
    ######################
    # COMPUTE P matrices #
    ######################
    "%XXX%" <- function(x, n) with(eigen(x), vectors %*% (values^n * t(vectors))) # function for projections 
    
    #################
    # Decompose ape::vcv #
    #################
    matrix.C <- ape::vcv(phy = handle.Phylogeny); matrix.U <- eigen(x = matrix.C)$vectors; matrix.W <- matrix(0, nrow = numeric.NumberOfSpecies, ncol = numeric.NumberOfSpecies); diag(matrix.W) <- eigen(x = matrix.C)$values
    matrix.P <- solve(matrix.U %*% (matrix.W %XXX% (0.5)) %*% t(matrix.U)); colnames(matrix.P) <- rownames(matrix.P) <- handle.Phylogeny$tip.label
    
    matrix.C <- ape::vcv(phy = handle.Phylogeny_DEPTH); matrix.U <- eigen(x = matrix.C)$vectors; matrix.W <- matrix(0, nrow = numeric.NumberOfSpecies, ncol = numeric.NumberOfSpecies); diag(matrix.W) <- eigen(x = matrix.C)$values
    matrix.P_DEPTH <- solve(matrix.U %*% (matrix.W %XXX% (0.5)) %*% t(matrix.U)); colnames(matrix.P_DEPTH) <- rownames(matrix.P_DEPTH) <- handle.Phylogeny$tip.label
    
    #######################
    # loop through traits #
    #######################
    for (p in 1:numeric.NumberOfTraits){
      ################
      # select trait #
      ################
      string.Search <- paste0("Trait", p)
      handle.SELECT <- t(matrix.RESULTS_TRAIT[,grep(pattern = string.Search, colnames(matrix.RESULTS_TRAIT))])
      
      ###############
      # RAW PROJECT #
      ###############
      handle.RESULTS_PROJECTIONS <- matrix.P %*% handle.SELECT
      rownames(handle.RESULTS_PROJECTIONS) <- colnames(handle.SELECT)
      list.RESULTS_PROJECTION[[p]] <- t(handle.RESULTS_PROJECTIONS)
      
      #################
      # DEPTH PROJECT #
      #################
      handle.RESULTS_PROJECTIONS_DEPTH <- matrix.P_DEPTH %*% handle.SELECT
      rownames(handle.RESULTS_PROJECTIONS_DEPTH) <- colnames(handle.SELECT)
      list.RESULTS_PROJECTION_DEPTH[[p]] <- t(handle.RESULTS_PROJECTIONS_DEPTH)
      
    }
    
    ###########
    # CONVERT #
    ###########
    matrix.RESULTS_PROJECTION <- do.call(cbind, list.RESULTS_PROJECTION)
    matrix.RESULTS_PROJECTION_DEPTH <- do.call(cbind, list.RESULTS_PROJECTION_DEPTH)
    matrix.RESULTS_PROJECTION <- cbind(matrix.RESULTS_PROJECTION, SimulationModelNumber = matrix.RESULTS_TRAIT[,ncol(matrix.RESULTS_TRAIT)])
    matrix.RESULTS_PROJECTION_DEPTH <- cbind(matrix.RESULTS_PROJECTION_DEPTH, SimulationModelNumber = matrix.RESULTS_TRAIT[,ncol(matrix.RESULTS_TRAIT)])
    colnames(matrix.RESULTS_PROJECTION) <- colnames(matrix.RESULTS_PROJECTION_DEPTH) <- colnames(matrix.RESULTS_TRAIT)
    
  }
  
  if (logical.PIC == F){matrix.RESULTS_PIC <- matrix.RESULTS_PIC_DEPTH <- NULL}
  return(list(RESULTS_TRAIT = data.frame(matrix.RESULTS_TRAIT), 
              RESULTS_PIC = data.frame(matrix.RESULTS_PIC), 
              RESULTS_PIC_DEPTH = data.frame(matrix.RESULTS_PIC_DEPTH), 
              RESULTS_PROJECT = data.frame(matrix.RESULTS_PROJECTION),
              RESULTS_PROJECT_DEPTH = data.frame(matrix.RESULTS_PROJECTION_DEPTH)))
}