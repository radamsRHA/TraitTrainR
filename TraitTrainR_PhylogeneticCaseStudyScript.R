#######################
# Load dependencies   #
#######################
library(TraitTrainR); library(ape); library(geiger); library(phytools); library(caret); library(corrplot)

##################
# DEFINE OPTIONS #
##################
numeric.ERROR <- 0; numeric.SE <- 0

#########################
# get example phylogeny #
#########################
handle.TargetTree <- read.tree(text = "((((Ladona_fulva:255.023267,Ephemera_danica:255.023267):68.976733,((Blatella_germanica:101.148321,Zootermopsis_nevadensis:101.148321):216.876761,(((((Copidosoma_floridanum:86.271212,Nasonia_vitripennis:86.271212):23.010504,Trichogramma_pretiosum:109.281716):126.902472,(((Orussus_abietinus:106.953420,Cephus_cinctus:106.953420):23.018917,(((Dufourea_novaeangliae:36.095130,Lasioglossum_albipes:36.095130):16.113401,((((((Bombus_terrestris:5.186677,Bombus_impatiens:5.186677):15.492116,Melipona_quadrifasciata:20.678794):5.021206,(Apis_florea:6.953128,Apis_mellifera:6.953128):18.746872):2.917824,Eufriesea_mexicana:28.617824):6.394660,Habropoda_laboriosa:35.012484):8.361841,Megachile_rotundata:43.374325):8.834206):39.641469,((((((Cardiocondyla_obscurior:21.650487,Solenopsis_invicta:21.650487):2.371290,(Acromyrmex_echinatior:7.059427,Atta_cephalotes:7.059427):16.962350):3.418271,Pogonomyrmex_barbatus:27.440048):7.575110,Camponotus_floridanus:35.015157):3.132369,Linepithema_humile:38.147527):12.620282,Harpegnathos_saltator:50.767809):41.082191):38.122337):98.027663,Athalia_rosae:228.000000):8.184188):72.965812,((((((Leptinotarsa_decemlineata:78.221283,Anoplophora_glabripennis:78.221283):20.022053,Dendroctonus_ponderosae:98.243336):15.244386,Tribolium_castaneum:113.487723):29.468798,Onthophagus_taurus:142.956520):15.584891,Agrilus_planipennis:158.541411):111.976371,(((((Heliconius_melpomene:56.157043,Danaus_plexippus:56.157043):27.358024,(Manduca_sexta:60.792195,Bombyx_mori:60.792195):22.722871):40.684933,Plutella_xylostella:124.200000):73.466744,Limnephilus_lunatus:197.666744):44.372315,(((Culex_quinquefasciatus:35.755285,Aedes_egypti:35.755285):24.173533,((Anopheles_gambiae:16.970563,Anopheles_funestus:16.970563):18.704446,Anopheles_albimanus:35.675010):24.253809):39.971182,(((((Drosophila_melanogaster:25.540985,Drosophila_pseudoobscura:25.540985):10.084320,Drosophila_grimshawi:35.625305):31.541797,(((Lucilia_cuprina:36.892316,Musca_domestica:36.892316):11.339494,Glossina_morsitans:48.231810):12.870694,Ceratitis_capitata:61.102504):6.064599):25.853183,Mayetiola_destructor:93.020286):4.469973,Lutzomyia_longipalpis:97.490259):2.409741):142.139059):28.478724):38.632218):5.932829,(((((((Halyomorpha_halys:50.856065,Oncopeltus_fasciatus:50.856065):25.007392,Cimex_lectularius:75.863457):11.837700,Gerris_buenoi:87.701157):10.147324,Homalodisca_vitripennis:97.848481):2.051519,(Acyrthosiphon_pisu:93.342828,Pachypsylla_venusta:93.342828):6.557172):122.100000,Frankliniella_occidentalis:222.000000):36.545614,Pediculus_humanus:258.545614):56.537216):2.942253):5.974917):36.590624,Catajapyx_aquilonaris:360.590624):12.545713,(((Daphnia_pulex:319.660629,Eurytemora_affinis:319.660629):26.684919,Hyalella_azteca:346.345548):16.516889,((((Metaseiulus_occidentalis:234.985824,Ixodes_scapularis:234.985824):44.729414,((((Parasteatoda_tepidariorum:68.827729,Latrodectus_hesperus:68.827729):45.631231,Stegodyphus_mimosarum:114.458960):67.290395,Loxosceles_reclusa:181.749355):67.650284,Centruroides_sculpturatus:249.399640):30.315599):11.530677,Tetranychus_urticae:291.245915):24.683949,Strigamia_maritima:315.929864):46.932572):10.273900);")
string.PATH_OUTPUT_FILE <- '~/Desktop/REVISIONS_9Aug24//Arthropod/RESULTS_Arthropod_XXX_YYY.txt'
string.PATH_OUTPUT_FILE <- gsub(pattern = "XXX", replacement = numeric.ERROR, x = string.PATH_OUTPUT_FILE)
string.PATH_OUTPUT_FILE <- gsub(pattern = "YYY", replacement = toString(numeric.SE), x = string.PATH_OUTPUT_FILE)


#############################
# Simulation Model Settings #
#############################
list.SimulationModelSettings <- list() # define an empty model list

##########################
# SET SIMULATION MODEL 1 #
##########################
numeric.NumberTrainingReps <- 10^4 # same number of replicates for all models in list.SimulationModelSettings
list.Rmatrix <- list(); for (i in 1:numeric.NumberTrainingReps){list.Rmatrix[[i]] <- matrix(1, nrow = 1, ncol = 1)} # three traits. Different rates for different traits can be specified here. 

######################
# First model is BM  #
######################
list.SimulationModelSettings[[1]] <- list(string.SimulationModel = "BM", 
                                          vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
                                          vector.AncestralState = rnorm(n = numeric.NumberTrainingReps),
                                          list.Rmatrix = list.Rmatrix)

#######################
# Second model is OU  #
#######################
list.SimulationModelSettings[[2]] <- list(string.SimulationModel = "OU", 
                                          vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
                                          vector.AncestralState = rnorm(n = numeric.NumberTrainingReps),
                                          vector.Alpha = runif(n = numeric.NumberTrainingReps, min = exp(-500), max = exp(1)),
                                          list.Rmatrix = list.Rmatrix)


######################
# third model is EB  #
######################
list.SimulationModelSettings[[3]] <- list(string.SimulationModel = "EB", 
                                          vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
                                          vector.AncestralState = rnorm(n = numeric.NumberTrainingReps), 
                                          vector.A = runif(n = numeric.NumberTrainingReps, min = log(10^-5)/373.1363, max = -0.000001),
                                          list.Rmatrix = list.Rmatrix)

##########################
# fourth model is kappa  #
##########################
list.SimulationModelSettings[[4]] <- list(string.SimulationModel = "kappa", 
                                          vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
                                          vector.AncestralState = rnorm(n = numeric.NumberTrainingReps),
                                          vector.kappa = runif(n = numeric.NumberTrainingReps, min = exp(-500), max = 1), 
                                          list.Rmatrix = list.Rmatrix)
##########################
# fifth model is lambda  #
##########################
list.SimulationModelSettings[[5]] <- list(string.SimulationModel = "lambda", 
                                          vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
                                          vector.AncestralState = rnorm(n = numeric.NumberTrainingReps), 
                                          vector.lambda = runif(n = numeric.NumberTrainingReps, min = exp(-500), max = 1), 
                                          list.Rmatrix = list.Rmatrix)

#######################
# six model is delta  #
#######################
list.SimulationModelSettings[[6]] <- list(string.SimulationModel = "delta", 
                                          vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
                                          vector.AncestralState = rnorm(n = numeric.NumberTrainingReps), 
                                          vector.delta = runif(n = numeric.NumberTrainingReps, min = exp(-500), max = 3), 
                                          list.Rmatrix = list.Rmatrix)


###########################
# seventh model is trend  #
###########################
list.SimulationModelSettings[[7]] <- list(string.SimulationModel = "trend", 
                                          vector.Sig2 = rexp(n = numeric.NumberTrainingReps, rate = 1), 
                                          vector.AncestralState = rnorm(n = numeric.NumberTrainingReps),
                                          vector.slope = runif(n = numeric.NumberTrainingReps, min = -100, max = 100),
                                          list.Rmatrix = list.Rmatrix)

####################
# SIMULATE TRAITS! #
####################
handle.RESULTS_TEST <- TraitTrain(handle.Phylogeny = handle.TargetTree, 
                                  list.SimulationModelSettings = list.SimulationModelSettings, 
                                  logical.PIC = F, logical.PROJECT = F, numeric.MeasurementError = numeric.ERROR)

##############################
# Compare with fitContinuous #
##############################
TEST_TRAIT <- handle.RESULTS_TEST$RESULTS_TRAIT
RESULTS <- matrix(nrow = nrow(TEST_TRAIT), ncol = 9)
colnames(RESULTS) <- c("AIC_BM", "AIC_OU", "AIC_EB", "AIC_Kappa", "AIC_Lambda", "AIC_Delta", "AIC_Trend", "TRUE_MODEL", "BEST_MODEL")

rep_count <- 0
for (testrep in 1:nrow(TEST_TRAIT)){
  
  print(rep_count)
  rep_count <- rep_count + 1
  vector.Trait <- as.numeric(TEST_TRAIT[testrep,1:length(handle.TargetTree$tip.label)])
  names(vector.Trait) <- handle.TargetTree$tip.label
  
  model.BM <- fitContinuous(phy = handle.TargetTree, dat = vector.Trait, SE = numeric.SE, model = "BM")
  model.OU <- fitContinuous(phy = handle.TargetTree, dat = vector.Trait, SE = numeric.SE, model = "OU")
  model.EB <- fitContinuous(phy = handle.TargetTree, dat = vector.Trait, SE = numeric.SE, model = "EB")
  model.kappa <- fitContinuous(phy = handle.TargetTree, dat = vector.Trait, SE = numeric.SE, model = "kappa")
  model.lambda <- fitContinuous(phy = handle.TargetTree, dat = vector.Trait, SE = numeric.SE, model = "lambda")
  model.delta <- fitContinuous(phy = handle.TargetTree, dat = vector.Trait, SE = numeric.SE, model = "delta")
  model.trend <- fitContinuous(phy = handle.TargetTree, dat = vector.Trait, SE = numeric.SE, model = "rate_trend")
  vector.AIC <- c(AIC(model.BM), AIC(model.OU), AIC(model.EB), AIC(model.kappa), AIC(model.lambda), AIC(model.delta), AIC(model.trend))
  names(vector.AIC) <- 1:7
  
  print(rep_count)
  RESULTS[rep_count,] <- c(vector.AIC, TEST_TRAIT$SimulationModelNumber[rep_count], paste0(x = names(vector.AIC)[vector.AIC==min(vector.AIC)], collapse = ""))
  
  write.table(x = RESULTS, file = string.PATH_OUTPUT_FILE, append = F, quote = F, sep = '\t', row.names = F)
  
}

