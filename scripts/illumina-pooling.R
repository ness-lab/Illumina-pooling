# Script for determining volume of each cleaned, post-PCR Illumina library
#   to be added pooled library for multiplexed sequencing

#### SETUP ####

# Load required packages
library(tidyverse)

# Load test data
testData <- read_csv("test-data/tesPrep_postPCR_cleaned_qubit.csv")

#### VARIABLE PARAMETERS ####

# Mean size of library fragments. Approximated by gel electrophoresis or
# Bioanalyser
mean_fragment_size <- 400

# Volume to be taken from individual, diluted pools to make up 
# the "final pool" that is sent for sequencing. This volume will
# be identical for all multiplexed samples. It may need to be increased
# if the initial volume to be taken from the undiluted libraries is < 2 uL,
# since pipetting error increases at small volumes
sample_pooling_vol <- 5

#### FIXED PARAMETERS ####

# Total number of samples to be multiplexed
total_samples <- nrow(testData)

# Maximum possible molarity of pool, rounded down to nearest 0.25
# Based on the lowest concentration among samples to be sequenced
# This approach will not allow the use of "undiluted", cleaned, post-PCR
# libraries. 
final_pool_molarity <- floor(min(testData_mod$nM_concentration) / 0.25) * 0.25

# Final volume of individual diluted pools
# Ensure enough of each pool is left for another round of 
# sequencing, if necessary. 
sample_final_volume <- 2 * sample_pooling_vol

# Volume of final pool
# This is the amount of library that will be sent for sequencing
final_pool_volume <- total_samples * sample_pooling_vol

#### FUNCTIONS ####

# Move to separate scripts?

#' Converts concentration from ng/uL to nM
#' 
#' @param ng_uL_concentration Concentration of sample in ng/uL
#' @param mean_fragment_size Mean fragment size of library
#' 
#' @return Molarity of sample in nanomolar (nM)
calculate_molarity <- function(ng_uL_concentration, mean_fragment_size){
  
  nM_concentration <- (ng_uL_concentration / (660 * mean_fragment_size)) * 10 ^ 6
  
  return(nM_concentration)
  
}

#' Determines initial volume taken from each undiluted sample
#' 
#' @param nM_concentration Molarity of sample in nanomolar (nM)
#' @param final_pool_molarity Molarity (in nM) of final pooled library to be sequenced
#' @param sample_final_volume Final lolume (in uL) of each sample to be added to the final
#'   library pool that is sent for sequencing
#'
#' @return The initial volume (in uL) to be removed from the undiluted library to create an
#'   individual sample diluted to the maximum "final_pool_molarity"
calculate_library_initial_volume <- function(nM_concentration, final_pool_molarity, sample_final_volume){
  
  initial_volume <- (final_pool_molarity * sample_final_volume) / nM_concentration
  
  return(initial_volume)
  
}

#### DETERMINE POOLING VOLUMES ####

data_out <- testData %>%
  
  # Calculate molarity of all libraries
  mutate(nM_concentration = calculate_molarity(Qubit, mean_fragment_size), 
         
         # Calculate initial library volume to create diluted libraries
         library_volume_uL = round(
           calculate_library_initial_volume(
             nM_concentration,
             final_pool_molarity, 
             final_pool_volume), 1),
    
         # Calculate volume of TE to be added to dilute libraries to 
         # final volume of "sample_pooling_vol"
         TE_vol_uL = round(sample_final_volume - library_volume_uL, 1))

