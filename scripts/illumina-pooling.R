# Script for determining volume of each cleaned, post-PCR Illumina library
#   to be added pooled library for multiplexed sequencing

#### SETUP ####

# Load required packages
library(tidyverse)

#### VARIABLE PARAMETERS ####

# Name of column with Qubit concentrations
qubit <- "max_library_qubit"

# Load test data
allData <- read_csv("~/github-repos/projects/glue-paper1/sequencing-prep/data/clean/deep3/deep3_lane1_libraryConcentrations.csv")
allData_noQubitRem <- read_csv("~/github-repos/projects/glue-paper1/sequencing-prep/data/clean/deep3/deep3_lane1_libraryConcentrations.csv") %>% 
  
  # Remove samples if they have no concentration
  filter(!(!!sym(qubit) == 0))

# Total volume of post-PCR libraries
available_vol <- 21

# Mean size of library fragments. Approximated by gel electrophoresis or
# Bioanalyser
mean_fragment_size <- 350

# Final volume of diluted pools. You should make sure to have at least
# 4 uL so that 2 uL can be taken twice in case additional sequencing
# is required downstream. It is not advisable to pipette less than 2uL
# due to greater pipetting error at small volumes. This value may need to be 
# increased if the initial volume to be taken from the undiluted libraries 
# is < 2 uL, again since pipetting error increases at small volumes
sample_pooling_vol <- 5

# Minimum molarity required by sequencing center
min_required_molarity <- 3

# Number of lanes across which to split samples
num_lanes <- 1

#### FIXED PARAMETERS ####

# Total number of samples to be multiplexed
total_samples <- nrow(allData)

# Final volume of individual diluted pools
# Ensure enough of each pool is left for at least one more
# round of sequencing, if necessary. 
sample_final_volume <- 2 * sample_pooling_vol

# Final volume of pools for serial dilutions (if required)
serial_dilution_final_vol <- 20

#### FUNCTIONS ####

#' Converts concentration from ng/uL to nM
#' 
#' @param ng_uL_concentration Concentration of sample in ng/uL
#' @param mean_fragment_size Mean fragment size of library in bp
#' 
#' @return Molarity of sample in nanomolar (nM)
calculate_molarity <- function(ng_uL_concentration, mean_fragment_size){
  
  nM_concentration <- (ng_uL_concentration / (660 * mean_fragment_size)) * 10 ^ 6
  nM_concentration <- round(nM_concentration, 1)
  
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

#' Creates dataframe with sample and TE volumes to create dilutions
#' 
#' @param allData Dataframe with all samples and their concentrations
#' @param qubit Column containing qubit concentrations
#' @param mean_fragment_size Mean fragment size of library in bp
#'
#' @return Dataframe with sample initial volume and TE volume appended as columns
output_data <- function(allData, qubit, mean_fragment_size){
  
  city <- allData %>% pull(city) %>% unique()
  print(city)
  # Calculate molarity of all samples in df
  data_out <- allData %>% 
    
    # Calculate concentration in nM
    # Need to conver 'qubit' character string to symbol and then evaluate to column
    mutate(nM_concentration = calculate_molarity((!!sym(qubit)), mean_fragment_size))
  
  data_out_lowMol <- data_out %>% filter(nM_concentration < min_required_molarity)
  data_out <- data_out %>% filter(nM_concentration >= min_required_molarity)
  print(sprintf("%s samples are below the minimum required molarity.", nrow(data_out_lowMol)))
  print(sprintf("%s samples meet the minimum molarity threshold", nrow(data_out)))
  print("Low molarity samples will be removed for now and added back at the end")
  
  # Figure out maximum molarity of final library pool. Based on minimum across
  # all samples so that each lane will have equal final molarity.
  # Molarity rounded down to nearest 0.25
  final_pool_molarity <- floor(min(data_out$nM_concentration) / 0.25) * 0.25
  final_pool_ng_uL <- data_out %>% pull(!!sym(qubit)) %>% min()
  
  print(sprintf("The final pool molarity will be %s nM", final_pool_molarity))
  print(sprintf("The final pool ng/uL will be %s ng/uL", final_pool_ng_uL))
  
  
  # Calculate minimum library concentration such that dilution pipettes >= 2 uL
  library_molarity_for_accuracy <- (sample_final_volume * final_pool_molarity) / 2
  print(sprintf("Libraries must be less than %s nM to pipette >= 2uL", 
                library_molarity_for_accuracy))
  
  
  # Add serial dilutions for samples that are too concentrated (i.e., require < 2uL library)
  data_out <- data_out %>% 
    mutate(serial_dilution_required = ifelse(nM_concentration > library_molarity_for_accuracy, 
                                             "Yes", "No"),
           dilution_factor = ifelse(serial_dilution_required == 'Yes', 
                                    ceiling(nM_concentration / library_molarity_for_accuracy),
                                    NA),
           conc_post_dilution =  round(nM_concentration / dilution_factor, 1),
           library_vol_forSerial = round((conc_post_dilution * 20) / nM_concentration, 1),
           TE_vol_forSerial = round(20 - library_vol_forSerial, 1))
  
  num_serials <- data_out %>% filter(serial_dilution_required == 'Yes') %>% nrow()
  print(sprintf("%s samples need to be serially diluted", num_serials))
  
  # Calculate initial library volume to create diluted libraries
  data_out <- data_out %>%
    mutate(
      library_volume_uL =
        ifelse(
          serial_dilution_required == 'Yes',
          calculate_library_initial_volume(conc_post_dilution,
                                           final_pool_molarity,
                                           sample_final_volume),
          calculate_library_initial_volume(nM_concentration,
                                           final_pool_molarity,
                                           sample_final_volume)
        ),
      
      library_volume_uL = round(library_volume_uL, 1),
      # Calculate volume of TE to be added to dilute libraries to
      # final volume of "sample_pooling_vol"
      TE_vol_uL = round(sample_final_volume - library_volume_uL, 1)
    )
  
  minVol <- min(data_out$library_volume_uL)
  maxVol <- max(data_out$library_volume_uL)
  if(minVol < 2){
    warning("Some samples have < 2 uL initial volume. Consider increasing
            'sample_pooling_vol' to reduce pipetting error ")
  }else if(maxVol > available_vol){
    warning(sprintf("Some samples require pipetting more than is available in the 
            post-PCR libraries (i.e., %s uL). Decrease 'sample_pooling_vol'", available_vol))
  }else{
    print(sprintf("The volumes used from the post-PCR libraries range from %s to %s uL", minVol, maxVol))
  }
  
  if(nrow(data_out_lowMol) > 0){
    data_out <- plyr::rbind.fill(data_out, data_out_lowMol)
  }
  
  return(data_out)

}

#### DETERMINE POOLING VOLUMES ####

# If pooling by group (e.g., city), run this code below.
data_out <- allData_noQubitRem %>% group_split(city) %>% map_dfr(., output_data, qubit, mean_fragment_size)
allData %>% filter(!!sym(qubit) == 0)
# If pooling all samples together, run this code below 
# data_out <- output_data(allData, qubit, mean_fragment_size)

# Split samples equally across lanes. Only necessary if sequencing across multiple lanes. 
# sequencing_lanes <- split(data_out, rep(1:num_lanes, each = total_samples))
# names(sequencing_lanes) <- paste0("Lane_", seq_along(sequencing_lanes))
# list2env(sequencing_lanes, envir = .GlobalEnv)

# Write Lanes to CSV
write_csv(data_out, file = "~/github-repos/projects/glue-paper1/sequencing-prep/data/clean/deep3/deep3_lane1_dilutions.csv")

test <- data_out %>% 
  group_by(city) %>% 
  summarise(n = n(),
            g = n * (800/686))
nrow(data_out)
sum(test$g)
