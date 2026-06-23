# R code for EA imidacloprid concentration data processing
# Author: Jiayi Liu
# Date: June 2026
# Database: [Accessed 16 Jan 2025] Water quality monitoring data GC-MS and LC-MS semi-quantitative screen https://environment.data.gov.uk/dataset/e85a7a52-7a75-4856-a0b3-8c6e4e303858
# Dataset used: LCMS_Target_and_Non_Targeted_Screening.csv



## R version 4.5.2
library(dplyr)  # v1.2.0
library(stringr)  # v1.6.0



### 1. Extract Records taken between 2014 and 2023 from Freshwater Sites -------
# i.e., positive detections
positives <- LCMS_Target_and_Non_Targeted_Screening %>%
  filter(year >= 2014 & year <= 2023) %>%
  filter(str_starts(SPT_DESC, "FRESHWATER")) %>%
  select(Sample_Site_ID, Sample_datetime, year, Concentration, SMPT_LONG_NAME, 
         SPT_DESC, SMC_DESC,
         SMPT_EASTING, SMPT_NORTHING, Longitude, Latitude,
         ARE_DESC, OPCAT_NAME, PURP_DESC, Compound_Name)



### 2. Compute Non-Detects -----------------------------------------------------
nondetects <- function(data = positives){
  IMI <- data %>% filter(Compound_Name == "Imidacloprid")
  non <- data %>% filter(!(Compound_Name == "Imidacloprid"))
  
  y1 <- non[!duplicated(non[c("Sample_datetime", "SMPT_EASTING", "SMPT_NORTHING")]),]
  y1$Concentration <- NA
  
  new_rows <- y1[!paste(y1$Sample_datetime, y1$SMPT_EASTING, y1$SMPT_NORTHING) %in% paste(IMI$Sample_datetime, IMI$SMPT_EASTING, IMI$SMPT_NORTHING), ]
  temp <- rbind(IMI, new_rows)
  
  RawAllVisits <- temp %>%
    arrange(Sample_Site_ID, Sample_datetime) %>%
    select(-Compound_Name) %>%
    mutate(Detection = ifelse(is.na(Concentration), 0, 1)) %>%
    mutate(SMPT_NORTHING = as.numeric(SMPT_NORTHING))
  
  return(RawAllVisits)
}
allSamples <- nondetects(data = positives)



### 3. Data Cleaning -----------------------------------------------------------
# Remove two sites that fall outside boundary of England (45401034 and 88000045)
# -- checked in ArcGIS Pro
allSamples <- allSamples %>%
  filter(!(Sample_Site_ID %in% c("45401034", "88000045"))) 
# 4,152 resulting samples

#saveRDS(allSamples, "./Raw_EA_2014OW.rds")


allSamples2019OW <- allSamples %>% filter(year >= 2019)
# 1,909 resulting samples
# River basin district (RBD) data produced from ArcGIS Pro was joined to create dataset "allSamples1909RBD.rds"



### 4. Running & Standing Waters -----------------------------------------------
# Standing Waters
standing_waters_names <- unique(allSamples$Sample_Site_ID[allSamples$SMC_DESC == "POND / LAKE / RESERVOIR WATER"])
standing_waters_samples_since2014 <- allSamples %>%
  filter(Sample_Site_ID %in% standing_waters_names)

# Running Waters
running_waters_samples_since2014 <- allSamples %>%
  filter(! Sample_Site_ID %in% standing_waters_names)



# END
