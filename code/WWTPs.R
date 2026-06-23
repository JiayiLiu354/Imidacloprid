# R code for statistical analyses --
# -- influence of the presence of WWTPs on imidacloprid detection 
# Author: Jiayi Liu
# Date: June 2026
# Database: [Accessed 19 May 2025] Wastewater Treatment in England https://www.data.gov.uk/dataset/d7e2c57b-110a-462b-97a0-9833e7d26cc2/wastewater-treatment-in-england
# Dataset used: UWWTR_Art15_24thOct2024.ods



## R version 4.5.2 
library(readODS)  # v2.3.2
library(dplyr)  # v1.2.0
library(binom)  # v1.1-1.1



### 1. Data Processing ---------------------------------------------------------
T_DischargePoints <- read_ods("UWWTR_Art15_24thOct2024.ods", sheet = 5)

df_dcp <- T_DischargePoints %>%
  filter(dcpState == 1,
         dcpSurfaceWaters == 1,
         dcpWaterBodyType == "FW") %>%
  select(uwwCode, dcpCode, dcpName, dcpBeginLife, dcpLatitude, dcpLongitude)

T_UWWTPS <- read_ods("UWWTR_Art15_24thOct2024.ods", sheet = 2)

df_dcp <- left_join(df_dcp, T_UWWTPS, by = "uwwCode") %>%
  select(uwwCode, uwwName, uwwLoadEnteringUWWTP, uwwCapacity, uwwLatitude, uwwLongitude,
         aggCode,
         dcpCode, dcpName, dcpBeginLife, dcpLatitude, dcpLongitude)


T_UWWTPAgglos <- read_ods("UWWTR_Art15_24thOct2024.ods", sheet = 3)

df_dcp <- left_join(df_dcp, T_UWWTPAgglos, by = "aggCode") 

df_UWWTPs <- df_dcp %>%
  select(uwwCode, uwwName, dcpBeginLife,
         uwwCapacity, uwwLoadEnteringUWWTP, aggGenerated,
         uwwLatitude, uwwLongitude,
         aggCode, aggLatitude, aggLongitude,
         dcpCode, dcpName, dcpLatitude, dcpLongitude)

duplicates <- df_UWWTPs[duplicated(df_UWWTPs$uwwCode), ]

UWWTPs <- df_UWWTPs %>%
  mutate(LAT = ifelse(uwwCode %in% duplicates$uwwCode, uwwLatitude, dcpLatitude),
         LON = ifelse(uwwCode %in% duplicates$uwwCode, uwwLongitude, dcpLongitude)) %>%
  distinct(uwwCode, .keep_all = TRUE) %>%
  select(uwwCode, uwwName, LAT, LON, uwwLoadEnteringUWWTP, aggGenerated, uwwCapacity)

# 1,217 resulting records



### 2. Logistic Regression -----------------------------------------------------
# Calculated numbers in each category produced from ArcGIS Pro
binom.confint(x = 49, n = 67, methods = "wilson")
binom.confint(x = 8, n = 18, methods = "wilson")

detection <- c(rep(1, 49), rep(0, 18), rep(1, 8), rep(0, 10))
wwtp <- c(rep(1, 67), rep(0, 18))
df <- data.frame(detection = detection, wwtp = factor(wwtp, labels = c("No", "Yes")))

model <- glm(detection ~ wwtp, data = df, family = binomial)
summary(model)

or_ci <- exp(cbind(OR = coef(model), confint(model)))
print(or_ci)



# END
