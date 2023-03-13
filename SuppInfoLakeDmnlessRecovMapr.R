# Lake TP dimensionless target-based recovery time calculator/mapper.
#
library(dplyr)
library(sf)
library(rgdal)
library(mapview)
library(rstudioapi)

# current_path <- getActiveDocumentContext()$path
# setwd(dirname(current_path))

## Dimensionless recovery target
target <- 0.75  #  This is the target (1-beta)/(1-alpha) value.

## Scenario: whether step or exponential influent drop
scenario <- 2 # Enter "1" for influent step drop, or "2" for exponential drop
eta <- 2.74E-4        # Decay coefficient for the influent [TP] (1/day)

## Select states to map
conus <- c("AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI",
  "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",
  "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")

lagos <- c("CT", "IA", "IL", "IN", "MA", "ME", "MI", "MN", "MO", "NH", "NJ", "NY", "OH", "PA", "RI", "VT", "WI")

# states2map <- conus
# states2map <- lagos
# states2map <- c('ME','NH','VT')
states2map <- "WI"

DF1 <- read.csv("SuppInfoLakesCONUS_NHDflow.csv")
DF1 <- filter(DF1, tau_d >= 7)
DF1 <- filter(DF1, tau_yr < 1000) # Removing apparent outliers
DF1 <- filter(DF1, ST %in% states2map)

theta <- 0.5606078
kr_over_kb <- 8.603011
phi <- -0.5039417

DF1$ks <- theta * DF1$tau_d^phi # ks is lake-specific
kr <- 2.31E-4
kb <- kr / kr_over_kb

DF1$A <- -(1 / DF1$tau_d + DF1$ks)
DF1$BE <- kr * DF1$ks
DF1$F <- -(kb + kr)

# eigenvalues
DF1$L1 <- (DF1$A + DF1$F + sqrt((DF1$A + DF1$F)^2 - 4 * (DF1$A * DF1$F - DF1$BE))) / 2
DF1$L2 <- (DF1$A + DF1$F - sqrt((DF1$A + DF1$F)^2 - 4 * (DF1$A * DF1$F - DF1$BE))) / 2

# coeffs for step influent drop
DF1$M1 <- (DF1$L2 - DF1$A + (DF1$BE / DF1$F)) / (DF1$L2 - DF1$L1)
DF1$M2 <- (DF1$A - DF1$L1 - (DF1$BE / DF1$F)) / (DF1$L2 - DF1$L1)

# coeffs for exponential influent drop
DF1$N1 <- ((DF1$F + eta) * (DF1$L2 - DF1$L1) * DF1$L1 / ((DF1$L1 + eta) * (DF1$L2 + eta) * (DF1$L2 - DF1$A))) * DF1$M1
DF1$N2 <- (eta / (DF1$L1 + eta)) * DF1$M1
DF1$N3 <- (eta / (DF1$L2 + eta)) * DF1$M2

DF1$Trecov_d <- 0

# Function(s) for calculating the time to lake recovery.
if (scenario == 1) { # step drop
  print("step drop")
  Pfun <- function(t, M_1, M_2, L_1, L_2) {
    M_1 * exp(L_1 * t) + M_2 * exp(L_2 * t) + target - 1
  }

  for (i in 1:length(DF1$COMID)) {
    u1 <- uniroot(function(x) Pfun(x, DF1$M1[i], DF1$M2[i], DF1$L1[i], DF1$L2[i]), lower = 0, upper = 10^8)
    DF1$Trecov_d[i] <- u1$root
  }

} else { # exponential drop
  print("exponential drop")
  Pfun <- function(t, eta, N_1, N_2, N_3, L_1, L_2) {
    N_1 * exp(-eta * t) + N_2 * exp(L_1 * t) + N_3 * exp(L_2 * t) + target - 1
  }

  for (i in 1:length(DF1$COMID)) {
    u1 <- uniroot(function(x) Pfun(x, eta, DF1$N1[i], DF1$N2[i], DF1$N3[i], DF1$L1[i], DF1$L2[i]), lower = 0, upper = 10^8)
    DF1$Trecov_d[i] <- u1$root
  }

}

DF1$Trecov_yr <- DF1$Trecov_d / 365.25
DF1$Trecov_yr <- round(DF1$Trecov_d / 365.25, 1)

## mapview map of recovery times
DF1$Target <- toString(target)
if (scenario == 1) {
  DF1$Scenario <- "Step influent drop"
} else {
  DF1$Scenario <- "Exponential influent drop"
}

Lakes <- select(DF1, Scenario, Target, COMID, LakeName, ST, HUC8, Latitude, Longitude, Q0001E, LakeVolume_m3, tau_d, tau_yr, Trecov_yr)

coordinates(Lakes) <- c("Longitude", "Latitude")
proj4string(Lakes) <- "+proj=longlat +datum=NAD83"

Qmin <- min(Lakes$Trecov_yr)
Q1 <- as.numeric(quantile(Lakes$Trecov_yr, probs = 0.25))
Q2 <- as.numeric(quantile(Lakes$Trecov_yr, probs = 0.5))
Q3 <- as.numeric(quantile(Lakes$Trecov_yr, probs = 0.75))
Qmax <- max(Lakes$Trecov_yr)

names(Lakes)[names(Lakes) == "Trecov_yr"] <- "Recovery Time yrs"

m2 <- mapview(Lakes, xcol = "Longitude", ycol = "Latitude", zcol = "Recovery Time yrs", at = c(Qmin, Q1, Q2, Q3, Qmax),
  legend = "TRUE", cex = 3, crs = 4326, grid = FALSE)
print(m2)
