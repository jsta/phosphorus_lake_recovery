# Lake TP concentration-based recovery time calculator/mapper.
#
library(dplyr)
library(sf)
library(rgdal)
library(mapview)
library(rstudioapi)

# current_path <- getActiveDocumentContext()$path
# setwd(dirname(current_path))

##
scenario <- 2 # Enter "1" for influent step drop, or "2" for exponential drop
eta <- 2.74E-4        # Decay coefficient for the influent [TP] (1/day).

## Drop magnitude and TP conc. target
alpha <- 0.2       # 1-alpha is the fractional decrease imposed on the influent concentration starting at t=0.
Cf <- 10           # This is the target for C1: the "acceptable" water column [TP], e.g., a WQ criterion.

## Select states to map
conus <- c("AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI",
  "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",
  "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")

lagos <- c("CT", "IA", "IL", "IN", "MA", "ME", "MI", "MN", "MO", "NH", "NJ", "NY", "OH", "PA", "RI", "VT", "WI")

# states2map <- conus
# states2map <- lagos
# states2map <- c('ME','NH','VT')
states2map <- "VT"

DF1 <- read.csv("SuppInfoLakesCONUS_SPARROWflow.csv")

DF1 <- filter(DF1, tau_d >= 7)
DF1 <- filter(DF1, tau_yr < 1000) # Removing three apparent outliers

DF1 <- filter(DF1, ST %in% states2map)

DF1 <- filter(DF1, tau_d >= 7) # Added 08/24/21

names(DF1)[names(DF1) == "LakeTP_ugL"] <- "Lake_TP_ppb"
DF1 <- filter(DF1, Lake_TP_ppb > 0)

DF1$Target_TP_ppb <- toString(Cf)
if (scenario == 1) {
  DF1$Scenario <- "Step influent drop"
} else {
  DF1$Scenario <- "Exponential influent drop"
}

#
theta <- 0.5606078
kr_over_kb <- 8.603011
phi <- -0.5039417

DF1$ks <- theta * DF1$tau_d^phi # ks is lake-specific
kr <- 2.31E-4      # This is the approximate annual mean kr from Chapra & Canale's calibration of Shagawa model.
kb <- kr / kr_over_kb

# Effective steady-state CSTR k
DF1$keff <- DF1$ks * kb / (kb + kr)

# Initial conditions
DF1$C10 <- DF1$Lake_TP_ppb
DF1$Ci0 <- (1 + DF1$keff * DF1$tau_d) * DF1$C10
DF1$Ci <- alpha * DF1$Ci0
DF1$beta <- Cf / DF1$C10

for (i in 1:length(DF1$COMID)) {
  if (DF1$beta[i] > 1) {
    DF1$beta[i] <- 1          # Any lake with beta >1 already has C1 below Cf, and its "recovery" time is therefore zero.
  }
}

# Logical constraint: must have beta > alpha (because C1 can't decrease more than Ci decreases)
DF1a <- filter(DF1, beta == 1)      # Lakes that already meet the goal
DF2 <- filter(DF1, beta < 1)      # Lakes that don't already meet the goal
DF2a <- filter(DF2, beta <= alpha) # Lakes that can't meet the goal, given a (1-alpha) decrease in influent TP
DF2 <- filter(DF2, beta > alpha)
###
DF2$A <- -(1 / DF2$tau_d + DF2$ks)
DF2$BE <- kr * DF2$ks
DF2$F <- -(kb + kr)

# eigenvalues
DF2$L1 <- (DF2$A + DF2$F + sqrt((DF2$A + DF2$F)^2 - 4 * (DF2$A * DF2$F - DF2$BE))) / 2
DF2$L2 <- (DF2$A + DF2$F - sqrt((DF2$A + DF2$F)^2 - 4 * (DF2$A * DF2$F - DF2$BE))) / 2

# coeffs for step influent drop
DF2$M1 <- (DF2$L2 - DF2$A + (DF2$BE / DF2$F)) * (1 - alpha) / ((DF2$L2 - DF2$L1) * (DF2$beta - alpha))
DF2$M2 <- (DF2$A - DF2$L1 - (DF2$BE / DF2$F)) * (1 - alpha) / ((DF2$L2 - DF2$L1) * (DF2$beta - alpha))

# coeffs for exponential influent drop
DF2$N1 <- ((DF2$F + eta) / (DF2$L2 - DF2$A)) * (DF2$L2 - DF2$A + (DF2$BE / DF2$F)) * DF2$L1 / ((DF2$L1 + eta) * (DF2$L2 + eta)) * (1 - alpha) / (DF2$beta - alpha)
DF2$N2 <- ((eta) / (DF2$L1 + eta)) * (DF2$L2 - DF2$A + (DF2$BE / DF2$F)) * (1 - alpha) / ((DF2$L2 - DF2$L1) * (DF2$beta - alpha))
DF2$N3 <- ((eta) / (DF2$L2 + eta)) * (DF2$A - DF2$L1 - (DF2$BE / DF2$F)) * (1 - alpha) / ((DF2$L2 - DF2$L1) * (DF2$beta - alpha))

DF2$Trecov_d <- 0

# Function(s) for calculating the time to lake recovery.
if (scenario == 1) { # step drop
  print("step drop")
  Pfun <- function(t, M_1, M_2, L_1, L_2) {
    M_1 * exp(L_1 * t) + M_2 * exp(L_2 * t) - 1
  }
  for (i in 1:length(DF2$COMID)) {
    u1 <- uniroot(function(x) Pfun(x, DF2$M1[i], DF2$M2[i], DF2$L1[i], DF2$L2[i]), lower = 0, upper = 10^8)
    DF2$Trecov_d[i] <- u1$root
  }

} else { # exponential drop
  print("exponential drop")
  Pfun <- function(t, eta, N_1, N_2, N_3, L_1, L_2) {
    N_1 * exp(-eta * t) + N_2 * exp(L_1 * t) + N_3 * exp(L_2 * t) - 1
  }
  for (i in 1:length(DF2$COMID)) {
    u1 <- uniroot(function(x) Pfun(x, eta, DF2$N1[i], DF2$N2[i], DF2$N3[i], DF2$L1[i], DF2$L2[i]), lower = 0, upper = 10^8)
    DF2$Trecov_d[i] <- u1$root
  }
}

DF2$Trecov_yr <- DF2$Trecov_d / 365.25
DF1a$Trecov_yr <- 0    # These are lakes that already meet or exceed (in the fall-beneath sense) the goal.
DF2a$Trecov_yr <- "NA"  # 10000 # These are lakes that can't meet the goal, given only (1-alpha) decrease in the influent.

DF1a <- select(DF1a, Scenario, Target_TP_ppb, COMID, LakeName, ST, Latitude, Longitude, LakeVolume_m3, FLOWcfs, tau_d, tau_yr, Lake_TP_ppb, Trecov_yr)
DF2 <- select(DF2, Scenario, Target_TP_ppb, COMID, LakeName, ST, Latitude, Longitude, LakeVolume_m3, FLOWcfs, tau_d, tau_yr, Lake_TP_ppb, Trecov_yr)
DF2a <- select(DF2a, Scenario, Target_TP_ppb, COMID, LakeName, ST, Latitude, Longitude, LakeVolume_m3, FLOWcfs, tau_d, tau_yr, Lake_TP_ppb, Trecov_yr)
DFx <- filter(DF2, Trecov_yr > 0)

Qmin <- min(DF1a$Trecov_yr)
Q1 <- min(DFx$Trecov_yr)
Q2 <- as.numeric(quantile(DFx$Trecov_yr, probs = 0.25))
Q3 <- as.numeric(quantile(DFx$Trecov_yr, probs = 0.5))
Q4 <- as.numeric(quantile(DFx$Trecov_yr, probs = 0.75))
Qmax <- max(DFx$Trecov_yr)

#
DF2 <- rbind(DF2, DF1a)
Lakes <- rbind(DF2, DF2a)

coordinates(Lakes) <- c("Longitude", "Latitude")
proj4string(Lakes) <- "+proj=longlat +datum=NAD83"

Lakes$Trecov_yr <- as.numeric(Lakes$Trecov_yr)
Lakes$Trecov_yr <- round(Lakes$Trecov_yr, 1)

names(Lakes)[names(Lakes) == "Trecov_yr"] <- "Recovery Time yrs"
m3 <- mapview(Lakes, xcol = "Longitude", ycol = "Latitude", zcol = "Recovery Time yrs", at = c(Qmin, Q1, Q2, Q3, Q4, Qmax),
  legend = "TRUE", cex = 3, crs = 4326, grid = FALSE)
print(m3)
