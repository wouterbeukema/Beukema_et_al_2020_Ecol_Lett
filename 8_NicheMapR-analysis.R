library(NicheMapR)
library(biogeo)
library(measurements)
library(data.table)

################################Constants################################

# General options;
writecsv <- 0
microdaily <- 1
runshade <- 1
runmoist <- 1
snowmodel <- 0
hourly <- 0
rainhourly <- 0
IR <- 0
message <- 0
fail <- 24

# Temporal and geographical parameters;
doynum <- 365
doy <- c(1:365)
idayst <- 1
ida <- 365
HEMIS <- 1
EC <- 0.0167238

# Ground texture and wind related parameters (for wind scaling);
RUF <- 0.004
Refhyt <- 2
Usrhyt <- 0.01
ZH <- 0
D0 <- 0
Z01 <- 0
Z02 <- 0
ZH1 <- 0
ZH2 <- 0

# Emissivity, reflectivity, atmospheric water and extinction;
SLE <- 0.96
REFL <- 0.10
CMH2O <- 1

# Aerosol extinction coefficient profile, values extracted from GADS for Madison WI;
TAI <-c(0.269904738,0.266147825, 0.262442906, 0.258789404, 0.255186744, 
       0.251634356, 0.248131676, 0.2412732, 0.234606887, 0.228128378, 0.221833385, 
       0.215717692, 0.20977715, 0.204007681, 0.198405272, 0.187685927, 0.177588357, 
       0.168082846, 0.159140695, 0.150734206, 0.142836655, 0.135422274, 0.128466227, 
       0.12194459, 0.115834329, 0.110113284, 0.104760141, 0.099754417, 0.09507644, 
       0.090707328, 0.086628967, 0.082823998, 0.07927579, 0.075968428, 0.072886691, 
       0.070016034, 0.067342571, 0.064853053, 0.062534858, 0.060375964, 0.058364941, 
       0.056490925, 0.054743609, 0.053113222, 0.051590514, 0.050166738, 0.046408775, 
       0.045302803, 0.044259051, 0.043271471, 0.042334415, 0.041442618, 0.040591184, 
       0.039775572, 0.038991583, 0.038235345, 0.037503301, 0.036792197, 0.036099067, 
       0.034101935, 0.033456388, 0.032817888, 0.032184949, 0.031556287, 0.030930816, 
       0.030307633, 0.029065372, 0.027825562, 0.027205981, 0.026586556, 0.025967391, 
       0.025348692, 0.024114005, 0.023498886, 0.021669152, 0.021066668, 0.019292088, 
       0.018144698, 0.016762709, 0.015451481, 0.014949794, 0.014224263, 0.013093462, 
       0.012670686, 0.012070223, 0.011164062, 0.010241734, 0.009731103, 0.009507687, 
       0.009212683, 0.008965785, 0.008827751, 0.008710756, 0.008574128, 0.008462605, 
       0.008446967, 0.008539475, 0.009015237, 0.009748444, 0.010586023, 0.011359647, 
       0.011901268, 0.012062153, 0.011735443, 0.010882215, 0.009561062, 0.007961182, 
       0.006438984, 0.005558204, 0.006133532, 0.009277754)


hori <- rep(0, 24)
VIEWF <- 1 - sum(sin(hori * pi / 180)) / length(hori)
solonly <- 0
lamb <- 0
IUV <- 0
minshade <- 90
maxshade <- 100
PCTWET <- 0


DEP <- c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200)
ERR <- 1.5


TIMINS <- c(0, 0, 1, 1)
TIMAXS <- c(1, 1, 0, 0)

TAIRhr <- rep(0, 24*doynum)
RHhr <- rep(0, 24*doynum)
WNhr <- rep(0, 24*doynum)
CLDhr <- rep(0, 24*doynum)
SOLRhr <- rep(0, 24*doynum)
RAINhr <- rep(0, 24*doynum)
ZENhr <- rep(-1, 24*doynum)

# Creating arrays of environmental variables that are assumed not to change with month for this simulation;
MAXSHADES <- rep(maxshade, doynum)
MINSHADES <- rep(minshade, doynum)
SLES <- rep(SLE, doynum)
REFLS <- rep(REFL, doynum)
PCTWET <- rep(PCTWET, doynum)

# Create a profile of soil properties with depth for each day to be run;
Numtyps <- 1 
Nodes <- matrix(data = 0, nrow = 10, ncol = doynum)
Nodes[1, 1:doynum] <- 10

# Soil thermal parameters;
Thcond <- 1.25
Density <- 2.560
SpecHeat <- 870
BulkDensity <- 2.56/2
SatWater <- 0.26
# SoilMoist<-rep(SoilMoist,timeinterval)

# Create depth-specific soil properties matrix;
soilprops <- matrix(data = 0, 
                    nrow = 10, 
                    ncol = 5)
soilprops[1, 1] <- BulkDensity
soilprops[1, 2] <- SatWater
soilprops[1, 3] <- Thcond
soilprops[1, 4] <- SpecHeat
soilprops[1, 5] <- Density

PE <- rep(0.7, 19)
KS <- rep(0.0058, 19)
BB <- rep(1.7, 19)
BD <- rep(1.3, 19)
DD <- rep(Density, 19)
L <- c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000
R1 <- 0.001
RW <- 2.5e+10
RL <- 2e+6
PC <- -1500
SP <- 10
IM <- 1e-06
MAXCOUNT <- 500
LAI <- rep(0.1, doynum)
rainmult <- 1
maxpool <- 10
evenrain <- 1
SoilMoist_Init <- rep(0.2, 10)
moists <- matrix(nrow = 10, ncol = doynum, data = 0)
moists[1:10, ] <- SoilMoist_Init

snowtemp <- 1.5
snowdens <- 0.375
densfun <- c(0.5979, 0.2178, 0.001, 0.0038)
snowmelt <- 1
undercatch <- 1
rainmelt <- 0.0125
snowcond <- 0
intercept <- 0
grasshade <- 0

tides <- matrix(data = 0, nrow = 24 * doynum, ncol = 3) 

####################Loop to run the model for all locations##############

lalat = read.table("fishnet0375degrees - corrected.csv", 
                 sep = ";",
                 dec = ",", 
                 header = T)
ALTTALL <- fread("altitude.csv", data.table = F)
slopeALL <- fread("slopedegr.csv", data.table = F)
azmuthALL <- fread("aspect.csv", data.table = F)

dum <- as.matrix(t(rep(NA,8760))) ## a NA vector to fill gaps
saltemp <- matrix(nrow = 0,ncol = 8760)
salrh <- matrix(nrow = 0,ncol = 8760)
saltempUd20 <- matrix(nrow = 0,ncol = 8760)
saltempUD30 <- matrix(nrow = 0,ncol = 8760)

for (i in 1:length(lalat$x))
{
longlat <- lalat[i,2:3] 
  
ALAT <- as.numeric(strsplit(conv_unit(longlat[2],"dec_deg","deg_dec_min")," ")[[1]])[1]
AMINUT <- as.numeric(strsplit(conv_unit(longlat[2],"dec_deg","deg_dec_min")," ")[[1]])[2]
ALONG <- as.numeric(strsplit(conv_unit(longlat[1],"dec_deg","deg_dec_min")," ")[[1]])[1]
ALMINT <- as.numeric(strsplit(conv_unit(longlat[1],"dec_deg","deg_dec_min")," ")[[1]])[2]
ALREF <- ALONG # reference longitude for time zone

ALTT <- ALTTALL[i,2]

slope <- slopeALL[i,2]

azmuth <- azmuthALL[i,2]

TMINN <- as.numeric(fread("temp-avg-min.csv", 
                          nrows = 1, 
                          skip = i-1))[2:366]
TMAXX <- as.numeric(fread("temp-avg-max.csv", 
                          nrows =1, 
                          skip = i-1))[2:366]
RHMINN <- as.numeric(fread("rh-avg-min.csv", 
                           nrows =1, 
                           skip = i-1))[2:366] 
RHMAXX <- as.numeric(fread("rh-avg-max.csv", 
                           nrows =1, 
                           skip = i-1))[2:366]
WNMINN <- as.numeric(fread("wind-avg-min.csv", 
                           nrows =1, 
                           skip = i-1))[2:366]
WNMAXX <- as.numeric(fread("wind-avg-max.csv", 
                           nrows = 1, 
                           skip = i-1))[2:366]
CCMINN <- as.numeric(fread("cloud-avg-min.csv", 
                           nrows = 1, 
                           skip = i-1))[2:366]*100
CCMAXX <- as.numeric(fread("cloud-avg-max.csv", 
                           nrows = 1, 
                           skip = i-1))[2:366]*100
RAINFALL <- as.numeric(fread("prep-avg.csv", 
                             nrows = 1, 
                             skip = i-1))[2:366]*24000   

tannul <- mean(c(TMINN, TMAXX)) 
tannulrun <- rep(tannul, doynum)
SoilMoist <- as.numeric(fread("soil-avg.csv", 
                              nrows = 1, 
                              skip = i-1))[2:366]
soilinit <- rep(tannul, 20)

# Fill holes in the data with NA (those holes are sea and water bodies);
if (is.na(sum(longlat)) | is.na(sum(ALTT)) | is.na(sum(TMINN)) | is.na(sum(TMAXX)) | is.na(sum(RHMINN)) | is.na(sum(RHMAXX)) | is.na(sum(WNMINN))
   | is.na(sum(WNMAXX)) | is.na(sum(CCMINN)) | is.na(sum(CCMAXX)) | is.na(sum(RAINFALL))
   | is.na(sum(slope)) | is.na(sum(azmuth)) | is.na(sum(tannul)) | is.na(sum(SoilMoist))) {
  write.table(dum, 
              file = "saltemp-avgdayC.csv", 
              sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE)
  write.table(dum, file = "salrhdayC.csv", 
              sep = ",", 
              append = TRUE)
  write.table(dum, file = "saltemp-UD-20-dayC.csv", 
              sep = ",", 
              append = TRUE)
  write.table(dum, file = "saltemp-UD-30-dayC.csv", 
              sep = ",", 
              append = TRUE)
}

  else {
# Input parameter vector;
microinput <- c(doynum, RUF, ERR, Usrhyt, Refhyt, Numtyps, Z01, Z02, ZH1, ZH2, idayst, ida, HEMIS, ALAT, AMINUT, ALONG, ALMINT, ALREF, slope, azmuth, ALTT, CMH2O, microdaily, tannul, EC, VIEWF, snowtemp, snowdens, snowmelt, undercatch, rainmult, runshade, runmoist, maxpool, evenrain, snowmodel, rainmelt, writecsv, densfun, hourly, rainhourly, lamb, IUV, RW, PC, RL, SP, R1, IM, MAXCOUNT, IR, message, fail, snowcond, intercept, grasshade, solonly, ZH, D0)

# Final input list - all these variables are expected by the input argument of the Fortran microclimate subroutine;
micro1 <- list(microinput = microinput, tides = tides, doy = doy, SLES = SLES, DEP = DEP, Nodes = Nodes, MAXSHADES = MAXSHADES, MINSHADES = MINSHADES, TIMAXS = TIMAXS, TIMINS = TIMINS, TMAXX = TMAXX, TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN, CCMAXX = CCMAXX, CCMINN = CCMINN, WNMAXX = WNMAXX, WNMINN = WNMINN, TAIRhr = TAIRhr, RHhr = RHhr, WNhr = WNhr, CLDhr = CLDhr, SOLRhr = SOLRhr, RAINhr = RAINhr, ZENhr = ZENhr, REFLS = REFLS, PCTWET = PCTWET, soilinit = soilinit, hori = hori, TAI = TAI, soilprops = soilprops, moists = moists, RAINFALL = RAINFALL, tannulrun = tannulrun, PE = PE, KS = KS, BB = BB, BD = BD, DD = DD, L = L, LAI = LAI)

# Run the microclimate model in Fortran;
micro <- microclimate(micro1) 
nyears <- 1
# Run the ectotherm model;
ectoCWfffShade<-ectotherm(rainfall = RAINFALL, 
                          elev = as.numeric(ALTT), 
                          alpha_sub = 1 - REFL, 
                          minshades = rep(minshade,
                          length(MAXSHADES)),
                          maxshades = MAXSHADES, 
                          longitude = longlat[1], 
                          latitude = longlat[2],
                          nyears = nyears, 
                          DEP = DEP, 
                          pct_wet = 100, 
                          nocturn = 1,
                          crepus = 1, 
                          diurn = 0,
                          maxdepth = 3) 

# Export salamander steady state temperature;
write.table(as.matrix(t(ectoCWfffShade$environ[,5])), 
            file = "saltemp-avgdayC.csv", 
            sep = ",", append = TRUE, 
            col.names = FALSE, 
            row.names = TRUE)
# Relative humidity;
write.table(as.matrix(t(ectoCWfffShade$environ[,14])), 
            file = "salrhdayC.csv", 
            sep = ",", 
            append = TRUE, col.names = FALSE) 
# Temperature at 20cm below surface;
write.table(as.matrix(t(ectoCWfffShade$soil[,8])), 
            file="saltemp-UD-20-dayC.csv", 
            sep = ",", 
            append = TRUE, col.names = FALSE) 
# Temperature at 30cm below surface (this is not use in the subsequent analysis);
write.table(as.matrix(t(ectoCWfffShade$environ[,9])), 
            file = "saltemp-UD-30-dayC.csv", 
            sep = ",", 
            append = TRUE, 
            col.names = FALSE) 
  } 
}


# Read the just-created S.sal temperature;
saltemp2<-read.table("saltemp-avgdayC.csv", sep=",", header = T)

#######################Estimating Bsal growth############################

# AMFP13/1 growth formula:  y=0.0004657988x^4  -0.0260187270x^3 + 0.4136326511x^2 + -1.6806661400x + 3.3541059987

# First, from 'raw' Te;
saltempraw <- as.matrix(saltemp2) # read.csv("saltemp_complete.csv"))
# Bsal growth as a function of S.sal temperature, adding limits of 0 growth;
salgrowraw <- 0.0004657988*saltempraw^4  - 0.0260187270*saltempraw^3 + 0.4136326511*saltempraw^2 - 1.6806661400*saltempraw + 3.3541059987 
salgrowraw<-ifelse(saltempraw < 2.654177,0,ifelse(saltempraw > 21.06894, 0, salgrowraw))
write.csv(salgrowraw, file = "salgrowth_completedayC.csv")
rm(salgrowraw)

gc();gc();gc();gc() # liberate some RAM

# Then, under humidity-constrained activity (assuming that salamanders retreat to a depth of 20cm when RH is below 85%);
relh <- read.csv("/mnt/e5c52b67-efad-4d60-aad6-3f45bd275b02/data0375deg/salrhdayC.csv", header = FALSE)
# Select only the rows with information;
relh <- subset.data.frame(relh, relh$V1 == 1)  
relh <- relh[2:7208,]
relh <- as.data.frame(sapply(relh,function(x) as.numeric(as.character(x))))
relh <- as.matrix(relh)

tem20 <- as.matrix(read.csv("/mnt/e5c52b67-efad-4d60-aad6-3f45bd275b02/data0375deg/saltemp-UD-20-dayC.csv"))
tem20 <- read.csv("/mnt/e5c52b67-efad-4d60-aad6-3f45bd275b02/data0375deg/saltemp-UD-20-dayC.csv", header = FALSE)
# Select only the rows with information;
tem20 <- subset.data.frame(tem20, tem20$V1 == 1)
tem20 <- tem20[2:7208,]
tem20 <- as.data.frame(sapply(tem20,function(x) as.numeric(as.character(x))))
tem20 <- as.matrix(tem20)

# Salamander temperature assuming that, if relative humidity falls below 85% the animal retreats 20cm below ground (and temperature becomes equal to that of surrounding soil);
saltemp85 <- ifelse(relh <= 0.85, tem20, saltempraw)
write.csv(saltemp85, file = "saltemp85rhdayC.csv")

gc();gc();gc();gc() # liberate some RAM

salgrowr85<-0.0004657988*saltemp85^4  -0.0260187270*saltemp85^3 + 0.4136326511*saltemp85^2 -1.6806661400*saltemp85 + 3.3541059987 ## Bsal growth as a function of humidity constrained S.sal temperature.
salgrowr85<-ifelse(saltemp85 < 2.654177,0,ifelse(saltemp85 > 21.06894, 0, salgrowr85)) # Adding limits of 0 growth.

# Bsal growth constrained by relative humidity;
write.csv(salgrowr85, file = "salgrowth_085dayC.csv")
