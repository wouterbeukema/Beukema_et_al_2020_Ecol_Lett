library(raster)
library(ncdf4)
library(humidity)
library(matrixStats)
library(rWind)

##############################ERA5 climate data########################

# Upload climate data layers;
# (ERA5 temperature data is in Kelvin instead of Celsius);
dewpoint <- brick("ERA5_dew_2010-2018_0300and1500hours.nc", varname = "d2m")
temp <- brick("ERA5_temp_2010-2018_0300and1500hours.nc", varname = "t2m")
cloud <- brick("ERA5_cloud_2010-2018_0300and1500hours.nc", varname = "tcc")
vwind <- brick("ERA5_vwind_2010-2018_0300and1500hours.nc", varname = "v10")
uwind <- brick("ERA5_uwind_2010-2018_0300and1500hours.nc", varname = "u10")
prep <- brick("ERA5_prep_2010-2018_0300and1500hours.nc", varname = "tp")
soil <- brick("ERA5_soil_2010-2018_0300and1500hours.nc", varname = "swvl1")

# ERA5 data longitude coordinates run from 0-360, so center on the Pacific and cut Europe at the Greenwich line. raster::rotate can be used to convert these to -180 - 180, but is very slow. I therefore added '360' to all the negative coordinates in the fishnet used to extract climate values, which does the trick. 
# Upload fishnet;
fishnet <- read.csv("fishnetERA5.csv", header = TRUE, row.names="Id", sep = ";", dec = ",")

# Extract site data;
dewpoint.sites <- data.frame(extract(dewpoint, fishnet, ncol = 2))
temp.sites <- data.frame(extract(temp, fishnet, ncol = 2))
cloud.sites <- data.frame(extract(cloud, fishnet, ncol = 2))
uwind.sites <- data.frame(extract(uwind, fishnet, ncol = 2))
vwind.sites <- data.frame(extract(vwind, fishnet, ncol = 2))
prep.sites <- data.frame(extract(prep, fishnet, ncol = 2))
soil.sites <- data.frame(extract(soil, fishnet, ncol = 2))

# First restrict column names to day-month-year. Determine maxima across these, then drop year, and average across day-month.
# Keep only day, month and year in column names;
names(dewpoint.sites) <- substring(names(dewpoint.sites),2,11)
names(temp.sites) <- substring(names(temp.sites),2,11)
names(cloud.sites) <- substring(names(cloud.sites),2,11)
names(uwind.sites) <- substring(names(uwind.sites),2,11)
names(vwind.sites) <- substring(names(vwind.sites),2,11)
names(prep.sites) <- substring(names(prep.sites),2,11)
names(soil.sites) <- substring(names(soil.sites),2,11)

# Calculate max and min dewpoint temp;
# Convert to matrix for matrixStats::rowMins, matrixStats::rowMaxs
dewpoint.matrix <- data.matrix(dewpoint.sites, rownames.force = NA)
# Determine maxima and minima for each day+month+year combination;
dewpoint.max <- sapply(unique(colnames(dewpoint.sites)), function(i)
  rowMaxs(dewpoint.matrix[,colnames(dewpoint.sites) == i]))
dewpoint.min <- sapply(unique(colnames(dewpoint.sites)), function(i)
  rowMins(dewpoint.matrix[,colnames(dewpoint.sites) == i]))
# Convert back to data.frame to be able to edit column names again;
dewpoint.max.df <- as.data.frame(dewpoint.max)
dewpoint.min.df <- as.data.frame(dewpoint.min)
# Drop year;
names(dewpoint.max.df) <- substring(names(dewpoint.max.df),6,10)
names(dewpoint.min.df) <- substring(names(dewpoint.min.df),6,10)
# Calculate average maximum and minimum for every month;
dewpoint.avg.max <- as.data.frame(sapply(unique(colnames(dewpoint.max.df)), function(i)
  rowMeans(dewpoint.max.df[,colnames(dewpoint.max.df) == i])))
dewpoint.avg.min <- as.data.frame(sapply(unique(colnames(dewpoint.min.df)), function(i)
  rowMeans(dewpoint.min.df[,colnames(dewpoint.min.df) == i])))
# Convert to celsius (no need to export output to file);
dewpoint.avg.max <- K2C(data.matrix(dewpoint.avg.max, rownames.force = NA))
dewpoint.avg.min <- K2C(data.matrix(dewpoint.avg.min, rownames.force = NA))

# Same for temp, incl. mean and conversion to celsius;
temp.matrix <- data.matrix(temp.sites, rownames.force = NA)
temp.max <- sapply(unique(colnames(temp.sites)), function(i)
  rowMaxs(temp.matrix[,colnames(temp.sites) == i]))
temp.min <- sapply(unique(colnames(temp.sites)), function(i)
  rowMins(temp.matrix[,colnames(temp.sites) == i]))
temp.max.df <- as.data.frame(temp.max)
temp.min.df <- as.data.frame(temp.min)
names(temp.max.df) <- substring(names(temp.max.df),6,10)
names(temp.min.df) <- substring(names(temp.min.df),6,10)
temp.avg.max <- as.data.frame(sapply(unique(colnames(temp.max.df)), function(i)
  rowMeans(temp.max.df[,colnames(temp.max.df) == i])))
temp.avg.min <- as.data.frame(sapply(unique(colnames(temp.min.df)), function(i)
  rowMeans(temp.min.df[,colnames(temp.min.df) == i])))
temp.avg.max <- K2C(data.matrix(temp.avg.max, rownames.force = NA))
temp.avg.min <- K2C(data.matrix(temp.avg.min, rownames.force = NA))
write.csv(temp.avg.max, file = "temp-avg-max.csv")
write.csv(temp.avg.min, file = "temp-avg-min.csv")
# Daily average;
names(temp.sites) <- substring(names(temp.sites),6,10)
temp.avg <- as.data.frame(sapply(unique(colnames(temp.sites)), function(i)
  rowMeans(temp.sites[,colnames(temp.sites) == i])))
temp.avg <- K2C(data.matrix(temp.avg, rownames.force = NA))
write.csv(temp.avg, file = "temp-avg.csv")

# Calculate Relative Humidity; The RH() function only runs on numerical data, so I convert the data.frames to data.matrix. rownames.force = NA makes sure that the dataframe structure is maintained. isK = TRUE means that these data are in Kelvin.
rh.max <- as.data.frame(RH((data.matrix(temp.avg.max, rownames.force = NA)), (data.matrix(dewpoint.avg.max, rownames.force = NA)), isK = FALSE))
rh.min <- as.data.frame(RH((data.matrix(temp.avg.min, rownames.force = NA)), (data.matrix(dewpoint.avg.min, rownames.force = NA)), isK = FALSE))
write.csv(rh.max, file="rh-avg-max.csv")
write.csv(rh.min, file="rh-avg-min.csv")

# Cloud cover;
cloud.matrix <- data.matrix(cloud.sites, rownames.force = NA)
cloud.max <- sapply(unique(colnames(cloud.sites)), function(i)
  rowMaxs(cloud.matrix[,colnames(cloud.sites) == i]))
cloud.min <- sapply(unique(colnames(cloud.sites)), function(i)
  rowMins(cloud.matrix[,colnames(cloud.sites) == i]))
cloud.max.df <- as.data.frame(cloud.max)
cloud.min.df <- as.data.frame(cloud.min)
names(cloud.max.df) <- substring(names(cloud.max.df),6,10)
names(cloud.min.df) <- substring(names(cloud.min.df),6,10)
cloud.avg.max <- as.data.frame(sapply(unique(colnames(cloud.max.df)), function(i)
  rowMeans(cloud.max.df[,colnames(cloud.max.df) == i])))
cloud.avg.min <- as.data.frame(sapply(unique(colnames(cloud.min.df)), function(i)
  rowMeans(cloud.min.df[,colnames(cloud.min.df) == i])))
write.csv(cloud.avg.max, file="cloud-avg-max.csv")
write.csv(cloud.avg.min, file="cloud-avg-min.csv")

# vwind;
vwind.matrix <- data.matrix(vwind.sites, rownames.force = NA)
vwind.max <- sapply(unique(colnames(vwind.sites)), function(i)
  rowMaxs(vwind.matrix[,colnames(vwind.sites) == i]))
vwind.min <- sapply(unique(colnames(vwind.sites)), function(i)
  rowMins(vwind.matrix[,colnames(vwind.sites) == i]))
vwind.max.df <- as.data.frame(vwind.max)
vwind.min.df <- as.data.frame(vwind.min)
names(vwind.max.df) <- substring(names(vwind.max.df),6,10)
names(vwind.min.df) <- substring(names(vwind.min.df),6,10)
vwind.avg.max <- as.data.frame(sapply(unique(colnames(vwind.max.df)), function(i)
  rowMeans(vwind.max.df[,colnames(vwind.max.df) == i])))
vwind.avg.min <- as.data.frame(sapply(unique(colnames(vwind.min.df)), function(i)
  rowMeans(vwind.min.df[,colnames(vwind.min.df) == i])))

# uwind;
uwind.matrix <- data.matrix(uwind.sites, rownames.force = NA)
uwind.max <- sapply(unique(colnames(uwind.sites)), function(i)
  rowMaxs(uwind.matrix[,colnames(uwind.sites) == i]))
uwind.min <- sapply(unique(colnames(uwind.sites)), function(i)
  rowMins(uwind.matrix[,colnames(uwind.sites) == i]))
uwind.max.df <- as.data.frame(uwind.max)
uwind.min.df <- as.data.frame(uwind.min)
names(uwind.max.df) <- substring(names(uwind.max.df),6,10)
names(uwind.min.df) <- substring(names(uwind.min.df),6,10)
uwind.avg.max <- as.data.frame(sapply(unique(colnames(uwind.max.df)), function(i)
  rowMeans(uwind.max.df[,colnames(uwind.max.df) == i])))
uwind.avg.min <- as.data.frame(sapply(unique(colnames(uwind.min.df)), function(i)
  rowMeans(uwind.min.df[,colnames(uwind.min.df) == i])))

# Calculate wind speed;
# Average maximum. First convert to matrix for later use in the uv2ds function, then get rid of NA values to be able to do calculations;
uwind.avg.max.matrix <- data.matrix(uwind.avg.max, rownames.force = NA)
vwind.avg.max.matrix <- data.matrix(vwind.avg.max, rownames.force = NA)
uwind.avg.max.matrix[is.na(uwind.avg.max.matrix)] <- 0
vwind.avg.max.matrix[is.na(vwind.avg.max.matrix)] <- 0
windspeed.max <- as.data.frame(uv2ds(uwind.avg.max.matrix, vwind.avg.max.matrix))
# Remove wind direction columns;
windspeed.max <- windspeed.max[, -c(1:366)]
# Average minimum;
uwind.avg.min.matrix <- data.matrix(uwind.avg.min, rownames.force = NA)
vwind.avg.min.matrix <- data.matrix(vwind.avg.min, rownames.force = NA)
uwind.avg.min.matrix[is.na(uwind.avg.min.matrix)] <- 0
vwind.avg.min.matrix[is.na(vwind.avg.min.matrix)] <- 0
windspeed.min <- as.data.frame(uv2ds(uwind.avg.min.matrix, vwind.avg.min.matrix))
windspeed.min <- windspeed.min[, -c(1:366)]
write.csv(windspeed.max, file="wind-avg-max.csv")
write.csv(windspeed.min, file="wind-avg-min.csv")

# Precipitation;
names(prep.sites) <- substring(names(prep.sites),6,10)
prep.avg <- as.data.frame(sapply(unique(colnames(prep.sites)), function(i)
  rowMeans(prep.sites[,colnames(prep.sites) == i])))
write.csv(prep.avg, file = "prep-avg.csv")

# Soil water content;
names(soil.sites) <- substring(names(soil.sites),6,10)
soil.avg <- as.data.frame(sapply(unique(colnames(soil.sites)), function(i)
  rowMeans(soil.sites[,colnames(soil.sites) == i])))
write.csv(soil.avg, file = "soil-avg.csv")

#######################Altitude, slope, aspect#########################
 
alt1 <- raster("b10g.bil")
alt2 <- raster("c10g.bil")
alt3 <- raster("f10g.bil")
alt4 <- raster("g10g.bil")
alt <- merge(alt1, alt2, alt3, alt4, tolerance = 0.1, filename = "srtm.img")

# Raster projected in ArcGIS.

alt <- raster("srtm_proj.img")
europe <- extent(-14.5, 42.5, 34, 62)
alt <- crop(alt, europe)

fishnet <- read.csv("fishnetTOPO.csv", header = TRUE, row.names="Id", sep = ";", dec = ",")
alt.sites <- data.frame(extract(alt, fishnet, ncol=2))
write.csv(alt.sites, file = "altitude.csv")

aspect <- raster("aspect.img")
aspect <- crop(aspect, europe)
aspect.sites <- data.frame(extract(aspect, fishnet, ncol=2))
write.csv(aspect.sites, file = "aspect.csv")

slopedegr <- raster("slopedegr.img")
slopedegr.sites <- data.frame(extract(slopedegr, fishnet, ncol=2))
write.csv(slopedegr.sites, file = "slopedegr.csv")
