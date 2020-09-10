library(scales)
library(matrixStats)

setwd("C:/Users/woute/Desktop/Mechanistic models/0.375 model April")

xy <- read.csv("xy.csv", header = TRUE, row.names="Id", sep = ";", dec = ",")

####################Salsal steady-state body temp#####################

# Upload data. Each column represents S. salamandra steady-state body temperature at one hour (24 in a day) for each day of the year; 
sal_85 <- read.csv("data/saltemp_085rh.csv", header = TRUE, sep = ",", dec = ".")

# Remove redundant columns at the beginning of the file;
sal_85 <- sal_85[-1:-3]

# Create function to obtain the average of each group of 24 columns (daily average);
byapply <- function(x, by, fun, ...)
{
  # Create index list
  if (length(by) == 1)
  {
    nc <- ncol(x)
    split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
  } else # 'by' is a vector of groups
  {
    nc <- length(by)
    split.index <- by
  }
  index.list <- split(seq(from = 1, to = nc), split.index)
  
  # Pass index list to fun using sapply() and return object
  sapply(index.list, function(i)
  {
    do.call(fun, list(x[, i], ...))
  })
}

# Run function;
sal_85.avg <- as.data.frame(byapply(sal_85, 24, rowMeans))

# Test to make sure it returns expected result;
sal_85.avg.test <- rowMeans(sal_85[, 25:48]) 
all.equal(sal_85.avg[, 2], sal_85.avg.test) 

# Identify streaks of five days above a desired temperature - in this case daily average of 20C;
# Convert to TRUE-FALSE to id values above 20 degrees;
sal_85.avg.above20 <- as.data.frame(sal_85.avg > 20)

# Function to detect a sequence of five same values;
detect_streak <- function(x) {
  streaks <- rle(x)$lengths
  return(any(streaks >= 5))
}

# Create a column named 'streaks' containing FALSE if no 5 consecutive values of the desired value were found, TRUE if they were;
sal_85.avg.above20[sal_85.avg.above20=="FALSE"]<-NA
sal_85.avg.above20$streaks <- vector(mode = "logical", length = nrow(sal_85.avg.above20))
for (i in 1:nrow(sal_85.avg.above20)) {
  sal_85.avg.above20$streaks[i] <- detect_streak(sal_85.avg.above20[i,2:11])
}
# Combine with xy data;
sal_85_20 <- cbind(xy, sal_85.avg.above20$streaks)
#write file;
write.csv(sal_85_20, file ="sal_85_20.csv")

# Determine mean Tb across the year;
sal_85_avg <- rowMeans(sal_85.avg, na.rm=TRUE)
sal_85_avg <- cbind(xy, sal_85_avg)
write.csv(sal_85_avg, file ="sal_85_avg.csv")
# Max;
sal_85_max <- rowMaxs(as.matrix(sal_85.avg, na.rm=TRUE))
sal_85_max <- cbind(xy, sal_85_max)
write.csv(sal_85_max, file ="sal_85_max.csv")
# Min;
sal_85_min <- rowMins(as.matrix(sal_85.avg, na.rm=TRUE))
sal_85_min <- cbind(xy, sal_85_min)
write.csv(sal_85_min, file ="sal_85_min.csv")

##########Bsal growth relative to Salsal steady-state body temp#######

# Upload data. Each column represents Bsal growth relative to Salsal steady-state body temp at one hour (24 in a day) for each day of the year
bsal_85 <- read.csv("data/salgrowth_085.csv", header = TRUE, sep = ",", dec = ".")

# Remove first four redundant rows;
bsal_85 <- bsal_85[-1:-3]

# Create daily average;
bsal_85.daily <- as.data.frame(byapply(bsal_85, 24, rowMeans))

# Create seasonal averages as preparation for map figure creation;
# Subset seasonal data;
winter1 <- subset(bsal_85.daily, select = 336:365) # 1 Dec to 31 Dec
winter2 <- subset(bsal_85.daily, select = 1:59) # 1 Jan to 28 Feb
spring <- subset(bsal_85.daily, select = 60:152) # 1 Mar to 31 May
summer <- subset(bsal_85.daily, select = 153:244) # 1 Jun to 31 Aug
autumn <- subset(bsal_85.daily, select = 245:335) # 1 Sept to 30 Nov
winter <- cbind(winter1, winter2)

# Determine mean for each season;
wintermean <- rowMeans(winter, na.rm = TRUE)
springmean <- rowMeans(spring, na.rm = TRUE)
summermean <- rowMeans(summer, na.rm = TRUE)
autumnmean <- rowMeans(autumn, na.rm = TRUE)

# Combine seasons, add xy data, write csv with mean Bsal growth per cell per season;
seasons <- cbind(wintermean, springmean, summermean, autumnmean)
seasonsxy <- cbind(xy, seasons)
write.csv(seasonsxy, file ="seasonsxy.csv")

# The code creates a vector that contains, for each cell, the day of the year during which growth is highest. When plotted across an 4-color algorithmic color map in for instance ArcGIS, one can show where Bsal growth peaks during which season. 
# Create a copy for reordering, as we want the file to start at spring with column number 1 to facilitate attributing a color map to the data;
springrescale <- spring
# Make a vector with new column names for spring;
springdays <- 1:93
colnames(springrescale) <- springdays
# Other seasons;
summerrescale <- summer
autumnrescale <- autumn
winterrescale <- winter
summerdays <- 94:185
autumndays <- 186:276
winterdays <- 277:365
colnames(summerrescale) <- summerdays
colnames(autumnrescale) <- autumndays
colnames(winterrescale) <- winterdays
daily.rescale <- cbind(springrescale, summerrescale, autumnrescale, winterrescale)

# The apply command below does not work on dataframes, list, etc containing NA. Therefore, convert NA to 0 first;
is.na.data.frame <- function(x)
  do.call(cbind, lapply(x, is.na))
daily.rescale[is.na(daily.rescale)] <- 0

# Select the colname which has the highest Bsal growth, and use that colname as cell value;
highestdailyres <- as.data.frame(colnames(daily.rescale)[apply(daily.rescale, 1, which.max)])
# Add column name;
colnames(highestdailyres) <- c("day")
# Add xy;
highestdailyresxy <- cbind(xy, highestdailyres)
# Write csv;
write.csv(highestdailyresxy, file ="highestdailyresxy.csv")
