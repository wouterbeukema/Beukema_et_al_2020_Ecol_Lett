library(lme4)
library(lmerTest)
library(MuMIn)
library(sjPlot)
library(ggplot2)

##############################Pre-processing#############################

surveydata <-read.table(file="surveydata.csv", quote="", sep = ";", header = TRUE, dec = ",")

surveydata <- surveydata[complete.cases(surveydata), ] # Remove NA;

surveydata$season <- as.factor(surveydata$season)
surveydata$sp <- as.factor(surveydata$sp)
surveydata$humidity <- as.numeric(surveydata$humidity)
surveydata$sector <- as.factor(surveydata$sector)

####################################lmer##################################

lmer1 <- lmer(C ~ sp +
                humidity + 
                minutessincesundown +
                sector + 
                (1|samplEvent/season),
              data = surveydata)

summary(lmer1)
# p-values are estimated via t-tests using the Satterthwaite approximations to degrees of freedom 

r.squaredGLMM(lmer1) 

plot(lmer1) # Check residuals

qqpoints <- qqnorm(resid(lmer1, which = -3)) # Check normality
qqline(resid(lmer1), col = "red") 

# Output table;
tab_model(lmer1, show.se = TRUE, show.stat = TRUE, show.df = TRUE)

