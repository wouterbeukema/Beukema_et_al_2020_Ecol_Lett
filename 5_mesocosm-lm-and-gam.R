library(sjPlot)
library(mgcv)
library(itsadug)

##################################LM#####################################

mesocosm <- read.table(file="mesocosm-data.csv", quote="", sep = ";", header = TRUE, dec = ",")

# transformed Tb of C3031 as it was right-skewed;
mesocosm.all$sqrtC3031 <- sqrt(mesocosm.all$C3031)

# Species-specific models;
C3031 <- lm(sqrtC3031 ~ C +
               Humid +
               Rain,
             data = mesocosm)
summary(C3031)
plot(C3031)
# Output table;
tab_model(C3031, show.se = TRUE, show.stat = TRUE, show.df = TRUE)

C3032 <- lm(C3032 ~ C +
              Humid +
              Rain,
            data = mesocosm)
summary(C3032)
plot(C3032)
# Output table;
tab_model(C3032, show.se = TRUE, show.stat = TRUE, show.df = TRUE)

C3015 <- lm(C3015 ~ C +
              Humid +
              Rain,
            data = mesocosm)
summary(C3015)
plot(C3015)
# Output table;
tab_model(C3015, show.se = TRUE, show.stat = TRUE, show.df = TRUE)

##################################GAM##################################

mesocosm-long <- read.table(file="mesocosm-data-longformat.csv", quote="", sep = ";", header = TRUE, dec = ",")

#set Env as reference
mesocosm-long$Type <- relevel(mesocosm-long$Type, ref = "Env")
mesocosm-long$OFType <- as.factor(mesocosm-long$Type)

#Transform Treatment into an ordered factor. This will permit testing for differences in the trend between the two treatments over time. 
mesocosm-long$OFType <- as.ordered(mesocosm-long$OFType)

contrasts(mesocosm-long$OFType) <- 'contr.treatment'

contrasts(mesocosm-long$OFType)

gam <- gam(C ~ OFType + s(Datetime) + s(Datetime, by = OFType) +
             s(Humid) +
             s(Rain),
           data = mesocosm-long)
par(mfrow = c(2,2))
gam.check(gam)
summary(gam)
plot(gam)

par(mfrow = c(1,2), cex = 1.1)
m <- plot_smooth(gam, 
                 view = "Datetime",
                 plot_all = "OFType", 
                 ylim=c(0,30), 
                 rm.ranef = TRUE)

plot_diff(gam, 
          view = "Datetime",
          comp = list(OFType = c("C3031","C3015")), 
          rm.ranef = TRUE) 
