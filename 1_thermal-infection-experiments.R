library(mgcv)
library(graphics)
library(itsadug)
library(ggplot2)
library(gamm4)
library(lattice)
library(lme4)

###########################Experiment 1################################

Exp1 <- read.table(file ="thermal-infection-experiment1.csv", sep = ";", header = TRUE)

ls.str(Exp1) 
Exp1$Treatment <- as.factor(Exp1$Treatment) 
Exp1$Gradient <- as.factor(Exp1$Gradient) 
Exp1$Indiv <- as.factor(Exp1$Indiv)

# Set non-inoculated animals as reference;
Exp1$Treatment <- relevel(Exp1$Treatment, ref = "0")

#Transform Treatment into an ordered factor. This permits testing for differences in the trend between the two treatments over time. 
Exp1$OFTreatment <- as.factor(Exp1$Treatment)
Exp1$OFTreatment <- as.ordered(Exp1$OFTreatment)
# Change contrast to treatment coding (difference curves);
contrasts(Exp1$OFTreatment) <- 'contr.treatment'
contrasts(Exp1$OFTreatment) # Inspect contrasts

#GAMM;
Exp1.orderedfactor <- bam(Avg_Tp ~ OFTreatment + s(Day) + s(Day, by = OFTreatment)
                          + s(Day, Indiv, bs = "fs", m = 1), 
                          data = Exp1, 
                          family = gaussian())
check_resid(Exp1.orderedfactor)
summary(Exp1.orderedfactor)

# The data need to be sorted to correct for autocorrelation;
Exp1 <- start_event(Exp1, column = "Day", event = c("Indiv", "Trial"))
# Calculation of the value of the autocorrelation parameter rho; 
(valRho <- acf(resid(Exp1.orderedfactor), plot = FALSE)$acf[2])

# GAMM incl. AR1 model;
Exp1.orderedfactor <- bam(Avg_Tp ~ OFTreatment + s(Day) + s(Day, by = OFTreatment)
                          + s(Day, Indiv, bs = "fs", m = 1),
                          AR.start = Exp1$start.event, rho = valRho,
                          family = gaussian(),
                          data = Exp1
)
check_resid(Exp1.orderedfactor)
summary(Exp1.orderedfactor)
# Note that the warning 'model has repeated 1-d smooths of same variable' can be ignored, as this simply indicates that you have smooths over time both in the random effects as in the fixed effects.

# Smooth and difference plots for the parametric term OFTreatment
par(mfrow = c(1,2), cex = 1.1)
plot_smooth(Exp1.orderedfactor, 
                 view = "Day",
                 plot_all = "OFTreatment", 
                 ylim=c(0,30), 
                 rm.ranef=TRUE)
plot_diff(Exp1.orderedfactor, 
          view = "Day",
          comp = list(OFTreatment = c("0","1")), 
          rm.ranef=TRUE)                   

# Plotting the random smooth (s(Day, Indiv, bs="fs", m=1)), first per Indiv, and subsequently averaged per Treatment. 
par(mfrow = c(1,2), cex = 1.1)
inspect_random(Exp1.orderedfactor, 
                    select = 3, 
                    main = 's(Day, Indiv, bs="fs", m=1)')             

inspect_random(Exp1.orderedfactor, 
                    select = 3, 
                    main = 'Averages', 
                    fun = mean, 
                    cond = list(Indiv = c("A1", "A2", "A4", "A5", "B1", "B2", "B3", "B4", "B5")))

inspect_random(Exp1.orderedfactor, 
                    select = 3, 
                    fun = mean, 
                    cond = list(Indiv = c("C1", "C2", "C3", "C4", "C5", "D1", "D2", "D3",                "D4", "D5")), 
                    add=TRUE, 
                    col = 'red', 
                    lty = 5)

legend('bottomleft',
       legend = c('non-exposed', 'exposed'),
       col = c('black', 'red'), 
       lty = c(1,5),
       bty = 'n')

###########################Experiment 2################################

Exp2 <-read.table(file = "thermal-infection-experiment2.csv", sep=";", header=TRUE, dec = ",")

ls.str(Exp2)
Exp2$Treatment <- as.factor(Exp2$Treatment)
Exp2$Indiv <- as.factor(Exp2$Indiv)

Exp2$Treatment <- relevel(Exp2$Treatment, ref = "1")

Exp2$OFTreatment <- as.factor(Exp2$Treatment)
Exp2$OFTreatment <- as.ordered(Exp2$OFTreatment)
contrasts(Exp2$OFTreatment) <- 'contr.treatment'
contrasts(Exp2$OFTreatment)

Exp2.orderedfactor <- bam(Avg_Tp ~ OFTreatment + s(Day) + s(Day, by = OFTreatment) 
                          + s(Day, Indiv, bs = "fs", m = 1), 
                          family = gaussian(),
                          data = Exp2)

Exp2 <- start_event(Exp2, column = "Day", event = c("Indiv", "Trial"))
(valRho <- acf(resid(Exp2.orderedfactor), plot = FALSE)$acf[2])

Exp2.orderedfactor <- bam(Avg_Tp ~ OFTreatment + s(Day) + s(Day, by = OFTreatment)
                          + s(Day, Indiv, bs = "fs", m = 1),
                          AR.start = Exp2$start.event, rho = valRho,
                          family = gaussian(),
                          data = Exp2
)
check_resid(Exp2.orderedfactor)
summary(Exp2.orderedfactor)

par(mfrow = c(1,2), cex = 1.1)
m <- plot_smooth(Exp2.orderedfactor, 
                 view = "Day",
                 plot_all = "OFTreatment", 
                 ylim=c(0,30), 
                 rm.ranef=TRUE)
plot_diff(Exp2.orderedfactor, 
          view = "Day",
          comp = list(OFTreatment = c("1","2")), 
          rm.ranef=TRUE)                   

par(mfrow = c(1,2), cex = 1.1)
inspect_random(Exp2.orderedfactor, 
                    select = 3, 
                    main = 's(Day, Indiv, bs="fs", m=1)')             

inspect_random(Exp2.orderedfactor, 
                    select = 3, 
                    main = 'Averages', 
                    fun = mean, 
                    cond = list(Indiv = c("A1", "A2", "A4", "A5", "B1", "B2", "B3", "B4", "B5")))

inspect_random(Exp2.orderedfactor, 
                    select = 3, 
                    fun = mean, 
                    cond = list(Indiv = c("C1", "C2", "C3", "C4", "C5", "D1", "D2", "D3", "D4", "D5")), 
                    add = TRUE, 
                    col = 'red', 
                    lty = 5)

legend('bottomleft',
       legend = c('non-infected', 'infected'),
       col = c('black', 'red'), 
       lty = c(1,5),
       bty = 'n')

#########################Experiment 2 Dry vs. wet#######################

Exp2dvsw <-read.table(file = "thermal-infection-experiment2-dryvswet.csv", sep=";", header=TRUE, dec = ",")

ls.str(Exp2dvsw)
Exp2dvsw$Treatment <- as.factor(Exp2dvsw$Treatment)
Exp2dvsw$Indiv <- as.factor(Exp2dvsw$Indiv)

Exp2dvsw$Treatment <- relevel(Exp2dvsw$Treatment, ref = "1")

Exp2dvsw$OFTreatment <- as.factor(Exp2dvsw$Treatment)
Exp2dvsw$OFTreatment <- as.ordered(Exp2dvsw$OFTreatment)
contrasts(Exp2dvsw$OFTreatment) <- 'contr.treatment'
contrasts(Exp2dvsw$OFTreatment)

Exp2dvsw.orderedfactor <- bam(Avg_Tp ~ OFTreatment + s(Day) + s(Day, by = OFTreatment) 
                              + s(Day, Indiv, bs = "fs", m = 1), 
                              family = gaussian(),
                              data = Exp2dvsw)

check_resid(Exp2dvsw.orderedfactor)
summary(Exp2dvsw.orderedfactor)

Exp2dvsw <- start_event(Exp2dvsw, column = "Day", event = c("Indiv", "Trial"))
(valRho <- acf(resid(Exp2dvsw.orderedfactor), plot = FALSE)$acf[2])

Exp2dvsw.orderedfactor <- bam(Avg_Tp ~ OFTreatment + s(Day) + s(Day, by = OFTreatment)
                              + s(Day, Indiv, bs = "fs", m = 1),
                              AR.start = Exp2dvsw$start.event, rho = valRho,
                              family = gaussian(),
                              data = Exp2dvsw
)
check_resid(Exp2dvsw.orderedfactor)
summary(Exp2dvsw.orderedfactor)

par(mfrow = c(1,2), cex = 1.1)
m <- plot_smooth(Exp2dvsw.orderedfactor, 
                 view = "Day",
                 plot_all = "OFTreatment", 
                 ylim=c(0,30), 
                 rm.ranef=TRUE)
plot_diff(Exp2dvsw.orderedfactor, 
          view = "Day",
          comp = list(OFTreatment = c("1","2")), 
          rm.ranef=TRUE)                   

par(mfrow = c(1,2), cex = 1.1)
inspect_random(Exp2dvsw.orderedfactor, 
                    select = 3, 
                    main = 's(Day, Indiv, bs="fs", m=1)')             

inspect_random(Exp2dvsw.orderedfactor, 
                    select = 3, 
                    main = 'Averages', 
                    fun=mean, 
                    cond=list(Indiv=c("A1", "A3", "A5", "B1", "B4", "C1", "C3", "C5", "D1", "D4")))

inspect_random(Exp2dvsw.orderedfactor, 
                    select=3, 
                    fun=mean, 
                    cond=list(Indiv=c("A2", "A4", "B2", "B3", "B5", "C2", "C4", "D2", "D3", "D5")), 
                    add=TRUE, 
                    col='red', 
                    lty=5)

legend('bottomleft',
       legend=c('wet', 'dry'),
       col=c('black', 'red'), 
       lty=c(1,5),
       bty='n')

###########################Experiment 3################################

Exp3 <- read.table(file = "thermal-infection-experiment3", sep=";", dec = ",", header = TRUE)

ls.str(Exp3)
Exp3$Treatment <- as.factor(Exp3$Treatment)
Exp3$Gradient <- as.factor(Exp3$Gradient)
Exp3$Indiv <- as.factor(Exp3$Indiv)

Exp3$Treatment <- relevel(Exp3$Treatment, ref = "1")

Exp3$OFTreatment <- as.factor(Exp3$Treatment)
Exp3$OFTreatment <- as.ordered(Exp3$OFTreatment)
contrasts(Exp3$OFTreatment) <- 'contr.treatment'
contrasts(Exp3$OFTreatment)

Exp3.orderedfactor <- bam(Avg_Tp ~ OFTreatment + s(Day) + s(Day, by = OFTreatment) 
                          + s(Day, Indiv, bs = "fs", m = 1), 
                          family = gaussian(),
                          data = Exp3)

Exp3 <- start_event(Exp3, column = "Day", event = c("Indiv", "Trial"))
(valRho <- acf(resid(Exp3.orderedfactor), plot = FALSE)$acf[2])

Exp3.orderedfactor <- bam(Avg_Tp ~ OFTreatment + s(Day) + s(Day, by = OFTreatment)
                          + s(Day, Indiv, bs = "fs", m = 1),
                          AR.start = Exp3$start.event, rho = valRho,
                          family = gaussian(),
                          data = Exp3
)

check_resid(Exp3.orderedfactor)
summary(Exp3.orderedfactor)

par(mfrow = c(1,2), cex = 1.1)
m <- plot_smooth(Exp3.orderedfactor, 
                 view = "Day",
                 plot_all = "OFTreatment", 
                 ylim=c(0,30), 
                 rm.ranef=TRUE)
plot_diff(Exp3.orderedfactor, 
          view = "Day",
          comp = list(OFTreatment = c("1","2")), 
          rm.ranef=TRUE)                   

par(mfrow = c(1,2), cex = 1.1)
inspect_random(Exp3.orderedfactor, 
                    select = 3, 
                    main = 's(Day, Indiv, bs="fs", m=1)')

inspect_random(Exp3.orderedfactor, 
                    select = 3, 
                    main = 'Averages', 
                    fun=mean, 
                    cond=list(Indiv=c("A1", "A2", "A4", "A5")))

inspect_random(Exp3.orderedfactor, 
                    select=3, 
                    fun=mean, 
                    cond=list(Indiv=c("C1", "C2", "C3", "C4", "C5", "D1", "D2", "D3",                "D4", "D5")), 
                    add=TRUE, 
                    col='red', 
                    lty=5)

legend('bottomleft',
       legend=c('non-infected', 'infected'),
       col=c('black', 'red'), 
       lty=c(1,5),
       bty='n')
