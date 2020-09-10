library(tidyverse)
library(growthcurver)
library(stringr)
library(splines)

##########################Load data and tidy############################

AMFP_counts<- read.csv("bsal_temperature_growth_counts.csv", header = TRUE, sep = ",")

head(AMFP_counts)

# Make proportional coverage and covert NAs to 0 as represent none counted;
AMFP_counts_plus_coverage<- AMFP_counts %>%
  mutate(Sporangia.count =ifelse(is.na(Sporangia_count), 0,Sporangia_count))%>%
  mutate(Sporangia.cov1 = sporan_mean * Sporangia.count) %>%
  mutate(Sporangia.cov = ifelse(is.na(Sporangia.cov1), 0,Sporangia.cov1)) %>%
  mutate(proportional_sporangia_coverage = Sporangia.cov/1000000)%>%
  select(Isolate, Plate, Well, Temp, Day, Spore_count, Sporangia.count, sporan_mean, 
         sporang_min, sporang_max, sporang_sd, Sporangia.cov, proportional_sporangia_coverage)

# Plot data-check for outliers
# Plot sporangia count and geom_smooth using method = 'loess' and formula 'y ~ x'; 
(plot_sporangia_count<- AMFP_counts_plus_coverage%>%
    ggplot(aes(x = Day)) +
    geom_jitter(aes(y=as.numeric(Sporangia.count), colour = as.factor(Temp)))+
    geom_smooth(aes(group=as.factor(Temp), color = as.factor(Temp),y=Sporangia.count))
)

# Plotting sporangia count and means;
(count_line_plots_sporangia<- AMFP_counts_plus_coverage %>%
    ggplot(aes(x = Day)) +
    geom_jitter(aes(y=Sporangia.count, colour = as.factor(Temp)))+
    stat_summary(fun.y=mean,geom="line",lwd=1,aes(group=as.factor(Temp), color = as.factor(Temp),y=Sporangia.count)))

head(AMFP_counts_plus_coverage)

##############(growth rate) plotting and equation extraction###########

# Data in growthcurver format;
# Make data into format as above: time first column, each well next column;
# Make data table for sporangia count;
AMFP_sporangia_counts_for_growthcurver <- AMFP_counts_plus_coverage %>%
  mutate(time = Day) %>%
  mutate(Well_id = paste("T",Temp,Plate,Well, sep = "_")) %>%
  select(time, Well_id, Sporangia.count)%>%
  group_by(time, Well_id) %>%
  summarise(mean_sporangia_count = mean(Sporangia.count)) %>%  # data requires measurements per well - used mean per well
  spread(key = Well_id, value = mean_sporangia_count)

# Run growthcurver - extracate r values;
# growthcurver growth curve for "per plate" - calculate all wells;
gc_AMFP_sporangia_count_out <- SummarizeGrowthByPlate(plate = AMFP_sporangia_counts_for_growthcurver,
                                 bg_correct = "none",
                                 plot_fit = TRUE, 
                                     plot_file = "gc_AMFP_cov_plots.pdf")


# Extract r values and make long data for plotting;
gc_AMFP_sporangia_count_r_values <- gc_AMFP_sporangia_count_out %>%
  select(sample,r) %>%
  mutate(Temperature = as.integer(str_sub(sample,3,-5)))

# Check for errors and outliers in growth curve fit;
# Check for error notes;
gc_AMFP_sporangia_count_out$vals$note ##NULL

# Check for outliers/poor fit- sigma is  representation of fit (residual sum of squares) - larger value poorer fit;
gc_out <- as_data_frame(gc_AMFP_sporangia_count_out)

# Plot a histogram of the sigma values in order to check for outliers;
hist(gc_out$sigma, main = "Histogram of sigma values", xlab = "sigma")

# PCA as another way to check for errors;
pca_gc_out <- as_data_frame(gc_out) 

# Prepare the gc_out data for the PCA;
rownames(pca_gc_out) <- pca_gc_out$sample

# Run PCA;
pca.res <- prcomp(pca_gc_out %>% select(k:sigma), center=TRUE, scale=TRUE)

# Plot the results;
as_data_frame(list(PC1=pca.res$x[,1],
                   PC2=pca.res$x[,2],
                   samples = rownames(pca.res$x))) %>% 
  ggplot(aes(x=PC1,y=PC2, label=samples)) + 
  geom_text(size = 3)

# Plot r values per temp and extract equation;
(r_value_plot <- gc_AMFP_sporangia_count_r_values %>%
    ggplot(aes(x=Temperature, y = r)) +
  geom_point())

# Finding fit curve-----
x<- gc_AMFP_sporangia_count_r_values$Temperature
y<- gc_AMFP_sporangia_count_r_values$r

gc_AMFP_sporangia_mean<- gc_AMFP_sporangia_count_r_values %>%
  group_by(Temperature) %>%
  summarise(r_mean = mean(r))
  
ymean<- gc_AMFP_sporangia_mean$r_mean
xmean <- gc_AMFP_sporangia_mean$Temperature

# Fit first degree polynomial equation;
fit  <- lm(y~x)
# second degree;
fit2 <- lm(y~poly(x,2,raw=TRUE))
# third degree;
fit3 <- lm(y~poly(x,3,raw=TRUE))
# fourth degree;
fit4 <- lm(y~poly(x,4,raw=TRUE))
# Generate range of 50 numbers starting from 0 and ending at 30;
xx <- seq(0,30, length=50)
plot(x,y,pch=19,ylim=c(0,20))
lines(xx, predict(fit, data.frame(x=xx)), col="red")
lines(xx, predict(fit2, data.frame(x=xx)), col="green")
lines(xx, predict(fit3, data.frame(x=xx)), col="blue")
lines(xx, predict(fit4, data.frame(x=xx)), col="purple")  

# Summaries of fit;
summary(fit)
summary(fit2)
summary(fit3)
summary(fit4)

# Multiple R-squared using different degrees of polynomials;
# 1st 0.02692
# 2nd 0.2851
# 3rd 0.2964
# 4th 0.6022

# Compare models using anova;
anova(fit, fit2)
# Model 1: y ~ x
# Model 2: y ~ poly(x, 2, raw = TRUE)
# Res.Df    RSS Df Sum of Sq      F  Pr(>F)  
# 1     20 222.12                              
# 2     19 163.20  1    58.928 6.8606 0.01688 * ##fit 2 is better
anova(fit2, fit3) # no significant improvement of fit
anova(fit2, fit4)
# Model 1: y ~ poly(x, 2, raw = TRUE)
# Model 2: y ~ poly(x, 4, raw = TRUE)
# Res.Df     RSS Df Sum of Sq      F   Pr(>F)   
# 1     19 163.196                                
# 2     17  90.797  2    72.399 6.7777 0.006848 ** ###fit4 is better

# Create a function that reflects the polynomial;
# y=ax^4 + bx^3 + cx^2 + dx + e

coef(fit4)
# (Intercept) poly(x, 4, raw = TRUE)1 poly(x, 4, raw = TRUE)2 poly(x, 4, raw = TRUE)3 
# 3.3541059987           -1.6806661400            0.4136326511           -0.0260187270 
# poly(x, 4, raw = TRUE)4 
# 0.0004657988 
#create function for the third order polynomial
fourth_order <- function(newdist, model) {
  coefs <- coef(model)
  #y = e + dx + cx^2 + bx^3 + ax^4
  res <- coefs[1] + (coefs[2] * newdist) + (coefs[3] * newdist^2) + (coefs[4] * newdist^3) + (coefs[5] * newdist^4)
  return(res)
}
yy <- fourth_order(xx,fit4)
plot(xx,yy,type="l", ylim = c(0,8))
lines(xmean,ymean,col="red")
s1 <- splinefun(xmean, ymean, method = "monoH.FC")
curve(s1(x), add = TRUE, col = "blue", n = 1001)
