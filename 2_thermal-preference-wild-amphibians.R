library(lme4)
library(lmerTest) 
library(MuMIn)
library(emmeans)
library(sjPlot)
library(ggplot2)
library(stats)

############################Thermal preference###########################

Tpref <- read.table(file ="amphibian-preference-data.csv", quote="", sep = ";", header = TRUE, dec = ",")

# Reorder so Salsal goes first;
Tpref$Sp <- factor(Tpref$Sp, levels=c('Salsal', 'Ichalp', 'Lishel', 'Bufbuf', 'Rantem'))

lmerTpref <- lmer(Tpref ~ 
                    Sp + 
                    (1|Indiv),
                  data = Tpref)

summary(lmerTpref)
r.squaredGLMM(lmerTpref) 
plot(lmerTpref)
qqpoints <- qqnorm(resid(lmerTpref, which = -3))
qqline(resid(lmerTpref), col = "red") 

tab_model(lmerTpref, show.se = TRUE, show.stat = TRUE, show.df = TRUE)

########################Accompanying field measurements###################

Tfield <- read.table(file ="amphibian-field-data-selected.csv", quote="", sep = ";", header = TRUE, dec = ",")

oneway.test(C ~ entity,
            data = Tfield,
            var.equal = FALSE)

# Make unique pairs for post hoc analyses;
allPairs <- expand.grid(levels(Tfield$entity), levels(Tfield$entity))
allPairs <- unique(t(apply(allPairs, 1, sort)))
allPairs <- allPairs[ allPairs[,1] != allPairs[,2], ]
allPairs

# Results;
allResults <- apply(allPairs, 1, function(p) {
  dat <- Tfield[ Tfield$entity %in% p, ]
  ret <- oneway.test(C ~ entity, data = dat, na.action = na.omit, var.equal = FALSE)
  ret$entity <- p
  ret
})
length(allResults)
## [1] 3
allResults[[15]]
