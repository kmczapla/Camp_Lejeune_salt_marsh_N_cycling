install.packages("lme4")
?lme4
library(lme4)

setwd("C://Users//Ken//Desktop//NRE data//R")
data1=read.csv("N_cycling.csv", header=TRUE)
head(data1)
#test data normality)
install.packages("e1071")
library(e1071)
#http://rstudio-pubs-static.s3.amazonaws.com/1563_1ae2544c0e324b9bb7f6e63cf8f9e098.html
?skewness
skewness(data1$Net.C.Budget, na.rm=TRUE)
skew.score<-function(c,x) (skewness(log(x+c)))^2
cval <- seq(0, 40, l = 101)
skew <- cval * 0
for (i in 1:length(cval)) skew[i] <- skewness(log(cval[i] + data1$sulfide..uM..5.cm))
plot(cval, skew, type = "l", ylab = expression(b[3](c)), xlab = expression(c))
abline(h = 0, lty = 3)
best.c<-optimise(skew.score, c(0,20), x=data1$Stem.Density)$minimum
best.c
sulfide.transformed<-log(data1$sulfide..uM..5.cm+best.c)
hist(log(data1$D14))
sulfide.transformed
skewness(, na.rm=TRUE)
shapiro.test(data1$Stem.Length)
hist(sqrt(data1$Methanotrophs))
boxplot(log(data1$Methane)~data1$Location+data1$Season)
SULFIDE
version
SULFIDE=data1$ab
qqPlot((data1$D14)) #sulfide log-normal?
qqPlot(data1$D14)
    
#test homogeneity of variance
install.packages("car")
library(car)
install.packages('haven')
library('haven')
leveneTest(SULFIDE~Location, data=data)
leveneTest(SULFIDE~Season, data=data)
leveneTest(SULFIDE~Treatment, data=data) 
?leveneTest
 #90 observations on the following variables:
    #block: a random variable with 6 levels nested in location, random
    #treatment: a factor with 2 levels (control and fertilized), fixed
    #location: a factor with 3 levels (edge, interior, Traps), fixed
    #observation: a numeric dependent variable
    #especially interested in interaction between treatment and location
    #treatment and season interaction?
    #repeated measures 
aov1=aov(log(data1$D15)~Location*Treatment*Season,data=data1)#two-way anova
na.omit(data1)
summary(aov1)
update
data$DNRAn <- as.numeric(data$DNRAn)

m1 <-lmer(data$DNRAn)~Location*Treatment*Location:Treatment+(1|Block:Location), data=data, REML=FALSE)
m2 <-lmer(log(data$DNRAn)~Location*Treatment*Season+(1|Block:Location), data=data, REML=TRUE)
#probably should use REML, but no REML has lower AIC
boxplot(data$DNRAn)~data$Location
summary(m1)
print(m1, correlation=TRUE)
summary(m2)
dim(vcov(m1))
print(m2, correlation=TRUE)
Anova(m1) 
data
TukeyHSD(m1, "Location")
library(lsmeans)
?lsmip
lsmip(m1, Treatment~Location)
aov1=aov(sqrt(sqrt(data1$Anaerolineae))~Location*Treatment*Season,data=data1)
summary(aov1)
TukeyHSD(aov1, "Season")
TukeyHSD(aov1, "Location")
TukeyHSD(aov1, "Treatment")
TukeyHSD(aov1, "Location:Treatment")
TukeyHSD(aov1, "Location:Season")
TukeyHSD(aov1, "Treatment:Season")
TukeyHSD(aov1, "Location:Treatment:Season")
install.packages("multcomp")
library(multcomp)
??multcomp
??glht
summary(glht(m1, linfct=mcp(Location="Tukey")))
summary(glht(m1, linfct=mcp(Treatment="Tukey")))
summary(glht(m1, linfct=mcp(Location:Treatment="Tukey")))
install.packages("agricolae")
library(agricolae)
tx <- with(m1, interaction(Location))
amod <- aov(Richness ~ tx, data=diversity)
library(agricolae)
interaction <-HSD.test(amod, "tx", group=TRUE)
interaction
install.packages("lsmeans")
library(lsmeans)
lsmeans(m1, pairwise~Location*Treatment, adjust="tukey")
?glht
m2 <-lmer(log(Above.ground.Biomass..g.m2.)~Location+Treatment+Season+Location*Treatment+(1|Block:Location), data=data, REML=FALSE)
AIC(m2)
AIC(m1)

summary(m2)
anova(m1)
Anova(m2, "II")
plot(m2)
hist(log(data1$Stem.Height))

# installing/loading the package:
if(!require(installr)) {
  install.packages("installr"); 
  require(installr)
} #load / install+load installr

updateR()