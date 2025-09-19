install.packages('piecewiseSEM')
library(piecewiseSEM)

setwd("C://Users//kencz//Desktop//R")
na.omit(SEMdata)

#DCERP SEMS
SEMdata=read.csv("N_Cycling_noNA.csv")
SEMdata

model <- psem(lm(sqrt(D14) ~ Temp + H2S + sqrt(NH4) + DOC + Fert + sqrt(DNRA14), SEMdata),
              lm(sqrt(DNRA14) ~ Temp + H2S + sqrt(NH4) + Fert + DOC, SEMdata))
summary(model, .progressBar = F)


