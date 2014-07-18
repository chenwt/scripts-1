rm(list=ls())

install.packages("languageR")
install.packages("Design")
library(party); library(languageR);

##load data---

mutfile = ''




data.controls <- cforest_unbiased(ntree=1000, mtry=3)
set.seed(47)
data.cforest <- cforest(Resp ~ x + y + z..., data = mydata, controls=data.controls)
data.cforest.varimp <- varimp(data.cforest, conditional = TRUE)