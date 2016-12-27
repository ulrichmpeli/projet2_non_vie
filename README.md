# projet2_non_vie



rm(list=ls())
setwd("C:/Users/Utilisateur/Documents/ENSAE 2016_2017/actuariat assurance non vie/DM2")
load("base_ensae_1.RData")

library(splines)
library(lsr)
library(MASS)
library(ROCR)
library(pROC)

library(tm)
library(utils)
library(SnowballC)
library(ggplot2)
library(rpart)
library(rpart.plot)
library(rattle)
library(tree)
library(hmeasure)
library(wordcloud)
library(texreg)
library(stargazer)
library(ade4)
library(doBy)
library(xlsx)
library(party)
library(ipred)
library(randomForest)
library(parallel)
library(ParallelPC)
library(ParallelForest)
library(doParallel)
library(dplyr)
library(class)
library(mgcv)
library(AER)
library(pscl)


# == Recherche et suppression des doublons ===========================================
data1 = base_ensae_1
descr = summary(data1)
table(duplicated(data1))
x = duplicated(data1)
data2 = cbind(data1,x)
data2 = data2[data2$x=="FALSE",]
table(duplicated(data2))
summary(data2)




#### ON CALCUL LA MOYENNE annuelle de la fréquence et sa variance empirique

vY<-data2$Nb1
vE<-data2$Exppdays/365
m<-sum(vY)/sum(vE)
v=sum((vY-m*vE)^2)/sum(vE)
cat("average=", m, "variance=",v, "phi=", v/m, "\n")

 #variance=phi*m

#si on considère la région
vX<-data2$Group2
for(i in 1:length(levels(vX))){
  vEi<-vE[vX==levels(vX)[i]]
  vYi<-vY[vX==levels(vX)[i]]
  mi<-sum(vYi)/sum(vEi)
  vi<-sum((vYi-mi*vEi)^2)/sum(vEi)
  cat("average=", mi, "variance=",vi, "phi=", vi/mi, "\n")
}

## sans variable explicative
Y<-data2$Nb1
E<-data2$Exppdays/365
(lambda<-sum(Y)/sum(E))
weighted.mean(Y/E,E)
dpois(0:3, lambda)*100



# découpage des variables

graph_freq_Nb(Category, type=1)
freq_sin(Nb1)


data2$Exppdays=data2$Exppdays/365

#data3=data2[,c(-2,-3,-18,-19,-20,-21)]
data3=data2

data3$CAge = cut(data3$Age, breaks = c(18, 23, 29, 33, 47,59,75), include.lowest = TRUE)
#data3$CAge = cut(data3$Age, breaks = c(18, 29, 47,75), include.lowest = TRUE)
#data3$CAge=data3$Age
table(data3$CAge)/nrow(data3)
data3$CBonus = cut(data3$Bonus, breaks = c(-50, -36, -16, 14, 64,150), include.lowest = TRUE)
#data3$CBonus=data3$Bonus
table(data3$CBonus)/nrow(data3)
data3$CGroup1 = cut(data3$Group1, breaks = c(1, 9, 15, 20), include.lowest = TRUE)
table(data3$CGroup1)/nrow(data3)
#data3$CPoldur=data3$Poldur
data3$CPoldur = cut(data3$Poldur, breaks = c(0, 0.5, 4,15) , include.lowest = TRUE)
table(data3$CPoldur)/nrow(data3)
data3$CValue = cut(data3$Value, breaks = c(1000, 12999,49996) , include.lowest = TRUE)
table(data3$CValue)/nrow(data3)
#data3$CDensity=data3$Density
data3$CDensity = cut(data3$Density, breaks = c(0,77.99,162.99,237.99,500) , include.lowest = TRUE)
table(data3$CDensity)/nrow(data3)

#data4 = data3[,c(-5,-6,-7,-8,-9,-13)]
#data5 = data4[,c(-1)]
data5=data3


#constitution de deux bases
alea = runif(nrow(data5))

test = subset(data5, (alea<quantile(alea,0.2)))
apprentissage = subset(data5, (alea>=quantile(alea,0.2)))



#premire regression
regp = glm (Nb1 ~offset(log(Exppdays))+CBonus+CAge+CGroup1+Occupation+CPoldur+Type+Group2+CDensity+Category+Adind+CValue, data =apprentissage , family = poisson (link ="log"))
summary(reg1)


freg<-formula(Nb1 ~offset(log(Exppdays))+CBonus+CAge+CGroup1+Occupation+CPoldur+Type+Group2+CDensity+Category+Adind+CValue)


pred1 = predict(reg1,newdata=apprentissage, type="response", se=TRUE)
pred2 = predict(reg1,newdata=test, type="response")


# NB2 variance form (wi=mui+alpha*mui^2)

regnb2<-glm.nb(freg, data=apprentissage)
summary(regnb2)




#negative binomiale regression
reg2<-glm(Nb1 ~offset(log(Exppdays))+CBonus+CAge+CGroup1+Occupation+CPoldur+Type+Group2+CDensity+Category+Adind+CValue, family = negative.binomial(1), data = data5)



#tester la surdispersion
dispersiontest(regp)


Y<-data2$Nb1
E<-data2$Exppdays/365
n=length(levels(X))
meani=rep(0,n)
variancei=rep(0,n)
Yi=rep(0,n)
X=as.factor(data2$Bonus)
for(i in 1:length(levels(X))){
  Ei=E[X==levels(X)[i]]
  Yi=Y[X==levels(X)[i]]
  meani[i]=weighted.mean(Yi/Ei,Ei)
  variancei[i]=sum((Yi-meani[i]*Ei)^2/sum(Ei))
  cat("Bonus", levels(X)[i], "average=", meani[i], "variance=", variancei[i], "\n")
 
}



plot(meani, variancei, cex=sqrt(Ei),col="grey", pch=19, xlab="Empirical average", ylab="Empirical variance")
points(meani, variancei, cex=sqrt(Ei))

#zero-Inflated Models

#seul la valeur influe sur le second paramètre 
reg<-zeroinfl(Nb1 ~offset(log(Exppdays))+CBonus+CAge+CGroup1+Occupation+CPoldur+Type+Group2+CDensity+Category+Adind+CValue
                |Bonus,
              data=data5, dist="poisson",
              link="logit")

summary(reg)

