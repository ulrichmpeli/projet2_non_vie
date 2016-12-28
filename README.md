
rm(list=ls())
setwd("C:/Users/Utilisateur/Documents/ENSAE 2016_2017/actuariat assurance non vie/DM2")
load("base_ensae_1.RData")

library(splines)
library(lsr)
library(MASS)
library(ROCR)
library(pROC)
library(car)
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




#===================== Surdispersion=============
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







#==== Stat desc
par(c(1,1))
barplot(table(data2$Nb1)/nrow(data2),ylab = "Proportion", xlab = "Nombre de sinistres", col = "orange", main = "Proportion du nombre de sinitres RC matériel" )
mean(data2$Nb1)

# == Traitement préliminaire des variables ===========================================

### Détection des effets non linéaires

## Age

reg1 = glm (Nb1 ~offset(log(Exppdays))+bs(Age,df=7) , data =data2 ,
            family = poisson (link ="log"))
summary(reg1)

reg2 = glm (Nb1 ~offset(log(Exppdays))+Age , data =data2 ,
            family = poisson (link ="log"))
summary(reg2)

reg3 = glm (Nb1 ~offset(log(Exppdays))+as.factor(Age) , data =data2 ,
            family = poisson (link ="log"))
summary(reg3)

u= seq (18 ,75 , by =1)
v=rep(365,length(u))
newd = data.frame(Age =u,Exppdays=v)

#en linéaire
y= predict(reg2 , newdata =newd , type ="response",se.fit = TRUE )

#discret
z= predict(reg3 , newdata =newd , type ="response",se.fit = TRUE )

#en spline
w=predict(reg1 , newdata =newd , type ="response",se.fit = TRUE )

#on trace pour montrer la non linéarité
plot (u,z$fit , col ="red",ylab = "Fréquence", xlab = "Age", 
      main="Non linéarité de l'âge")
lines(u,y$fit, col="blue")
lines(u,w$fit, col="darkgreen")
legend( 40,0.35,legend=c("Catégorielle","Linéaire","Splines linéaires M=7",
                         "intervalle confiance à 95% splines " ),
        col=c("red", "blue","darkgreen","yellowgreen"),lty=1:1, cex=0.8)
polygon(c(u,rev(u)),c(w$fit+2*w$se.fit,rev(w$fit-2*w$se.fit)), col='yellowgreen', border=NA)



#un échantillon aléatoire
base1=data2[sample(nrow(data2),80000),]


#les 20 000 restant
base2<-data2[!( data2$PolNum %in%  base1$PolNum),]

deg=seq(1,10,by=1)
vec_var_app=rep(0,10)
vec_var_test=rep(0,10)

#nombre de degré de liberté

for(i in 1:10){
  reg_app<- glm (Nb1 ~offset(log(Exppdays))+bs(Age,df=deg[i]) , data =base1 ,
                 family = poisson (link ="log"))
  vec_var_app [i]=var(residuals(reg_app,type = "response"))
  reg_test<-glm (Nb1 ~offset(log(Exppdays))+bs(Age,df=deg[i]) , data =base2,
                 family = poisson (link ="log"))
  vec_var_test[i]=var(residuals(reg_test,type="response"))
}


plot(vec_var_app,ylab="Variance des résidus", xlab = "Degré du polynôme",
     col="blue",main="Variance des résidus en utilisant l'échantillon d'apprentissage")

lines(vec_var_app,ylab="Variance des résidus", xlab = "Degré du polynôme",
      col="blue",main="Variance des résidus en utilisant l'échantillon d'apprentissage")


plot(vec_var_test,ylab="Variance des résidus", xlab = "Degré du polynôme",
     col="blue",main="Variance des résidus en utilisant l'échantillon test")

lines(vec_var_test,ylab="Variance des résidus", xlab = "Degré du polynôme",
      col="blue",main="Variance des résidus en utilisant l'échantillon test")


##############Value



reg5 <- glm (Nb1 ~Value+offset(log(Exppdays)) , data =data2 ,
             family = poisson (link ="log"))
summary(reg5)


####???n découpe la valeur tous les 1000
max(data2$Value)
min(data2$Value)
point=seq(1000,50000,by=1000)
X=cut(data2$Value,breaks = point )

reg3 <- glm (Nb1 ~offset(log(Exppdays))+ cut(Value,breaks = point ) ,
             data =data2 , family = poisson (link ="log"))
summary(reg3)

u= seq (1000 ,50000 , by =1000)
v=rep(365,length(u))
newd = data.frame(Value =u, Exppdays=v)
y= predict(reg3 , newdata =newd , type ="response",se.fit = TRUE )

z= predict(reg5 , newdata =newd , type ="response",se.fit = TRUE )

plot (u,y$fit , col =" red ", xlab = "Valeur du véhcilue", 
      ylab = "Fréquence", main="Linéarité de Value")
lines(u,z$fit,col="blue")
legend( 10,0.2,legend=c("facteur","continue"),
        col=c("red", "blue"),lty=1:2, cex=0.8)



#============== Bonus =======================# 

##### pour le bonus
for(i in 1:10){
  reg_app<- glm (Nb1 ~offset(log(Exppdays))+bs(Bonus,df=deg[i]) ,
                 data =base1 , family = poisson (link ="log"))
  vec_var_app [i]=var(residuals(reg_app,type = "response"))
  reg_test<-glm (Nb1 ~offset(log(Exppdays))+bs(Bonus,df=deg[i]) ,
                 data =base2, family = poisson (link ="log"))
  vec_var_test[i]=var(residuals(reg_test,type="response"))
}
plot(vec_var_app,ylab="Variance des résidus", xlab = "Degré du polynôme", 
     col="blue",main="Variance des résidus en utilisant l'échantillon d'apprentissage")
lines(vec_var_app,ylab="Variance des résidus", xlab = "Degré du polynôme",
      col="blue",main="Variance des résidus en utilisant l'échantillon d'apprentissage")


plot(vec_var_test,ylab="Variance des résidus", xlab = "Degré du polynôme", 
     col="blue",main="Variance des résidus en utilisant l'échantillon test")
lines(vec_var_test,ylab="Variance des résidus", xlab = "Degré du polynôme", 
      col="blue",main="Variance des résidus en utilisant l'échantillon test")


reg3 <- glm (Nb1 ~factor(Bonus)+ offset (log (Exppdays)) , data =data2 , 
             family = poisson (link ="log"))
summary(reg3)


reg4 <- glm (Nb1 ~Bonus+ offset (log (Exppdays)) , data =data2 , 
             family = poisson (link ="log"))
summary(reg4)


reg6 <- glm (Nb1 ~bs(Bonus,df=6)+ offset (log (Exppdays)) , data =data2 ,
             family = poisson (link ="log"))
summary(reg6)


u= seq(-50,150,by=10)
v=rep(365,21)

newd = data.frame(Bonus =u,Exppdays=v)


y= predict(reg4 , newdata =newd , type ="response",se.fit = TRUE )
z= predict(reg3 , newdata =newd , type ="response",se.fit = TRUE )
w= predict(reg6 , newdata =newd , type ="response",se.fit = TRUE )

plot (u,z$fit , col ="red",ylab = "Fréquence", xlab = "Bonus",
      main="Non linéarité de Bonus")
lines(u,y$fit, col="blue")
lines(u,w$fit, col="darkgreen")


legend( -50,0.6,legend=c("facteur","continue","splines M=6",
                          "intervalle confiance 95% splines" ),
        col=c("red", "blue","darkgreen", "yellowgreen"),lty=1:1, cex=0.8)
polygon(c(u,rev(u)),c(w$fit+2*w$se.fit,rev(w$fit-2*w$se.fit)),
        col='yellowgreen', border=NA)




#on s'intéresse à Poldur
reg3 <- glm (Nb1 ~Poldur+offset(log(Exppdays)) , data =data2 ,
             family = poisson (link ="log"))
summary(reg3)

reg4 <- glm (Nb1 ~factor(Poldur)+offset(log(Exppdays)) , data =data2 ,
             family = poisson (link ="log"))
summary(reg4)

u= seq (0 ,15 , by =1)
v=rep(365,16)

newd = data.frame(Poldur =u, Exppdays=v)

z= predict(reg4 , newdata =newd , type ="response",se.fit = TRUE )
y= predict(reg3 , newdata =newd , type ="response",se.fit = TRUE )


plot (u,z$fit , col ="red",xlab = "Ancienneté du contrat", ylab = "Fréquence", main="Linéarité de Poldur")
lines(u,y$fit,col="blue")
legend( 10,0.2,legend=c("facteur","continue"),col=c("red", "blue"),
        lty=1:1, cex=0.8)


#effet linéaire OK

#on s'intéresse à la densité
reg3 <- glm (Nb1 ~Density +offset(log(Exppdays)), data =data2 ,
             family = poisson (link ="log"))
summary(reg3)


point=seq(10,300,by=10)

reg4 <- glm (Nb1 ~offset(log(Exppdays))+cut(Density,breaks = point ) ,
             data =data2 ,  family = poisson (link ="log"))
summary(reg4)

u= seq(10,300,by=10)
v=rep(356,30)
newd = data.frame(Density =u, Exppdays=v)
y= predict(reg4 , newdata =newd , type ="response",se.fit = TRUE )
plot (u,y$fit , col =" red ", xlab="Densité de population", ylab="Fréquence",
      main="Linéarité de Density" )

z= predict(reg3 , newdata =newd , type ="response",se.fit = TRUE)
lines(u,z$fit,col="blue")
legend( 50,0.5,legend=c("facteur","continue"),col=c("red", "blue"),
        lty=1:1, cex=0.8)



#== Transformation des variables continues

par(mfrow=c(1,1))
arbre=rpart(Nb1~Age, data=data2,method="poisson",cp=4e-4, minsplit=15000) 
prp(arbre,type=2,extra=1 , box.col=5, main="Age")

arbre=rpart(Nb1~Group1, data=data2,method="poisson",cp=4e-4, minsplit=50000) 
prp(arbre,type=2,extra=1 , box.col=5, main="Group1")

arbre=rpart(Nb1~Bonus, data=data2,method="poisson",cp=4e-4, minsplit=30000) 
prp(arbre,type=2,extra=1, box.col=5, main="Bonus" )

arbre=rpart(Nb1~Poldur, data=data2,method="poisson",cp=4e-4, minsplit=15000) 
prp(arbre,type=2,extra=1, box.col=5, main="Poldur" )

arbre=rpart(Nb1~Value, data=data2,method="poisson",cp=4e-4, minsplit=15000) 
prp(arbre,type=2,extra=1, box.col=5, main="Value" )

arbre=rpart(Nb1~Density, data=data2,method="poisson",cp=4e-4, minsplit=15000) 
prp(arbre,type=2,extra=1, box.col=5, main="Density" )





#=== découpage des variables




data2$Exppdays=data2$Exppdays/365

#data3=data2[,c(-2,-3,-18,-19,-20,-21)]
data3=data2

data3$CAge = cut(data3$Age, breaks = c(18, 21, 25, 29, 33,47,59,100), include.lowest = TRUE)

table(data3$CAge)/nrow(data3)
data3$CBonus = cut(data3$Bonus, breaks = c(-50, -36, -16, 14, 54,150), include.lowest = TRUE)

table(data3$CBonus)/nrow(data3)
data3$CGroup1 = cut(data3$Group1, breaks = c(1, 9, 15, 20), include.lowest = TRUE)
table(data3$CGroup1)/nrow(data3)

data3$CPoldur = cut(data3$Poldur, breaks = c(0, 0.5, 4,9,15) , include.lowest = TRUE)
table(data3$CPoldur)/nrow(data3)
data3$CValue = cut(data3$Value, breaks = c(1000, 14999,49996) , include.lowest = TRUE)
table(data3$CValue)/nrow(data3)

data3$CDensity = cut(data3$Density, breaks = c(0,83.99,162.99,237.99,500) , include.lowest = TRUE)
table(data3$CDensity)/nrow(data3)

#data4 = data3[,c(-5,-6,-7,-8,-9,-13)]
#data5 = data4[,c(-1)]
data5=data3



#======= Regroupement des modalités des variables

graph_freq_Nb(Group2, type=1)
freq_sin(Nb1)


#Type
reg =  glm (Nb1 ~offset(log(Exppdays))+Type, data = data2 ,
            family = poisson )
summary(reg)
data5$Type=relevel(data5$Type,"B")
reg = glm (Nb1 ~offset(log(Exppdays))+Type, data = data5 ,
           family = poisson (link ="log"))
summary(reg)
data5$Type=relevel(data5$Type,"C")
reg = glm (Nb1 ~offset(log(Exppdays))+Type, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Type=relevel(data5$Type,"D")
reg = glm (Nb1 ~offset(log(Exppdays))+Type, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Type=relevel(data5$Type,"E")
reg = glm (Nb1 ~offset(log(Exppdays))+Type, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Type=relevel(data5$Type,"F")
reg = glm (Nb1 ~offset(log(Exppdays))+Type, data = data5 , family = poisson (link ="log"))
summary(reg) # On peut regrouper les modalités E et F


linearHypothesis ( reg,c("TypeE = TypeF "))

# Occupation

reg = glm (Nb1 ~offset(log(Exppdays))+Occupation, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Occupation=relevel(data5$Occupation,"Housewife")
reg = glm (Nb1 ~offset(log(Exppdays))+Occupation, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Occupation=relevel(data5$Occupation,"Retired")
reg = glm (Nb1 ~offset(log(Exppdays))+Occupation, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Occupation=relevel(data5$Occupation,"Self-employed")
reg = glm (Nb1 ~offset(log(Exppdays))+Occupation, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Occupation=relevel(data5$Occupation,"Unemployed")
reg = glm (Nb1 ~offset(log(Exppdays))+Occupation, data = data5 , family = poisson (link ="log"))
summary(reg)

# Group2

reg = glm (Nb1 ~offset(log(log(Exppdays)))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Group2=relevel(data5$Group2,"M")
reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Group2=relevel(data5$Group2,"N")
reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Group2=relevel(data5$Group2,"O")
reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Group2=relevel(data5$Group2,"P")
reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Group2=relevel(data5$Group2,"Q")
reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Group2=relevel(data5$Group2,"R")
reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Group2=relevel(data5$Group2,"S")
reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Group2=relevel(data5$Group2,"T")
reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)
data5$Group2=relevel(data5$Group2,"U")
reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)


data5$Group2=relevel(data5$Group2,"U")
reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
summary(reg)



linearHypothesis (reg ,c( "Group2O=Group2P", "Group2O=Group2T","Group2O=Group2L",
                          
                           "Group2Q=Group2N"))



linearHypothesis (reg ,c( "Group2T=Group2P", "Group2O=Group2T", "Group2T=Group2L",
                           
                           "Group2N=Group2Q"))









#constitution de deux bases
alea = runif(nrow(data5))

test = subset(data5, (alea<quantile(alea,0.2)))
apprentissage = subset(data5, (alea>=quantile(alea,0.2)))



#premire regression
regp = glm (Nb1 ~offset(log(Exppdays))+CBonus+CAge+CGroup1+Occupation+CPoldur
            +Type+Group2+CDensity+Category+Adind+CValue, data =apprentissage
            , family = poisson (link ="log"))
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

