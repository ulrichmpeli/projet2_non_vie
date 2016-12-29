#=====On importe les base

rm(list=ls())
setwd("C:/Users/Utilisateur/Documents/ENSAE 2016_2017/actuariat assurance non vie/DM2")
load("base_ensae_1.RData")

# Les packages

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
require(boot)


# == Recherche et suppression des doublons ===========================================
data1 = base_ensae_1
descr = summary(data1)
table(duplicated(data1))
x = duplicated(data1)
data2 = cbind(data1,x)
data2 = data2[data2$x=="FALSE",]
table(duplicated(data2))
summary(data2)
data2$Exppdays=data2$Exppdays/365 #on ramène l'exposition entre 0 et 1




#===================== Surdispersion=============
#### ON CALCUL LA MOYENNE annuelle de la fréquence et sa variance empirique

vY<-data2$Nb1
vE<-data2$Exppdays
m<-sum(vY)/sum(vE)
v=sum((vY-m*vE)^2)/sum(vE)
cat("average=", m, "variance=",v, "phi=", v/m, "\n")

 #variance=phi*m
#phi >1 donc surdispersio,

#phi avec Group2
vX<-data2$Group2
for(i in 1:length(levels(vX))){
  vEi<-vE[vX==levels(vX)[i]]
  vYi<-vY[vX==levels(vX)[i]]
  mi<-sum(vYi)/sum(vEi)
  vi<-sum((vYi-mi*vEi)^2)/sum(vEi)
  cat("average=", mi, "variance=",vi, "phi=", vi/mi, "\n")
}

## lambda sans variable explicative
Y<-data2$Nb1
E<-data2$Exppdays
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

reg1 = glm (Nb1 ~offset(log(Exppdays))+bs(Age,df=5) , data =data2 ,
            family = poisson (link ="log"))
summary(reg1)

reg2 = glm (Nb1 ~offset(log(Exppdays))+Age , data =data2 ,
            family = poisson (link ="log"))
summary(reg2)

reg3 = glm (Nb1 ~offset(log(Exppdays))+as.factor(Age) , data =data2 ,
            family = poisson (link ="log"))
summary(reg3)

u= seq (18 ,75 , by =1)
v=rep(1,length(u))
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
legend( 40,0.35,legend=c("Catégorielle","Linéaire","Splines linéaires M=5",
                         "intervalle confiance à 95% splines " ),
        col=c("red", "blue","darkgreen","yellowgreen"),lty=1:1, cex=0.8)
polygon(c(u,rev(u)),c(w$fit+2*w$se.fit,rev(w$fit-2*w$se.fit)), col='yellowgreen', border=NA)


###Recherche du nombre de degr à utiliser pour les splines
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


#========Value


reg5 <- glm (Nb1 ~Value+offset(log(Exppdays)) , data =data2 ,
             family = poisson (link ="log"))
summary(reg5)


####on découpe la valeur tous les 1000
max(data2$Value)
min(data2$Value)
point=seq(1000,50000,by=1000)
X=cut(data2$Value,breaks = point )

reg3 <- glm (Nb1 ~offset(log(Exppdays))+ cut(Value,breaks = point ) ,
             data =data2 , family = poisson (link ="log"))
summary(reg3)

u= seq (1000 ,50000 , by =1000)
v=rep(1,length(u))
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
v=rep(1,21)

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




#======== à Poldur
reg3 <- glm (Nb1 ~Poldur+offset(log(Exppdays)) , data =data2 ,
             family = poisson (link ="log"))
summary(reg3)

reg4 <- glm (Nb1 ~factor(Poldur)+offset(log(Exppdays)) , data =data2 ,
             family = poisson (link ="log"))
summary(reg4)

u= seq (0 ,15 , by =1)
v=rep(1,16)

newd = data.frame(Poldur =u, Exppdays=v)

z= predict(reg4 , newdata =newd , type ="response",se.fit = TRUE )
y= predict(reg3 , newdata =newd , type ="response",se.fit = TRUE )


plot (u,z$fit , col ="red",xlab = "Ancienneté du contrat", ylab = "Fréquence", main="Linéarité de Poldur")
lines(u,y$fit,col="blue")
legend( 10,0.2,legend=c("facteur","continue"),col=c("red", "blue"),
        lty=1:1, cex=0.8)


#effet linéaire OK

#=============Density
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



#== Transformation des variables continues en facteurs

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


#fonction pour voir les différences entres modalités

graph_freq_Nb=function(var,type){
  X=data10$Group2
  E=data10$Exppdays
  Y=data10$Nb1
  FREQ=levels(X)
  moyenne=variance=n=rep(NA,length(FREQ))
  for(k in 1:length(FREQ)){
    moyenne[k] =weighted.mean(Y[X==FREQ[k]]/E[X==FREQ[k]],
                              E[X==FREQ[k]])
    variance[k]=weighted.mean((Y[X==FREQ[k]]/E[X==FREQ[k]]-
                                 moyenne[k])^2,E[X==FREQ[k]])
    n[k]=sum(E[X==FREQ[k]])
  }
  w=barplot(n,names.arg=FREQ,col="yellow",axes=FALSE)
  mid=w[,1]
  axis(4,)
  par(new=TRUE)
  IC1=moyenne+1.96/sqrt(n)*sqrt(variance)
  IC2=moyenne-1.96/sqrt(n)*sqrt(variance)
  moyenneglobale=sum(Y)/sum(E)
  
  if(type==1){
    plot(mid,moyenne,ylim=c(min(c(IC1,IC2)-diff(range(c(IC1,IC2)))/4),max(c(IC1,IC2))),type="l",col="#FF00FF",axes=FALSE,xlab="",ylab="",xlim=c(0,1.2*length(FREQ)+.5))
    for(i in (-10):20) segments(min(mid)-.8,i/40,max(mid)+.8,i/40,lty=2,col="grey")
    segments(mid,IC1,mid,IC2,col="darkorchid4")
    segments(mid-.1,IC1,mid+.1,IC1,col="darkorchid4")
    segments(mid-.1,IC2,mid+.1,IC2,col="darkorchid4")
    points(mid,moyenne,pch=15,col="#FF00FF")
    axis(2,at=seq(0,.3,by=.05))
    abline(h=moyenneglobale,lty=2,col="darkorchid3")}
  
  if(type==2){
    cic=c(100*(IC1/moyenneglobale),100*(IC2/moyenneglobale))
    YL=c(min(cic)-diff(range(cic))/4,max(cic))
    plot(mid,(moyenne/moyenneglobale)*100,ylim=YL,
         type="l",col="#FF00FF",axes=FALSE,xlab="",ylab="",
         xlim=c(0,1.2*length(FREQ)+.5))
    for(i in (-10):20) segments(min(mid)-.8,i*25,max(mid)+.8,i*25,lty=2,col="grey")
    segments(mid,100*(IC1/moyenneglobale),mid,
             (IC2/moyenneglobale)*100,col="darkorchid4")
    segments(mid-.1,100*(IC1/moyenneglobale),mid+.1,
             (IC1/moyenneglobale)*100,col="darkorchid4")
    segments(mid-.1,100*(IC2/moyenneglobale),mid+.1,
             (IC2/moyenneglobale)*100,col="darkorchid4")
    points(mid,100*(moyenne/moyenneglobale),pch=15,col="#FF00FF")
    axis(2,at=seq(0,300,by=50))
    abline(h=100,lty=2,col="darkorchid3")}
  
  mtext("Exposure", 4, line=-2, cex=.8)
  if(type==1){mtext("Annualized Frequency", 
                    2, line=+2, cex=.8)}
  if(type==2){mtext("Annualized Frequency (multiplier, %)", 
                    2, line=+2, cex=.8)}
}

### Graphique des modalités
graph_freq_Nb(Group2, type=1)




#fonction pour faire des tableaux de contigence
freq_sin = function(nom="Nb1"){
  N = data2[,"Nb1"]
  E= data2$Exppdays
  X1=data2$Group2
  X2=cut(data2$Value,c(0,1000,10000,50000),right=FALSE)
  N_polices = table(X1,X2)
  E_agg=aggregate(E, by = list(X1 = X1, X2 = X2), sum)
  N_exposition=N_polices
  N_exposition[1:nrow(N_exposition),1:ncol(N_exposition)]=
    matrix(E_agg$x,nrow(N_exposition),ncol(N_exposition))
  N_agg=aggregate(N, by = list(X1 = X1, X2 = X2), sum)
  N_sinistres=N_polices
  N_sinistres[1:nrow(N_sinistres),1:ncol(N_sinistres)]=
    matrix(N_agg$x,nrow(N_sinistres),ncol(N_sinistres))
  N_sinistres
  Freq_sinistres = N_sinistres/N_exposition
  return(list(Exp=N_exposition,Sin=N_sinistres,Freq=Freq_sinistres))}


#tableau de contingence Freq / Group2
freq_sin(Nb1)

#=== On veut repérer les modalités non significatives

###Type
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

reg = glm (Nb1 ~offset(log(Exppdays))+Group2, data = data5 , family = poisson (link ="log"))
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




# test de ficher
linearHypothesis (reg ,c( "Group2O=Group2P","Group2O=Group2L",
                           "Group2S=Group2T",
                          "Group2Q=Group2N"))

# On regroupe les modalités
data10 = data5

data10$Group21[data10$Group2 == "L"] = "LOP"
data10$Group21[data10$Group2 == "M"] = "M"
data10$Group21[data10$Group2 == "N"] = "NQ"
data10$Group21[data10$Group2 == "O"] = "LOP"
data10$Group21[data10$Group2 == "P"] = "LOP"
data10$Group21[data10$Group2 == "Q"] = "NQ"
data10$Group21[data10$Group2 == "R"] = "R"
data10$Group21[data10$Group2 == "S"] = "ST"
data10$Group21[data10$Group2 == "T"] = "ST"
data10$Group21[data10$Group2 == "U"] = "U"
data10$Group21 = as.factor(data10$Group21)



reg = glm (Nb1 ~offset(log(Exppdays))+Group21, data = data10 , family = poisson (link ="log"))
summary(reg)


#on regarde si c'est bien significativement différent
graph_freq_Nb(Group21,type = 1)

data10$Group2=data10$Group21


#on fusionne modalité pour type
levels(data10$Type)=c("A", "B", "C", "D", "EF", "EF")

graph_freq_Nb(Type,type = 1)

data10$Type=relevel(data10$Type,"C")
reg = glm (Nb1 ~offset(log(Exppdays))+Type, data = data10 , family = poisson (link ="log"))
summary(reg)

# A n'est pas significatif par rapport à B
linearHypothesis(reg,c("TypeA=TypeB"))

#on fusionne à nouveau les modalités
levels(data10$Type)=c("AB", "AB", "C", "D", "EF")
reg = glm (Nb1 ~offset(log(Exppdays))+Type, data = data10 , family = poisson (link ="log"))
summary(reg)

graph_freq_Nb(Type,type = 1)


#==========

#==Effets Croisés
lst_p= levels (data10$Group2)
lst_m= levels ( data10$Type )
REG1 = glm (Nb1~ Type + Group2 + offset (log(Exppdays)), data =data10 , family =
    poisson )
nd= data.frame (Group2 = rep ( lst_p, length ( lst_m)),Type = rep ( lst_m, each = length ( lst_p)),
 Exppdays =1)
y1= predict(REG1 , newdata =nd , type ="response")
my1 = matrix ( y1,length (lst_p),length ( lst_m))


persp(my1, theta = 135, phi = 30, col = "green3", scale = TRUE,
                ltheta = -120, shade = 0.75, border = TRUE, box = TRUE, ticktype = "detailed",xlab = "X", ylab = "Y", zlab = "Sinc( r )")





#constitution de deux bases
alea = runif(nrow(data10))

test = subset(data10, (alea<quantile(alea,0.2)))
apprentissage = subset(data10, (alea>=quantile(alea,0.2)))



#premire regression
regp = glm (Nb1 ~offset(log(Exppdays))+CBonus+CAge+CGroup1+Occupation+CPoldur
            +Type+Group2+CDensity+Category+Adind+CValue, data =apprentissage
            , family = poisson (link ="log"))
summary(regp)

BIC(regp)


#on prend le critère BIC pour savoir quels variables conserver
selection = stepAIC(regp, k=log(nrow(apprentissage)))


#latex
stargazer(selection)

#prediction
pred1 = predict(selection,newdata=apprentissage, type="response")
pred2 = predict(selection,newdata=test, type="response")

#Mean square error
MSE<-sum((test$Nb1-pred2)^2)/nrow(test)
MSE

##On fait du quasipoisson pour voir le paramètre phi
regqp = glm (Nb1 ~offset(log(Exppdays))+CBonus+CAge+CGroup1+Occupation+CPoldur
            +Type+Group2+CDensity, data =apprentissage
            , family = quasipoisson (link ="log"))
summary(regqp)

#phi proche de 1 

##===== test de surdispersion
dispersiontest(regp)




#freg<-formula(Nb1 ~offset(log(Exppdays))+CBonus+CAge+CGroup1+Occupation+CPoldur+Type+Group2+CDensity+Category+Adind+CValue)




#=== modèle NB2 variance form (wi=mui+alpha*mui^2)

regnb2<-glm.nb(freg, data=apprentissage)
summary(regnb2)


#====negative binomiale regression
reg2<-glm(Nb1 ~offset(log(Exppdays))+CBonus+CAge+CGroup1+Occupation+CPoldur+Type+Group2+CDensity+Category+Adind+CValue, family = negative.binomial(1), data = data5)



#tester la surdispersion
dispersiontest(regp)


#======graphique pour la surdiperion ( ne marche pas)
Y<-data2$Nb1
E<-data2$Exppdays
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

#=============================zero-Inflated Models

#Bonus influe sur le second paramètre 
regzi<-zeroinfl(Nb1 ~offset(log(Exppdays))+CBonus+CAge+CGroup1+Occupation+CPoldur+Type+Group2+CDensity+Category+Adind+CValue
                |bs(Bonus),
              data=data5, dist="poisson",
              link="logit")

summary(regzi)


#===effet de bonus sur la proba de ne pas déclarer
regZIbm <- zeroinfl ( Nb1 ~1 |
                       + as.factor(Bonus),offset = log (Exppdays),
                        data = data10 , dist = "poisson",link ="logit")

regZIbm2 <- zeroinfl ( Nb1 ~1 |
                        Bonus,offset = log (Exppdays),
                      data = data10 , dist = "poisson",link ="logit")

regZIbm3 <- zeroinfl ( Nb1 ~1 |
                         bs(Bonus,df=3),offset = log (Exppdays),
                       data = data10 , dist = "poisson",link ="logit")

u=seq(-50,150, by=10)
v=rep(1,length(u))
B <- data.frame ( Bonus =u, Exppdays=v)
z <- predict ( regZIbm2 , newdata =B, type ="zero")
y <- predict ( regZIbm , newdata =B, type ="zero")
plot (u,y , col ="red", type ="l",lwd =2,
       ylim =c(0 ,1),xlab = "Bonus",ylab =" Probabilite de ne pas declarer un sinistre ")
lines(u,z,col="blue")
lines(u,x, col="darkgreen")

legend( 50,0.8,legend=c("Facteur","Linéaire", "Splines"),col=c("red", "blue", "darkgreen"),lty=1:1, cex=0.8)



#Age
regZIbm <- zeroinfl ( Nb1 ~bs(Age,df=5) |
                        bs(Age,df=5), offset = log (Exppdays),
                      data = data10 , dist = "poisson",link ="logit")


u= seq (18 ,75 , by =1)
v=rep(1,length(u))
B <- data.frame ( Age =u, Exppdays=v)
y <- predict ( regZIbm , newdata =B, type ="zero")
plot (u,y , col ="red", type ="l",lwd =2,
      ylim =c(0 ,1),xlab = "Age",ylab =" Probabilite de ne pas declarer un sinistre ")



#=======Modèle

regzi<-zeroinfl(Nb1 ~CBonus+CAge+CGroup1+Occupation+CPoldur+Type+Group2+CDensity+Category+Adind+CValue
                |bs(Bonus),offset = log (Exppdays),
                data=apprentissage, dist="poisson",
                link="logit")

summary(regzi)

# Procédure BIC
selectionzi = stepAIC(regzi, k=log(nrow(apprentissage)))


stargazer(selectionzi)

# Prévision
pred1 = predict(selectionzi,newdata=apprentissage, type="response")
pred2 = predict(selectionzi,newdata=test, type="response")

#MSE
MSE<-sum((test$Nb1-pred2)^2)/nrow(test)
MSE

## comparaison des deux modèles : test du Vuong
vuong(selectionzi, selection)
