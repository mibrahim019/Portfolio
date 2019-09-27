##BIOCHAR UNIVARIATE ANALYSIS
ni<-read.table("Trans_uni3_musa18-02-17.txt",header=T)#Includes envt. variables
View(ni)

b<-read.table("bio_musa2017.txt",header=T)# without log_abundance
View(b)
str(b)

if(!require(installr)){install.packages("installr");require(installr)}
updateR()

##Factorise indendent variables
Treat<-as.factor(b$treat)
sit<-as.factor(b$Site)

hist(b$count)
#GLMS
m1<-glm(count~(sit+order+Treat+pH+WHC+WC+SOM)^2,family=quasipoisson(link="log"),data=b)
require(car)
Anova(m1,test="F")
summary(m1)

m2<-update(m1,~.-order:pH)
Anova(m2,test="F")
m3<-update(m2,~.-pH:SOM)
Anova(m3,test="F")
m4<-update(m3,~.-sit:Treat)
Anova(m4,test="F")
m5<-update(m4,~.-sit:WC)
Anova(m5,test="F")
m6<-update(m5,~.-order:Treat)
Anova(m6,test="F")
m7<-update(m6,~.-order:WHC)
Anova(m7,test="F")
m8<-update(m7,~.-order:SOM)
Anova(m8,test="F")
m9<-update(m8,~.-Treat:pH)
Anova(m9,test="F")
m10<-update(m9,~.-Treat:SOM)
Anova(m10,test="F")
m11<-update(m10,~.-sit:SOM)
Anova(m11,test="F")
m12<-update(m11,~.-sit:WHC)
Anova(m12,test="F")
m13<-update(m12,~.-Treat:WC)
Anova(m13,test="F")
m14<-update(m13,~.-Treat:WHC)
Anova(m14,test="F")
m15<-update(m14,~.-WHC:SOM)
Anova(m15)
summary(m15)###P-VALUE 0.001

par(mfrow=c(2,2))
plot(m15)

(9171.5-3797.4)/9171.5#r2=59%

#Negative binomial
require(MASS)
m01<-glm.nb(count~Treat+order+sit+WHC+pH+SOM,maxit=100,data=b)
Anova(m01)
m02<-update(m01,~.-pH)
anova(m02)
m03<-update(m02,~.-WHC)
anova(m03)
anova(m03,m02)
m04<-update(m03,~.-SOM)
anova(m04)
summary(m04)
m05<-update(m04,~.-sit)
anova(m05)
anova(m04,m05)
summary(m05)##Best model(MAM&AIC)
exp(1.06)
hist(residuals(m05))
par(mfrow=c(1,1))

plot(fitted(m04),residuals(m04), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(m04), residuals(m04)))

plot(fitted(m05),residuals(m05), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(m05), residuals(m05)))

plot(residuals(m05)~predict(m05,type="response"))

####plot(b$WHC,residuals(m01))

##GEE model
require(geepack)
M.gee1<-geeglm(count~Treat+order+WHC+pH+SOM+WC,family=poisson,id=sit,corstr = "ar1",data =b)
anova(M.gee1)
summary(M.gee1)

mgee<-update(M.gee1,~.-pH)
anova(mgee)
summary(mgee)

mgee1<-update(mgee,~.-WC)
anova(mgee1)
summary(mgee1)

mgee2<-update(mgee1,~.-SOM)
anova(mgee2)
summary(mgee2)

mgee3<-update(mgee2,~.-WHC)
anova(mgee3)
summary(mgee3)##Best Model due to scatter in "Residuals v Fitted"

mgee4<-geeglm(count ~ Treat*order, family = poisson, data = b, 
       id = sit, corstr = "ar1")
anova(mgee4)
summary(mgee4)

#Residual v fitted for GEEmodel
plot(fitted(mgee3),residuals(mgee3), main= "GEEglm Model",xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(mgee3), residuals(mgee3)))

plot(fitted(mgee4),residuals(mgee4), main= "GEEglm Model",xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(mgee4), residuals(mgee4)))

##MIXED MODEL GLMM
library(nlme)
Mlme1 <- lme(count~Treat+order+Treat:order,random= ~1|sit,method="ML", data = b)
anova(Mlme1)
summary(Mlme1)##Best model because of lower AIC Value
mlme2<-update(Mlme1,~.-Treat:order)
anova(mlme2)
anova(mlme2,Mlme1)

#Residual v fitted GLMM
plot(fitted(Mlme1),residuals(Mlme1),main="GLMM or mixed model",xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(Mlme1), residuals(Mlme1)))
plot(residuals(Mlme1)~predict(Mlme1,type="response"))

#LM and log transformation of the count data due to the variance
newda<-log(b$abundance+1)
hist(newda)
seconddata<-cbind(b,newda)
View(seconddata)

lm1<-lm(newda~(Treat+sit+order+pH+WHC+WC+SOM)^2,data=seconddata)
require(car)
Anova(lm1)
summary(lm1)
#Model diagnostic tools
par(mfrow=c(2,2))
plot(lm1)

#update lm1 reduce Forder:Ftreat
lm2<-update(lm1,~.-Treat:order)
Anova(lm2)

lm3<-update(lm2,~.-order:WHC)
Anova(lm3)

lm4<-update(lm3,~.-Treat:sit)
Anova(lm4)

lm5<-update(lm4,~.-order:WC)
Anova(lm5)
summary(lm5)

lm6<-update(lm5,~.-order:pH)
Anova(lm6)

lm7<-update(lm6,~.-order:SOM)
Anova(lm7)

lm8<-update(lm7,~.-sit:SOM)
Anova(lm8)
summary(lm8)
plot(lm8)

lm9<-update(lm8,~.-Treat:WHC)
Anova(lm9)

lm10<-update(lm9,~.-pH:SOM)
Anova(lm10)

lm11<-update(lm10,~.-sit:WHC)
Anova(lm11)

lm12<-update(lm11,~.-sit:pH )
Anova(lm12)

lm13<-update(lm12,~.-Treat:pH)
Anova(lm13)

lm14<-update(lm13,~.-Treat:WC)
Anova(lm14)
summary(lm14)

lm15<-update(lm14,~.-Treat:SOM)
Anova(lm15)

lm16<-update(lm15,~.-WHC:SOM)#MAM
Anova(lm16)
par(mfrow=c(2,2))

plot(lm16)
summary(lm16)##R2 = 0.5517=55%

par(mfrow=c(1,2))
##Bargraph for my model##
require(sciplot)
bargraph.CI(sit,count,group=order,data=b)
##sites difference based on treatment
bargraph.CI(sit,count,group=Treat,data=b)
##orders based on treatment differece
bargraph.CI(order,count,group=Treat,data=b)
##Orders differences based on site
bargraph.CI(order,count,group=sit,data=b)
##Bargraph model treatment efficincy based on Order
bargraph.CI(Treat,count,group=order,legend=T,data=b)

##Power analysis
require(pwr)
pwr.chisq.test(w = NULL, N = 120, df =(2-1)*(4-1), sig.level = 0.05, power = .80)

pwr.chisq.test(w = 0.3, N = 120, df =(2-1)*(4-1), sig.level =0.001, power =NULL)##achieved power=0.3194425



bargraph.CI(Site,abundance,group=treat,data=ni,legend=T,ylab="individual(m¯²)",xlab="Sites")

##Bargraph model based on Order
bargraph.CI(treat,abundance,group=order,data=ni,legend=T,ylab="individual(m¯²)",xlab="Substrate treatments")

##shapiro.test
shapiro.test(ni$abundance)#significance means it is different from a normal distribution

#Kruskal-Wallis rank sum test (Post-hoc)
kw<-kruskal.test(abun,Ftreat)#non significance means the treatments are identical
kw

#DCA BIOCHAR
biodca<-decorana(biochar)
plot(biodca)
ordiplot(biodca)
ordirgl(biodca)

ordiplot(biodca, display = 'sp', type = 'n')
orditorp(biodca, display = 'si')

biodca
x<-cbind(biochar$pod_lo,biochar$symphy_lo,biochar$ent_lo)
cor(x)

#PCA BIOCHAR best
biopca<-rda(biochar)
plot(biopca)
biplot(biopca)
biopca

ordiplot(biopca, type = 'n')
points(biopca,display='sites')


bix<-princomp(biochar,scores=TRUE,cor=TRUE)
summary(bix)
plot(bix)

screeplot(biopca,type="l")
screeplot(biopca)
biplot(biopca)

#ENVIRONMENTAL Variables (Soil parameters)
envbio<-read.table("soil_para_env.txt",header=T)
View(envbio)

bcharmod<-rda(biochar~pH+WHC+WC+SOM,envbio)
bcharmod
ordiplot(bcharmod)
plot(bcharmod)

anova(bcharmod)
anova(bcharmod,by="term",step=2000)
anova(bcharmod,by="axis",step=2000)

##Correlation test
cor.test(envbio$SOM,envbio$WHC)
plot(envbio$SOM,envbio$WHC)


plot(bcharmod,xlim=c(-1,1.5),ylim=c(-1,1.5),type="n")
points(bcharmod,display="bp")
text(bcharmod,display="spec",col="blue",cex=0.9)
text(bcharmod,display="bp",col="red")

plot(bcharmod,xlim=c(-1,1.5),ylim=c(-6,4),type="n")
points(bcharmod,display="bp")
text(bcharmod,display="sites",col="black",cex=0.9)
text(bcharmod,display="bp",col="red")


#NMDS:Load Package: MASS!Not important
bcnm<-vegdist(biochar)
bcnmds<-isoMDS(bcnm)
ordiplot(bcnmds,type="t")

#single linkage "jaccard"
bclass<-vegdist(biochar,method="jacc")
bclus<-hclust(bclass,"complete")
bclus
plot(bclus)
rect.hclust(bclus,3)


#Overlay of cluster analysis and ordination
grp<-cutree(bclus,3)
grp
bchargrp<-rda(biochar)
plot(bchargrp,display="sites")
ordihull(bchargrp,grp,lty=2,col="red")

