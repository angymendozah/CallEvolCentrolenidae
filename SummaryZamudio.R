devAskNewPage(ask = FALSE)
###Summary analisys between species###

setwd("E:/2021CallsPhD/2021Data/")
install.packages("nortest")#
#library(sjPlot)
library(phytools)
library(caper)#
library(nlme)
library(geiger)#
library(corrplot)
library(picante)#
library(colorRamps)
library(nortest)
library(ggplot2)#
library(factoextra)
library(dplyr)
library(letsR)
library(reshape)
library(nadiv)
library(MCMCglmm)
library(metricTester)
library(PhyloMeasures)
library(lme4)
?MCMCglmm
#Leer Arbol
CentroTree<-read.nexus("outjan2021.tre")
plot(CentroTree, no.margin=T, cex=0.5)

##Leer completo
cantosprom0<-read.table("CantosandExternal.txt", header=T)
names(cantosprom0)[1]<-"Species"
rownames(cantosprom0)<-cantosprom0$Species

####Species diversity
#Range Overlap
atlas<-read.table(file="PAM2019-2.txt", header=T)
atlas<-read.table(file="PAM2019Diverses.txt", header=T)
head(atlas)
colnames(atlas) #Species Names

#SUbclados 
#[13] "Centroleninae"      "Hyalinobatrachinae" "Hyalinobatrachium" 
#[16] "Hyalino17sp"        "CentroNymph"        "Centrolene"        
#[19] "X8genera"           "X5genera"
cantospromAll<-subset(cantosprom0, cantosprom0$X8genera==1) #Centroleninae
dim(cantospromAll)
cantospromAll<-cantosprom0#[,1:14]
#cantospromAll<-cantospromAll[,1:14]

###########

#unir PPCATemp2D a cantospromAll 
PhyPCATemp<-PPCATemp2D[match(rownames(cantospromAll), rownames(PPCATemp2D))]
cantospromAll$PC1Temp<-PhyPCATemp
write.table(cantospromAll, file="cantosPhyPCATemp.txt")
View(cantospromAll)

###############################################################
###################Models of evolution#########################
#######################
aic.w<-function(aic){
  d.aic<-aic-min(aic)
  exp(-1/2*d.aic)/sum(exp(-1/2*d.aic))
}
CantospromModels<-cantospromAll
CantospromModels[3]<-CantospromModels[3]/1000
CantospromModels[4]<-CantospromModels[4]/1000
CantospromModels[13]<-CantospromModels[13]/1000


Fit.Models<-function(dataframe,trait) {
  x<-setNames(dataframe[,trait], rownames(dataframe))
  x <- x[!is.na(x)]
  tree<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, names(x)))
  laks<-phylosig(tree,x,method="lambda",test=TRUE)
  fitBM<-fitContinuous(tree,x,bounds=list(lambda=c(0,1.1)))
  fitLam<-fitContinuous(tree,x,model="lambda", bounds=list(lambda=c(0,1.1)))
  fitK<-fitContinuous(tree,x,model="kappa", bounds=list(a=c(-1, 9999)))
  fitOU<-fitContinuous(tree,x,model="OU", bounds=list(alpha=c(0, 99999999999999)))
  fitWN<-fitContinuous(tree,x,model="white")
  aic.vals<-setNames(c(laks$lambda, laks$P, fitWN$opt$aic, fitBM$opt$aic, fitLam$opt$aic, fitOU$opt$aic, fitK$opt$aic),c("WN","BM", "Lambda","OU", "K"))
  return (aic.vals)
}
?phylosig
?fitContinuous
names(cantospromAll)
Mytable<-matrix(rep(NA,42),6,7)
colnames(Mytable)<-c("lambda","Pval","White","BM","lambda","OU", "Kappa")
rownames(Mytable)<-c("Frec","BandW","Notes","SylD","NRate","PRare")
Mytable[1,]<-Fit.Models(CantospromModels,3) #Peak Freq sin SVL Correction
Mytable[2,]<-Fit.Models(CantospromModels,4)
Mytable[3,]<-Fit.Models(cantospromAll,5)
Mytable[4,]<-Fit.Models(cantospromAll,6)
Mytable[5,]<-Fit.Models(cantospromAll,7)
Mytable[6,]<-Fit.Models(residuales$resid,1)

Mytable[1,]<-Fit.Models(x,13)
View(Mytable)
#Ver parametros particulares
x<-setNames(residuales$resid[,1], rownames(cantospromAll))
x <- x[!is.na(x)]
x<-log(x+1)#Peak frequency and bandwidth
treefit<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, names(x)))
fitContinuous(treefit,x,model="lambda", bounds=list(lambda=c(0,1.1)))
fitContinuous(treefit,x,model="OU", bounds=list(alpha=c(0, 999999)))
fitContinuous(treefit,x,model="delta")
residuales
hist(x)


###Statistical test

lam1<-phylosig(treefit,x,method="lambda")
fitBrownian<-brownie.lite(paintSubTree(treefit,treefit$edge[1,1],state="1"),x)
LR<-2*(lam1$logL-fitBrownian$logL1)
P.lr<-pchisq(LR,df=1,lower.tail=FALSE)
P.lr

##############################################################
#######################PGLS and MCMCglmm######################
ggplot(cantospromAll, aes(y=noterate, x=EVI))+geom_point()
ggplot(data3, aes(y=PeakFreq, x=EVI, size=SVL))+
  geom_point(size=data3$SVL/10, color= "darkcyan") +
  xlab("Enhanced Vegetatio Index (EVI)") +
  ylab("Peak Frequency (Hz)")+
  geom_smooth(method="lm", color='#2C3E50')+
  theme_classic()
hist(data3$noterate)
names(cantospromAll)
#[1] "Species"     "SVL"         "PeakFreq"    "ResidSVL"     "Bandwidht"  
#[6] "MaxNotes"    "SyllDur"     "noterate"    "pulserate"  "PC1Temp"    
#[11] "Genus"       "Diversity"   "EVI"         "Temperature
data3<-na.omit(subset(cantospromAll, 
                      select=c("Species","SVL", "EVI","PeakFreq", "Temperature")))

x<-setNames(data3$SVL, rownames(data3))
y<-setNames(data3$PeakFreq, rownames(data3))
plot(valor, y)
residuales<-phyl.resid(tree =tree3, x, y)
summary(residuales)

data3$PC1Temp<-log(data3$PC1Temp+abs(min(data3$PC1Temp))+1)
data3$noterate<-log(data3$noterate+1)
data3$SyllDur<-log(data3$SyllDur)
data3$MaxNotes<-log(data3$MaxNotes)
tree3<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, row.names(data3)))
tree33<-ladderize(tree3, right = TRUE)

is.ultrametric(tree3)
valor<-setNames(data3$PeakFreq,   rownames(data3))
obj3<-phytools::contMap(tree33, valor, type="fan", plot=FALSE,res=200)
plot(peaktree2, type="fan")
n<-length(obj3$cols)
pal<-c("#FF0000FF", "#FF9900FF", "#CCFF00FF", "#33FF00FF",
       "#00FFFFFF" ,"#3300FFFF")
obj3$cols[1:n]<-colorRampPalette(pal, space="Lab")(n)
obj3$cols[1:n]<-rev(obj3$col[1:n])
plot(obj3, type="fan", lwd=5,fsize=c(0.6,0.8),outline=F)
rainbow(10)
dev.print(pdf,"PeakResid.pdf")


#############
envir<-cor.test(data3$EVI, data3$Temperature)
summary(envir)
#lambda PF SYLL 
comp.data2<-comparative.data(tree3, data3, names.col="Species", vcv.dim=2, warn.dropped=TRUE)
#qplot(EVI,Temperature, data = data3) +geom_point()
m2<-pgls(Temperature~EVI, data=comp.data2, lambda="ML")
lillie.test(m2$residuals) # SVL+EVI OK
summary(m2)
m3<-pgls(Bandwidht~Temperature, data=comp.data2, lambda="ML")
lillie.test(m3$residuals) # SVL+EVI OK
summary(m3)
m4<-pgls(Bandwidht~EVI+Temperature, data=comp.data2, lambda="ML")
lillie.test(m4$residuals) # SVL+EVI OK
summary(m4)

#ResidSVL #Bandwidht #pulserate noterate PC1Temp
data3<-na.omit(subset(cantospromAll, select=c("Species","noterate", "EVI", "Temperature")))
tree3<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, row.names(data3)))

pglsAngy <- gls(SVL~Temperature, correlation = corPagel
                (1, phy = tree3), data = data3, method = "ML")
lillie.test(pglsAngy$residuals) 
as.data.frame(summary(pglsAngy)$tTable)

pglsAngy2 <- gls(pulserate~Temperature, correlation = corMartins(98, phy = tree3), data = data3, method = "ML")
lillie.test(pglsAngy2$residuals) 
as.data.frame(summary(pglsAngy2)$tTable)

pglsAngy3 <- gls(pulserate~EVI+Temperature, correlation = corMartins(1, phy = tree3), data = data3, method = "ML")
lillie.test(pglsAngy3$residuals) 
as.data.frame(summary(pglsAngy3)$tTable)

######MCMCglmm cuando los residules no son normales
library(MCMCglmm)
a<-1000
prior1<-list(R=list(V=1, nu=0.002), 
             G=list(G1=list(V=1, nu=0.002)))
prior2 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000)))
prior3 <- list(R = list(V = diag(1), nu = 0.002, fix = TRUE),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))
prior4 <- list(R = list(V = diag(1), nu = 0.002, fix = TRUE),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))

prior5 <- list(R = list(V = diag(1), nu = 0.002, fix = TRUE),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))

#####################
#SUbclados 
#[13] "Centroleninae"      "Hyalinobatrachinae" "Hyalinobatrachium" 
#[16] "Hyalino17sp"        "CentroNymph"        "Centrolene"        
#[19] "X8genera"           "X5genera"
cantospromAll<-cantosprom0
cantospromAll<-subset(cantosprom0, cantosprom0$Centroleninae==1) #Centroleninae

rownames(cantospromAll)
#[1] "Species"     "SVL"         "PeakFreq"    "ResidSVL"     "Bandwidht"  
#[6] "MaxNotes"    "SyllDur"     "noterate"    "pulserate"  "PC1Temp"    
#[11] "Genus"       "Diversity"   "EVI"         "Temperature
data3<-na.omit(subset(cantospromAll, select=c("Species", "noterate","SVL",
                                              "EVI", "Temperature")))
#data3$pulserate<-log(data3$pulserate+1)

MyMCMC<-function(grupo,rasgo){
  cantospromAll<-subset(cantosprom0, cantosprom0[,grupo]==1)
  data3<-na.omit(subset(cantospromAll, select=c("Species", as.character(rasgo),"SVL",
                                                "EVI", "Temperature")))
  if (rasgo=="puserate"){
    data3$pulserate<-log(data3$pulserate+1)
  } 
  else 
    
    tree3<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, row.names(data3)))
  Ainv4<-inverseA(tree3, scale=F)$Ainv
  modelo1<-as.formula(paste(rasgo, "~EVI"))
  model1<-MCMCglmm(fixed=modelo1, random=~Species+Temperature, data=data3, 
                   ginverse=list(Species=Ainv4), prior = prior3,
                   nitt=100000, thin=10, burnin=300, verbose=FALSE)
  Vegetation<-summary(model1)$solutions[2,]
  
  modelo2<-as.formula(paste(rasgo, "~Temperature"))
  model2<-MCMCglmm(fixed=modelo2, random=~Species+EVI, data=data3, 
                   ginverse=list(Species=Ainv4), prior = prior3,
                   nitt=100000, thin=10, burnin=300, verbose=FALSE)
  Bio8<-summary(model2)$solutions[2,]
  
  unidos<-rbind(Vegetation,Bio8)
  return(unidos)
}

#16	Centroleninae
#17	Hyalinobatrachinae
#18	Hyalinobatrachium
#19	Hyalino17sp
#20	CentroNymph
#21	Centrolene
#22	8genera
#23	5genera

#Centroleninae
MyMCMC(17, "PeakFreq")
MyMCMC(22, "SyllDur")
##################

#SUbclados 
#[13] "Centroleninae"      "Hyalinobatrachinae" "Hyalinobatrachium" 
#[16] "Hyalino17sp"        "CentroNymph"        "Centrolene"        
#[19] "X8genera"           "X5genera"
cantospromAll<-cantosprom0
cantospromAll<-subset(cantosprom0, cantosprom0$X8genera==1) #Centroleninae

#rownames(cantospromAll)
#[1] "Species"     "SVL"         "PeakFreq"    "ResidSVL"     "Bandwidht"  
#[6] "MaxNotes"    "SyllDur"     "noterate"    "pulserate"  "PC1Temp"    
#[11] "Genus"       "Diversity"   "EVI"         "Temperature
data3<-na.omit(subset(cantospromAll, select=c("Species", "SyllDur","pulserate",
                                              "noterate", "MaxNotes")))
#data3$pulserate<-log(data3$pulserate+1)
tree3<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, row.names(data3)))
###

Ainv4<-inverseA(tree3, scale=F)$Ainv
model2.1<-MCMCglmm(fixed=PeakFreq~EVI, random=~Species+Temperature+SVL, data=data3, 
                   ginverse=list(Species=Ainv4), prior = prior4,
                   nitt=100000, thin=10, burnin=300, verbose=FALSE)
summary(model2.1)
model2.2<-MCMCglmm(fixed=PeakFreq~Temperature, random=~Species+EVI, data=data3, 
                   ginverse=list(Species=Ainv4), prior = prior3,
                   nitt=100000, thin=10, burnin=300, verbose=FALSE)
summary(model2.2)

a<-1000
prior1<-list(R=list(V=1, nu=0.002), 
             G=list(G1=list(V=1, nu=0.002)))
prior2 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000)))
prior3 <- list(R = list(V = diag(1), nu = 0.002, fix = TRUE),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))
prior4 <- list(R = list(V = diag(1), nu = 0.002, fix = TRUE),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))

prior5 <- list(R = list(V = diag(1), nu = 0.002, fix = TRUE),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))



model2.1<-MCMCglmm(fixed=noterate~EVI, random=~Species+Temperature, data=data3, 
                   ginverse=list(Species=Ainv4), prior = prior3,
                   nitt=100000, thin=10, burnin=300, verbose=FALSE)



model3.2<-MCMCglmm(fixed=PeakFreq~EVI, random=~Species+Temperature+SVL, data=data3, 
                   ginverse=list(Species=Ainv4), prior = prior4,
                   nitt=500000, thin=10, burnin=300, verbose=FALSE)
summary(model3.2)
model3.3<-MCMCglmm(fixed=PeakFreq~Temperature, random=~Species+EVI+SVL, data=data3, 
                   ginverse=list(Species=Ainv4), prior = prior4,
                   nitt=500000, thin=10, burnin=300, verbose=FALSE)
summary(model3.3)
#in noterate/pulserate~ Mixed model equations singular: use a (stronger) prior

plot(model2.1)

dev.off()
##check traces
traceplot(model2.1$VCV[,2])
#check autocorrelation
autocorr.diag(model2.1$Sol)
autocorr.diag(model2.1$VCV)
# acf plot for the first random term in our model (the animal term)
acf(model2.1$VCV[,1], lag.max = 20)
# Checking convergence for our fixed factors (against an equal second model)
gelman.diag(mcmc.list(model2.1$Sol, model2.2$Sol))

#Median of the fixed factor
median(model2.1$Sol[,1])

################################
#####Community comparitons######
################################
#############################

###Construir archivos para el LMM
####Species diversity
#Range Overlap
atlas<-read.table(file="PAM2019-2.txt", header=T)
atlas<-read.table(file="PAM2019Diverses.txt", header=T)
head(atlas)
colnames(atlas) #Species Names

atlas1<-as.matrix(atlas)
overlap<-lets.overlap(atlas1, method = "Chesser&Zink", xy = F)


#SUbclados 
#[13] "Centroleninae"      "Hyalinobatrachinae" "Hyalinobatrachium" 
#[16] "Hyalino17sp"        "CentroNymph"        "Centrolene"        
#[19] "X8genera"           "X5genera"
############funcion para pares de especies###########
#####################################################
matrixLMM<-function (Dcall, overlap, Tree) {
  matrixcall  <- melt(as.matrix(Dcall))
  matrixSpOverlap <- melt(as.matrix(overlap))
  matrixphylo<-cophenetic(Tree)                    
  
  matrixcall2<-cbind(matrixcall, 
                     c(paste(matrixcall$X1,matrixcall$X2,sep="_")),
                     c(paste(matrixcall$X2,matrixcall$X1,sep="_"))) 
  colnames(matrixcall2)[4]<-"paste"
  colnames(matrixcall2)[5]<-"duplicate"
  colnames(matrixcall2)[3]<-"Call"
  
  matrixSpOverlap2<-cbind(matrixSpOverlap, c(paste(matrixSpOverlap$X1,matrixSpOverlap$X2,sep="_")))
  colnames(matrixSpOverlap2)[4]<-"paste"
  colnames(matrixSpOverlap2)[3]<-"Coexistence"
  
  matrixphylo  <- melt(as.matrix(matrixphylo))
  matrixphylo2<-cbind(matrixphylo, c(paste(matrixphylo$X1,matrixphylo$X2,sep="_")))
  colnames(matrixphylo2)[4]<-"paste"
  colnames(matrixphylo2)[3]<-"phylodist"
  #head(matrixcall2)
  #head(matrixSpOverlap2)
  myMatrix<-inner_join(matrixcall2, matrixSpOverlap2, by="paste")
  myMatrix<-inner_join(myMatrix, matrixphylo2, by="paste")
  
  suppressWarnings(
    for(i in 1:nrow(myMatrix)){
      if(sum(myMatrix$paste %in% myMatrix$duplicate[i])>=1){
        myMatrix$duplicate[myMatrix$duplicate %in% myMatrix$paste[i]]<-"rep"
      }
    }) 
  myMatrix1<-na.omit(myMatrix)
  rownames(myMatrix1)<-myMatrix1$paste
  myMatrix2<-myMatrix1[,c("X1.x", "X2.x", "Call", "Coexistence", "phylodist")]
  colnames(myMatrix2)[1:2]<-c("Species1", "Species2")
  rownames(myMatrix2)<-rownames(myMatrix1)
  myMatrix2<-filter(myMatrix2, Call>0)
  return (myMatrix2)
}
##############################################################

componenst<-c(1,2,3,49)
colnames(cantosprom0)

##############################
coexistencia<-function(grupo){
  if(grupo<16) {
    cantospromAll<-cantosprom0
  } else 
  {cantospromAll<-subset(cantosprom0, cantosprom0[,grupo]==1)}
  
  if(sd(na.omit(cantospromAll$noterate))>0) {
    components<-na.omit(cantospromAll[,c(3:8)])
    TempCall<-na.omit(cantospromAll[,c(5:8)])#Sólo variables temporales
  } else  
  {components<-na.omit(cantospromAll[,c(3,4,6,8)])
  TempCall<-na.omit(cantospromAll[,c(6,8)])#Sólo variables temporales
  }
  
#names(cantospromAll)
species2<-rownames(components)
to.dropW<-setdiff(CentroTree$tip.label, rownames(components))
Tree0<-drop.tip(CentroTree, to.dropW)
ScaleCall<-scale(components)

cantonoPhy2<-princomp(ScaleCall)
cantonoPhy$loadings<-canto$L
cantonoPhy$scores<-canto$S


cantonoPhy[]
fviz_pca_biplot(cantonoPhy, col.ind = "#00AFBB",
             repel = TRUE)

canto<-phyl.pca(Tree0, ScaleCall)
biplot(canto)
library(factoextra)
fviz_pca_biplot(canto)
canto2<-phyl.pca(Tree0, ScaleCall, method="lambda")
call2D<-canto$S[,1:2]
Dcall2D<-dist(call2D)
SpecCall<-na.omit(cantospromAll[,c(3:4)])#Sólo variables spectrales
speciesTemp<-rownames(TempCall) #95 spp
speciesSpec<-rownames(SpecCall) #91
TempCall<-scale(TempCall)
SpecCall<-scale(SpecCall)
to.drop3<-setdiff(CentroTree$tip.label, speciesTemp) # speciesSpec
to.drop4<-setdiff(CentroTree$tip.label, speciesSpec) # speciesSpec
TreeTemp<-drop.tip(CentroTree, to.drop3) #95
TreeSpec<-drop.tip(CentroTree, to.drop4) #56
PPCATemp<-phyl.pca(TreeTemp, TempCall, method="lambda")
PPCASpec<-phyl.pca(TreeSpec, SpecCall, method="lambda")
PPCATemp2D<-PPCATemp$S[,1:2]
PPCASpec2D<-PPCASpec$S[,1:2]
PPCATempDist<-dist(PPCATemp2D)
PPCASpecDist<-dist(PPCASpec2D)

#Depurando la PAM
overlap1<-subset(overlap, rownames(overlap) %in% species2)
overlap2<-subset(t(overlap1), colnames(overlap1) %in% species2)

try4<-matrixLMM(Dcall2D, overlap2, Tree0)

Ainv4<-inverseA(Tree0, scale=F)$Ainv
prior1<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
model2.0<-MCMCglmm(Call~Coexistence, random=~Species1, data=try4, 
                   ginverse=list(Species1=Ainv4), prior = prior1,
                   nitt=100000, thin=10, burnin=300, verbose=FALSE)
todos<-summary(model2.0)$solutions[2,]

##################################
#solo variables temporales/spectrales
overlap1t<-subset(overlap, rownames(overlap) %in% speciesTemp) 
overlap2t<-subset(t(overlap1t), colnames(overlap1t) %in% speciesTemp) #speciesTemp

overlap1s<-subset(overlap, rownames(overlap) %in% speciesSpec) 
overlap2s<-subset(t(overlap1s), colnames(overlap1s) %in% speciesSpec) #speciesTemp

tryTemp<-matrixLMM(PPCATempDist, overlap2t, TreeTemp) #/ (PPCATempDist, overlap2t, Tree4)
trySpec<-matrixLMM(PPCASpecDist, overlap2s, TreeSpec) #/ (PPCATempDist, overlap2t, Tree4)

##con MCMCGLMM Temporals
Ainv4<-inverseA(TreeTemp, scale=F)$Ainv
prior1<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
model2.1<-MCMCglmm(Call~Coexistence, random=~Species1, data=tryTemp, 
                   ginverse=list(Species1=Ainv4), prior = prior1,
                   nitt=100000, thin=10, burnin=300, verbose=FALSE)
temporales<-summary(model2.1)$solutions[2,]

##Spectrals
Ainv4<-inverseA(TreeSpec, scale=F)$Ainv
prior1<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
model2.2<-MCMCglmm(Call~Coexistence, random=~Species1, data=trySpec, 
                   ginverse=list(Species1=Ainv4), prior = prior1,
                   nitt=100000, thin=10, burnin=300, verbose=FALSE)
espectrales<-summary(model2.2)$solutions[2,]

completo<-rbind(todos, espectrales, temporales)
return(completo)
}

grupo15<-coexistencia(15) #Centrolenidae
grupo16<-coexistencia(16) #Centroleninae
grupo17<-coexistencia(17) #Hyalinobatrachinae
grupo20<-coexistencia(20) #CentroNymph
grupo21<-coexistencia(21) #Centrolene
grupo22<-coexistencia(22) #Cochranellini

warnings()


################################
################################
round(sort(effectiveSize(model2.2$Sol)))

songGLMS<-glmer(Call~Coexistence+(1|phylodist)+(1|Species1)+(1|Species2),
               data=trySpec, family=gaussian(link=log))
summary(songGLMS)

######
tryTempshort<-subset(tryTemp, phylodist<0.05)
write.table(try4, file="DistancesCall.txt", sep=";")
write.table(tryTemp, file="DistancesCallTemp.txt", sep=";")
write.table(trySpec, file="DistancesCallSpec.txt", sep=";")
hist(tryTemp$phylodist)
ggplot(data=trySpec, aes(y=Call, x=Coexistence, color=phylodist))+
  geom_point()+scale_color_gradient(low="green", high="red")
View(try4)

head(mtcars)

############################
######## null model#########
############################
library(lme4)
###establecer condiciones de atributos simulados
Mytable<-matrix(rep(NA,14),2,7)
colnames(Mytable)<-c("lambda","Pval","White","BM","lambda","OU", "Kappa")
rownames(Mytable)<-c("PC1","PC2")
Mytable[1,]<-Fit.Models(call2D,1) #/call2D / PPCATemp2D / PPCASpec2D
Mytable[2,]<-Fit.Models(call2D,2) #/call2D / PPCASpec2D / PPCATemp2D
Mytable

#Ver parametros particulares #PPCATemp2D/ call2D / PPCASpec2D
x<-setNames(call2D[,1], rownames(call2D))
treefit<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, names(x)))
fitContinuous(treefit,x,model="lambda") #, bounds=list(alpha=c(0,999)))
fitContinuous(treefit,x,model="lambda")

x<-setNames(call2D[,2], rownames(call2D))
treefit<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, names(x)))
fitContinuous(treefit,x,model="OU", bounds=list(alpha=c(0,999)))
fitContinuous(treefit,x,model="lambda")

###########################
set.seed(42)
?fastBM
rownames(modelsPC1s)
modelsPC1<-fastBM(Tree0, alpha=63.23, nsim=1000, sig2=162.94)
modelsPC2<-fastBM(rescale(Tree0,model="lambda",0.43), nsim=1000, sig2=13.193)

modelsPC1t<-fastBM(TreeTemp, nsim=100, alpha=31.92, sig2=65.38)
modelsPC2t<-fastBM(TreeTemp, nsim=100, sig2=18.98)

modelsPC1s<-fastBM(TreeSpec, nsim=100, alpha=7.49, sig2=28.93)
modelsPC2s<-fastBM(TreeSpec, nsim=100, alpha=26.16, sig2=24.1)

modelsPC1a<-subset(modelsPC1, subset=rownames(modelsPC1) %in% Tree0$tip.label)
modelsPC2a<-subset(modelsPC1, subset=rownames(modelsPC1) %in% Tree0$tip.label)

modelsPC1ta<-subset(modelsPC1t, subset=rownames(modelsPC1t) %in% TreeTemp$tip.label)
modelsPC2ta<-subset(modelsPC1t, subset=rownames(modelsPC1t) %in% TreeTemp$tip.label)

modelsPC1sa<-subset(modelsPC1s, subset=rownames(modelsPC1s) %in% TreeSpec$tip.label)
modelsPC2sa<-subset(modelsPC1s, subset=rownames(modelsPC1s) %in% TreeSpec$tip.label)

################################################
rownames(myMatrix3)
hist(myMatrix3$Call)
myMatrixNull<-try4[,]
myMatrixNull<-tryTemp[,]
myMatrixNull<-trySpec[,]
Estimate<-rep(NA, 100)

suppressWarnings(
  for(i in 1:length(Estimate)){
    myMatrixNull[,3]<-NA
    pairs<-cbind(modelsPC1sa[,i], modelsPC2sa[,i])
    Dcall<-dist(pairs)
    matrixcall<-melt(as.matrix(Dcall))
    matrixcall2<-cbind(matrixcall,
                       c(paste(matrixcall$X1,matrixcall$X2,sep="_")),
                       c(paste(matrixcall$X2,matrixcall$X1,sep="_")))
    colnames(matrixcall2)[4]<-"paste"
    colnames(matrixcall2)[5]<-"duplicate"
    colnames(matrixcall2)[3]<-"Call"
    for(j in 1:nrow(matrixcall2)){
      if(sum(matrixcall2$paste %in% matrixcall2$duplicate[j])>=1){
        matrixcall2$duplicate[matrixcall2$duplicate %in% matrixcall2$paste[j]]<-"rep"
      }
    }
    matrixcall2<-na.omit(matrixcall2)
    rownames(matrixcall2)<-matrixcall2$paste
    myMatrixNull[,3]<-matrixcall2$Call
    ###
    songGLM<-lmer(Call~Coexistence+(1|phylodist)+(1|Species1)+(1|Species2),
                        data=myMatrixNull)#, family=gaussian(link=log))
    EstimateA<-coef(summary(songGLM))[2,1]
    Estimate[i]<-EstimateA
    print(i)
  }
)
summary(Estimate)
hist(Estimate)
View(myMatrixNull)
length(speciesTemp)
length(speciesSpec)
setdiff(TreeSpec$tip.label, rownames(call2D))
plot(Call~Coexistence, data=tryTemp)
setdiff(colnames(atlassong2), TreeSpec$tip.label)
dim(PPCATemp2D)
names(tryTemp)

########################################
#######Community approach###############
########################################


#################################################
###################################################
#species2, Tree0, call2D
Elliottest<-function (Species, PCA){
atlassong<-subset(atlas, select=colnames(atlas) %in% Species) # /speciesTemp /SpeciesSpec
atlassong1<-subset(atlassong, subset=rowSums(atlassong) > 1)
atlassong2<-subset(atlassong1, select=colSums(atlassong1) > 0)

#Functional Dispersion Laliberte and Legendre (2010)
callNew<-subset(PCA, rownames(PCA)%in%colnames(atlassong2))
functionalDisp<-FDis(ordination.results=callNew, road.map=atlassong2)
SpeciesRichness<- as.matrix(rowSums(atlassong2))
RangeSize<- as.matrix(colSums(atlassong2))
atlasfunct<-functionalDisp*atlassong2
#columnas
AcousticField<-colSums(atlasfunct)/RangeSize
AcousticField<-AcousticField[,1]

treeA<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, colnames(atlassong2)))
weights = runif(length(treeA$tip.label))
names(weights) = treeA$tip.label
busqueda<-mpd.query(treeA,atlassong2,TRUE,null.model="sequential",
                    abundance.weights=weights, reps=1000)
atlasNRI<-busqueda*atlassong2
meanNRI<-colSums(atlasNRI)/RangeSize
meanNRI<-meanNRI[,1]

comundiadX<-data.frame(AcousticField, meanNRI)

#ggplot(comundiadX, aes(y=AcousticField, x=meanNRI))+
#  geom_point(color= "darkcyan") +
#  xlab("NRI") +   ylab("Functional disp")+
#  geom_smooth(method="lm", color='#2C3E50')+
#  theme_classic()

#Original
#prepped <- prepFieldData(tree=treeA, picante.cdm=atlassong2)
#results <- calcField(prepped, metrics="NAW_MPD")
#MPD<-setNames(results$NAW_MPD, results$species) 
#species<-names(MPD)

Elliot<-data.frame(AcousticField, meanNRI)
Elliot<-subset(Elliot, subset=Elliot$AcousticField>0)
Elliot<-cbind(Elliot, rownames(Elliot))
colnames(Elliot)[3]<-"Species"
return (Elliot)
}

ComunTest<-function (grupo){
  if(grupo<16) {
    cantospromAll<-cantosprom0
  } else 
  {cantospromAll<-subset(cantosprom0, cantosprom0[,grupo]==1)}
  
  if(sd(na.omit(cantospromAll$noterate))>0) {
    components<-na.omit(cantospromAll[,c(3:8)])
    TempCall<-na.omit(cantospromAll[,c(5:8)])#Sólo variables temporales
  } else  
  {components<-na.omit(cantospromAll[,c(3,4,6,8)])
  TempCall<-na.omit(cantospromAll[,c(6,8)])#Sólo variables temporales
  }

  ambiente<-cantosprom0[,9:10]
  species2<-rownames(components)
  to.dropW<-setdiff(CentroTree$tip.label, rownames(components))
  Tree0<-drop.tip(CentroTree, to.dropW)
  ScaleCall<-scale(components)
  canto<-phyl.pca(Tree0, ScaleCall)
  call2D<-canto$S[,1:2]
  Dcall2D<-dist(call2D)
  SpecCall<-na.omit(cantospromAll[,c(3:4)])#Sólo variables spectrales
  speciesTemp<-rownames(TempCall) #95 spp
  speciesSpec<-rownames(SpecCall) #91
  TempCall<-scale(TempCall)
  SpecCall<-scale(SpecCall)
  to.drop3<-setdiff(CentroTree$tip.label, speciesTemp) # speciesSpec
  to.drop4<-setdiff(CentroTree$tip.label, speciesSpec) # speciesSpec
  TreeTemp<-drop.tip(CentroTree, to.drop3) #95
  TreeSpec<-drop.tip(CentroTree, to.drop4) #56
  PPCATemp<-phyl.pca(TreeTemp, TempCall, method="lambda")
  PPCASpec<-phyl.pca(TreeSpec, SpecCall, method="lambda")
  PPCATemp2D<-PPCATemp$S[,1:2]
  PPCASpec2D<-PPCASpec$S[,1:2]
  PPCATempDist<-dist(PPCATemp2D)
  PPCASpecDist<-dist(PPCASpec2D)
  
  #Depurando la PAM
  overlap1<-subset(overlap, rownames(overlap) %in% species2)
  overlap2<-subset(t(overlap1), colnames(overlap1) %in% species2)  
  
  Elliot<-Elliottest(species2, call2D)
  Elliott<-Elliottest(speciesTemp, PPCATemp2D)
  Elliots<-Elliottest(speciesSpec, PPCASpec2D)
  
  Tree1<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, rownames(Elliot)))
  Treet<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, rownames(Elliott)))
  Trees<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, rownames(Elliots)))
  
  Comunidad0<-merge(Elliot, ambiente, all=F, by="row.names")
  Comunidadt<-merge(Elliott, ambiente, all=F, by="row.names")
  Comunidads<-merge(Elliots, ambiente, all=F, by="row.names")
  
  Ainv4<-inverseA(Tree1, scale=F)$Ainv
  Ainvt<-inverseA(Treet, scale=F)$Ainv
  Ainvs<-inverseA(Trees, scale=F)$Ainv
  
  model0<-MCMCglmm(AcousticField~MPD, data=Comunidad0, random=~Temperature+EVI,
                   ginverse=list(Species=Ainv4), prior = prior3,
                   nitt=1000000, thin=10, burnin=300, verbose=FALSE)
  model.t<-MCMCglmm(AcousticField~MPD, data=Comunidadt, random=~Temperature+EVI,
                    ginverse=list(Species=Ainvt), prior = prior3,
                    nitt=500000, thin=10, burnin=300, verbose=FALSE)
  model.s<-MCMCglmm(AcousticField~MPD, data=Comunidads, random=~Temperature+EVI,
                    ginverse=list(Species=Ainvs), prior = prior3,
                    nitt=500000, thin=10, burnin=300, verbose=FALSE)
  
  todos<-summary(model0)$solutions[2,]
  temporales<-summary(model.t)$solutions[2,]
  espectrales<-summary(model.s)$solutions[2,]
  completo<-rbind(todos, espectrales, temporales)
  
  return(completo)
}

ComunTest15<-ComunTest(15) #Centrolenidae
ComunTest16<-ComunTest(16) #Centroleninae
ComunTest17<-ComunTest(17) #Hyalinobatrachinae
ComunTest20<-ComunTest(20) #CentroNymph
ComunTest21<-ComunTest(21) #Centrolene
ComunTest22<-ComunTest(22) #Cochranellini

 
plot(SpeciesRichness,functionalDisp)
plot(RangeRichness, Varianza_aFNew)
plot(RangeRichness, AcousticFieldNew)
plot(AcousticFieldOld, AcousticFieldNew2)

##############################
#[13] "Centroleninae"      "Hyalinobatrachinae" "Hyalinobatrachium" 
#[16] "Hyalino17sp"        "CentroNymph"        "Centrolene"        
#[19] "X8genera"           "X5genera"
cantospromAll<-cantosprom0
cantospromAll<-subset(cantosprom0, cantosprom0$Hyalino17sp==1)
cantospromAll<-subset(cantospromAll, subset=rownames(cantospromAll) %in% colnames(atlas))
components<-na.omit(cantospromAll[,c(3:8)])
#components<-na.omit(cantospromAll[,c(3,4,6,8)])#Hyalinos
species2<-rownames(components)
to.dropW<-setdiff(CentroTree$tip.label, rownames(components))
Tree0<-drop.tip(CentroTree, to.dropW)
ScaleCall<-scale(components)
canto<-phyl.pca(Tree0, ScaleCall)
call2D<-canto$S[,1:2]
Dcall2D<-dist(call2D)

TempCall<-na.omit(cantospromAll[,c(5:8)])#Sólo variables temporales
#TempCall<-na.omit(cantospromAll[,c(6,8)]) #Hyalinos
SpecCall<-na.omit(cantospromAll[,c(3:4)])#Sólo variables spectrales

speciesTemp<-rownames(TempCall) #95 spp
speciesSpec<-rownames(SpecCall) #91

TempCall<-scale(TempCall)
SpecCall<-scale(SpecCall)

to.drop3<-setdiff(CentroTree$tip.label, speciesTemp) # speciesSpec
to.drop4<-setdiff(CentroTree$tip.label, speciesSpec) # speciesSpec

TreeTemp<-drop.tip(CentroTree, to.drop3) #95
TreeSpec<-drop.tip(CentroTree, to.drop4) #56

PPCATemp<-phyl.pca(TreeTemp, TempCall, method="lambda")
PPCASpec<-phyl.pca(TreeSpec, SpecCall, method="lambda")

PPCATemp2D<-PPCATemp$S[,1:2]
PPCASpec2D<-PPCASpec$S[,1:2]
#############################

#################################################
#SPecies = speciesSpec / speciesTemp /species2
#Tree= Tree0 / TreeSpec / TreeTemp
#PCA= call2D / PPCATemp2D / PPCASpec2D
#dim(ElliotNull)
new.lambda<-0.5

ambiente<-cantosprom0[,9:10]
a<-1000
prior3 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))



Elliot<-Elliottest(species2, call2D)
Elliott<-Elliottest(speciesTemp, PPCATemp2D)
Elliots<-Elliottest(speciesSpec, PPCASpec2D)

dim(Elliot)
plot(Call2D)

Tree1<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, rownames(Elliot)))
Treet<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, rownames(Elliott)))
Trees<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, rownames(Elliots)))

Comunidad0<-merge(Elliot, ambiente, all=F, by="row.names")
Comunidadt<-merge(Elliott, ambiente, all=F, by="row.names")
Comunidads<-merge(Elliots, ambiente, all=F, by="row.names")

Ainv4<-inverseA(Tree1, scale=F)$Ainv
Ainvt<-inverseA(Treet, scale=F)$Ainv
Ainvs<-inverseA(Trees, scale=F)$Ainv

model0<-MCMCglmm(AcousticField~MPD, data=Comunidad0, random=~Temperature+EVI,
                 ginverse=list(Species=Ainv4), prior = prior3,
                 nitt=1000000, thin=10, burnin=300, verbose=FALSE)
summary(model0)
model.t<-MCMCglmm(AcousticField~MPD, data=Comunidadt, random=~Temperature+EVI,
                  ginverse=list(Species=Ainvt), prior = prior3,
                  nitt=500000, thin=10, burnin=300, verbose=FALSE)
summary(model.t)
model.s<-MCMCglmm(AcousticField~MPD, data=Comunidads, random=~Temperature+EVI,
                  ginverse=list(Species=Ainvs), prior = prior3,
                  nitt=500000, thin=10, burnin=300, verbose=FALSE)
summary(model.s)

class(model.s)
summary(model.s)[]
#
###########################
resultados<-read.table(file="mcmcglmmresultsjun.txt", header=T, sep="\t")
names(resultados)
dim(resultados)

resultados2<- resultados[order(resultados$Metric),]
#####
plot(resultados2$Mean[7:1], 1:7, yaxt="n", ylab="", xlim=c(-20,30),
     pch=19, main="Overall Call")
grid()
columnas<-resultados2$Clade[1:7]
axis(2, at=7:1, as.graphicsAnnot(columnas), las=2, lwd=0.5)
arrows(resultados2$under[7:1], 1:7, resultados2$upper[7:1], 1:7,  lwd=2.5, code=0)#, col="blue")
abline(v=0, lty=2, col="red")

plot(resultados2$Mean[14:8], 1:7,yaxt="n", ylab="", xlim=c(-30,20),
     pch=19, main="Spectral")
grid()
columnas<-resultados2$Clade[8:14]
axis(2, at=7:1, as.graphicsAnnot(columnas), las=2, lwd=0.5)
arrows(resultados2$under[14:8], 1:7, resultados2$upper[14:8], 1:7,  lwd=2.5, code=0)#, col="blue")
abline(v=0, lty=2, col="red")

plot(resultados2$Mean[21:15], 1:7,
     yaxt="n", ylab="", xlim=c(-40,20), pch=19,
     main="Temporal")
grid()
columnas<-resultados2$Clade[15:21]
axis(2, at=7:1, as.graphicsAnnot(columnas), las=2, lwd=0.5)
arrows(resultados2$under[21:15], 1:7, resultados2$upper[21:15], 1:7,  lwd=2.5, code=0)#, col="blue")
abline(v=0, lty=2, col="red")




####################

pglsAngy <- gls(AcousticField~MPD, 
                correlation = corPagel(0.5, phy = Tree1, fixed =TRUE), 
                data=Elliot, method = "ML")
pglsAngyt <- gls(AcousticField~MPD,
                 correlation = corPagel(0.5, phy = Treet, fixed =TRUE), 
                 data=Elliott, method = "ML")
pglsAngys <- gls(AcousticField~MPD, 
                 correlation = corPagel(0.5, phy = Trees, fixed =TRUE), 
                 data=Elliots, method = "ML")

lillie.test(pglsAngy$residuals) #
lillie.test(pglsAngyt$residuals) #
lillie.test(pglsAngys$residuals) #
dev.off()
summary(pglsAngy)
summary(pglsAngyt)
summary(pglsAngys)


pglsAngy$coefficients
########
View(Comunidad2)
View(cantosprom0)
ambiente<-cantosprom0[,9:10]
plot(ambiente)

#################
df <- data.frame()
ggplot(Elliot, aes(MPD,AcousticField))+geom_point()+theme_classic()+
#ggplot(df) + geom_point() +
  xlim(0,0.2)+ylim(-2, 10)+
  geom_abline(intercept = 0.5861, slope=3.8859, color='gray')+
  geom_abline(intercept =1.2040, slope=-21.1443, color='blue')+
  geom_abline(intercept = 7.9267, slope=20.7118, color='red')
#  geom_text(aes(label=Species),vjust=0) +
  geom_abline(intercept=coef(pglsAngys)[1], slope=coef(pglsAngys)[2])

#1.7786	-5.5883	/ -0.1138	1.5249	-4.1564	10.7014 #COchranellini temporal
#0.4733	4.4978/ -2.1561	-0.1187	8.8851	50.1192 #Hyalinosubclade spectral
#0.5861	3.8859 / 1.2040	7.9267	-21.1443	20.7118 #Centrolenidae Overall >3


#Intercept	MPD	Intercept	MPD	Intercept	MPD	2.5 Intercept	9.75 intercept	2.5 MPD	9.75 MPD
#1.0794	-0.2142	1.3120	-1.2801	0.3099	3.3701	2.3937	6.1615	-26.0449	-2.4788
abline(fit)
plot(lm(Y~0.2142X+1.0794))


test2<-comparative.data(data=Elliot, Tree1, names.col="Species", vcv.dim=2, warn.dropped=TRUE)
pglsAngyt<-pgls(AcousticField~MPD, data=test2, lambda="ML")
summary(pglsAngyt)

install.packages("adephylo")
library(adephylo)
biplot(canto)
canto$Evec
canto3<-ppca(Tree0, ScaleCall)

summary(canto)
tempcol <- rep("grey",7)
barplot(canto$eig, main='pPCA eigenvalues',cex.main=1.8,col=tempcol)


call2Da<-as.data.frame(call2D)
Genus<-subset(cantospromAll$Genus, subset=rownames(cantospromAll) %in% rownames(call2Da))
dim(cantospromAll)

ggplot(call2Da, aes(PC1,PC2, color=Genus))+geom_point(size=3)+ theme_minimal()+
  scale_color_brewer(palette="Paired") + theme(legend.position="bottom")
  geom_text(aes(label=rownames(call2D)),hjust=0, vjust=0)+theme_minimal()

scale_color_jcolors(palette = "pal8")
  theme_classic()

ggplot(Elliot, aes(MPD,AcousticField))+geom_point()+
  geom_text(aes(label=Species),hjust=0, vjust=0) +
  geom_abline(intercept=coef(pglsAngy)[1], slope=coef(pglsAngy)[2])



AcousticField0<-(Elliot$MPD)
names(AcousticField0)<-row.names(Elliot)
treeplot<-drop.tip(CentroTree, setdiff(CentroTree$tip.label, names(AcousticField0)))
obj<-contMap(treeplot, AcousticField0,fsize=c(0.5,1),outline=FALSE)
plot(obj, type="fan", fsize=c(0.5,1),outline=FALSE)
  
#summary(lm(AcousticField~MPD, data=Elliot))

#################
prior1<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
prior2 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000)))
a <- 1000
prior3 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))
prior4<- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))

prior = list(R = list(V = diag(2)/3, n = 2),
                G = list(G1 = list(V = diag(2)/3, n = 2),
                G2 = list(V = diag(2)/3, n = 2)))
prior.m1a.3 <- list(R = list(V = 1, nu = 0))
Ainv4<-inverseA(Tree1, scale=F)$Ainv
Ainv4<-inverseA(Treet, scale=F)$Ainv
Ainv4<-inverseA(Trees, scale=F)$Ainv

####
model2.1<-MCMCglmm(AcousticField~MPD, data=Elliot, 
                   ginverse=list(Species=Ainv4), prior = prior.m1a.3,
                   nitt=100000, thin=10, burnin=300, verbose=FALSE)
summary(model2.1)

plot(model2.1)
round(sort(effectiveSize(model2.1$Sol)))
geweke.diag(model2.1$Sol)
rownames(modelsPC1t)



#######################
######Null Model#######
#######################
modelsPC1a<-subset(modelsPC1, subset=rownames(modelsPC1) %in% rownames(Elliot))
modelsPC2a<-subset(modelsPC1, subset=rownames(modelsPC1) %in% rownames(Elliot))

modelsPC1ta<-subset(modelsPC1t, subset=rownames(modelsPC1t) %in% rownames(Elliott))
modelsPC2ta<-subset(modelsPC1t, subset=rownames(modelsPC1t) %in% rownames(Elliott))

modelsPC1sa<-subset(modelsPC1s, subset=rownames(modelsPC1s) %in% rownames(Elliots))
modelsPC2sa<-subset(modelsPC1s, subset=rownames(modelsPC1s) %in% rownames(Elliots))

ElliotNull<-Elliot
atlassong<-subset(atlas, select=colnames(atlas) %in% species2) 
atlassongNull<-subset(atlassong, select=colnames(atlassong) %in% rownames(ElliotNull))
atlassongNull<-subset(atlassongNull, subset=rowSums(atlassongNull) > 0)
SpeciesRichnessNull<- as.matrix(rowSums(atlassongNull))
RangeSizeNull<- as.matrix(colSums(atlassongNull))
Estimate<-rep(NA, 100)
new.lambda<-0.5
system.time(
for(i in 1:length(Estimate)){
    matrixNull<-as.matrix(cbind(modelsPC1a[,i], modelsPC2a[,i]))
    functionalDispNull<-FDis(ordination.results=matrixNull, road.map=atlassongNull)
    atlasfunctNull<-functionalDispNull*atlassongNull
    D_VolumeNull<- t(atlasfunctNull)%*%SpeciesRichnessNull 
    AcousticFieldNull<- D_VolumeNull/RangeSizeNull #Dispersión acustica por especie
    ElliotNull[,1]<-AcousticFieldNull[,1]
    pglsNull<-gls(AcousticField~MPD, correlation = corPagel(new.lambda, phy = Tree1, fixed = TRUE), data=ElliotNull)
    EstimateA<-pglsNull$coefficients[1]
    Estimate[i]<-EstimateA
    print(i)
})
summary(Estimate)
summary(pglsAngyt)

contMap(Tree1, AcousticFieldNull)

####################
ElliotNull<-Elliott #Elliott /Elliots
atlassong<-subset(atlas, select=colnames(atlas) %in% speciesTemp) 
atlassongNull<-subset(atlassong, select=colnames(atlassong) %in% rownames(ElliotNull))
atlassongNull<-subset(atlassongNull, subset=rowSums(atlassongNull) > 0)
atlassongNull<-subset(atlassongNull, select=colSums(atlassongNull) > 0)
SpeciesRichnessNull<- as.matrix(rowSums(atlassongNull))
RangeSizeNull<- as.matrix(colSums(atlassongNull))
Estimate<-matrix(rep(NA, 200), 100,2)

system.time(
  for(i in 1:nrow(Estimate)){
    matrixNull<-as.matrix(cbind(modelsPC1ta[,i], modelsPC2ta[,i]))
    functionalDispNull<-FDis(ordination.results=matrixNull, road.map=atlassongNull)
    atlasfunctNull<-functionalDispNull*atlassongNull
    D_VolumeNull<- t(atlasfunctNull)%*%SpeciesRichnessNull 
    AcousticFieldNull<- D_VolumeNull/RangeSizeNull #Dispersión acustica por especie
    ElliotNull[,1]<-AcousticFieldNull[,1]
    pglsNull<-gls(AcousticField~MPD, 
                  correlation = corPagel(new.lambda, phy = Treet, fixed = TRUE), data=ElliotNull)
    Estimate[i,]<-pglsNull$coefficients
    print(i)
  })
summary(Estimate)

summary(pglsAngy)
summary(pglsAngyt)
summary(pglsNull)

dim(atlassong)
dim(atlassongNull)
dim(matrixNull)

setdiff(rownames(matrixNull),colnames(atlassongNull))
View(matrixNull)

x<-Estimate[,2]
intervalo<-c(quantile(x,.025), quantile(x,.975))
###################
prior2 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000)))
a <- 1000
prior3 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))
prior4<- list(R = list(V = diag(1), nu = 0.002),
              G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))

prior = list(R = list(V = diag(2)/3, n = 2),
             G = list(G1 = list(V = diag(2)/3, n = 2),
                      G2 = list(V = diag(2)/3, n = 2)))

prior.m1a.3 <- list(R = list(V = 1, nu = 0))



#summary(pglsNull)
Estimate1<-data.frame(Estimate)
summary(Estimate)

summary(songGLMS)
real<-1.079
real<-pglsAngyt$coefficients[1]
p<-ggplot(Estimate1, aes(x=Estimate[,1])) #+ geom_histogram(aes(y=..density..), binwidth=1, fill="white",colour="gray")
p+geom_density(alpha=0.3)+ geom_vline(aes(xintercept=real),
 color="blue", linetype="dashed", size=1)

#########
install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R",
                 type="source")
library(coefplot2)


