###Summary analisys between species###

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

#Read Tree
CentroTree<-read.nexus("outjan2021.tre")
plot(CentroTree, no.margin=T, cex=0.5)

##Read traits file
cantosprom0<-read.table("CantosandExternal.txt", header=T)
names(cantosprom0)[1]<-"Species"
rownames(cantosprom0)<-cantosprom0$Species

#Subclados 
#[13] "Centroleninae"      "Hyalinobatrachinae" "Hyalinobatrachium" 
#[16] "Hyalino17sp"        "CentroNymph"        "Centrolene"        
#[19] "X8genera"           "X5genera"
cantospromAll<-subset(cantosprom0, cantosprom0$X8genera==1) #Centroleninae
dim(cantospromAll)

#If the analysis is for all species in the family
cantospromAll<-cantosprom0#[,1:14]

###checking trait changes along the tree##

data3<-na.omit(subset(cantospromAll, 
                      select=c("Species","SVL", "EVI","PeakFreq", "Temperature")))

x<-setNames(data3$SVL, rownames(data3))
y<-setNames(data3$PeakFreq, rownames(data3))
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
#Subclades 
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
model3.3<-MCMCglmm(fixed=PeakFreq~Temperature, random=~Species+EVI+SVL, data=data3, 
                   ginverse=list(Species=Ainv4), prior = prior4,
                   nitt=500000, thin=10, burnin=300, verbose=FALSE)
summary(model3.3)

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
#####Community comparisons######
################################

#Species diversity
#Range Overlap
atlas<-read.table(file="PAM2019-2.txt", header=T)
head(atlas)
colnames(atlas) #Species Names

atlas1<-as.matrix(atlas)
overlap<-lets.overlap(atlas1, method = "Chesser&Zink", xy = F)

##############################################
##########Pairwise species analysis###########
##############################################

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
#subset only spectral or temporal parameters
overlap1t<-subset(overlap, rownames(overlap) %in% speciesTemp) 
overlap2t<-subset(t(overlap1t), colnames(overlap1t) %in% speciesTemp) #speciesTemp

overlap1s<-subset(overlap, rownames(overlap) %in% speciesSpec) 
overlap2s<-subset(t(overlap1s), colnames(overlap1s) %in% speciesSpec) #speciesTemp

tryTemp<-matrixLMM(PPCATempDist, overlap2t, TreeTemp) #/ (PPCATempDist, overlap2t, Tree4)
trySpec<-matrixLMM(PPCASpecDist, overlap2s, TreeSpec) #/ (PPCATempDist, overlap2t, Tree4)

##MCMCGLMM Temporals
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

########################################
#######Community approach###############
########################################

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

#####
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
