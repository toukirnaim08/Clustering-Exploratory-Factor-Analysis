# A general summary :
#   1)  Read csv from file as dataframe.
#   2)  firstly conduct EFA on the data set  
#   3)  secondly conduct cluster analysis on the data set


#load external sources
source('mosthighlycorrelated.r')
source('circle_cor.R')
source('kmo.r')
source('residual_stats.R')

library(readr)
library(gridExtra)
library(grid)
library(gridBase)
library(dplyr)
library(ggplot2)
library(rebus)
library(psych)
library(GPArotation)
library(cluster)
library(MASS)


#read file from the same directory of this R script
currentDirectory <-dirname(rstudioapi::getSourceEditorContext()$path)
filepathDir<- paste(currentDirectory,"/chocolate.csv",sep = "")
chocolateRawData <- read.csv(filepathDir, row.names=NULL)

View(chocolateRawData)
dim(chocolateRawData)

#split dataset into two based on expert and amateur.
expertRawData<-chocolateRawData[chocolateRawData$Role == "expert", ]
dim(expertRawData)
View(expertRawData[,2:15])
amateurRawData<-chocolateRawData[chocolateRawData$Role == "amateur", ]
rownames(amateurRawData) <- NULL
dim(amateurRawData)
View(amateurRawData)

#function for readable float value
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

#function for making table
createTable<- function(title1,Desc1,Desc2){
  grid.newpage()
  tempSummary<-c(title1,Desc1,Desc2)
  tbl <- tableGrob(tempSummary)
  grid.arrange(tbl)
}


#check whether large amount of variable ratio to conduct analysis
expertRawDataDim<-dim(expertRawData)
ratio1<- expertRawDataDim[1]/expertRawDataDim[2]

amateurRawDataDim<-dim(amateurRawData)
ratio2<-amateurRawDataDim[1]/amateurRawDataDim[2]
#create ratio table
createTable("Amount of variable ratio to conduct analysis",
            paste("Experts Responses Ratio",ratio1,sep = ": "),
            paste("Amateurs Responses Ratio",ratio2,sep = ": "))

# ################################################ Cronbach Alpha #############################
#check cronbach alpha to check the information is valid and we can trust
#since dataset has catagorical and numerical data so we should check cronbach alpha 
psych::alpha(expertRawData[,2:15]) #raw alpha 0.19 there might be ambiguity
psych::alpha(amateurRawData[,2:15])  #raw alpha 0.38 not relyable

# ################################################ Correlation #############################
#Correlation
#top corelation (multi coliniarity have blue positive green neagtive)
expertRawData.cas <- cor(expertRawData[,2:15], use="complete.obs")
colnames(expertRawData.cas) <- c("Chocolate.Aroma","Milk.Aroma","Sweetness","Acidity","Bitterness","Chocolate.Flavour","Milk.Flavour","Caramel.Flavour","Vanilla.Flavour","Astringency","Crispy.Texture","Melting.Texture","Sticky.Texture","Granular.Texture")
#top 10 correaltions
mosthighlycorrelated(expertRawData.cas,10)
#correaltions plot
circle_cor(expertRawData.cas)

amateurRawData.cas <- cor(amateurRawData[,2:15], use="complete.obs")
colnames(amateurRawData.cas) <- c("Chocolate.Aroma","Milk.Aroma","Sweetness","Acidity","Bitterness","Chocolate.Flavour","Milk.Flavour","Caramel.Flavour","Vanilla.Flavour","Astringency","Crispy.Texture","Melting.Texture","Sticky.Texture","Granular.Texture")
#top 10 correaltions
mosthighlycorrelated(amateurRawData.cas,10)
#correaltions plot
circle_cor(amateurRawData.cas)

# ################################################ Determinant test #############################
#determinant test not too close to zero (multicolniarity even for factor analysis)
det(expertRawData.cas) > 1e-05
det(amateurRawData.cas) > 1e-05

# ################################################ Bartlett’s test #############################
#cor matrix is not identity matrix(Bartlett’s test)
cortest.bartlett(expertRawData.cas,n=expertRawDataDim[1])#(significantly different from identity matrix)
cortest.bartlett(amateurRawData.cas,n=amateurRawDataDim[1])

# ################################################ kmo statistics ############################# 
#kmo how we sample our data appropriantly
expertRawDataKmo <- kmo(expertRawData.cas)
cbind(round(expertRawDataKmo$overall,2), expertRawDataKmo$report)
sort(round(expertRawDataKmo$individual,2))

amateurRawDataKmo <- kmo(amateurRawData.cas)
cbind(round(amateurRawDataKmo$overall,2), amateurRawDataKmo$report)
sort(round(amateurRawDataKmo$individual,2))

# ################################################ number of factors to extract ############################# 
#pca on standertised data
expertRawData.pca <- prcomp(expertRawData.cas, scale=TRUE) # >75% rules we should pick 2
# Summary of PCA results
summary(expertRawData.pca)
# How many components to retain?
quartz()
screeplot(expertRawData.pca, type="lines",pch=10,col='blue',main="Screeplot of expert responses")
#eigenvalues
round((expertRawData.pca$sdev)^2 , 3)

#pca on standertised data
amateurRawData.pca <- prcomp(amateurRawData.cas, scale=TRUE) # >75% rules we should pick 2
# Summary of PCA results
summary(amateurRawData.pca)
# How many components to retain?
quartz()
screeplot(amateurRawData.pca, type="lines",pch=10,col='blue',main="Screeplot of amateur responses")
#eigenvalues
round((amateurRawData.pca$sdev)^2 , 3)


# ################################################ factor analysis no rotation ############################# 
#pca
nf <- 3  #performing 2 , 3 and 4 i have decided to choose 3 factors as it gives me better outputs
#h2 comnality u2 uniqness some split loading
expertRawData.pca.none <- principal(expertRawData.cas, nfactors=nf, rotate="none")
print.psych(expertRawData.pca.none, cut = 0.4, sort = TRUE)

amateurRawData.pca.none <- principal(amateurRawData.cas, nfactors=nf, rotate="none")
print.psych(amateurRawData.pca.none, cut = 0.4, sort = TRUE)

#ML extraction
expertRawData.ml.none <- factanal(cov=expertRawData.cas,factors=nf,rotation="none")
print(expertRawData.ml.none, digits = 3, cutoff = 0.4, sort = TRUE)

amateurRawData.ml.none <- factanal(cov=amateurRawData.cas,factors=nf,rotation="none")
print(amateurRawData.ml.none, digits = 3, cutoff = 0.4, sort = TRUE)

#PA Factoring
expertRawData.paf.none <- fa(expertRawData.cas,nfactors=nf,fm="pa")
print.psych(expertRawData.paf.none, cut = 0.4, sort = TRUE)

amateurRawData.paf.none <- fa(amateurRawData.cas,nfactors=nf,fm="pa")
print.psych(amateurRawData.paf.none, cut = 0.4, sort = TRUE)

  
# ################################################ factor analysis with rotation ############################# 
#pca with varimax rotation
expertRawData.pca.var <- principal(expertRawData.cas, nfactors=nf, rotate="varimax")
print.psych(expertRawData.pca.var, cut = 0.4, sort = TRUE)
sort(round(expertRawData.pca.var$communality,2)) #very few communalities satisfied

amateurRawData.pca.var <- principal(amateurRawData.cas, nfactors=nf, rotate="varimax")
print.psych(amateurRawData.pca.var, cut = 0.4, sort = TRUE)
sort(round(amateurRawData.pca.var$communality,2)) #very few communalities satisfied
  
# check residuals. 
expertRawData_resids<-factor.residuals(expertRawData.cas, expertRawData.pca.var$loadings)
residual_stats(expertRawData_resids)

amateurRawData_resids<-factor.residuals(amateurRawData.cas, amateurRawData.pca.var$loadings)
residual_stats(amateurRawData_resids)

#plot of rotated expert factors
par(mfrow=c(1,3))
plot(loadings(expertRawData.pca.var),xlim=c(-3,3),ylim=c(-0.5,1.2),cex=c(1,1),pch=16,col='skyblue',main="PCA with Varimax") 
text(loadings(expertRawData.pca.var)-c(0.03,0.0,0.03,0.0,0.03,0.0,0.03,0.03,0.03,0.03,0.03,0.03), labels=dimnames(expertRawData.cas)[[1]],cex=1.2) 
abline(h=0,lty=2,col='gray')
abline(v=0,lty=2,col='gray')
  #on to that first factor, some 2nd some in between them
  
#plot of rotated amateur factors
par(mfrow=c(1,3))
plot(loadings(amateurRawData.pca.var),xlim=c(-3,3),ylim=c(-0.5,1.2),cex=c(1,1),pch=16,col='skyblue',main="PCA with Varimax") 
text(loadings(amateurRawData.pca.var)-c(0.03,0.0,0.03,0.0,0.03,0.0,0.03,0.03,0.03,0.03,0.03,0.03), labels=dimnames(amateurRawData.cas)[[1]],cex=1.2) 
abline(h=0,lty=2,col='gray')
abline(v=0,lty=2,col='gray')
  #on to that first factor, some 2nd some in between them
  
# Maximum Likelihood Extraction with Varimax rotation
# for expert data
expertRawData.ml.var <- factanal(cov=expertRawData.cas,factors=nf,rotation="varimax")
print(expertRawData.ml.var , digits = 3, cutoff = .4, sort = TRUE)
# check communalities 
sort(round(1 - expertRawData.ml.var$uniqueness, 2))
# check residuals
expertRawData_resids<-factor.residuals(expertRawData.cas, expertRawData.ml.var$loadings)
residual_stats(expertRawData_resids)

# check factor scores (An orthogonal rotation is appropriate as correlations are close to zero)
round(expertRawData.pca.var$r.scores,3) #shoudint think oblique rotation
#plot of rotated expert factors
plot(loadings(expertRawData.ml.var),xlim=c(-3,3),ylim=c(-0.5,1.2),cex=c(1,1),pch=16,col='lightgreen',main="MLE with Varimax") 
text(loadings(expertRawData.ml.var)-c(0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03), labels=dimnames(expertRawData.cas)[[2]],cex=1.2) 
abline(h=0,lty=2,col='gray')
abline(v=0,lty=2,col='gray')  

# for amateur data
amateurRawData.ml.var <- factanal(cov=amateurRawData.cas,factors=nf,rotation="varimax")
print(amateurRawData.ml.var , digits = 3, cutoff = .4, sort = TRUE)
# check communalities
sort(round(1 - amateurRawData.ml.var$uniqueness, 2))
# check residuals
amateurRawData_resids<-factor.residuals(amateurRawData.cas, amateurRawData.ml.var$loadings)
residual_stats(amateurRawData_resids)
# check factor scores (An orthogonal rotation is appropriate as correlations are close to zero)
round(amateurRawData.pca.var$r.scores,3) #shoudint think oblique rotation
#plot of rotated amateur factors
plot(loadings(amateurRawData.ml.var),xlim=c(-3,3),ylim=c(-0.5,1.2),cex=c(1,1),pch=16,col='lightgreen',main="MLE with Varimax") 
text(loadings(amateurRawData.ml.var)-c(0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03), labels=dimnames(amateurRawData.cas)[[2]],cex=1.2) 
abline(h=0,lty=2,col='gray')
abline(v=0,lty=2,col='gray')  
  
# Principal Axis factoring with rotation.
# for expert data
expertRawData.paf.var <- fa(expertRawData.cas,nfactors=nf,rotate="varimax",fm="pa")
print.psych(expertRawData.paf.var, cut = 0.4, sort = TRUE)
# check communalities 
sort(round(expertRawData.paf.var$communality,2))
# check residuals.
expertRawData_resids<-factor.residuals(expertRawData.cas, expertRawData.paf.var$loadings)
residual_stats(expertRawData_resids)
# check factor scores
round(expertRawData.paf.var$r.scores,2)
#loadings plot
plot(loadings(expertRawData.paf.var),xlim=c(-3,3),ylim=c(-0.5,1.2),cex=c(1,1),pch=16,col='coral',main="PAF with Varimax") 
text(loadings(expertRawData.paf.var)-c(0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03), labels=dimnames(expertRawData.cas)[[2]],cex=1.2) 
abline(h=0,lty=2,col='gray')
abline(v=0,lty=2,col='gray')

# for amateur data
amateurRawData.paf.var <- fa(amateurRawData.cas,nfactors=nf,rotate="varimax",fm="pa")
print.psych(amateurRawData.paf.var, cut = 0.4, sort = TRUE)
# check communalities 
sort(round(amateurRawData.paf.var$communality,2))
# check residuals.
amateurRawData_resids<-factor.residuals(amateurRawData.cas, amateurRawData.paf.var$loadings)
residual_stats(amateurRawData_resids)
# check factor scores
round(amateurRawData.paf.var$r.scores,2)
#loadings plot
plot(loadings(amateurRawData.paf.var),xlim=c(-3,3),ylim=c(-0.5,1.2),cex=c(1,1),pch=16,col='coral',main="PAF with Varimax") 
text(loadings(amateurRawData.paf.var)-c(0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03), labels=dimnames(amateurRawData.cas)[[2]],cex=1.2) 
abline(h=0,lty=2,col='gray')
abline(v=0,lty=2,col='gray')
  

# ################################################ clustering #############################   
#modify dataset for clustering
View(chocolateRawData[2:15])
modifiedChocolateData <- as.data.frame(t(chocolateRawData[2:15]))
View(modifiedChocolateData)
# Calculate the distance matrix, using euclidean distance
d = dist(modifiedChocolateData,method="euclidean")
d

######################## Hierarchical Clustering ###############
# single linkage 
hs <- hclust(d,method="single")
max(hs$height)
heightseq <- seq(0, 80, 10)
quartz()
plot(hs, col="darkblue", col.main="#45ADA8", col.lab="#7C8071",col.axis="#F38630", lwd=1.5, lty=3,cex=.75, sub='', hang=-1, axes=FALSE,main='Chocolates: AGNES with Single Linkage')
axis(side=2, at=heightseq, col="#F38630",labels=FALSE, lwd=1)
mtext(heightseq, side=2, at=heightseq,line=1, col="#A38630", las=2)

#complete linkage
hc <- hclust(d,method="complete")
max(hc$height)
heightseq <- seq(0, 150, 10)
quartz()
plot(hc, col="darkblue", col.main="#45ADA8", col.lab="#7C8071",col.axis="#F38630", lwd=1.5, lty=3,cex=.75, sub='', hang=-1, axes=FALSE,main='Chocolates: AGNES with Complete Linkage')
axis(side=2, at=heightseq, col="#F38630",labels=FALSE, lwd=1)
mtext(heightseq, side=2, at=heightseq,line=1, col="#A38630", las=2)

# Average Linkage
ha = hclust(d,method="average")
max(ha$height)
heightseq <- seq(0, 110, 10)
quartz()
plot(ha, col="darkblue", col.main="#45ADA8", col.lab="#7C8071",col.axis="#F38630", lwd=1.5, lty=3,cex=.75, sub='', hang=-1, axes=FALSE,main='Chocolates: AGNES with Average Linkage')
axis(side=2, at=heightseq, col="#F38630",labels=FALSE, lwd=1)
mtext(heightseq, side=2, at=heightseq,line=1, col="#A38630", las=2)

#wards mathod
hw <- hclust(d,method="ward.D2")
max(hw$height)
heightseq <- seq(0, 220, 20)
quartz()
plot(hw, col="darkblue", col.main="#45ADA8", col.lab="#7C8071",col.axis="#F38630", lwd=1.5, lty=3,cex=.75, sub='', hang=-1, axes=FALSE,main='Chocolates: AGNES with Wards Method')
axis(side=2, at=heightseq , col="#F38630",labels=FALSE, lwd=1)
mtext(heightseq , side=2, at=heightseq ,line=1, col="#A38630", las=2)

######################## Partitional Clustering ###############
#partitional clustering ## k mean
# Determine the sample size n and the number of variables we have p
n = dim(modifiedChocolateData)[1]
p = dim(modifiedChocolateData)[2]
# Compute variances
SSE <- (n - 1) * sum(apply(modifiedChocolateData,2,var)) 
# We compute the SSE for 1 to 10 clusters (ie k=1,...,10) and plot the 10 values
for (i in 2:10) {
  SSE[i] <- sum(kmeans(modifiedChocolateData,centers=i,nstart=25)$withinss)
}
dev.off()
plot(1:10, SSE, type="b", xlab="Number of Clusters", ylab="Sum of squares within groups",pch=19, col="blue")

pc.km <- kmeans(modifiedChocolateData, centers = 3)
dev.off()
clusplot(d, pc.km$cluster, diss = TRUE, cex=0.7,col.p='midnightblue', col.clus='seagreen3',col.txt='blue', labels=3, main="Clusters")

######################## Validation ###############
# cycle through plots and choose using the 'best' silhouette
asw <- numeric(10)
for (k in 2:10) {
  km <- kmeans(modifiedChocolateData, centers = k)
  si <- silhouette(km$cluster,d)
  asw[k] <- mean(si[,3])
}
k.best <- which.max(asw)
cat("silhouette-optimal number of clusters:", k.best, "\n")
plot(1:10, asw, type= "h", main = "k-means clustering assessment",
     xlab= "k  (# clusters)", ylab = "average silhouette width")
axis(1, k.best, paste("best",k.best,sep="\n"), col = "seagreen3", col.axis = "seagreen3")

# Produce the silhouette plot
pc.km <- kmeans(modifiedChocolateData, centers = 3)
si <- silhouette(pc.km$cluster,d)
ssi <- summary(si)
plot(si, col = c("green", "blue"),main='silhouette plot')


##create new table group by role with mean ratings
meanChocolateRawData<-aggregate(chocolateRawData[,2:15], by=list(chocolateRawData$Role), FUN=mean)
View(meanChocolateRawData)

barplot(as.matrix(meanChocolateRawData[,2:14]),cex.names=0.6,col = c(1,2),ylim=c(0,12),las=2)
legend("topright", fill = 1:2, legend = c('Amateur','Expert'), 
       horiz = TRUE, inset = c(0,-0.1), xpd = TRUE)


