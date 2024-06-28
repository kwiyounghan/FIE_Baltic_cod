##############################
###### PCA analysis 
###### for 115 random samples with and without inversions and only inversions 
##############################
## 00. Load Libraries
#update.packages(checkBuilt = TRUE, ask = FALSE)
#install.packages("dplyr")
#install.packages("magrittr")
#install.packages("ggplot2")
#install.packages("reshape2")
library(dplyr)
#install.packages("pcadapt")
library(pcadapt)

## 00.1 Set colors 
hex_codes <- c("#ff595e","#ffca3a","#8AC926","#1982c4","#6A4C93")
show_col(hex_codes)

## 01. Load metadata
setwd("/home/khan/cod3/11_pca")
## 01.1 load population map
popmap <- read.csv("all_sample_metadata_233.csv", h=T,sep=",")
popmap$'SeqSample' #no need to reorder samples 

#exclude 2 samples of low sequencing depth
grep("J32540",popmap$'SeqSample')
popmap[c(140,141),]
popmap <- popmap[-c(140,141),]
dim(popmap)
##exclude samples that are not tempEBC = only "random" samples
popmap <- popmap[popmap$class=="temp",]

## 01.2 set colors and shapes for each group 
head(popmap)
pop <- popmap[,'EBC_class']
pop
pop_col <- ifelse(pop=="1996",hex_codes[1],
                  ifelse(pop=="2002",hex_codes[2],
                         ifelse(pop=="2008",hex_codes[3],
                                ifelse(pop=="2014",hex_codes[4],
                                       ifelse(pop=="2019",hex_codes[5],"black"
                                       )))))
##################
###### PCA with inversions 
##################
## 02.1 load geno data
list.files(pattern=".bed")
f <- "cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.01_tempEBC_115_clean"
## 02.2 Run pcadapt 
geno.012 <- read.pcadapt(input=paste0(f,".bed"), type="bed")
x <- pcadapt(input = geno.012, K=7)
dim(x$scores)
dim(x$loadings)
attributes(x)
#rm(geno.012)
#save(x, file="cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.01_tempEBC_115_clean_pcadapt.Rdata")
#load(file="cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.01_tempEBC_115_clean_pcadapt.Rdata")

## 03. Plot
## 03.1 Scree plot
png(filename=paste0(f,".screeplot.png",sep=""),width=1000,height=1000,res=300, units='px')
plot(x,option="screeplot")
dev.off()

## 03.2 Plot PCA, PC1 and PC2, pc3 and pc4
png(filename=paste0(f,".PC12.png",sep=""),width=4000,height=2000,res=300, units='px')
par(mfrow=c(1,2),oma=c(1,1,2,1))
for (i in c(1,3)){
  var1 <- round(cumsum((x$singular.values[i])^2)*100,digits=2)
  var2 <- round(cumsum((x$singular.values[i+1])^2)*100,digits=2)
  
  plot(x$scores[,i],x$scores[,i+1], 
       col=pop_col, pch=19,
       cex=1.8,xlab=paste("PC",i," ",var1,"%",sep=""),ylab=paste("PC",i+1," ",var2,"%",sep=""),cex.lab=1.5)
  legend("topleft",title=c("A","","B")[i],title.font=2,legend="",bty="n")
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x="top",
       legend=c("1996","2002","2008","2014","2019"),
       col=c(hex_codes[1],hex_codes[2],hex_codes[3],hex_codes[4],hex_codes[5]),
       pch=c(19,19,19,19,19),bty='n',pt.cex=3, cex=1.4,horiz=TRUE)
dev.off()


## 03.3 Plot loadings
x.loadings <- as.data.frame(x$loadings)
x.loadings <- abs(x.loadings)

png(filename=paste0(f,".pc12_loadings.png",sep=""),width=4000,height=3000,res=300, units='px')
par(mfrow=c(2,1))
for(i in 1:2){
  plot(x.loadings[,i],col=alpha("black",0.4),cex=1,pch=20, xlab="chromosomal positions",ylab=paste0("PC",i," loadings"))
}
dev.off()

png(filename=paste0(f,".pc34_loadings.png",sep=""),width=4000,height=3000,res=300, units='px')
par(mfrow=c(2,1))
for(i in 3:4){
  plot(x.loadings[,i],col=alpha("black",0.4),cex=1,pch=20, xlab="chromosomal positions",ylab=paste0("PC",i," loadings"))
}
dev.off()


##################
###### PCA without inversions 
##################
## 02.1 load geno data
list.files(pattern=".bed")
f <- "cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.01_tempEBC_115_clean_noINV_50_10_0.5"

## 02.2 Run pcadapt 
geno.012 <- read.pcadapt(input=paste0(f,".bed"), type="bed")
x <- pcadapt(input = geno.012, K=7)
dim(x$scores)
attributes(x)
#rm(geno.012)
#save(x, file="cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.01_tempEBC_115_clean_noINV_10_0.5_pcadapt.Rdata")
#load(file="cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.01_tempEBC_115_clean_noINV_10_0.5_pcadapt.Rdata")
#dim(x$loadings)## 2030929       7

## 03. Plot
## 03.1 Scree plot
png(filename=paste0(f,".screeplot.png",sep=""),width=1000,height=1000,res=300, units='px')
plot(x,option="screeplot")
dev.off()

## 03.2 Plot PCA, PC1 and PC2
png(filename=paste0(f,".PC1.png",sep=""),width=2000,height=2000,res=300, units='px')
par(mfrow=c(1,1),oma=c(1,1,1,1),mar=c(5,5,3,2))
for (i in c(1)){
  var1 <- round(cumsum((x$singular.values[i])^2)*100,digits=2)
  var2 <- round(cumsum((x$singular.values[i+1])^2)*100,digits=2)
  
  plot(x$scores[,i],x$scores[,i+1], 
       col=pop_col, pch=19,
       cex=1.8,xlab=paste("PC",i," ",var1,"%",sep=""),ylab=paste("PC",i+1," ",var2,"%",sep=""),cex.lab=1.5)
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x="top",
       legend=c("1996","2002","2008","2014","2019"),
       col=c(hex_codes[1],hex_codes[2],hex_codes[3],hex_codes[4],hex_codes[5]),
       pch=c(19,19,19,19,19),bty='n',pt.cex=3, cex=1.4,horiz=TRUE)
dev.off()

# identify(x$scores[,1],x$scores[,2],tolerance = 1)
# popmap[62,]
# popmap[63,]
# popmap[91,]
## 03.3 Plot loadings
x.loadings <- as.data.frame(x$loadings)
x.loadings <- abs(x.loadings)

png(filename=paste0(f,".pc12_loadings.png",sep=""),width=4000,height=3000,res=300, units='px')
par(mfrow=c(2,1))
for(i in 1:2){
  plot(x.loadings[,i],col=alpha("black",0.4),cex=1,pch=20, xlab="chromosomal positions",ylab=paste0("PC",i," loadings"))
}
dev.off()

png(filename=paste0(f,".pc34_loadings.png",sep=""),width=4000,height=3000,res=300, units='px')
par(mfrow=c(2,1))
for(i in 3:4){
  plot(x.loadings[,i],col=alpha("black",0.4),cex=1,pch=20, xlab="chromosomal positions",ylab=paste0("PC",i," loadings"))
}
dev.off()



##################
###### PCA on inversions
##################
# load geno data
list.files(pattern=".bed")
#which files are the target files?
huh <- c(28,29,30)
tit <- c("LG02","LG07","LG12")
#pcadapt
f <- "cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.01_tempEBC_115_clean"
png(filename=paste0(f,"_INV.PC1_2.png",sep=""),width=4000,height=1500,res=300, units='px')
par(mfrow=c(1,3),oma=c(1,1,3,1),mar=c(4,5,5,1))
for (i in huh){
  f <- list.files(pattern=".bed")[i]
  geno.012 <- read.pcadapt(input=f, type="bed")
  x <- pcadapt(input = geno.012, K=2)
  var1 <- round(cumsum((x$singular.values[1])^2)*100,digits=2)
  var2 <- round(cumsum((x$singular.values[2])^2)*100,digits=2)
  
  plot(x$scores[,1],x$scores[,2], col=pop_col,pch=19, cex=2,xlab=paste("PC1 ",var1,"%",sep=""),ylab=paste("PC2 ",var2,"%",sep=""),cex.lab=1.6,cex.axis=1.5)
  title(tit[i-27], cex.main=2)
  #legend("topleft",title=c("A","B","C")[i-27],title.font=2,legend="",bty="n",cex=2)
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 1, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x="top",legend=c("1996","2002","2008","2014","2019"),
       col=c(hex_codes[1],hex_codes[2],hex_codes[3],hex_codes[4],hex_codes[5]),
       pch=c(19,19,19,19,19),bty='n',pt.cex=3, cex=1.4,horiz=T)

dev.off()