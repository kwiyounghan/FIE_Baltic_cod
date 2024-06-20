#!/usr/bin/env Rscript --vanilla
## all the calculations/functions in this script are replicated from cvtkpy

## feed in a directory with five afreq files from plink2
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1){
  stop("I think you forgot arguments\n first arg=directory name, second arg=output prefix")
}
directory <- as.character(args[1])
out_pre <- as.character(args[2])

setwd(directory)
afreqFiles <- list.files(pattern=".afreq$")
# read in af matrix
t1_afreq <- read.table(afreqFiles[1],header=F)
t2_afreq <- read.table(afreqFiles[2],header=F)
t3_afreq <- read.table(afreqFiles[3],header=F)
t4_afreq <- read.table(afreqFiles[4],header=F)
t5_afreq <- read.table(afreqFiles[5],header=F)

# plot histograms of allele frequencies
png(paste0(out_pre,"_af_histogram.png"),width=800,height=700,res=72)
par(mfrow=c(2,3))
hist(t1_afreq$V5, freq=FALSE,breaks=100)
hist(t2_afreq$V5, freq=FALSE,breaks=100)
hist(t3_afreq$V5, freq=FALSE,breaks=100)
hist(t4_afreq$V5, freq=FALSE,breaks=100)
hist(t5_afreq$V5, freq=FALSE,breaks=100)
dev.off()

# combine af of each time points
df <- rbind(t1_afreq$V5,t2_afreq$V5,t3_afreq$V5,t4_afreq$V5,t5_afreq$V5)	#column: snp sites row: different time points
# AF change among two consequtive time points
df_diff <- diff(df); dim(df_diff)
# covariance of the subtracts (AF changes)
df_cov <- cov(t(df_diff))

# plot raw covariance values
png(paste0(out_pre,"_raw_covariance.png"),width=800,height=700,res=72)
plot(c(1:3),df_cov[1,2:4],col=rainbow(3)[1],type="b",ylim=c(min(df_cov),max(df_cov)))
points(c(2:3),df_cov[2,3:4],col=rainbow(3)[2],type="b")
points(c(3:3),df_cov[3,4:4],col=rainbow(3)[3],type="b")
abline(h=0,lty=2,lwd=0.5)
dev.off()

## bias correction
## need frequency matrix as well as number of individuals for each population
## calculate hets and correct for bias , need freq , depth and ndiploids matrices
freqs <- rbind(t1_afreq$V5,t2_afreq$V5,t3_afreq$V5,t4_afreq$V5,t5_afreq$V5)	#column: snp sites row: different time points

dep <- rbind(t1_afreq$V6,t2_afreq$V6,t3_afreq$V6,t4_afreq$V6,t5_afreq$V6)	#column: snp sites row: different time points

ndip = c(20,20,20,20,20) # number of samples per time point, simulation takes 20 iindividuals per time points

## caculate heterozygosity denominator for covariance values
hets <- 2* freqs * (1-freqs)
hets <- hets * dep / (dep-1) # depth correction
hets <- hets * 2 * ndip / (2*ndip - 1) # ndiploids correction
mean_het <- rowMeans(hets,na.rm=TRUE)
het_denom <- matrix(nrow=4,ncol=4)
for (i in 1:4){
	for (j in 1:4){
		if (i==j){het_denom[i,j] <- mean_het[i]}
		if (i < j){het_denom[i,j] <- mean_het[i]}
		if (i > j){het_denom[i,j] <- mean_het[j]}
			}}
het_denom <- het_denom / 2

## bias correction
ave_bias <- c(rep(0,dim(freqs)[1]))
dep_cor <- 1/dep
dip_cor <- 1/(2*ndip)
dip_cor <- dip_cor + 1/(2*dep*ndip) #ndip is muliplied for each column
ave_bias <- ave_bias + rowMeans( 0.5 * hets * (dip_cor + dep_cor), na.rm=TRUE )

# variance correction -- correcting variance = diagonal elements of covariance matrix
var_cor <- c(rep(0,dim(freqs)[1]-1))
var_cor = var_cor - ave_bias[1:length(ave_bias)-1] - ave_bias[2:length(ave_bias)]

# covariance correction
cov_cor <- c(rep(0,dim(freqs)[1]-2)) # the dimention of this correction is number of covariance values of k=1 in a covariance matrix = number of off-diagonal cells  , in my case 3, R*T-1 (when T+1 is number of your time points, R is number of replicates )
cov_cor <- cov_cor + ave_bias[2:(dim(freqs)[1]-1)]

# define df_cov
df_diff <- diff(freqs)
df_cov <- cov(t(df_diff))
# correct bias,
for (i in 1:nrow(df_cov)){
	for(j in 1:nrow(df_cov)){
		if(i == j){
			df_cov[i,j] <- df_cov[i,j] + var_cor[i]}
		if(i-j == -1){
			df_cov[i,j] <- df_cov[i,j] + cov_cor[i]}
		if(i-j == 1){
			df_cov[i,j] <- df_cov[i,j] + cov_cor[j]}
	}}
df_cov <- df_cov/ het_denom
write.table(df_cov,file=paste0(out_pre,"_covMatrix.txt"),quote=FALSE, sep="	")
## plot bias corrected covariance
png(paste0(out_pre,"_covariance.png"),width=800,height=700,res=72)
plot(c(1:3),df_cov[1,2:4],col=rainbow(3)[1],type="b",ylim=c(min(df_cov),max(df_cov)))
points(c(2:3),df_cov[2,3:4],col=rainbow(3)[2],type="b")
points(c(3:3),df_cov[3,4:4],col=rainbow(3)[3],type="b")
abline(h=0,lty=2,lwd=0.5)
dev.off()
