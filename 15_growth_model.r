##modified script from Hoehne et al. 2019
setwd("~/GrowthCurve")
remove(list=ls())
# Use the same colors in all graphs!
colors<- c("#ff595e","#ffca3a","#8AC926","#1982c4","#6A4C93")

#----------------------------
# 1. Install and load libraries
#----------------------------
#install.packages("R2jags")  # If using JAGS for the first time, you need to download JAGS 4.3.0, install and then use this line to load into R
#I installed the JAGS then had to restart R to load the library after installing R2jags
library(R2jags)
library(bayesplot)

#----------------------------
# 2. Prepare input data
#----------------------------
## Load data
df <- read.csv("table1_FIE_metadata.csv", sep="\t")

## remove samples that are outlier in phenotype (measurement error) and WBC individuals
rm_samples <- read.csv(file="~/cod3/11_pca/outlier_sample_WBC_disguising_EBC.csv",sep="\t")
dim(rm_samples)
df <- df[!(df$SeqSample %in% rm_samples$SeqSample),]
## 154 samples because of the samples with low sequencing depth are included. I will include these samples for phenotype analysis, but exclude in genotype analysis.
dim(df) # 154

# Routine to prepare data
n.fish <- dim(df)[1]

# Transposing the "horizontal" input structure of annulli readings to a "vertical" nested structure that can be read by the JAGS model
length <- NULL
fish.length <- NULL
ID.fish <- NULL
age <- NULL
ID.bin <- rep(NA,n.fish)
otolith.length <- NULL
bin.list <- as.character(unique(df$bins))
n.bin <- length(bin.list)
# transform data
df[,8:16]
R <- 8 #column at which the reading starts
RE <- 16 #column at which the reading ends

for (i in 1:n.fish){
  temp <- df[i,R:RE]
  range <- range(which(!is.na(temp)))
  from <- (R+range[1]-1)
  till <- from+range[2]-1
  temp2 <- as.numeric(df[i,from:till])
  length <- c(length,temp2)
  age <- c(age,c(1:length(temp2)))
  ID.fish <- c(ID.fish,rep(i,length(temp2)))
  temp3 <- which(bin.list==as.character(df$bins[i]))
  ID.bin[i] <- temp3
  fish.length <- c(fish.length,df$length[i])
  otolith.length <- c(otolith.length,df[i,'otoLength_micron'])
}

# Packing data and z-transformation
data.jags = list(
  n.obs = length(length),     # number of observations = number of measurements of otolith rings
  age = age,                  # age from otolith readings
  y = length,                 # response variable # here distance to otolith rings
  ID.fish = ID.fish,
  n.fish = n.fish,
  ID.bin = ID.bin,            # bin = grouping by catch year, catch year 1996-1998 are grouped in bin1 together
  n.bin = length(unique(df$bins)),
  fish.length = fish.length,  # fish body length at catch
  otolith.length = otolith.length # otolith radii
)



#----------------------------
# 3.von Bertalanffy growth model (hierarchical)
#----------------------------

sink("m4_norm_final.txt")
cat("model {# first line
    for (i in 1:n.obs) {
      y[i] ~ dnorm(OT_radius_hat[i],tol) #here tol being the observation error
              ## this means observed y values fall on the normal distribution with 'true' value as a mean the observation error
      OT_radius_hat[i] <- linf.i[ID.fish[i]]*(1-exp(-k.i[ID.fish[i]]*(age[i]-t0.i[ID.fish[i]])))
    }
    # fish level
    for (i in 1:n.fish){
      linf.i[i] ~ dnorm(linf.bin[ID.bin[i]] ,1/(sd.linf.bin[ID.bin[i]])^2) T(0,)
      k.i[i] ~ dnorm(k.bin[ID.bin[i]] , 1/(sd.k.bin[ID.bin[i]])^2) T(0,)
      t0.i[i] ~ dnorm(t0.bin[ID.bin[i]], 1/(sd.t0.bin[ID.bin[i]])^2) T(-1,1)
    }

    for (i in 1:n.bin){
      linf.bin[i] ~ dgamma(linf.all*sd.linf,sd.linf)
      sd.linf.bin[i] ~ dt(0,1,1)T(0,)

      k.bin[i] ~ dgamma(k.all*sd.k,sd.k)
      sd.k.bin[i] ~  dt(0,1,1)T(0,)

      t0.bin[i] ~ dgamma(t0.all*sd.t0,sd.t0)
      sd.t0.bin[i] ~ dt(0,1,1)T(0,)
     }

    # PRIOR PARAMETERS
    k.all ~ dgamma(0.001,0.001)
    linf.all ~ dgamma(0.001,0.001)
    t0.all ~ dgamma(0.001,0.001)
    sd.t0 ~ dgamma(0.001,0.001)
    sd.linf ~ dgamma(0.001,0.001)
    sd.k ~ dgamma(0.001,0.001)

    tol ~ dgamma(0.001,0.001)


    # # DERIVED for individual fish otolith
    for (j in 1:n.fish){ # size at age (1 to 8 years) for each bin (time points)
    for (i in 1:8){
    y.pre[i,j] <- linf.i[j] * (1-exp(-k.i[j] * (i-t0.i[j]))) ## level on otolith
    }
    }

    }",fill = TRUE)
sink()



# 3.2: Initial values
inits <- list(
  list(k.all=0.2,linf.all=15),
  list(k.all=0.04,linf.all=10),
  list(k.all=0.1,linf.all=30)
)

# 3.3: Settings
# Parameters monitored
params <- c("k.bin","sd.k.bin",  "k.i"  ,  "k.all", "sd.k",              # k
            "linf.bin","sd.linf.bin","linf.i",  "linf.all","sd.linf" ,   # Linf
            "t0.bin","sd.t0.bin","t0.i", "t0.all","sd.t0",               # t0
            "y.pre","tol")                                               # 1/var at observation level

# 3.4: MCMC settings, initialization
ni <- 100000   #number of total iterations per chain (including burn in; default: 2000)
nt <- 10       #thinning rate
nb <- 10000    #length of burn in, i.e. number of iterations to discard at the beginning.
nc <- 3        #number of Markov chains (default: 3)

# 3.5: Running and updating
out <- jags(data.jags, inits, params, "m4_norm_final.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb)

#out=update(out,n.iter = 100000, n.thin = 10)  # Run a few iterations just to make sure the model is working
save(out,file="results_m4_norm_allSamples154_final.RData") # saves the BUGS output

#load("results_m4_norm_allSamples154_final.RData")
#out$model

#----------------------------
# 4.Check and plot results
#----------------------------

# 4.1: visual examination of model convergence
out.mcmc <- as.mcmc(out)
mcmc_combo(out.mcmc,combo=c("dens","trace"),pars = c("k.all","linf.all","t0.all","tol"))
mcmc_intervals(out.mcmc,regex_pars = c("k.i*"))
mcmc_intervals(out.mcmc,regex_pars = c("linf.i*"))

# 4.2 get a summary of statistics (e.g., to also inspect R.hat values).
# Note that the DIC estimate is dependent on the surveyed parameters
print(out$BUGSoutput$summary ,digits=3)
hist(out$BUGSoutput$summary[,'Rhat'])
out$BUGSoutput$DIC       # get DIC estimate (Note that the DIC is never a steady point estimate. To ensure robustness of model selection, we recommend to run the analysis 3 times and record the lowest DIC estimate for a model where all relevant parameters have converged. If available, use additional model evaluation citeria (e.g., out-of-sample prediciton))

out$BUGSoutput$median$k.all
out$BUGSoutput$median$linf.all
out$BUGSoutput$median$t0.all
out$BUGSoutput$median$k.bin
out$BUGSoutput$median$linf.bin
out$BUGSoutput$median$t0.bin

# 4.3 : Boxplots of estimated individual Linf and K based on catch year
# Linf
boxplot(out$BUGSoutput$median$linf.i~data.jags$ID.bin,ylab=expression(paste("L"[infinity]," otolith", " (mm)")), col=colors, cex.lab=1.3, cex.axis=1.07, font = 2, font.lab = 2, xaxt="n", mgp =c(1,1,0),xlab="")
title(xlab="Bins", mgp =c(1,1,0), font=2, cex.lab=1.4)
legend("topright",legend=c("1996","2002","2008","2014","2019"),col=colors,pch=15,cex=1,bg="transparent",box.lty=0,y.intersp=0.8)

# K
boxplot(out$BUGSoutput$median$k.i~data.jags$ID.bin,ylab=expression(paste("Growth coefficient",italic(" K"))),  col=colors, cex.lab=1.15, cex.axis=1.15, font = 2, font.lab = 2, xaxt="n", mgp =c(1,1,0),xlab="")
title(xlab="Bins", mgp =c(1,1,0), font=2, cex.lab=1.4)
legend("topleft",legend=c("1996","2002","2008","2014","2019"),col=colors,pch=15,cex=1,bg="transparent",box.lty=0,y.intersp=0.8)

# Phi
phi <- log(out$BUGSoutput$median$k.i,base=10) + 2 * log(out$BUGSoutput$median$linf.i,base=10)
phi_bin <- log(out$BUGSoutput$median$k.bin,base=10) + 2 * log(out$BUGSoutput$median$linf.bin,base=10)

boxplot(phi~data.jags$ID.bin,ylab=expression(paste("Growth performance index",italic(" phi"))),  col=colors, cex.lab=1.15, cex.axis=1.15, font = 2, font.lab = 2, names=c("1996","2002","2008","2014","2019"), mgp =c(1,1,0),xlab="Catch year")
legend("topright",legend=c("1996","2002","2008","2014","2019"),col=colors,pch=15,cex=1,bg="transparent",box.lty=0,y.intersp=0.8)


# 4.4 : plotting predicted VBGMs
gf <- function(l,k,x,t) {l*(1-exp(-k*(x-t)))}
age.sim=c(1:8)
plot(0,0,xlim=c(1,8),ylim=c(0,7),ylab="Otolith radius length (mm)",xlab="Age (years)",mgp = c(2, 0.5, 0),cex.lab=1.2)
for (i in 1:data.jags$n.bin){
  temp_k=quantile(out$BUGSoutput$sims.list$k.bin[,i],c(0.025,0.5,0.975))
  temp_l=quantile(out$BUGSoutput$sims.list$linf.bin[,i],c(0.025,0.5,0.975))
  temp_t=quantile(out$BUGSoutput$sims.list$t0.bin[,i],c(0.025,0.5,0.975))
  lines(age.sim,gf(temp_l[2],temp_k[2],age.sim,temp_t[2]),col=colors[i],lwd=3)
  #lines(age.sim,gf(temp_l[1],temp_k[1],age.sim,temp_t[1]),col=colors[i],lty=3)
  #lines(age.sim,gf(temp_l[3],temp_k[3],age.sim,temp_t[3]),col=colors[i],lty=3)
}
points(data.jags$age,data.jags$y, pch= 19, col = colors[data.jags$ID.bin[data.jags$ID.fish]])
legend("topleft",legend=c("1996","2002","2008","2014","2019"),col=colors,lwd=2,cex=1,bg="transparent",box.lty=0,y.intersp=0.8)


# ## Examine credible interval for each catch year
# par(mfrow=c(3,2))
# for (i in 1:n.bin){
#   plot(0,0,xlim=c(1,7),ylim=c(0,7),ylab="Otolith radius length (mm)",xlab="Age (years)",mgp = c(2, 0.5, 0),cex.lab=1.2)
#   temp_k=quantile(out$BUGSoutput$sims.list$k.bin[,i],c(0.025,0.5,0.975))
#   temp_l=quantile(out$BUGSoutput$sims.list$linf.bin[,i],c(0.025,0.5,0.975))
#   temp_t=quantile(out$BUGSoutput$sims.list$t0.bin[,i],c(0.025,0.5,0.975))
#   lines(age.sim,gf(temp_l[2],temp_k[2],age.sim,temp_t[2]),col=colors[i],lwd=3)
#   lines(age.sim,gf(temp_l[1],temp_k[1],age.sim,temp_t[1]),col=colors[i],lty=3)
#   lines(age.sim,gf(temp_l[3],temp_k[3],age.sim,temp_t[3]),col=colors[i],lty=3)
#
#   points(data.jags$age[data.jags$ID.bin[data.jags$ID.fish]==i],data.jags$y[data.jags$ID.bin[data.jags$ID.fish]==i], pch= 1, col = colors[i])
#   legend("topleft",legend=c("1996","2002","2008","2014","2019")[i],col=colors[i],lwd=2,cex=2,bg="transparent",box.lty=0,y.intersp=0.6)
# }

all_param <- rbind(c(out$BUGSoutput$median$linf.all,out$BUGSoutput$median$k.all,out$BUGSoutput$median$t0.all),cbind(out$BUGSoutput$median$linf.bin,out$BUGSoutput$median$k.bin,out$BUGSoutput$median$t0.bin),cbind(out$BUGSoutput$median$linf.i,out$BUGSoutput$median$k.i,out$BUGSoutput$median$t0.i))
dim(all_param)
colnames(all_param) <- c("linf","k","t0")
rownames(all_param) <- c("all","bin1","bin2","bin3","bin4","bin5",df$fishID)
head(all_param)
write.table(all_param, file="these_table3_growth_parameters.csv",sep="\t",quote=FALSE)

#----------------------------
# 5.Fish length at otolith relationship, Back-calculation (Hussy et al 2018)
#----------------------------
Oa <- length## otolith radii
Lc <- fish.length[ID.fish]##fish body length at catch
Oc <- otolith.length[ID.fish]##total length of otolith (or the radius?) at catch
L0 <- 4.3 # in mm at age0
O0 <- 0.01 # in mm at age0

La <- Lc + (Oa - Oc)*(Lc-L0)*(Oc-O0)^(-1)
La

# 5.1 Visualize
plot(Oa,La,xlab="otolith length (mm)",ylab="fish length (mm)")
points(Oc,Lc,col="red")
legend("topleft",legend=c("at age","at catch"),col=c("black","red"),pch=c(1,1))

# 5.2 Linf - back-calculated fish length
Lc <- fish.length
Oa <- out$BUGSoutput$median$linf.i
Oc <- otolith.length
La <- Lc + (Oa - Oc)*(Lc-L0)*(Oc-O0)^(-1)

boxplot(La~data.jags$ID.bin,ylab=expression(paste("L"[infinity]," otolith", " (mm)")), col=colors, cex.lab=1.3, cex.axis=1.07, font = 2, font.lab = 2, xaxt="n", mgp =c(1,1,0),xlab="")
title(xlab="Bins", mgp =c(1,1,0), font=2, cex.lab=1.4)
legend("topright",legend=c("1996","2002","2008","2014","2019"),col=colors,pch=15,cex=1,bg="transparent",box.lty=0,y.intersp=0.8)

la_val <- boxplot(La~data.jags$ID.bin)
la_val$stats[3,]
#[1] 1150.0207  943.9672  911.5783  622.6696  539.8754

#----------------------------
# 5. Juvenile growth, check 1-year olds
#----------------------------
# 5.1 Calculate means and sd
aggregate(df$distance_min_1, list(df$bins), FUN=sd)
####
# Group.1         x
# 1    bin1 0.1948053
# 2    bin2 0.1621913
# 3    bin3 0.1745603
# 4    bin4 0.1097144
# 5    bin5 0.1276450
####
aggregate(df$distance_min_1, list(df$bins), FUN=mean)
####
# Group.1         x
# 1    bin1 0.6267397
# 2    bin2 0.6284471
# 3    bin3 0.5735613
# 4    bin4 0.6710356
# 5    bin5 0.6396841
####

# 5.2 Test the variances of each bin
bartlett.test(df$distance_min_1~df$bins)
####
## Bartlett's K-squared = 10.818, df = 4, p-value = 0.02868
####

# 5.3 Back-calculate the fish length for more interpretable results and plot
Oa <- df$distance_min_1## otolith radii
Lc <- fish.length##fish body length at catch
Oc <- otolith.length##total length of otolith (or the radius?) at catch
L0 <- 4.3 # in mm at age0
O0 <- 0.01 # in mm at age0
La <- Lc + (Oa - Oc)*(Lc-L0)*(Oc-O0)^(-1)

boxplot(La~df$bins,ylab="back-calculated fish length",xlab="catch year",names=c("1996-1998","2002","2008","2014","2019"),col=colors)
legend("topright",legend=c("1996-1998","2002","2008","2014","2019"),col=colors,pch=15,cex=1,bg="transparent",box.lty=0,y.intersp=0.8)
mtext("Bartlett's test p-value = 0.02868")



#----------------------------
# 6. Plot four parameters in one figure
#----------------------------
png(filename="growth_param_all.png",width=4000,height=4000,res=300, units='px')
par(mfrow=c(2,2),oma=c(2,2,5,2),mar=c(2,5,2,2))

# Linf
boxplot(out$BUGSoutput$median$linf.i~data.jags$ID.bin,ylab=expression(paste("L"[infinity]," otolith", " (mm)")), col=colors, cex.lab=1.5, cex.axis=1.5, font = 2, font.lab = 2, xaxt="n", mgp =c(1,1,0),xlab="")
title(main="(a)",line=1,adj=0,cex=4)
boxplot(out$BUGSoutput$median$k.i~data.jags$ID.bin,ylab=expression(paste("Growth coefficient",italic(" k"))),  col=colors, cex.lab=1.5, cex.axis=1.5, font = 2, font.lab = 2, xaxt="n", mgp =c(1,1,0),xlab="")
title(main="(b)",line=1,adj=0,cex=4)
boxplot(phi~data.jags$ID.bin,ylab=expression(paste("Growth performance index ",phi)),  col=colors, cex.lab=1.5, cex.axis=1.5, font = 2, font.lab = 2, xaxt="n", mgp =c(1,1,0),xlab="")
title(main="(c)",line=1,adj=0,cex=4)
boxplot(La~df$bins,ylab="back-calculated fish length (mm)",col=colors,cex.lab=1.5, cex.axis=1.5, font = 2, font.lab = 2, xaxt="n", mgp =c(1,1,0),xlab="")
title(main="(d)",line=1,adj=0,cex=4)
par(fig = c(0, 1, 0, 1), oma = c(0, 1, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("top",horiz=T,cex=2,xpd=T,legend=c("1996-1998","2002","2008","2014","2019"),col=colors,pch=15,bg="transparent",box.lty=1,y.intersp=0.8,pt.cex=3)
mtext(text="catch year",side=1, cex=2,line=-2)
dev.off()



#----------------------------
# 7. Check Residuals
#----------------------------
n_pre=length(1:7)
residuals=array(NA,dim=c(n.fish,n_pre))
dim(out$BUGSoutput$sims.list$y.pre)# 27000     8   154

out$BUGSoutput$sims.list$y.pre[,8,1]

for (i in 1:n.fish){
  print(i)
  length_hat=array(NA,c(n_pre,3))
  for (j in 1:n_pre){
    temp=out$BUGSoutput$sims.list$y.pre[,j,i]
    length_hat[j,]=quantile(temp,c(0.025,0.5,0.975))
  }
  temp=which(data.jags$ID.fish == i)
  for (j in 1:length(temp)){
    residuals[i,data.jags$age[temp][j]]=(length_hat[j,2]-data.jags$y[temp][j])/data.jags$y[temp][j]
  }
}
dim(residuals)## n.fish = 154 age=7

png(filename="growth_residuals.png",width=2000,height=2000,res=300, units='px')
par(mar=c(4,4,3,2))
boxplot(residuals, xlab="age (years)", names=c(1:7),ylab="Residuals (otolith length)",outline=F)
abline(h=0)
dev.off()

sum(data.jags$age==1);sum(data.jags$age==2);sum(data.jags$age==3);sum(data.jags$age==4);sum(data.jags$age==5);sum(data.jags$age==6);sum(data.jags$age==7)
# [1] 154
# [1] 154
# [1] 142
# [1] 83
# [1] 38
# [1] 14
# [1] 1
