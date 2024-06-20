setwd("~/Desktop/WORK/Project_Fisheries_induced_evol/12_pi_pixy")
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
list.files(pattern="invariant_cod3")
#dxy fst pi

#convert pixy output to long format
source("pixy_to_long.r")
pixy_to_long()

#this put all output files in "output" folder into a long format
pixy_folder <- getwd()
#for all
#pixy_files <- list.files(pixy_folder,pattern="invariant_cod3", full.names = TRUE)[c(1,4)]
# exclude samples from pilot study with 100 bp read lengths
pixy_files <- list.files(pixy_folder,pattern="invariant_cod3", full.names = TRUE)[c(2,3)]


pixy_df <- pixy_to_long(pixy_files)
head(pixy_df,10)
nrow(pixy_df)
pixy_df <- pixy_df %>% filter(chromosome != "NC_002081.1") %>% filter(!grepl('NW',chromosome))
## plot pi from different years in one
## plot dxy and fst 1996 vs. in. two
## then combine these three plots into one
win=50000
quartz()
par(mfrow=c(2,1),oma=c(3,0,2,1),mar=c(1,4,1,1))
#assign values
years <- c(1996,2002,2008,2014,2019)
 for (i in 1:5){
	df <- pixy_df %>%
		filter(statistic == "avg_pi") %>%
		filter(pop1 == years[i]) %>% as.data.frame(.)
	df$cum_pos <- seq(from=win/2, by=win, length.out = nrow(df))/1000000
	assign(paste("los",i,sep="_"),loess.smooth(x= as.numeric(df$cum_pos), y=as.numeric(df$value),span=0.005,degree=0,cex=0.2, family="gaussian",evaluation=5000))
 }
 #this sets up the whole plot
# type='n' gives nothing on the plot area
plot(x= as.numeric(df$cum_pos), y=as.numeric(df$value),
     xaxs="i", type="n",cex=0.2,
     xlim=c(0,max(as.numeric(df$cum_pos))),
     ylim=c(0.000,0.015),
     xaxt="n",xlab="",ylab="")

# set background color to distinguish chromosomes
lim <- par("usr")
lg_names <- unique(df$chromosome)
for(j in 1:length(lg_names)){
  rect(min(df[df$chromosome==lg_names[j],'cum_pos']),lim[3]-1,
  max(df[df$chromosome==lg_names[j],'cum_pos']),lim[4]+1,
       col=ifelse(as.numeric(j)%%2==0,"snow2","snow1"),border="transparent")
}
title(line=2,ylab="Nucleotide Diversity", xlab="", main="")
#color
col <- c("#ff595e","#ffca3a","#8AC926","#1982c4","#6A4C93")
#plot each year
lines(x=los_1$x,y=los_1$y, lwd=1,col=col[1])
lines(x=los_2$x,y=los_2$y, lwd=1,col=col[2])
lines(x=los_3$x,y=los_3$y, lwd=1,col=col[3])
lines(x=los_4$x,y=los_4$y, lwd=1,col=col[4])
lines(x=los_5$x,y=los_5$y, lwd=1,col=col[5])

#legend
legend("topleft",xpd=T,
       legend=c("1996","2002","2008", "2014", "2019"),
       col=col,
       lty=c(1,1,1,1,1),
       lwd=c(4,4,4,4,4),
       horiz=T,bty='n',inset=c(0,-0.14))


##dxy
head(pixy_df)
unique(pixy_df$pop2)
df_dxy <- pixy_df %>%
	filter(statistic == "avg_dxy")

df <- df_dxy %>% filter(pop1== years[1]) %>% filter(pop2== years[5]) %>% as.data.frame(.)

head(df_dxy)
unique(df_dxy$pop2)

years <- c(2002,2008,2014,2019)
for (i in c(1,3,4){
	df <- df_dxy %>% filter(pop1== "1996") %>% filter(pop2== years[i]) %>% as.data.frame(.)
	df$cum_pos <- seq(from=win/2, by=win, length.out = nrow(df))/1000000
	assign(paste("los",i,sep="_"),loess.smooth(x= as.numeric(df$cum_pos) , y=as.numeric(df$value),span=0.005,degree=0,cex=0.2, family="gaussian",evaluation=5000))
}

df <- df_dxy %>% filter(pop1== "2008") %>% filter(pop2== "1996") %>% as.data.frame(.)
df$cum_pos <- seq(from=win/2, by=win, length.out = nrow(df))/1000000
assign(paste("los",2,sep="_"),loess.smooth(x= as.numeric(df$cum_pos) , y=as.numeric(df$value),span=0.005,degree=0,cex=0.2, family="gaussian",evaluation=5000))

#this sets up the whole plot
# type='n' gives nothing on the plot area
plot(x= as.numeric(df$cum_pos), y=as.numeric(df$value),
     xaxs="i", type="n",cex=0.2,
     xlim=c(0,max(as.numeric(df$cum_pos))),
     ylim=c(0.00,0.015),
     xaxt="n",xlab="",ylab="")

# set background color to distinguish chromosomes
lim <- par("usr")
lg_names <- unique(df$chromosome)
for(j in 1:length(lg_names)){
  rect(min(df[df$chromosome==lg_names[j],'cum_pos']),lim[3]-1,
  max(df[df$chromosome==lg_names[j],'cum_pos']),lim[4]+1,
       col=ifelse(as.numeric(j)%%2==0,"snow2","snow1"),border="transparent")
}
#titles and labels
title(line=2,ylab="Dxy", xlab="", main="")
col <- c("#ffca3a","#8AC926","#1982c4","#6A4C93")
lines(x=los_1$x,y=los_1$y, lwd=1,col=col[1])
lines(x=los_2$x,y=los_2$y, lwd=1,col=col[2])
lines(x=los_3$x,y=los_3$y, lwd=1,col=col[3])
lines(x=los_4$x,y=los_4$y, lwd=1,col=col[4])
#legend
legend("topleft",xpd=T,
       legend=c("1996 vs 2002","1996 vs 2008", "1996 vs 2014", "1996 vs 2019"),
       col=col,
       lty=c(1,1,1,1),
       lwd=c(4,4,4,4),
       horiz=T,bty='n',inset=c(0,-0.14))

#for ticks on x axis
lg_names <- paste0("LG",seq(1,23))

xtick <- c()
for (j in unique(df[,3])){
  xtick <- c(xtick,mean(df[df$chromosome==j,8]))}
axis(side=1,padj=0, at=xtick,labels =lg_names, tick=FALSE,
     cex.axis=0.7,font=2, mgp=c(3,0.3,0.1))

mtext("Chromosome",side=1,line=1.8,cex=1.4,outer=T)

quartz.save("pi_dxy_cod3.png",dpi=2000)


## calculate whole genome nucleotide diversity
## equation based on pixy website https://pixy.readthedocs.io/en/latest/output.html
## (window 1 count_diffs + window 2 count_diffs) / (window 1 comparisons + window 2 comparisons)

##calculate genomewide dxy pi
## pi
inp <- read.table(list.files(pattern=".txt")[3], sep="\t", header=T)#pi no pilot
head(inp)

wg_pi <- c()
for (i in c(1996,2002,2008,2014,2019)){
	df <- inp %>% filter(pop == i) %>% as.data.frame(.)
wg_pi <- c(wg_pi,sum(df$count_diffs,na.rm=T)/sum(df$count_comparisons,na.rm=T))

}
wg_pi # [1] 0.007174294 0.007402100 0.007779840 0.007602320 0.007634114

# dxy
inp <- read.table(list.files(pattern=".txt")[2], sep="\t", header=T)#dxy no pilot
head(inp)

wg_dxy <- c()
i=2008
for (i in c(2002,2008,2014,2019)){
	df <- inp %>% filter(pop1 == 1996) %>% filter(pop2== i) %>% as.data.frame(.)
wg_dxy <- c(wg_dxy,sum(df$count_diffs,na.rm=T)/sum(df$count_comparisons,na.rm=T))

}

wg_dxy # [1] 0.007242150         NaN 0.007343077 0.007344347
## the order of pop1 and pop2 not consistent, so
df <- inp %>% filter(pop1 == 2008) %>% filter(pop2== 1996) %>% as.data.frame(.)
wg_dxy <- c(wg_dxy,sum(df$count_diffs,na.rm=T)/sum(df$count_comparisons,na.rm=T))
wg_dxy #[1] 0.007242150         NaN 0.007343077 0.007344347 0.007429637
# 1996vs2002	1996vs2008	1996vs2014	1996vs2019
#0.007242150	0.007429637	0.007343077 0.007344347
