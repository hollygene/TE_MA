setwd("~/Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/BedGraphs/D20") #sets your working directory 
#!/usr/bin/env Rscript
library(data.table) 

samples <- c("~/Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/BedGraphs/D20/HM-D20-10.depth") #here is where you put your filenames if you want to run a loop

depthtab <- fread("/Users/hollymcqueary/Dropbox/McQueary/Paradoxus_MA/BedGraphs/D20/HM-D20-10.depth", col.names = c("Chromosome", "Position", "Depth")) # this reads in your files into an array called args 
c = 0 # variable for chromosome
k = 0 
cstart <- 0 #variable for chromosome start
cend <- 0 #variable for chromosome end 
g = sub(".*depth_","",args[1]) 
#pdf(paste0(g,"_SW_SP.pdf"))
med = 2*median(depthtab$Depth, na.rm=TRUE) # this returns a value that is twice the median depth of

a <- unique(depthtab$Chromosome) # this returns a vector with any duplicate rows removed  
a[5:8] <- unique(depthtab$Chromosome)[7:10] # this is where you rearrange your chromosomes in case they aren't in the correct order
a[9] <- unique(depthtab$Chromosome)[5]
a[10:16] <- unique(depthtab$Chromosome)[11:17]
a[17] <- unique(depthtab$Chromosome)[6]
op <- par(mar=c(1,1,4,1)+0.1, oma = c(5,4,4,1)+0.1, mfcol=c(1,16)) # this sets the parameters for your plot 


for(i in a){
  
	if(grepl("chrM", i)){next} 
	chrm <- depthtab[depthtab$Chromosome == i,2:3]
  W <- 1:as.integer(max(depthtab$Position)/1000)*1000
  
  W <- append(W,max(W)+1000)
  meanDepth <- 0
  for(j in 1:length(W)){meanDepth[j] <- mean(depthtab$Depth[depthtab$Position>(W[j]-1000)&depthtab$Position<=W[j]])}
  
	par(mar=c(0,0,2,0))
	if(c%%4==0){color="black"}else if(c%%4==1){color="#E44848"}else if(c%%4==2){color="#9FCDF0"}else if(c%%4==3){color="#E0F09F"}
  
# col.chrm = as.factor(depthtab$Chromosome)	

# pdf("test.chrms.pdf")
	plot(W[1:max(W)], meanDepth[1:max(W)], col=color, ylab="Coverage", main="Coverage By Chromosome")
# gplot <- 
	
	# dev.off()
	
	# plot(W[1:max(W)], meanDepth[1:max(W)],pch=19,cex=0.2, ylim=c(0,med), axes = FALSE, col= color)
	title(main = paste0(sub(".*chr", "", i), "\n", median(depthtab$Depth, na.rm = TRUE), "x"), cex.main = 1)
	k = k + max(W)
	axis(side = 1, at = max(W), labels = k/1000, cex.axis = 1, pos=0)
	if(c==0){axis(side = 2, at = seq(0,med,(med/5)),pos=0)}
  #abline(h=mean(Depth),col="red")
  abline(h=median(depthtab$Depth, na.rm = TRUE), col="blue")
	abline(h=0, col = "black")
		# detach(chrm) # don't want this 
	c=c+1
}



title(main = g, xlab="genomic position (kbp)", ylab="mean depth", outer = TRUE, line = 2)
par(op)
#dev.off()
c=0
k=0
rm(list=ls())