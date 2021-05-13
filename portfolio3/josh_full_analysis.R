set.seed(1982)
library(pacman)
p_load(R2jags)
p_load(polspline)

groupSize <- 4
ntrials <- 10
pi <- 1.4
ntokens <- 20
vals <- seq(1,21,1) #possible values to contribute - from 0 to 20 tokens



rawDat <- read.csv("HerrmannThoeniGaechterDATA.csv", skip = 3) # Public goods game

# create variable for GINI. Data from 
# http://hdr.undp.org/sites/default/files/reports/269/hdr_2009_en_complete.pdf,

rawDat$gini <- c()
rawDat$gini[rawDat$city=="Melbourne"]=34.3
rawDat$gini[rawDat$city=="Minsk"]=25.3
rawDat$gini[rawDat$city=="Chengdu"]=38.5
rawDat$gini[rawDat$city=="Copenhagen"]=28.7
rawDat$gini[rawDat$city=="Bonn"]=31.9
rawDat$gini[rawDat$city=="Athens"]=34.4
rawDat$gini[rawDat$city=="Seoul"]=31.6
rawDat$gini[rawDat$city=="Samara"]=37.5
rawDat$gini[rawDat$city=="Zurich"]=32.7
rawDat$gini[rawDat$city=="St. Gallen"]=32.7
rawDat$gini[rawDat$city=="Istanbul"]=41.9
rawDat$gini[rawDat$city=="Nottingham"]=34.8
rawDat$gini[rawDat$city=="Dnipropetrovs'k"]=26.1
rawDat$gini[rawDat$city=="Boston"]=41.1

# extract every third line - data file has lines representing others responses and we don't need that
redDat <- rawDat[seq(1,length(rawDat$sessionid),3),]

group_names <- unique(redDat$groupid)
ngroups <- length(group_names)

# THIS WILL REMOVE SUBJECTS WITH MISSING DATA IN NO PUNISHMENT CONDITION
ngroups <- 269

subject_names <- unique(redDat$subjectid)
nsubjects <- length(subject_names)

# data for no punishment condition #
c_no_punish <- array(0,c(groupSize,ntrials,ngroups))
Ga_no_punish <- array(0,c(ntrials,ngroups))
Gc_no_punish <- array(0,c(groupSize,ntrials,ngroups))
missing <- array(0,ngroups)



for (g in 1:ngroups) {
  c_no_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][1:10],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][11:20],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][21:30],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][31:40])
  
  Ga_no_punish[,g] <- colMeans(c_no_punish[,,g])
  
  missing[g] <- is.na(c_no_punish[1,1,g])
  
  for (s in 1:groupSize) {
    Gc_no_punish[,,g] <- colSums(c_no_punish[-s,,g])
  }
}

# data for punishment condition #
c_punish <- array(0,c(groupSize,ntrials,ngroups))
Ga_punish <- array(0,c(ntrials,ngroups))
Gc_punish <- array(0,c(groupSize,ntrials,ngroups))

for (g in 1:ngroups) {
  c_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][1:10],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][11:20],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][21:30],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][31:40])
  
  Ga_punish[,g] <- colMeans(c_punish[,,g])
  
  for (s in 1:groupSize) {
    Gc_punish[,,g] <- colSums(c_punish[-s,,g])
  }
}

# compile data from each condition into 4D matrix
c <- array(0,c(groupSize,ntrials,ngroups,2))
c[,,,1] <- c_no_punish
c[,,,2] <- c_punish

Ga <- array(0,c(ntrials,ngroups,2))
Ga[,,1] <- Ga_no_punish
Ga[,,2] <- Ga_punish

Gc <- array(0,c(groupSize,ntrials,ngroups,2))
Gc[,,,1] <- Gc_no_punish
Gc[,,,2] <- Gc_punish

c_choice_index <- c

Gini <- array(0,c(ngroups))
for (g in 1:ngroups) {
  Gini[g] <- mean(redDat$gini[redDat$groupid==group_names[g]&redDat$p=="P-experiment"])
}

Ga_punish <- Ga_punish[,!is.na(Gini)]
Ga_no_punish <- Ga_no_punish[,!is.na(Gini)]

c <- c[,,!is.na(Gini),]
Ga <- Ga[,!is.na(Gini),]
Gc <- Gc[,,!is.na(Gini),]
Gini <- Gini[!is.na(Gini)]

#redefine number of groups after removing those without civic scores
ngroups <- length(Gini)

#-------------------  Apply decay model to data ------------------------------------------------------

data <- list("ntrials","ngroups","groupSize","Gini","c") #data inputted into jags
params <- c("c0","gamma","beta0_c0","betaGini_c0","beta0_gamma","betaGini_gamma")

setwd("C:/Users/au199986/Dropbox/Inequality")

# - run jags code
decay.samples <- jags(data, inits=NULL, params,
                      model.file ="decay.txt",
                      n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)

# Bayes factors and credible intervals for effect of Gini on initial contribution
prior <- dnorm(0,1)
fit.posterior <- logspline(decay.samples$BUGSoutput$sims.list$betaGini_c0)
posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
decay.betaGini_c0.BF <- prior/posterior

decay.betaGini_c0.CIL <- qlogspline(0.025,fit.posterior) #2.5% CI
decay.betaGini_c0.CIH <- qlogspline(0.975,fit.posterior) #97.5% CI

# Bayes factors and credible intervals for effect of Gini on decay in contribution
prior <- dnorm(0,1)
fit.posterior <- logspline(decay.samples$BUGSoutput$sims.list$betaGini_gamma)
posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point = 0
decay.betaGini_gamma.BF <- prior/posterior

decay.betaGini_gamma.CIL <- qlogspline(0.025,fit.posterior) #2.5% CI
decay.betaGini_gamma.CIH <- qlogspline(0.975,fit.posterior) #97.5% CI

#-------------------  Apply CC model to data ------------------------------------------------------

#standardise Gini score for analysis - we need to do this because model is more complex and won't converge otherwise
Gini_unst <- Gini #used for plotting unstandardised
Gini <- (Gini-mean(Gini))/sd(Gini)

data <- list("groupSize", "ngroups", "ntrials", "ntokens","vals","c","Ga","Gini") #data inputted into jags
params <- c("beta0_pbeta","betaGini_pbeta",
            "mu_pbeta_probit") #parameters we'll track in jags

# - run jags code
CC.samples <- jags(data, inits=NULL, params,
                   model.file ="CC_Gini.txt",
                   n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)

prior <- dnorm(0,1)
fit.posterior <- logspline(CC.samples$BUGSoutput$sims.list$beta0)
posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
CC.beta0.BF <- prior/posterior

CC.beta0.CIL <- qlogspline(0.025,fit.posterior) #2.5% CI
CC.beta0.HIL <- qlogspline(0.975,fit.posterior) #97.5% CI

prior <- dnorm(0,1)
fit.posterior <- logspline(CC.samples$BUGSoutput$sims.list$betaGini_pbeta)
posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
CC.betaGini.BF <- prior/posterior

CC.betaGini.CIL <- qlogspline(0.025,fit.posterior) #2.5% CI
CC.betaGini.CIH <- qlogspline(0.975,fit.posterior) #97.5% CI

#---------------- plots ---------------------------------------

bmp(file="decay_plot.bmp",
    width=13, height=13, units="in", res=100)

layout(rbind(c(1,1,3,3,4,4), 
             c(1,1,3,3,5,5),
             c(2,2,6,6,7,7), 
             c(2,2,6,6,8,8)))

par(cex=1.1)

###############################################################
# Plot of data - actual contributions by Gini coefficient
# all contributions overplotted with low alpha - shows decay trend

contrib <- cbind(c[1,,,1],c[2,,,1],c[3,,,1],c[4,,,1]) # treat all subjects within groups independently
gini_array <- rep(Gini_unst,4) # reshape gini array so it matches groups

plot(contrib[,gini_array<quantile(gini_array)[2]][,1],
     col=rgb(0,0,0,.02),pch=15,xlim = c(0,10),ylim = c(-2,22),
     cex=2,xlab = "Trial",ylab = "Contribution",main="A: Lower Quartile: Gini < 31.6",frame=FALSE)
for (i in 1:sum(gini_array<quantile(gini_array)[2])){
  points(contrib[,gini_array<quantile(gini_array)[2]][,i],col=rgb(0,0,0,.02),pch=15,cex=2)
}

plot(contrib[,gini_array>quantile(gini_array)[4]][,1],
     col=rgb(0,0,0,.02),pch=15,xlim = c(0,10),ylim = c(-2,22),
     cex=2,xlab = "Trial",ylab = "Contribution",main="B: Upper Quartile: Gini > 37.5",frame=FALSE)
for (i in 1:sum(gini_array>quantile(gini_array)[4])){
  points(contrib[,gini_array>quantile(gini_array)[4]][,i],col=rgb(0,0,0,.02),pch=15,cex=2)
}

###############################################################
# Plot of posterior decay models - lowest and highest Gini scores
t <- seq(1,10,1)

c0 <- mean(decay.samples$BUGSoutput$sims.list$beta0_c0)+(min(Gini_unst)*mean(decay.samples$BUGSoutput$sims.list$betaGini_c0))
gamma <- mean(decay.samples$BUGSoutput$sims.list$beta0_gamma)+min(Gini_unst)*mean(decay.samples$BUGSoutput$sims.list$betaGini_gamma)
Y <- c0 * exp(-gamma*t)
plot(Y,type='l',frame=FALSE,lwd=2,ylim=c(0,16),
     cex=2,col="dark blue",xlab = "Trial",ylab = "Contribution",main="C: Posterior Decay Models")
legend (x = 4, y = 15, legend = c("Lowest Gini","Highest Gini"), 
        col = c("dark blue","red"), bty = "n", lwd = 2)

c0 <- mean(decay.samples$BUGSoutput$sims.list$beta0_c0)+(max(Gini_unst)*mean(decay.samples$BUGSoutput$sims.list$betaGini_c0))
gamma <- mean(decay.samples$BUGSoutput$sims.list$beta0_gamma)+max(Gini_unst)*mean(decay.samples$BUGSoutput$sims.list$betaGini_gamma)
Y <- c0 * exp(-gamma*t)
lines(Y,type='l',lwd=2,col="red")

plot(density(decay.samples$BUGSoutput$sims.list$betaGini_c0),frame=FALSE,lwd=2,ylim=c(0,20),
     cex=2,xlab = "Parameter Estimate",ylab = "Density",main="D: Gini Effect: Initial Contribution")

plot(density(decay.samples$BUGSoutput$sims.list$betaGini_gamma),frame=FALSE,lwd=2,
     cex=2,xlab = "Parameter Estimate",ylab = "Density",main="E: Gini Effect: Contribution Decay")

###############################################################
# Plot of posterior contribution preference model - lowest and highest Gini scores
X <- seq(1,20,1)

beta0_pbeta <- CC.samples$BUGSoutput$sims.list$beta0_pbeta  
betaGini_pbeta <- CC.samples$BUGSoutput$sims.list$betaGini_pbeta

std_pbeta_low <- pnorm(beta0_pbeta+(min(Gini)*betaGini_pbeta))
std_pbeta_high <- pnorm(beta0_pbeta+(max(Gini)*betaGini_pbeta))

Y <- X*mean(std_pbeta_low)
plot(Y,type='l',frame=FALSE,lwd=2,ylim=c(0,20),xlim=c(-1,20),
     cex=2,col="dark blue",xlab = "Believed Contribution",ylab = "Preferred Contribution",main="F: Posterior Preference Models")
legend (x = 0, y = 19, legend = c("Lowest Gini","Highest Gini"), 
        col = c("dark blue","red"), bty = "n", lwd = 2)

Y <- X*mean(std_pbeta_high)
lines(Y,type='l',lwd=2,col="red")

plot(density(std_pbeta_low),frame=FALSE,lwd=2,
     cex=2,xlab = "Parameter Estimate",ylab = "Density",main="G: Slope Lowest Gini")

plot(density(std_pbeta_high),frame=FALSE,lwd=2,
     cex=2,xlab = "Parameter Estimate",ylab = "Density",main="H: Slope Highest Gini")

dev.off()

#--------------- winnings -------------------------------------

winnings <- array(0,c(ngroups))
for (g in 1:ngroups) {
  winnings[g] <- sum(colSums(c[,,g,1])*pi)
}

cor.test(Gini_unst,winnings)