#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Reproducible Results for:
##  Anderegg et al. 
## "Aridity drives coordinated trait shifts but not decreased trait variance across the geographic range of eight Australian trees"
## In review: New Phytologist
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This document reproduces the analyses and figures presented in Anderegg et al. using the data currently available in the 'trait_data' directory in this github repository
# For questions or to report bugs, please contact Leander Anderegg: leanderegg@gmail.com
# last updated: 03 June 2020

  
# set working directory if needed:
#setwd("")
# create a directory for generated results
results_dirname <- "Results_20200603" # "EucTraits_Github_Repo/Results_20200603
dir.create(results_dirname)




require(lme4)
require(lmerTest)
require(nlme)
require(RColorBrewer)
require(dplyr)
require(reshape2)
require(tidyr)
require(lmodel2)
require(car)
require(MuMIn)
#source("ggplot_helpers.R") # ggplot helper functions to get rid of pesky background

# set a nice palette
pal <- brewer.pal(9, name = "Set1")
pallight <- paste0(pal, "77")
palette(pallight)

## make a standard error function
se <- function(x){
  se <- sd(x, na.rm=T)/(length(na.omit(x)) ^ .5)
  return(se)
}

#### VIF function from lmer models: (from http://jonlefcheck.net/2012/12/28/dealing-with-multicollinearity-using-variance-inflation-factors/)
vif.mer <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

# create a function for plotting SMA regression lines (Type II regressions) to data, fitted using lmodel2
plot.MAR <- function(xvar, yvar, data, method="SMA", linecol, lwd=1, lty=1) {
  if(method=="SMA") meth = 3
  if(method=="MA") meth = 2
  if(method=="OLS") meth =1
  if(length(t(as.vector(data[!(is.na(data[,xvar]) | is.na(data[,yvar])), xvar])))<3){
    #return(rep(NA, times=7))
    break()
  }
  else{
    if(var(data[,yvar], na.rm=T)==0){
      break()
    }
    else{
      tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), data))
      intercept <- tmp.mod$regression.results$Intercept[meth]
      slope <- tmp.mod$regression.results$Slope[meth]
      yhat <- tmp.mod$x * slope + intercept
      lines(yhat[order(tmp.mod$x)]~tmp.mod$x[order(tmp.mod$x)], col=linecol, lwd=lwd, lty=lty)
      
    }
  }
}


#Mypairs
#Make fancy pair plots (modified from Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.)
Mypairs <- function(Z, col, pt.cex=0.8) {
  MyVarx <- colnames(Z)
  pairs(Z, labels = MyVarx,
        cex.labels =  2,
        lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
          panel.cor(x, y, digits, prefix, cex.cor)}, 
        upper.panel =  function(x, y) { 
          points(x, y, 
                 pch = 16, cex = pt.cex, 
                 col = as.numeric(factor(col)))
          abline(lm(y~x))})
  #print(P)
}


panel.cor <- function(x, y, digits=1, prefix="", cex.cor = .5)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- ifelse(abs(cor(x, y,use="pairwise.complete.obs"))<.5, .2, .3) 
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  #if(missing(cex.cor)) { cex <- 0.9/strwidth(txt) } else {
  #  cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex.cor * r)
}



##### ++++++++++++++++++++++++++++++++++++++ #############
################### BEGIN: LOAD DATA #####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# import soil data from the Soil and Landscape Grid of Australia (Grundy et al. 2015)
soilsums <- read.csv("data/Soils_summaries0-60cm_clim_20200106.csv") [,-1]
# this was created with the 'Cd-SoilsExtraction.R' code, also in this repo

# import raw branch-level traits
traits.b0 <- read.csv("data/TraitsAll_branch_20200512.csv", header=T, row.names = 1 )
# NOTE: these are raw data with problematic points identified by the 'Flag_XXX' columns (outliers, young leaves, lost or damaged samples, etc identified with flag >0)
# kill bad trait values so I don't always have to filter them.
traits.b0$LDMC[which(traits.b0$Flag_LDMC>0)] <- NA
traits.b0$SLA[which(traits.b0$Flag_LDMC>0)] <- NA
traits.b0$LMA[which(traits.b0$Flag_LDMC>0)] <- NA
traits.b0$Al_As[which(traits.b0$Flag_Area>0)] <- NA
traits.b0$SLA[which(traits.b0$Flag_Area>0)] <- NA
traits.b0$LMA[which(traits.b0$Flag_Area>0)] <- NA

# remove three additional outliers identified during analysis. Results are qualitatively robust to removal of these outliers.
#WD of ACAC-PER-B-5-1 is a weird low outlier that needs removed
traits.b0$WD[which(traits.b0$Branchtag=="ACAC-PER-B-5-1")] <- NA
#LMA  "VIMI-FREY-C-1-c" is also weird low outlier
traits.b0$LMA[which(traits.b0$Branchtag == "VIMI-FREY-C-1-c")] <- NA
# bad ESAL tree for LMA, but LDMC and Al_As seem to be OK
traits.b0$LMA[which(traits.b0$Treetag=="ESAL-BEN-C-1")] <- NA
# calculate Huber value instead of Al_As
traits.b0$hub <- 1/traits.b0$Al_As
# add in scaled climate for use with trait-climate models
traits.b0 <- traits.b0 %>% group_by(Species) %>% mutate(MDc_scaled = scale(MDc), PPTc_scaled = scale(PPTc), PETc_scaled=scale(PETc))

# merge soil data
soils.b <- soilsums[match(traits.b0$Plottag, soilsums$Plottag),-1]
traits.b <- data.frame(traits.b0, soils.b)
# renaming soil for residplots function later on
traits.b$pPC1fert <- traits.b$PC1fert
traits.b$pPC2depth <- traits.b$PC2depth
# ACAC has no DBHs (multistemed in general), so making a dummy value of 0 so my trait-climate model fitting function works for ACAC
traits.b$TreeDBH[which(traits.b$Species=="ACAC")] <- 0 # fill with dummy variable so my function works for fitting trait-env models



# aggregate traits to tree average 
# (need this for most plotting, which has too much overplotting at branch level)
traits.t0 <- data.frame(traits.b  %>% group_by (Species,Site,Plot,Tree,Plottag,Sitetag,Treetag,TreeDBH) %>% summarise(
  AI=mean(AI), pAI = mean(pAI),  MD = mean(MD), pMD = mean(pMD),
  PET=mean(PET),pPET=mean(pPET), PPT=mean(PPT), pPPT=mean(pPPT),
  MDc = mean(MDc), MATc = mean(MATc),
  PETc=mean(PETc), PPTc=mean(PPTc), Tminc = mean(Tminc), Tmaxc = mean(Tmaxc),
  Pminqc = mean(Pminqc), rhc = mean(rhc),
  Lat = mean(Lat), Lon= mean(Lon), StandBA = mean(StandBA, na.rm=T),
  mVolume = mean(Volume, na.rm=T),  sdWD = sd(WD, na.rm=T), WD = mean(WD, na.rm=T),
  sdLDMC = sd(LDMC, na.rm=T), LDMC = mean(LDMC, na.rm=T),
  sdLMA = sd(LMA, na.rm=T), LMA=mean(LMA, na.rm=T),
  avg_nleaves = mean(nleaves, na.rm=T), avg_tot_area = mean(tot_area, na.rm=T),
  avg_area = mean(avg_area, na.rm=T), median_area = mean(median_area, na.rm=T),
  max_area = max(max_area, na.rm=T), min_area = min(min_area, na.rm=T),
  sdAl_As = sd(Al_As, na.rm=T), Al_As = mean(Al_As, na.rm=T),
  sdhub = sd(hub, na.rm=T), hub = mean(hub, na.mr=T)
)
)
traits.t0$pchs <- rep(16, times=length(traits.t0$Species))
traits.t0$pchs[which(traits.t0$Species=="ACAC")]<- 3
# add in soils data
soils.t <- soilsums[match(traits.t0$Plottag, soilsums$Plottag),-1]
traits.t <- data.frame(traits.t0, soils.t)



# aggregate to plot level
traits.p0 <- data.frame(traits.t %>% group_by (Species,Site,Plot,Plottag,Sitetag) %>% summarise(
  AI=mean(AI,na.rm=T), pAI = mean(pAI,na.rm=T),  MD = mean(MD,na.rm=T), pMD = mean(pMD,na.rm=T), PET=mean(PET,na.rm=T),pPET=mean(pPET), PPT=mean(PPT), pPPT=mean(pPPT),
  MDc = mean(MDc), MATc = mean(MATc),
  PETc=mean(PETc), PPTc=mean(PPTc), Tminc = mean(Tminc), Tmaxc = mean(Tmaxc),
  Pminqc = mean(Pminqc), rhc = mean(rhc),
  Lat = mean(Lat), Lon= mean(Lon), StandBA = mean(StandBA, na.rm=T), meanDBH = mean(TreeDBH, na.rm=T),
  avg_nleaves = mean(avg_nleaves, na.rm=T), avg_area = mean(avg_area, na.rm=T), max_area = max(max_area, na.rm=T), min_area=min(min_area, na.rm=T),
  sdWD = sd(WD, na.rm=T), medWD = median(WD, na.rm=T), WD = mean(WD, na.rm=T),
  sdLMA = sd(LMA, na.rm=T), medLMA = median(LMA, na.rm=T), LMA = mean(LMA, na.rm=T),
  sdLDMC = sd(LDMC, na.rm=T), medLDMC = median(LDMC, na.rm=T), LDMC = mean(LDMC, na.rm=T),
  sdAl_As = sd(Al_As, na.rm=T), medAl_As = median(Al_As, na.rm=T), Al_As=mean(Al_As, na.rm=T),
  sdhub = sd(hub, na.rm=T), medhub = median(hub, na.rm=T), hub = mean(hub, na.rm=T))
  
)
# add in soils data
soils.p <- soilsums[match(traits.p0$Plottag, soilsums$Plottag),-1]
traits.p <- data.frame(traits.p0, soils.p)


# aggregate to site level
traits.s <- data.frame(traits.p %>% group_by (Species,Site,Sitetag) %>% summarise(
  AI=mean(AI, na.rm=T),  MD = mean(MD, na.rm=T), 
  PET=mean(PET, na.rm=T), PPT=mean(PPT, na.rm=T),
  MDc = mean(MDc), MATc = mean(MATc),
  PETc=mean(PETc), PPTc=mean(PPTc), Tminc = mean(Tminc), Tmaxc = mean(Tmaxc),
  Pminqc = mean(Pminqc), rhc = mean(rhc),
  Lat = mean(Lat), Lon= mean(Lon), StandBA = mean(StandBA, na.rm=T), meanDBH = mean(meanDBH, na.rm=T),
  CLY=mean(CLY),SND=mean(SND), SLT = mean(SLT),
  PTO=mean(PTO), NTO=mean(NTO), AWC=mean(AWC), BDW=mean(BDW),
  ECE=mean(ECE), DES=mean(DES), DER=mean(DER), PC1fert = mean(PC1fert),
  PC2depth=mean(PC2depth), Tmin=mean(Tmin), Elev=mean(Elev), MAT=mean(MAT),
  avg_nleaves = mean(avg_nleaves, na.rm=T), avg_area = mean(avg_area, na.rm=T), max_area = max(max_area, na.rm=T),
  min_area=min(min_area, na.rm=T),
  sdWD = sd(WD, na.rm=T), medWD = median(WD, na.rm=T), WD = mean(WD, na.rm=T),
  sdLMA = sd(LMA, na.rm=T), medLMA = median(LMA, na.rm=T), LMA = mean(LMA, na.rm=T),
  sdLDMC = sd(LDMC, na.rm=T), medLDMC = median(LDMC, na.rm=T), LDMC = mean(LDMC, na.rm=T),
  sdAl_As = sd(Al_As, na.rm=T), medAl_As = median(Al_As, na.rm=T), Al_As=mean(Al_As, na.rm=T),
  sdhub = sd(hub, na.rm=T), medhub = median(hub, na.rm=T), hub = mean(hub, na.rm=T)
)

)
traits.s$AIc <- traits.s$PPTc/traits.s$PETc

# adding log.Al_As and huber value
traits.b$log.Al_As <- log(traits.b$Al_As, base=10)
traits.t$log.Al_As <- log(traits.t$Al_As, base=10)
traits.p$log.Al_As <- log(traits.p$Al_As, base=10)
traits.s$log.Al_As <- log(traits.s$Al_As, base=10)


traits.b$log.hub <- log(traits.b$hub, base=10)
traits.t$log.hub <- log(traits.t$hub, base=10)
traits.p$log.hub <- log(traits.p$hub, base=10)
traits.s$log.hub <- log(traits.s$hub, base=10)



# climate quantiles for the whole species distributions (for code used to create this, email LDL Anderegg)
quants <- read.csv("data/Climate_Quantiles_allspp_20200516.csv")[,-1]
quants$Species <- quants$spp
levels(quants$Species) <- list(ACAC="acac",ESAL="esal",EMARG="emarg",COCA="coca",OVAT="ovat",VIMI="vimi",AMYG="amyg",OBLI="obli")

########### Species Means ################################
traits.sp <- traits.b %>% group_by(Species) %>% summarise(WD = mean(WD, na.rm=T), LMA=mean(LMA, na.rm=T), LDMC=mean(LDMC, na.rm=T), hub = mean(hub, na.rm=T))
traits.sp$log.hub <- log(traits.sp$hub, base=10)

traits.sp$MD.50 <- quants$quantCMDcALA.50[match(traits.sp$Species, quants$Species)]
traits.sp$MD.90 <- quants$quantCMDcALA.90[match(traits.sp$Species, quants$Species)]


### load summary of intraspecific residual patterns ########
# respats <- read.csv("Intraspecific_Results/Resid_Patterns_20180811.csv")
#respats <- read.csv("Intraspecific_Results/Resid_Patterns_20200116.csv") # with hub and updated trait-climate relationships with soil and stand
#respats <- read.csv("Intraspecific_Results/Resid_Patterns_20200209.csv") # updated when cleaning code. droped Al_As columns
respats <- read.csv("data/Resid_Patterns_20200526.csv") # updated with NP R2R
  # ignore the Warning message. not sure why it throws it.

################# END: LOAD DATA #######################################





##### ++++++++++++++++++++++++++++++++++++++ #############
################### BEGIN: VISUALIZE TRAIT-CLIMATE #####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## . FIG 1: Trait-Env relationships W.HUBER & MD Chelsa ####################
palette(pallight)
#quartz(width=3.4, height=6)
jpeg(file=paste0("./",results_dirname,"/Fig1_TraitClimate_ChelsaMD.jpg"), width=3.4, height=6, units = "in",res = 600)
par(mfrow=c(4,1), mar=c(0,4,0,1), mgp=c(2.2,1,0), oma=c(3.5,0,3.5,0), cex.lab=1.3)
for( pp in c(1:3)){
  k <- c("WD","LMA","LDMC")[pp]
  plot(get(k)~MDc, traits.t, pch=pchs, cex=.9, col=Species, ylab="", xaxt="n")
  # points(get(k)~MD.50, traits.sp, pch=16, cex=2)
  if(k=="WD") {
    legend(x = -1200, y=1, legend =c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl")
           , text.font=3, pch=c(3,rep(16, times=7)), lty=1,lwd=1.5, col=pal[c(1,5,4,3,7,8,2,6)], xpd=NA, ncol=4, bty="n")
    mtext(expression(paste("WD ",(g/cm^3))), side=2, line=2.2, font=2)
  }
  if(k == "LMA"){
    mtext(expression(paste("LMA ",(g/cm^2))), side=2, line=2.2, font=2)
  }
  if(k == "LDMC"){
    mtext(expression(paste("LDMC ",(g/g))), side=2, line=2.2, font=2)
  }
  
  for (j in 1:length(levels(traits.p$Species))){
    i <- levels(traits.p$Species)[j] 
    tmp <- traits.t[which(traits.t$Species==i),]
    tmpmod <- lm(get(k)~MDc, tmp)
    tmp1 <- arrange(data.frame(x=tmpmod$model$MDc, y=tmpmod$fitted.values),x)
    lines(tmp1$y[c(1,nrow(tmp1))]~tmp1$x[c(1,nrow(tmp1))], col=pal[j], lwd=3)
    #plot.MAR(xvar = "MDc", yvar = k,data= traits.p[which(traits.p$Species==i),], linecol = pal[j], lwd=3)
  }
}

# add in log.huber panel
plot(log.hub~MDc, traits.t, pch=pchs, cex=.9, ylab="", col=Species, xaxt="n", yaxt="n")
k <- "log.hub"
# points(get(k)~MD.50, traits.sp, pch=16, cex=2)
for (j in 1:length(levels(traits.p$Species))){
  i <- levels(traits.p$Species)[j] 
  tmp <- traits.t[which(traits.t$Species==i),]
  tmpmod <- lm(get(k)~MDc, tmp)
  tmp1 <- arrange(data.frame(x=tmpmod$model$MDc, y=tmpmod$fitted.values),x)
  lines(tmp1$y[c(1,nrow(tmp1))]~tmp1$x[c(1,nrow(tmp1))], col=pal[j], lwd=3)
  #plot.MAR(xvar = "MDc", yvar = k,data= traits.p[which(traits.p$Species==i),], linecol = pal[j], lwd=3)
}
axis(2, labels = c("0.01","0.05","0.10","0.15"), at=log(c(0.01,0.05,0.10,0.15), base=10))
axis(1)
mtext(text = "MD (PET-PPT, mm)", side=1, line=2.5)
mtext(expression(paste("HV ",(mm^2/cm^2))), side=2, line=2.2, font=2)
#quartz.save(file=paste0("./",results_dirname,"/Fig1_TraitClimate_ChelsaMD.pdf"),type = "pdf" )
dev.off()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## . FIG S4: Trait-Env Across Species####################
pallighter <- paste0(pal, "22")
palette(pallighter)
quartz(width=3.4, height=6)
par(mfrow=c(4,1), mar=c(0,4,0,1), mgp=c(2.2,1,0), oma=c(3.5,0,3.5,0), cex.lab=1.3)
for( pp in c(1:3)){
  k <- c("WD","LMA","LDMC")[pp]
  plot(get(k)~MDc, traits.t, pch=pchs, cex=.9, col=Species, ylab="", xaxt="n")
  
  if(k=="WD") {
    legend(x = -1200, y=1, legend =c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl")
           , text.font=3, pch=c(3,rep(16, times=7)), lty=1,lwd=1.5, col=pal[c(1,5,4,3,7,8,2,6)], xpd=NA, ncol=4, bty="n")
    mtext(expression(paste("WD ",(g/cm^3))), side=2, line=2.2, font=2)
    mod <- lm(WD~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])
    abline(mod, lwd=2, lty=2, col="grey")
    mtext(text=paste0("p=",round(summary(mod)$coefficients[2,4], 3)), side=3, line=-1.5, adj=.1)
  }
  if(k == "LMA"){
    mtext(expression(paste("LMA ",(g/cm^2))), side=2, line=2.2, font=2)
    mod <- lm(LMA~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])
    abline(mod, lwd=2)
    mtext(text=paste0("p=",round(summary(mod)$coefficients[2,4], 3)), side=3, line=-1.5, adj=.1)
    
    
  }
  if(k == "LDMC"){
    mtext(expression(paste("LDMC ",(g/g))), side=2, line=2.2, font=2)
    mod <- lm(LDMC~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])
    abline(mod, lwd=2)
    mtext(text=paste0("p=",round(summary(mod)$coefficients[2,4], 3)), side=3, line=-1.5, adj=.1)
    
    
  }
  
  for (j in 1:length(levels(traits.p$Species))){
    i <- levels(traits.p$Species)[j] 
    tmp <- traits.t[which(traits.t$Species==i),]
    tmpmod <- lm(get(k)~MDc, tmp)
    tmp1 <- arrange(data.frame(x=tmpmod$model$MDc, y=tmpmod$fitted.values),x)
    lines(tmp1$y[c(1,nrow(tmp1))]~tmp1$x[c(1,nrow(tmp1))], col=pallight[j], lwd=3)
    points(get(k)~MD.50, traits.sp[which(traits.sp$Species==i),], pch=ifelse(i=="ACAC",3,16), cex=2.5, col=pal[j])
    #plot.MAR(xvar = "MDc", yvar = k,data= traits.p[which(traits.p$Species==i),], linecol = pal[j], lwd=3)
  }
}

# add in log.huber panel
plot(log.hub~MDc, traits.t, pch=pchs, cex=.9, ylab="", col=Species, xaxt="n", yaxt="n")
k <- "log.hub"
# points(get(k)~MD.50, traits.sp, pch=16, cex=2)
for (j in 1:length(levels(traits.p$Species))){
  i <- levels(traits.p$Species)[j] 
  tmp <- traits.t[which(traits.t$Species==i),]
  tmpmod <- lm(get(k)~MDc, tmp)
  tmp1 <- arrange(data.frame(x=tmpmod$model$MDc, y=tmpmod$fitted.values),x)
  lines(tmp1$y[c(1,nrow(tmp1))]~tmp1$x[c(1,nrow(tmp1))], col=pallight[j], lwd=3)
  points(get(k)~MD.50, traits.sp[which(traits.sp$Species==i),], pch=ifelse(i=="ACAC",3,16), cex=2.5, col=pal[j])
  #plot.MAR(xvar = "MDc", yvar = k,data= traits.p[which(traits.p$Species==i),], linecol = pal[j], lwd=3)
}
mod <- lm(log.hub~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])
abline(mod, lwd=2)
mtext(text=paste0("p=",round(summary(mod)$coefficients[2,4], 3)), side=3, line=-1.5, adj=.1)

axis(2, labels = c("0.01","0.05","0.10","0.15"), at=log(c(0.01,0.05,0.10,0.15), base=10))
axis(1)
mtext(text = "MD (PET-PPT, mm)", side=1, line=2.5)
mtext(expression(paste("HV ",(mm^2/cm^2))), side=2, line=2.2, font=2)
quartz.save(file=paste0("./",results_dirname,"/FigS4_TraitClimate_ChelsaMD_SpeciesMean.pdf"),type = "pdf" )


### significance of among-species trait relationships
summary(lm(WD~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])) # P=0.05235 .
summary(lm(LMA~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])) # P=0.0309 *
summary(lm(LDMC~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])) # P=0.008513 **
summary(lm(log.hub~MD.50, traits.sp[which(traits.sp$Species!="ACAC"),])) # P=0.003756 **








#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############. FIG S1: Colinearity among variables #####
palette(pal)


quartz(width=7, height=6)
Mypairs(traits.s %>% dplyr::select(Tmin=Tminc,MAT=MATc,PPT=PPTc,PET=PETc,MD=MDc,rh=rhc, AI=AIc,PC1=PC1fert,PC2=PC2depth), col=traits.s$Species, pt.cex=1.3)

legend("bottom",legend=levels(traits.s$Species), col=pal[1:8], pch=16, ncol=4, xpd=NA, bty="n", cex=.8, inset=1.05)
quartz.save(file=paste0("./",results_dirname,"/FigS2_Pairplot_climatesoil_chelsa.pdf"),type = "pdf" )













#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###### BEGIN: TRAIT-CLIMATE Model Selection ########
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# procedure:
# 1) fit host of candidate models (all single climate and soil variables, with and without DBH and Stand Basal Area tree covariates if measured)
# 2) Select best model via AICc
# 3) 


# useful model diagnostics plots. Not all that are suggested by icebreakR, but a reasonable sample
residplots <- function(model, dataz, climvar, res.type="normalized", Species = "",Trait="", mixed = F, path=getwd(), vers = "XXXXX", write.jpeg=T ){
  if(write.jpeg==T){
    jpeg(filename = paste(path,"/",Species,Trait,"_Resids_",vers,".jpeg", sep=""),width = 7, height=6, units = "in", res=600)
  }
  else{quartz(width=7, height=6, title=paste(dataz$Species[1], Trait))}
  
  par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(0,0,2,0))
  if(mixed==T){
    ref.group <- ranef(model)
    ref.var.group1 <- tapply(residuals(model, type="pearson", level=1),
                             dataz[,names(ref.group[1])], var)
    ref.var.group2 <- tapply(residuals(model, type="pearson", level=1),
                             dataz[,names(ref.group[2])], var)
    qqp(ref.group[[1]][[1]], main="Tree Random Effects")
    qqp(ref.group[[2]][[1]], main="Plot Random Effects")
    # qqnorm(ref.group[[1]][[1]], main="Q-Q Normal - group Random Effects")
    # qqline(ref.group[[1]][[1]], col="red")
  }
  else{
    # qqnorm(residuals(model), main = "raw"); qqline(residuals(model))
    qqp(residuals(model), main="raw resids")
  }
  # qqnorm(residuals(model, type=res.type), main="normalized"); qqline(residuals(model, type=res.type))
  qqp(residuals(model, type=res.type), main="normalized")
  mod.resids <- residuals(model, type=res.type)
  scatter.smooth(mod.resids~dataz[,climvar],xlab=climvar); abline(h=0, col="red")
  scatter.smooth(mod.resids~dataz$Diameter); abline(h=0, col="red")
  # scatter.smooth(mod.resids~dataz$BAI10yr); abline(h=0, col="red")
  boxplot(mod.resids~factor(dataz$Plottag), notch=T, varwidth=T); abline(h=0)
  #quartz(width=9, height=4)
  if(write.jpeg==T){
    dev.off()
    jpeg(filename = paste(path,"/",Species,Trait,"_Vars_",vers,".jpeg", sep=""),width = 9, height=4, units = "in", res=600)
  }
  else{ quartz(width=9, height=4)}
  par(mfrow=c(1,4), mar=c(4,4,1,1), oma=c(0,0,2,0))
  ref.group[[2]]$MDc <- traits.p$pMD[match(rownames(ref.group[[2]]),traits.p$Plottag)]
  ref.group[[2]]$PPTc <- traits.p$pPPT[match(rownames(ref.group[[2]]),traits.p$Plottag)]
  ref.group[[2]]$PETc <- traits.p$pPET[match(rownames(ref.group[[2]]),traits.p$Plottag)]
  ref.group[[2]]$PC1fert <- traits.p$PC1fert[match(rownames(ref.group[[2]]),traits.p$Plottag)]
  ref.group[[2]]$PC2depth <- traits.p$PC2depth[match(rownames(ref.group[[2]]),traits.p$Plottag)]
  ref.group[[2]]$var <- ref.var.group2[match(rownames(ref.group[[2]]), names(ref.var.group2))]
  
  ref.group[[1]]$MDc <- traits.t$pMD[match(rownames(ref.group[[1]]),traits.t$Treetag)]
  ref.group[[1]]$PETc <- traits.t$pPET[match(rownames(ref.group[[1]]),traits.t$Treetag)]
  ref.group[[1]]$PPTc <- traits.t$pPPT[match(rownames(ref.group[[1]]),traits.t$Treetag)]
  ref.group[[1]]$PC1fert <- traits.t$PC1fert[match(rownames(ref.group[[1]]),traits.t$Treetag)]
  ref.group[[1]]$PC2depth <- traits.t$PC2depth[match(rownames(ref.group[[1]]),traits.t$Treetag)]
  ref.group[[1]]$var <- ref.var.group1[match(rownames(ref.group[[1]]), names(ref.var.group1))]
  
  dataz$TreeEff <- ref.group[[1]]$`(Intercept)`[match(dataz$Treetag,rownames(ref.group[[1]]))]
  plot(ref.group[[2]]$`(Intercept)`~ref.group[[2]][,climvar],xlab=climvar, main="Plot Ran Effs")
  abline(h=0)
  # plot(TreeEff~pMD, dataz)
  # abline(h=0)
  plot(ref.group[[1]]$`(Intercept)`~ref.group[[1]][,climvar], dataz,xlab=climvar, main="Tree Ran Effs")
  abline(h=0)
  scatter.smooth(ref.group[[2]]$var~ref.group[[2]][,climvar], main="within Plot var",xlab=climvar)
  scatter.smooth(ref.group[[1]]$var~ref.group[[1]][,climvar], main="within Tree var",xlab=climvar)
  if(write.jpeg==T){dev.off()}
  #   plot(TreeEff~pPPT, dataz)
  #   abline(h=0)
}


fit.traits <- function(dataz, trait, species, diam = F, standBA = F){
  mods <- list()
  if(diam==T){
    mods$null <- lmer(get(trait)~ Diameter + (1|Plottag) + (1|Treetag), data = dataz, subset=Species==species, REML=F)
    mods$PET <- lmer(get(trait)~Diameter +PETc_scaled + (1|Plottag) + (1|Treetag), data = dataz, subset=Species==species, REML=F)
    mods$PPT <- lmer(get(trait)~Diameter +PPTc_scaled + (1|Plottag) + (1|Treetag), data = dataz, subset=Species==species, REML=F)
    mods$MD <- lmer(get(trait)~Diameter +MDc_scaled + (1|Plottag) + (1|Treetag), data = dataz, subset=Species==species, REML=F)
    mods$fert <- lmer(get(trait)~Diameter + PC1fert + (1|Plottag) + (1|Treetag), data = dataz, subset=Species==species, REML=F)
    mods$depth <- lmer(get(trait)~Diameter + PC2depth + (1|Plottag) + (1|Treetag), data = dataz, subset=Species==species, REML=F)
    mods$DBH <- lmer(get(trait)~Diameter + TreeDBH + (1|Plottag) + (1|Treetag), data = dataz, subset=Species==species, REML=F)
    mods$PETdbh <- update(mods$PET,.~TreeDBH + . )
    mods$PPTdbh <- update(mods$PPT,.~TreeDBH + . )
    mods$MDdbh <- update(mods$MD,.~TreeDBH + . )
    mods$fertdbh <- update(mods$fert,.~TreeDBH + . )
    mods$depthdbh <- update(mods$depth,.~TreeDBH + . )
    if(standBA == T){
    mods$StandBA <- lmer(get(trait)~Diameter + StandBA + (1|Plottag) + (1|Treetag), data = dataz, subset=Species==species, REML=F)
    mods$PETba <- update(mods$PET,.~ StandBA + . )
    mods$PPTba <- update(mods$PPT,.~ StandBA + . )
    mods$MDba <- update(mods$MD,.~ StandBA + . )
    mods$fertba <- update(mods$fert,.~ StandBA + . )
    mods$depthba <- update(mods$depth,.~ StandBA + . )
    mods$dbhba <- update(mods$DBH,.~ StandBA + . )
    mods$PETdbhba <- update(mods$PET,.~ StandBA + TreeDBH + . )
    mods$PPTdbhba <- update(mods$PPT,.~ StandBA + TreeDBH + . )
    mods$MDdbhba <- update(mods$MD,.~ StandBA + TreeDBH + . )
    mods$fertdbhba <- update(mods$fert,.~ StandBA + TreeDBH + . )
    mods$depthdbhba <- update(mods$depth,.~ StandBA + TreeDBH + . )
    
    
    print(AICc(mods$null,mods$PET,mods$PPT,mods$MD,mods$fert,mods$depth,mods$DBH,mods$StandBA, 
              mods$PETdbh,mods$PPTdbh,mods$MDdbh,mods$fertdbh,mods$depthdbh,
              mods$PETba,mods$PPTba,mods$MDba,mods$fertba,mods$depthba,
              mods$PETdbhba,mods$PPTdbhba,mods$MDdbhba,mods$fertdbhba,mods$depthdbhba, mods$dbhba)[order(AICc(mods$null,mods$PET,mods$PPT,mods$MD,mods$fert,mods$depth,mods$DBH,mods$StandBA, 
              mods$PETdbh,mods$PPTdbh,mods$MDdbh,mods$fertdbh,mods$depthdbh,
              mods$PETba,mods$PPTba,mods$MDba,mods$fertba,mods$depthba,
              mods$PETdbhba,mods$PPTdbhba,mods$MDdbhba,mods$fertdbhba,mods$depthdbhba, mods$dbhba)$AICc,decreasing = F),])
    }
    if(standBA==F){
      print(AICc(mods$null,mods$PET,mods$PPT,mods$MD,mods$fert,mods$depth,mods$DBH, mods$PETdbh,mods$PPTdbh,mods$MDdbh,mods$fertdbh,mods$depthdbh)[order(AICc(mods$null,mods$PET,mods$PPT,mods$MD,mods$fert,mods$depth,mods$DBH, mods$PETdbh,mods$PPTdbh,mods$MDdbh,mods$fertdbh,mods$depthdbh)$AICc,decreasing = F),])
    }
  }
  if(diam==F){
    mods$null <- lmer(get(trait)~ 1 + (1|Plottag) + (1|Treetag),
                      data = dataz, subset=Species==species, REML=F)
    mods$PET <- lmer(get(trait)~PETc_scaled + (1|Plottag) + (1|Treetag),
                     data = dataz, subset=Species==species, REML=F)
    mods$PPT <- lmer(get(trait)~PPTc_scaled + (1|Plottag) + (1|Treetag),
                     data = dataz, subset=Species==species, REML=F)
    mods$MD <- lmer(get(trait)~MDc_scaled + (1|Plottag) + (1|Treetag),
                    data = dataz, subset=Species==species, REML=F)
    mods$fert <- lmer(get(trait)~ PC1fert + (1|Plottag) + (1|Treetag),
                      data = dataz, subset=Species==species, REML=F)
    mods$depth <- lmer(get(trait)~ PC2depth + (1|Plottag) + (1|Treetag),
                       data = dataz, subset=Species==species, REML=F)
    mods$DBH <- lmer(get(trait)~ TreeDBH + (1|Plottag) + (1|Treetag),
                     data = dataz, subset=Species==species, REML=F)
    mods$PETdbh <- update(mods$PET,.~TreeDBH + . )
    mods$PPTdbh <- update(mods$PPT,.~TreeDBH + . )
    mods$MDdbh <- update(mods$MD,.~TreeDBH + . )
    mods$fertdbh <- update(mods$fert,.~TreeDBH + . )
    mods$depthdbh <- update(mods$depth,.~TreeDBH + . )
    if(standBA == T){
      mods$StandBA <- lmer(get(trait)~ StandBA + (1|Plottag) + (1|Treetag),
                           data = dataz, subset=Species==species, REML=F)
      mods$PETba <- update(mods$PET,.~ StandBA + . )
    mods$PPTba <- update(mods$PPT,.~ StandBA + . )
    mods$MDba <- update(mods$MD,.~ StandBA + . )
    mods$fertba <- update(mods$fert,.~ StandBA + . )
    mods$depthba <- update(mods$depth,.~ StandBA + . )
    mods$dbhba <- update(mods$DBH,.~ StandBA + . )
    mods$PETdbhba <- update(mods$PET,.~ StandBA + TreeDBH + . )
    mods$PPTdbhba <- update(mods$PPT,.~ StandBA + TreeDBH + . )
    mods$MDdbhba <- update(mods$MD,.~ StandBA + TreeDBH + . )
    mods$fertdbhba <- update(mods$fert,.~ StandBA + TreeDBH + . )
    mods$depthdbhba <- update(mods$depth,.~ StandBA + TreeDBH + . )
      print(AICc(mods$null,mods$PET,mods$PPT,mods$MD,mods$fert,mods$depth,mods$DBH,mods$StandBA, 
              mods$PETdbh,mods$PPTdbh,mods$MDdbh,mods$fertdbh,mods$depthdbh,
              mods$PETba,mods$PPTba,mods$MDba,mods$fertba,mods$depthba,
              mods$PETdbhba,mods$PPTdbhba,mods$MDdbhba,mods$fertdbhba,mods$depthdbhba, mods$dbhba)[order(AICc(mods$null,mods$PET,mods$PPT,mods$MD,mods$fert,mods$depth,mods$DBH,mods$StandBA, 
              mods$PETdbh,mods$PPTdbh,mods$MDdbh,mods$fertdbh,mods$depthdbh,
              mods$PETba,mods$PPTba,mods$MDba,mods$fertba,mods$depthba,
              mods$PETdbhba,mods$PPTdbhba,mods$MDdbhba,mods$fertdbhba,mods$depthdbhba, mods$dbhba)$AICc,decreasing = F),])
    }
    if(standBA==F){
      print(AICc(mods$null,mods$PET,mods$PPT,mods$MD,mods$fert,mods$depth,mods$DBH, mods$PETdbh,mods$PPTdbh,mods$MDdbh,mods$fertdbh,mods$depthdbh)[order(AICc(mods$null,mods$PET,mods$PPT,mods$MD,mods$fert,mods$depth,mods$DBH, mods$PETdbh,mods$PPTdbh,mods$MDdbh,mods$fertdbh,mods$depthdbh)$AICc,decreasing = F),])
    }
  }
  return(mods)
}

### ** WD modeling ##########
## first make a directory to store things in:
#dir.create("/Users/leeanderegg/Dropbox/Trade-offs project/Intraspecific_Results/Trait-Clim_ModCrit/")
# dir.create("/Users/leeanderegg/Dropbox/Trade-offs project/Intraspecific_Results/Trait-Clim_ModCrit/WoodDensity")
dir.create(paste0("./",results_dirname,"/Trait-Clim_ModCrit"))

#current_version <- "20180812"
# current_version <- "20200110"
current_version <- "20200521"
current_path <- paste0(results_dirname,"/Trait-Clim_ModCrit")

# before modeling, should note that:
testAMYG <- lmer(WD~Site + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="AMYG", REML=F)
testAMYGd <- lmer(WD~Site + Diameter + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="AMYG", REML=F)
anova(testAMYG, testAMYGd)
# extremely significant diameter effect for AMYG
testOBLI <- lmer(WD~Site + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="OBLI", REML=F)
testOBLId <- lmer(WD~Site + Diameter + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="OBLI", REML=F)
anova(testOBLI, testOBLId)
# extremely significant diameter effect for OBLI
testOVAT <- lmer(WD~Site + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="OVAT", REML=F)
testOVATd <- lmer(WD~Site + Diameter + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="OVAT", REML=F)
anova(testOVAT, testOVATd)
# marginally signicant diameter effect for OVAT (p=0.047, dAIC=1.95)
testVIMI <- lmer(WD~Site + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="VIMI", REML=F)
testVIMId <- lmer(WD~Site + Diameter + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="VIMI", REML=F)
anova(testVIMI, testVIMId)
# no effect for VIMI
testEMARG <- lmer(WD~Site + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG", REML=F)
testEMARGd <- lmer(WD~Site + Diameter + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG", REML=F)
anova(testEMARG, testEMARGd)
# extremely significant diameter effect for EMARG
testCOCA <- lmer(WD~Site + (1|Treetag), data = traits.b, subset=Species=="COCA", REML=F)
testCOCAd <- lmer(WD~Site + Diameter + (1|Treetag), data = traits.b, subset=Species=="COCA", REML=F)
anova(testCOCA,testCOCAd)
# marginal effect for COCA p=0.1075 (note: removed Plottag because estimated near 0 so throwing isSingular errors)
testESAL <- lmer(WD~Site + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL", REML=F)
testESALd <- lmer(WD~Site + Diameter + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL", REML=F)
anova(testESAL, testESALd)
# no effect for ESAL
testACAC <- lmer(WD~Site + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC", REML=F)
testACACd <- lmer(WD~Site + Diameter + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC", REML=F)
anova(testACAC,testACACd)
# marginal effect for ACAC p=0.06764

# need diameter for:

# AMYG, OBLI, OVAT, EMARG, COCA, ACAC (all sig and marginal effects)

# no need for diameter for:
#ESAL, VIMI

##
ACACwd <- fit.traits(dataz=traits.b, trait="WD",species="ACAC",standBA=F, diam = T)
  # PET best (MD close)

residplots(model=ACACwd$PET, dataz=traits.b[which(traits.b$Species=="ACAC" & traits.b$WD>0),]
           ,res.type = "pearson",mixed = T
           , Species = "ACAC", Trait="WD", vers=current_version
           , climvar="PETc", path=current_path)

## Variance pattern?
scatter.smooth(abs(resid(ACACwd$PET))~round(ACACwd$PET@frame$PETc_scaled,digits = 3))
boxplot(resid(ACACwd$PET)~round(traits.b$PETc[which(traits.b$Species=="ACAC" & traits.b$WD>0)], digits=0));abline(h=0, col="red")

  # variance appears to peak in the middle
ACACsite <- lme(WD~Diameter + PETc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="ACAC")
ACACsite1 <- lme(WD~Diameter + PETc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="ACAC", weights = varIdent (form = ~1|Sitetag))
ACACsite2 <- lme(WD~Diameter + PETc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="ACAC", weights = varIdent(~PETc_scaled))
ACACsite3 <- lme(WD~Diameter + PETc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="ACAC", weights = varPower (form = ~ PETc_scaled))
ACACsite4 <- lme(WD~Diameter + PETc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="ACAC", weights = varPower (form = ~ PETc))
AICc(ACACsite, ACACsite1, ACACsite2, ACACsite3, ACACsite4) 
# so there's a difference bewteen sites, but it's not related to PET -> it's hump shaped at mid values of PET




ESALwd <- fit.traits(dataz=traits.b, species="ESAL",trait="WD",diam=F, standBA=F)
# Missing DBH from a number of ESAL. So have to apply AIC outside function
unlist(lapply(ESALwd, FUN=AICc))
# MD best followed closely by PET
anova(ESALwd$MD, ESALwd$null) # p =0.00012

residplots(model=ESALwd$MD, dataz=traits.b[which(traits.b$Species=="ESAL" & traits.b$WD>0),]
           ,res.type = "pearson",mixed = T
           , Species = "ESAL", Trait="WD", vers=current_version
           , climvar = "MDc",write.jpeg = T, path=current_path)


### Variance Constant?
scatter.smooth(abs(resid(ESALwd$MD))~round(ESALwd$MD@frame$MDc_scaled,digits = 3))
boxplot(resid(ESALsite)~round(traits.b$pMD[which(traits.b$Species=="ESAL" & traits.b$WD>0)], digits=0))

# Seems so
ESALsite <- lme(WD~MDc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="ESAL")
ESALsite1 <- lme(WD~MDc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="ESAL", weights = varIdent (form = ~1|Sitetag))
ESALsite2 <- lme(WD~MDc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="ESAL", weights = varFixed (~MDc_scaled))

AICc(ESALsite, ESALsite1, ESALsite2)
# I need go no further. there is clearly no variance pattern here



EMARGwd <- fit.traits(dataz=traits.b, species="EMARG",trait="WD",diam=T, standBA=F)
anova(EMARGwd$depthdbh, EMARGwd$depth, EMARGwd$null) # depth p =0.003996, dbh p=0.086
anova(EMARGwd$PPT, EMARGwd$null) # p =0.01397 *



residplots(model=EMARGwd$depth, dataz=traits.b[which(traits.b$Species=="EMARG" & traits.b$WD>0),]
           ,res.type = "pearson",mixed = T
           , Species = "EMARG", Trait="WD", vers=current_version
           , climvar="PC2depth", path=current_path)

### Variance
scatter.smooth(abs(resid(EMARGwd$depth))~round(EMARGwd$depth@frame$PC2depth,digits = 3))
scatter.smooth(abs(resid(EMARGwd$PPTdbh))~round(EMARGwd$PPTdbh@frame$PPTc_scaled,digits = 3))
boxplot(resid(EMARGwd$PPT)~round(EMARGwd$PPT@frame$PPTc_scaled,digits = 1))


EMARGsite <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="EMARG")
EMARGsite1 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="EMARG", weights = varIdent (form = ~1|Sitetag))
EMARGsite2 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="ESAL", weights = varFixed(~PPTc_scaled))
EMARGsite3 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="EMARG", weights = varPower (form = ~PPTc_scaled))
EMARGsite4 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="EMARG", weights = varPower (form = ~PPTc))
AIC(EMARGsite, EMARGsite1, EMARGsite2, EMARGsite3, EMARGsite4) 
anova(EMARGsite, EMARGsite4) # no need for different variances
# PPT best, but not strong
# no residual trend


COCAwd <- fit.traits(dataz=traits.b, species="COCA",trait="WD",diam=T, standBA=F)
anova(COCAwd$PET, COCAwd$null) # p =0.1663
anova(COCAwd$DBH, COCAwd$null) # p =0.1514
anova(COCAwd$PETdbh, COCAwd$null) #p = 0.1055

residplots(model=COCAwd$PET, dataz=traits.b[which(traits.b$Species=="COCA" & traits.b$WD>0),]
           ,res.type = "pearson",mixed = T
           , Species = "COCA-PET-", Trait="WD", vers=current_version
           , climvar="PETc", path=current_path)
residplots(model=COCAwd$null, dataz=traits.b[which(traits.b$Species=="COCA" & traits.b$WD>0),]
           ,res.type = "pearson",mixed = T
           , Species = "COCA-null-", Trait="WD", vers=current_version
           , climvar="PETc", path=current_path)


# Plot Ranefs, Tree ranefs, within plot var, and within tree var seem to increase with MD. as do resids?
scatter.smooth(abs(resid(COCAwd$PET))~round(COCAwd$PET@frame$PETc_scaled,digits = 3))

COCAsite <- lme(WD~Diameter , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="COCA")
COCAsite1 <- lme(WD~Diameter , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="COCA", weights = varIdent (form = ~1|Sitetag))
COCAsite2 <- lme(WD~Diameter , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="COCA", weights = varFixed(~ PETc))
COCAsite3 <- lme(WD~Diameter , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="COCA", weights = varPower (form = ~ PETc_scaled))
COCAsite4 <- lme(WD~Diameter , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="COCA", weights = varPower (form = ~ PETc))
AIC(COCAsite, COCAsite1, COCAsite2, COCAsite3, COCAsite4) # 
anova(COCAsite, COCAsite2) # var increases with PET
# true with null model and PET model
#*** does appear to be an increase in variance at drier places        



### Stats for Tasssie Species ###

OVATwd <- fit.traits(dataz=traits.b, species="OVAT",trait="WD",diam=T, standBA=T)
anova(OVATwd$fertdbh, OVATwd$fert, OVATwd$null) # p =0.0009721 ***
anova(OVATwd$depth, OVATwd$null) # p =0.00111 **
anova(OVATwd$PPT, OVATwd$null) # p=0.002936 **
anova(OVATwd$DBH, OVATwd$null) # p=0.58
anova(OVATwd$StandBA, OVATwd$null) # p=0.87


residplots(model=OVATwd$fert, dataz=traits.b[which(traits.b$Species=="OVAT" & traits.b$WD>0),]
           ,res.type = "pearson",mixed = T
           , Species = "OVAT-fert", Trait="WD", vers=current_version, climvar="PC1fert", path=current_path)

#looks like the PET is really needed to yank one site's plot random effects back in. otherwise, they're strong negative resids
scatter.smooth(abs(resid(OVATwd$PPT))~round(OVATwd$PPT@frame$PPTc_scaled,digits = 3))
# not much residual pattern

OVATsite <- lme(WD~Diameter + PC1fert  , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0 & traits.b$Species=="OVAT"),], subset=Species=="OVAT")
OVATsite1 <- lme(WD~Diameter + PC1fert  , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="OVAT"),], subset=Species=="OVAT", weights = varIdent (form = ~1|Sitetag))
OVATsite2 <- lme(WD~Diameter + PC1fert  , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="OVAT"),], subset=Species=="OVAT", weights = varFixed( ~ PPTc_scaled))
OVATsite3 <- lme(WD~Diameter + PC1fert  , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="OVAT"),], subset=Species=="OVAT", weights = varPower (form = ~ PPTc_scaled))
OVATsite4 <- lme(WD~Diameter + PC1fert  , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="OVAT"),], subset=Species=="OVAT", weights = varPower (form = ~PPTc))
AIC(OVATsite, OVATsite1, OVATsite2, OVATsite3) # OVATsite is best with low df
# no evidence with any fixed effect structure (fert, PPT, fert + PPT)



VIMIwd <- fit.traits(dataz=traits.b, species="VIMI",trait="WD",diam=F, standBA=T)
anova(VIMIwd$fert, VIMIwd$null) # p =0.0001391 ***
anova(VIMIwd$PPT, VIMIwd$null) # p=0.001789 **
anova(VIMIwd$DBH, VIMIwd$null) # p=0.24
anova(VIMIwd$StandBA, VIMIwd$null) # p=0.175

residplots(model=VIMIwd$fert, dataz=traits.b[which(traits.b$Species=="VIMI" & traits.b$WD>0),]
           ,res.type = "pearson",mixed = T
           , Species = "VIMI", Trait="WD", vers=current_version
           , climvar= "PC1fert", path=current_path)

scatter.smooth(abs(resid(VIMIwd$PPT))~round(VIMIwd$PPT@frame$PPTc_scaled,digits = 3))

VIMIsite <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0 & traits.b$Species=="VIMI"),], subset=Species=="VIMI")
VIMIsite1 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="VIMI"),], subset=Species=="VIMI", weights = varIdent (form = ~1|Sitetag))
VIMIsite2 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="VIMI"),], subset=Species=="VIMI", weights = varFixed ( ~ PPTc))
VIMIsite3 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="VIMI"),], subset=Species=="VIMI", weights = varPower (form = ~ PPTc_scaled))
VIMIsite4 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="VIMI"),], subset=Species=="VIMI", weights = varPower (form = ~ PPTc))
AIC(VIMIsite, VIMIsite1, VIMIsite2, VIMIsite3, VIMIsite4) # VIMIsite is best with low df






AMYGwd <- fit.traits(dataz=traits.b, species="AMYG",trait="WD",diam=T, standBA=T)
anova(AMYGwd$PPT, AMYGwd$null) # p=0.0001 ***
anova(AMYGwd$DBH, AMYGwd$null) # p=0.872
anova(AMYGwd$StandBA, AMYGwd$null) # p=0.781


residplots(model=AMYGwd$PPT, dataz=traits.b[which(traits.b$Species=="AMYG" & traits.b$WD>0),]
           ,res.type = "pearson",mixed = T
           , Species = "AMYG", Trait="WD", vers=current_version
           , climvar = "PPTc", path=current_path)
scatter.smooth(abs(resid(AMYGwd$PPT))~round(AMYGwd$PPT@frame$PPTc_scaled,digits = 3))
# not much evidence for changing variance with drought

AMYGsite <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="AMYG")
AMYGsite1 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="AMYG", weights = varIdent (form = ~1|Sitetag))
AMYGsite2 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="AMYG", weights = varFixed( ~ PPTc))
AMYGsite3 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="AMYG", weights = varPower (form = ~ PPTc_scaled))
AMYGsite4 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0),], subset=Species=="AMYG", weights = varPower (form = ~ PPTc))

AIC(AMYGsite, AMYGsite1, AMYGsite2, AMYGsite3, AMYGsite4) # AMYGsite is best with low df
anova(AMYGsite, AMYGsite4) # p = 0.0797 , marginaly significant decrease



OBLIwd <- fit.traits(dataz=traits.b, species="OBLI",trait="WD",diam=T, standBA=T)
anova(OBLIwd$fert, OBLIwd$null) # p=0.043 *
anova(OBLIwd$PPT, OBLIwd$null) # p=0.1 .
anova(OBLIwd$DBH, OBLIwd$null) # p=0.183
anova(OBLIwd$StandBA, OBLIwd$null) # p=0.566

residplots(model=OBLIwd$fert, dataz=traits.b[which(traits.b$Species=="OBLI" & traits.b$WD>0),]
           ,res.type = "pearson",mixed = T
           , Species = "OBLI", Trait="WD", vers=current_version
           , climvar="PC1fert", path=current_path)

scatter.smooth(abs(resid(OBLIwd$PPT))~round(OBLIwd$PPT@frame$PPTc_scaled,digits = 3))
# huge withinplot and maybe within tree variance in lowest PET plot. But not much other evidence

OBLIsite <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0 & traits.b$Species=="OBLI"),], subset=Species=="OBLI")
OBLIsite1 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="OBLI"),], subset=Species=="OBLI", weights = varIdent (form = ~1|Sitetag))
OBLIsite2 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="OBLI"),], subset=Species=="OBLI", weights = varFixed (~ PPTc))
OBLIsite3 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="OBLI"),], subset=Species=="OBLI", weights = varPower (form = ~ PPTc_scaled))
OBLIsite4 <- lme(WD~Diameter + PPTc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$WD>0& traits.b$Species=="OBLI"),], subset=Species=="OBLI", weights = varPower (form = ~ PPTc))
AICc(OBLIsite, OBLIsite1, OBLIsite2, OBLIsite3, OBLIsite4) # OBLIsite is best with low df
anova(OBLIsite, OBLIsite1, OBLIsite3) # SUPER significant effect! p=1e-04
scatter.smooth(abs(resid(OBLIsite3,type = "normalized"))~round(OBLIsite3$data$PPTc_scaled,digits = 3))
## hump shaped increase at mid PPT best

#**** Very interesting ****








########### ** LDMC modeling ##########



ACACldmc <- fit.traits(dataz=traits.b, species="ACAC",trait="LDMC",diam=F, standBA=F)
anova(ACACldmc$MD, ACACldmc$null) # p=9.4e-06 ***


residplots(model=ACACldmc$MD, climvar= "MDc", dataz=traits.b[which(traits.b$Species=="ACAC" & traits.b$LDMC>0),]
           ,res.type = "pearson",mixed = T
           , Species = "ACAC", Trait="LDMC", vers=current_version, path=current_path)

scatter.smooth(abs(resid(ACACldmc$MD))~round(ACACldmc$MD@frame$MDc_scaled,digits = 3))

ACACsite <- lme(LDMC~ MDc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="ACAC")
ACACsite1 <- lme(LDMC~ MDc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="ACAC", weights = varIdent (form = ~1|Sitetag))
ACACsite2 <- lme(LDMC~ MDc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="ACAC", weights = varFixed(~ MDc))
ACACsite3 <- lme(LDMC~ MDc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="ACAC", weights = varPower (form = ~MDc_scaled))
ACACsite4 <- lme(LDMC~ MDc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="ACAC", weights = varPower (form = ~ MDc))
AIC(ACACsite, ACACsite1, ACACsite2, ACACsite3, ACACsite4) # so there's a difference bewteen sites, hump shaped
anova(ACACsite, ACACsite3) # 

ESALldmc <- fit.traits(dataz=traits.b, species="ESAL",trait="LDMC",diam=F, standBA=F)
anova(ESALldmc$PET, ESALldmc$null) # p=0.00036 ***
anova(ESALldmc$fert, ESALldmc$null) # p=0.07697 .

residplots(model=ESALldmc$PET,climvar="PETc", dataz=traits.b[which(traits.b$Species=="ESAL" & traits.b$LDMC>0),]
           ,res.type = "pearson",mixed = T
           , Species = "ESAL", Trait="LDMC", vers=current_version, path=current_path)


## PET 
ESALsite <- lme(LDMC~PETc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="ESAL")
ESALsite1 <- lme(LDMC~PETc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="ESAL", weights = varIdent (form = ~1|Sitetag))
ESALsite2 <- lme(LDMC~ PETc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0 & traits.b$Species=="ESAL"),], weights = varFixed (~PETc))
ESALsite3 <- lme(LDMC~ PETc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="ESAL", weights = varPower (form = ~ PETc_scaled))
ESALsite4 <- lme(LDMC~ PETc_scaled , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="ESAL", weights = varPower (form = ~ PETc))
AIC(ESALsite, ESALsite1, ESALsite2, ESALsite3,ESALsite4)
anova(ESALsite, ESALsite1) # sites different, symmetrical
anova(ESALsite, ESALsite2) # p=0.0257 for trough shaped, similar for increasing (increasing slightly better given fewer parameters)
scatter.smooth(abs(resid(ESALldmc$PET))~ESALldmc$PET@frame$PETc_scale)
plot(abs(resid(ESALsite2))~ESALsite2$data$pPETscale)



EMARGldmc <- fit.traits(dataz=traits.b, species="EMARG",trait="LDMC",diam=F, standBA=F)
anova(EMARGldmc$PET, EMARGldmc$null) # p=0.0457 *
anova(EMARGldmc$fert, EMARGldmc$null) # p=0.28
anova(EMARGldmc$DBH, EMARGldmc$null) # p= 0.78

bestvar <- "PETc"
residplots(model=EMARGldmc$PET,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="EMARG" & traits.b$LDMC>0),]
           ,res.type = "pearson",mixed = T
           , Species = "EMARG", Trait="LDMC", vers=current_version, path=current_path)
scatter.smooth(abs(resid(EMARGwd$PET))~round(EMARGwd$PET@frame[,paste0(bestvar,"_scaled")],digits = 3))

EMARGsite <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="EMARG")
EMARGsite1 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="EMARG", weights = varIdent (form = ~1|Sitetag))
EMARGsite2 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="EMARG", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
EMARGsite3 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="EMARG", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
EMARGsite4 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="EMARG", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
EMARGsite5 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="EMARG", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(EMARGsite, EMARGsite1, EMARGsite2, EMARGsite3,EMARGsite4, EMARGsite5)
anova(EMARGsite, EMARGsite3) # strong differences between sites, not climate related



COCAldmc <- fit.traits(dataz=traits.b, species="COCA",trait="LDMC",diam=F, standBA=F)
anova(COCAldmc$fert, COCAldmc$null) # p=0.008252 **
anova(COCAldmc$PPT, COCAldmc$null) # p=0.0307
anova(COCAldmc$DBH, COCAldmc$null) # p= 0.66

bestvar <- "PC1fert" # no climate relationship
residplots(model=COCAldmc$fert,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="COCA" & traits.b$LDMC>0),]
           ,res.type = "pearson",mixed = T
           , Species = "COCA", Trait="LDMC", vers=current_version, write.jpeg = T, path=current_path)


scatter.smooth(abs(resid(COCAldmc$PPT))~round(COCAldmc$PPT@frame$PPTc_scaled,digits = 3))

bestvar <- "PPTc" # best climate var rather than fert
COCAsite <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="COCA")
COCAsite1 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="COCA", weights = varIdent (form = ~1|Sitetag))
COCAsite2 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="COCA", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
COCAsite3 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
COCAsite4 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
COCAsite5 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(COCAsite, COCAsite1, COCAsite2, COCAsite3,COCAsite4, COCAsite5)
anova(COCAsite, COCAsite3) # strong differences between sites, not climate related







### Stats for Tasssie Species ###

OVATldmc <- fit.traits(dataz=traits.b, species="OVAT",trait="LDMC",diam=F, standBA=T)
anova(OVATldmc$fert, OVATldmc$null) # p=4.401e-06 ***
anova(OVATldmc$PPT, OVATldmc$null) # p=0.00017 ***
anova(OVATldmc$DBH, OVATldmc$null) # p= 0.416
anova(OVATldmc$StandBA, OVATldmc$null) # p= 0.30


# NOTE: a number of these models reduce the Plot random effect to 0, resulting in 'Singular' models.
# this could, in theory, influence thinks like LRT of these models. I refit the top models without the plot random effect just in case.
OVATnull <- lmer(LDMC~  1 + (1|Treetag), data = traits.b, subset=Species=="OVAT" , REML=F)
OVAT1 <- lmer(LDMC~PPTc_scaled  + (1|Treetag), data = traits.b, subset=Species=="OVAT", REML=F)
OVAT2 <- lmer(LDMC~PC1fert  + (1|Treetag), data = traits.b, subset=Species=="OVAT", REML=F)
AICc(OVATnull, OVAT1, OVAT2)
anova(OVATnull, OVAT2) # p= 7.294e-08 ***
# if anything, LRT with singular model UNDERESTIMATES significance.

bestvar <- "PC1fert" # no climate relationship
residplots(model=OVATldmc$fert,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="OVAT" & traits.b$LDMC>0),]
           ,res.type = "pearson",mixed = T
           , Species = "OVAT", Trait="LDMC", vers=current_version, write.jpeg = T, path=current_path)


scatter.smooth(abs(resid(OVATldmc$PPT))~round(OVATldmc$PPT@frame[,"PPTc_scaled"],digits = 3))
# strong visual decrease in resids with lower PPT
bestvar <- "PPTc"
OVATsite <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT")
OVATsite1 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT", weights = varIdent (form = ~1|Sitetag))
OVATsite2 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
OVATsite3 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
OVATsite4 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
OVATsite5 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(OVATsite, OVATsite1, OVATsite2, OVATsite3,OVATsite4, OVATsite5)
anova(OVATsite, OVATsite4) # strong differences between sites, decreasing at low PPT

# repeating with PC1fert
bestvar <- "PC1fert"
OVATsite <- lme(as.formula (paste("LDMC~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT")
OVATsite1 <- lme(as.formula (paste("LDMC~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT", weights = varIdent (form = ~1|Sitetag))
OVATsite2 <- lme(as.formula (paste("LDMC~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
OVATsite3 <- lme(as.formula (paste("LDMC~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~ I(",bestvar," + 5)", sep=""))))
OVATsite4 <- lme(as.formula (paste("LDMC~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OVAT", weights = varConstPower (form = as.formula(paste("~",bestvar, sep=""))))
AIC(OVATsite, OVATsite1, OVATsite2, OVATsite3,OVATsite4)
anova(OVATsite, OVATsite3) # strong differences, <.0001
# strong differences with fert as well, but varIdent best, followed by varConstPower



VIMIldmc <- fit.traits(dataz=traits.b, species="VIMI",trait="LDMC",diam=F, standBA=T)
anova(VIMIldmc$PPT, VIMIldmc$null) # p=1.996e-05 ***
anova(VIMIldmc$fert, VIMIldmc$null) # p=0.005791 **
anova(VIMIldmc$DBH, VIMIldmc$null) # p= 0.427
anova(VIMIldmc$StandBA, VIMIldmc$null) # p= 0.235

bestvar <- "PPTc" 
residplots(model=VIMIldmc$PPT,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="VIMI" & traits.b$LDMC>0),]
           ,res.type = "pearson",mixed = T
           , Species = "VIMI", Trait="LDMC", vers=current_version, path=current_path)


scatter.smooth(abs(resid(VIMIldmc$PPT))~round(VIMIldmc$PPT@frame$PPTc_scaled,digits = 3))

# maybe increasing
VIMIsite <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="VIMI")
VIMIsite1 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="VIMI", weights = varIdent (form = ~1|Sitetag))
VIMIsite2 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="VIMI", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
VIMIsite3 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="VIMI", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
VIMIsite4 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="VIMI", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
VIMIsite5 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="VIMI", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(VIMIsite, VIMIsite1, VIMIsite2, VIMIsite3,VIMIsite4, VIMIsite5)
anova(VIMIsite, VIMIsite4) # p=0.0086 decreasing at low PPT






AMYGldmc <- fit.traits(dataz=traits.b, species="AMYG",trait="LDMC",diam=F, standBA=T)
anova(AMYGldmc$PPT, AMYGldmc$null) # p=0.007 ***
anova(AMYGldmc$fert, AMYGldmc$null) # p=0.00186 **
anova(AMYGldmc$DBH, AMYGldmc$null) # p= 0.73
anova(AMYGldmc$StandBA, AMYGldmc$null) # p= 0.57

bestvar <- "PC1fert" 
residplots(model=AMYGldmc$PPT,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="AMYG" & traits.b$LDMC>0),]
           ,res.type = "pearson",mixed = T
           , Species = "AMYG", Trait="LDMC", vers=current_version, path=current_path)


scatter.smooth(abs(resid(AMYGldmc$PPT))~round(AMYGldmc$PPT@frame$PPTc_scaled,digits = 3))

bestvar <- "PPTc"
AMYGsite <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="AMYG")
AMYGsite1 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="AMYG", weights = varIdent (form = ~1|Sitetag))
AMYGsite2 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="AMYG", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
AMYGsite3 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="AMYG", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
AMYGsite4 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="AMYG", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
AMYGsite5 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="AMYG", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(AMYGsite, AMYGsite1, AMYGsite2, AMYGsite3,AMYGsite4, AMYGsite5)
anova(AMYGsite, AMYGsite4) # n.s.




OBLIldmc <- fit.traits(dataz=traits.b, species="OBLI",trait="LDMC",diam=F, standBA=T)
anova(OBLIldmc$MD, OBLIldmc$null) # p=2.11e-05 ***
anova(OBLIldmc$DBH, OBLIldmc$null) # p= 0.641
anova(OBLIldmc$StandBA, OBLIldmc$null) # p= 0.84

bestvar <- "MDc"
residplots(model=OBLIldmc$MD,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="OBLI" & traits.b$LDMC>0),]
           ,res.type = "pearson",mixed = T
           , Species = "OBLI", Trait="LDMC", vers=current_version, write.jpeg = T, path=current_path)

scatter.smooth(abs(resid(OBLIldmc$MD))~round(OBLIldmc$MD@frame$MDc_scaled,digits = 3))


OBLIsite <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OBLI")
OBLIsite1 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OBLI", weights = varIdent (form = ~1|Sitetag))
OBLIsite2 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OBLI", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
OBLIsite3 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OBLI", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
OBLIsite4 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OBLI", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
OBLIsite5 <- lme(as.formula (paste("LDMC~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LDMC>0),], subset=Species=="OBLI", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(OBLIsite, OBLIsite1, OBLIsite2, OBLIsite3,OBLIsite4, OBLIsite5)
anova(OBLIsite, OBLIsite1) # strong differences between sites, not climate related
anova(OBLIsite, OBLIsite4) # 0.0047, so not beautifully described by varPower, but decreasing
# decreasing with increasing MD





########### ** LMA modeling ##########



ACAClma <- fit.traits(dataz=traits.b, species="ACAC",trait="LMA",diam=F, standBA=F)
anova(ACAClma$MD, ACAClma$null) # 3.059e-08 ***
anova(ACAClma$fert, ACAClma$null) # 0.04121 *

bestvar <- "MDc"
residplots(model=ACAClma$MD,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="ACAC" & traits.b$LMA>0),]
           ,res.type = "pearson",mixed = T
           , Species = "ACACnew", Trait="LMA", vers=current_version, write.jpeg = T, path=current_path)

scatter.smooth(abs(resid(ACAClma$MD))~round(ACAClma$MD@frame$MDc_scaled,digits = 3))

ACACsite <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ACAC")
ACACsite1 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ACAC", weights = varIdent (form = ~1|Sitetag))
ACACsite2 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ACAC", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
ACACsite3 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ACAC", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
ACACsite4 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ACAC", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
ACACsite5 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ACAC", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(ACACsite, ACACsite1, ACACsite2, ACACsite3,ACACsite4, ACACsite5)
anova(ACACsite, ACACsite5) # n.s.



ESALlma <- fit.traits(dataz=traits.b, species="ESAL",trait="LMA",diam=F, standBA=F)
anova(ESALlma$MD, ESALlma$null) # 0.0156 *
anova(ESALlma$fert, ESALlma$null) # 0.3046


bestvar <- "MDc"
residplots(model=ESALlma$MD,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="ESAL" & traits.b$LMA>0),]
           ,res.type = "pearson",mixed = T
           , Species = "ESAL", Trait="LMA", vers=current_version, write.jpeg = T, path=current_path)


scatter.smooth(abs(resid(ESALlma$MD))~round(ESALlma$MD@frame$MDc_scaled,digits = 3))

ESALsite <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ESAL")
ESALsite1 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ESAL", weights = varIdent (form = ~1|Sitetag))
ESALsite2 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ESAL", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
ESALsite3 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ESAL", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
ESALsite4 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ESAL", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
ESALsite5 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="ESAL", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(ESALsite, ESALsite1, ESALsite2, ESALsite3,ESALsite4, ESALsite5)
anova(ESALsite, ESALsite5) # 0.0026, increasing with MD



EMARGlma <- fit.traits(dataz=traits.b, species="EMARG",trait="LMA",diam=F, standBA=F)
anova(EMARGlma$PPT, EMARGlma$null) # 0.00186 **
anova(EMARGlma$fert, EMARGlma$null) # 0.04848 *
anova(EMARGlma$DBH, EMARGlma$null) # p= 0.667

bestvar <- "PPTc"
residplots(model=EMARGlma$PPT,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="EMARG" & traits.b$LMA>0),]
           ,res.type = "pearson",mixed = T
           , Species = "EMARGnew", Trait="LMA", vers=current_version, write.jpeg = T, path=current_path)


scatter.smooth(abs(resid(EMARGlma$PPT))~round(EMARGlma$PPT@frame$PPTc_scaled,digits = 3))


EMARGsite <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="EMARG")
EMARGsite1 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="EMARG", weights = varIdent (form = ~1|Sitetag))
EMARGsite2 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="EMARG", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
EMARGsite3 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="EMARG", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
EMARGsite4 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="EMARG", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
EMARGsite5 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="EMARG", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(EMARGsite, EMARGsite1, EMARGsite2, EMARGsite3,EMARGsite4, EMARGsite5)
anova(EMARGsite, EMARGsite1) # 0.0034
anova(EMARGsite, EMARGsite3) # 0.0091 strong differences, hump shaped or non-climatic 




COCAlma <- fit.traits(dataz=traits.b, species="COCA",trait="LMA",diam=F, standBA=F)
anova(COCAlma$PPT, COCAlma$null) # 0.01138 *
anova(COCAlma$fert, COCAlma$null) # 0.002987 **
anova(COCAlma$DBH, COCAlma$null) # p= 0.41

bestvar <- "PC1fert" 
residplots(model=COCAlma$fert,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="COCA" & traits.b$LMA>0),]
           ,res.type = "pearson",mixed = T
           , Species = "COCAnew", Trait="LMA", vers=current_version, write.jpeg = F)

scatter.smooth(abs(resid(COCAlma$PPT))~round(COCAlma$PPT@frame$PPTc_scaled,digits = 3))

bestvar <- "PPTc" # 
COCAsite <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA")
COCAsite1 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA", weights = varIdent (form = ~1|Sitetag))
COCAsite2 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
COCAsite3 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
COCAsite4 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
COCAsite5 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(COCAsite, COCAsite1, COCAsite2, COCAsite3,COCAsite4, COCAsite5)
anova(COCAsite, COCAsite4) # n.s.
# same is true with fert:
# bestvar <- "PC1fert" # 
# COCAsite <- lme(as.formula (paste("LMA~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA")
# COCAsite1 <- lme(as.formula (paste("LMA~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA", weights = varIdent (form = ~1|Sitetag))
# COCAsite2 <- lme(as.formula (paste("LMA~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~ ",bestvar, sep=""))))
# COCAsite3 <- lme(as.formula (paste("LMA~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~ I(",bestvar," + 5)", sep=""))))
# COCAsite4 <- lme(as.formula (paste("LMA~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="COCA", weights = varConstPower (form = as.formula(paste("~ I(",bestvar," + 5)", sep=""))))
# AICc(COCAsite, COCAsite1, COCAsite2, COCAsite3,COCAsite4)
# anova(COCAsite, COCAsite2,COCAsite4) # no diffs?










### Stats for Tasssie Species ###

OVATlma <- fit.traits(dataz=traits.b, species="OVAT",trait="LMA",diam=F, standBA=T)
anova(OVATlma$PPT, OVATlma$null) # 0.00055 ***
anova(OVATlma$fert, OVATlma$null) # 4.703e-05 ***
anova(OVATlma$DBH, OVATlma$null) # p= 0.74
anova(OVATlma$StandBA, OVATlma$null) # p= 0.27


bestvar <- "PPTc" 
residplots(model=OVATlma$PPT,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="OVAT" & traits.b$LMA>0),]
           ,res.type = "pearson",mixed = T
           , Species = "OVAT", Trait="LMA", vers=current_version,write.jpeg = F)
bestvar <- "PC1fert" 
residplots(model=OVATlma$PPT,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="OVAT" & traits.b$LMA>0),]
           ,res.type = "pearson",mixed = T
           , Species = "OVAT", Trait="LMA", vers=current_version,write.jpeg = T, path=current_path)


scatter.smooth(abs(resid(OVATlma$PPT))~round(OVATlma$PPT@frame$PPTc_scaled,digits = 3))
# strong visual decrease in resids with lower PPT

bestvar <- "PPTc"
OVATsite <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OVAT")
OVATsite1 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OVAT", weights = varIdent (form = ~1|Sitetag))
OVATsite2 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OVAT", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
OVATsite3 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
OVATsite4 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
OVATsite5 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
OVATsite6 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OVAT", weights = varExp (form = as.formula(paste("~ I(",bestvar,"_scaled + 3)", sep=""))))
AICc(OVATsite, OVATsite1, OVATsite2, OVATsite3,OVATsite4, OVATsite5, OVATsite6)
anova(OVATsite, OVATsite5) # 0.0846 # n.s. for PC1fert







VIMIlma <- fit.traits(dataz=traits.b, species="VIMI",trait="LMA",diam=F, standBA=T)
anova(VIMIlma$MD, VIMIlma$null) # 0.05232
anova(VIMIlma$fert, VIMIlma$null) # 0.368
anova(VIMIlma$DBH, VIMIlma$null) # p= 0.38
anova(VIMIlma$StandBA, VIMIlma$null) # p= 0.727

bestvar <- "MDc" 
residplots(model=VIMIlma$MD,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="VIMI" & traits.b$LMA>0& traits.b$Branchtag != "VIMI-FREY-C-1-c"),]
           ,res.type = "pearson",mixed = T
           , Species = "VIMInew2", Trait="LMA", vers=current_version, write.jpeg = T, path=current_path)

scatter.smooth(abs(resid(VIMIlma$MD))~round(VIMIlma$MD@frame$MDc_scaled,digits = 3))

# maybe hump shaped?
VIMIsite <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="VIMI")
VIMIsite1 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="VIMI", weights = varIdent (form = ~1|Sitetag))
VIMIsite2 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="VIMI", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
VIMIsite3 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="VIMI", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
VIMIsite4 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="VIMI", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
VIMIsite5 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="VIMI", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(VIMIsite, VIMIsite1, VIMIsite2, VIMIsite3,VIMIsite4, VIMIsite5)
anova(VIMIsite, VIMIsite4) # n.s.






AMYGlma <- fit.traits(dataz=traits.b, species="AMYG",trait="LMA",diam=F, standBA=T)
anova(AMYGlma$PPTba,AMYGlma$PPT, AMYGlma$null) # 0.0003065 ***, + ba p=0.013757
anova(AMYGlma$fert, AMYGlma$null) # 0.009048 **
anova(AMYGlma$DBH, AMYGlma$null) # p= 0.93
anova(AMYGlma$StandBA, AMYGlma$null) # p= 0.02752 *


bestvar <- "PPTc" 
residplots(model=AMYGlma$PPT,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="AMYG" & traits.b$LMA>0),]
           ,res.type = "pearson",mixed = T
           , Species = "AMYG", Trait="LMA", vers=current_version,write.jpeg = T, path=current_path)

scatter.smooth(abs(resid(AMYGlma$PPT))~round(AMYGlma$PPT@frame$PPTc_scaled,digits = 3))

# hump shape
AMYGsite <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="AMYG")
AMYGsite1 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="AMYG", weights = varIdent (form = ~1|Sitetag))
AMYGsite2 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="AMYG", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
AMYGsite3 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="AMYG", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
AMYGsite4 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="AMYG", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
AMYGsite5 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="AMYG", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(AMYGsite, AMYGsite1, AMYGsite2, AMYGsite3,AMYGsite4, AMYGsite5)
anova(AMYGsite, AMYGsite1) # 0.0283 





OBLIlma <- fit.traits(dataz=traits.b, species="OBLI",trait="LMA",diam=F, standBA=T)
anova(OBLIlma$MDba,OBLIlma$MD, OBLIlma$null) # 0.003479 ** + 0.006 adding BA
anova(OBLIlma$fert, OBLIlma$null) # 0.007351 **
anova(OBLIlma$DBH, OBLIlma$null) # p= 0.15
anova(OBLIlma$StandBA, OBLIlma$null) # p= 0.0195 * ALSO decreasing LMA with stand BA

bestvar <- "MDc"
residplots(model=OBLIlma$MDba,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="OBLI" & traits.b$LMA>0),]
           ,res.type = "pearson",mixed = T
           , Species = "OBLI", Trait="LMA", vers=current_version, write.jpeg = T, path=current_path)

scatter.smooth(abs(resid(OBLIlma$MDba))~round(OBLIlma$MDba@frame$MDc_scaled,digits = 3))

OBLIsite <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OBLI")
OBLIsite1 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OBLI", weights = varIdent (form = ~1|Sitetag))
OBLIsite2 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OBLI", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
OBLIsite3 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OBLI", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
OBLIsite4 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OBLI", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
OBLIsite5 <- lme(as.formula (paste("LMA~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$LMA>0),], subset=Species=="OBLI", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(OBLIsite, OBLIsite1, OBLIsite2, OBLIsite3,OBLIsite4, OBLIsite5)
anova(OBLIsite, OBLIsite5) # 0.039, decrease at higher MD







########### ** log.hub modeling ##########



ACAChub <- fit.traits(dataz=traits.b, species="ACAC",trait="log.hub",diam=F, standBA=F)
anova(ACAChub$MD, ACAChub$null) # p=0.0259 *
anova(ACAChub$fert, ACAChub$null) # p=0.12


bestvar <- "MDc"
residplots(model=ACAChub$MD,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="ACAC" & traits.b$hub>0),]
           ,res.type = "pearson",mixed = T
           , Species = "ACAC", Trait="loghub", vers=current_version, write.jpeg = T, path=current_path)

scatter.smooth(abs(resid(ACAChub$MD))~round(ACAChub$MD@frame$MDc_scaled,digits = 3))

ACACsite <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ACAC")
ACACsite1 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ACAC", weights = varIdent (form = ~1|Sitetag))
ACACsite2 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ACAC", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
ACACsite3 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ACAC", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
ACACsite4 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ACAC", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
ACACsite5 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ACAC", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(ACACsite, ACACsite1, ACACsite2, ACACsite3,ACACsite4, ACACsite5)
anova(ACACsite, ACACsite3) # 0.0028, trough shaped with MD





ESALhub <- fit.traits(dataz=traits.b, species="ESAL",trait="log.hub",diam=F, standBA=F)
anova(ESALhub$MD, ESALhub$null) # p=0.02821 *
anova(ESALhub$fert, ESALhub$null) # p=0.30


bestvar <- "MDc"
residplots(model=ESALhub$MD,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="ESAL" & traits.b$hub>0),]
           ,res.type = "pearson",mixed = T
           , Species = "ESAL", Trait="loghub", vers=current_version, write.jpeg = T, path=current_path)


scatter.smooth(abs(resid(ESALhub$MD))~round(ESALhub$MD@frame$MDc_scaled,digits = 3))

ESALsite <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ESAL")
ESALsite1 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ESAL", weights = varIdent (form = ~1|Sitetag))
ESALsite2 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ESAL", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
ESALsite3 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ESAL", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
ESALsite4 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ESAL", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
ESALsite5 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="ESAL", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(ESALsite, ESALsite1, ESALsite2, ESALsite3,ESALsite4, ESALsite5)
anova(ESALsite, ESALsite3) #n.s.





EMARGhub <- fit.traits(dataz=traits.b, species="EMARG",trait="log.hub",diam=F, standBA=F)
anova(EMARGhub$PET, EMARGhub$null) # p= 0.2616
anova(EMARGhub$fert, EMARGhub$null) # p= 0.69
anova(EMARGhub$DBH, EMARGhub$null) # p= 0.42

bestvar <- "PETc"
residplots(model=EMARGhub$PET,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="EMARG" & traits.b$hub>0),]
           ,res.type = "pearson",mixed = T
           , Species = "EMARG", Trait="loghub", vers=current_version
           , write.jpeg = T, path=current_path)

scatter.smooth(abs(resid(EMARGhub$PET))~round(EMARGhub$PET@frame$PETc_scaled,digits = 3))

EMARGsite <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="EMARG")
EMARGsite1 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="EMARG", weights = varIdent (form = ~1|Sitetag))
EMARGsite2 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="EMARG", weights = varFixed ( as.formula(paste("~",bestvar, sep=""))))
EMARGsite3 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="EMARG", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
EMARGsite4 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="EMARG", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
EMARGsite5 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="EMARG", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(EMARGsite, EMARGsite1, EMARGsite2, EMARGsite3,EMARGsite4, EMARGsite5)
anova(EMARGsite, EMARGsite4) #0.06
anova(EMARGsite, EMARGsite1) # strong differences, hump shaped?




COCAhub <- fit.traits(dataz=traits.b, species="COCA",trait="log.hub",diam=F, standBA=F)
anova(COCAhub$depth, COCAhub$null) # p=0.003959 **
anova(COCAhub$PPT, COCAhub$null) # p=0.002484 **
anova(COCAhub$DBH, COCAhub$null) # p= 0.50

bestvar <- "PPTc" # 
residplots(model=COCAhub$depth,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="COCA" & traits.b$hub>0),]
           ,res.type = "pearson",mixed = T
           , Species = "COCA2", Trait="loghub", vers=current_version
           , write.jpeg = T, path=current_path)


scatter.smooth(abs(resid(COCAhub$PPT))~round(COCAhub$PPT@frame$PPTc_scaled,digits = 3))

bestvar <- "PPTc" # 
COCAsite <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="COCA")
COCAsite1 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="COCA", weights = varIdent (form = ~1|Sitetag))
COCAsite2 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="COCA", weights = varFixed ( as.formula(paste("~",bestvar, sep=""))))
COCAsite3 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
COCAsite4 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
COCAsite5 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="COCA", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(COCAsite, COCAsite1, COCAsite2, COCAsite3,COCAsite4, COCAsite5)
anova(COCAsite, COCAsite3) #n.s.








### Stats for Tasssie Species ###

OVAThub <- fit.traits(dataz=traits.b, species="OVAT",trait="log.hub",diam=F, standBA=T)
anova(OVAThub$PPT, OVAThub$null) # p=0.002893 **
anova(OVAThub$fert, OVAThub$null) # p=0.005059 **
anova(OVAThub$DBH, OVAThub$null) # p= 0.4689
anova(OVAThub$StandBA, OVAThub$null) # p= 0.81

bestvar <- "PPTc"
residplots(model=OVAThub$PPT,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="OVAT" & traits.b$hub>0),]
           ,res.type = "pearson",mixed = T
           , Species = "OVAT", Trait="loghub", vers=current_version
           , write.jpeg = T, path=current_path)

scatter.smooth(abs(resid(OVAThub$PPT))~round(OVAThub$PPT@frame$PPTc_scaled,digits = 3))

# strong visual decrease in resids with lower PPT
OVATsite <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OVAT")
OVATsite1 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OVAT", weights = varIdent (form = ~1|Sitetag))
OVATsite2 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OVAT", weights = varFixed ( as.formula(paste("~",bestvar, sep=""))))
OVATsite3 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
OVATsite4 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
OVATsite5 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OVAT", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(OVATsite, OVATsite1, OVATsite2, OVATsite3,OVATsite4, OVATsite5)
# anova(OVATsite, OVATsite3) #n.s.





VIMIhub <- fit.traits(dataz=traits.b, species="VIMI",trait="log.hub",diam=F, standBA=T)
anova(VIMIhub$fert, VIMIhub$null) # p=6.81e-05 ***
anova(VIMIhub$PPT, VIMIhub$null) # p=0.00099 **
anova(VIMIhub$DBH, VIMIhub$null) # p= 0.1297
anova(VIMIhub$StandBA, VIMIhub$null) # p= 0.911

bestvar <- "PC1fert" # 
residplots(model=VIMIhub$fert,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="VIMI" & traits.b$hub>0),]
           ,res.type = "pearson",mixed = T
           , Species = "VIMI", Trait="loghub", vers=current_version
           , write.jpeg = T, path=current_path)

scatter.smooth(abs(resid(VIMIhub$PPT))~round(VIMIhub$PPT@frame$PPTc_scaled,digits = 3))

# maybe hump shaped?
VIMIsite <- lme(as.formula (paste("log.hub~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$hub>0),], subset=Species=="VIMI")
VIMIsite1 <- lme(as.formula (paste("log.hub~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$hub>0),], subset=Species=="VIMI", weights = varIdent (form = ~1|Sitetag))
VIMIsite2 <- lme(as.formula (paste("log.hub~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$hub>0),], subset=Species=="VIMI", weights = varPower (form = as.formula(paste("~ ",bestvar, sep=""))))
VIMIsite4 <- lme(as.formula (paste("log.hub~ ",bestvar, sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[which(traits.b$hub>0),], subset=Species=="VIMI", weights = varConstPower (form = as.formula(paste("~ ",bestvar, sep=""))))
AIC(VIMIsite, VIMIsite1, VIMIsite2,VIMIsite4)
anova(VIMIsite, VIMIsite2) # no evidnece

bestvar <- "PPTc"
VIMIsite <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="VIMI")
VIMIsite1 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="VIMI", weights = varIdent (form = ~1|Sitetag))
VIMIsite2 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="VIMI", weights = varFixed ( as.formula(paste("~",bestvar, sep=""))))
VIMIsite3 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="VIMI", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
VIMIsite4 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="VIMI", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
VIMIsite5 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="VIMI", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(VIMIsite, VIMIsite1, VIMIsite2, VIMIsite3,VIMIsite4, VIMIsite5)
anova(VIMIsite, VIMIsite3) #n.s.




AMYGhub <- fit.traits(dataz=traits.b, species="AMYG",trait="log.hub",diam=F, standBA=T)
anova(AMYGhub$PPT, AMYGhub$null) # p=5.1e-06 ***
anova(AMYGhub$fert, AMYGhub$null) # p=0.0002024 ***
anova(AMYGhub$DBH, AMYGhub$null) # p= 0.60
anova(AMYGhub$StandBA, AMYGhub$null) # p= 0.74

bestvar <- "PPTc"
residplots(model=AMYGhub$PPT,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="AMYG" & traits.b$hub>0),]
           ,res.type = "pearson",mixed = T
           , Species = "AMYG", Trait="loghub", vers=current_version
           , write.jpeg = T, path=current_path)


scatter.smooth(abs(resid(AMYGhub$PPT))~round(AMYGhub$PPT@frame$PPTc_scaled,digits = 3))

# strong visual decrease in resids with lower PPT
AMYGsite <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="AMYG")
AMYGsite1 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="AMYG", weights = varIdent (form = ~1|Sitetag))
AMYGsite2 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="AMYG", weights = varFixed ( as.formula(paste("~",bestvar,"_scaled", sep=""))))
AMYGsite3 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="AMYG", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
AMYGsite4 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="AMYG", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
AMYGsite5 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="AMYG", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(AMYGsite, AMYGsite1, AMYGsite2, AMYGsite3,AMYGsite4, AMYGsite5)
anova(AMYGsite, AMYGsite3) # 0.0066 trough shaped





OBLIhub <- fit.traits(dataz=traits.b, species="OBLI",trait="log.hub",diam=F, standBA=T)
anova(OBLIhub$depth, OBLIhub$null) # p=0.001977 **
anova(OBLIhub$PPT, OBLIhub$null) # p=0.04235 *
anova(OBLIhub$DBH, OBLIhub$null) # p= 0.369
anova(OBLIhub$StandBA, OBLIhub$null) # p= 0.5388
# lots of singular errors because tree random effect estimated as 0. refitting just to check
OBLIhub1 <- lmer(log.hub ~ PC2depth + (1 | Plottag), traits.b, subset=Species=="OBLI", REML=F)
OBLIhub2 <- lmer(log.hub ~ MDc + (1 | Plottag), traits.b, subset=Species=="OBLI", REML=F)
OBLIhub3 <- lmer(log.hub ~ 1 + (1 | Plottag), traits.b, subset=Species=="OBLI", REML=F)
AICc(OBLIhub1, OBLIhub2, OBLIhub3)
# still holds

bestvar <- "PC2depth"
residplots(model=OBLIhub$depth,climvar=bestvar, dataz=traits.b[which(traits.b$Species=="OBLI" & traits.b$hub>0),]
           ,res.type = "pearson",mixed = T
           , Species = "OBLI2", Trait="loghub", vers=current_version
           , write.jpeg = T, path=current_path)

scatter.smooth(abs(resid(OBLIhub$depth))~round(OBLIhub$depth@frame$PC2depth,digits = 3))

bestvar <- "PPTc"
OBLIsite <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OBLI")
OBLIsite1 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OBLI", weights = varIdent (form = ~1|Sitetag))
OBLIsite2 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OBLI", weights = varFixed ( as.formula(paste("~",bestvar, sep=""))))
OBLIsite3 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OBLI", weights = varPower (form = as.formula(paste("~",bestvar,"_scaled", sep=""))))
OBLIsite4 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OBLI", weights = varPower (form = as.formula(paste("~ I(",bestvar,"_scaled + 5)", sep=""))))
OBLIsite5 <- lme(as.formula (paste("log.hub~ ",bestvar,"_scaled", sep="")) , random = ~ 1 | Plottag/Treetag, data = traits.b[-which(is.na(traits.b$log.hub)),], subset=Species=="OBLI", weights = varPower (form = as.formula(paste("~",bestvar, sep=""))))
AICc(OBLIsite, OBLIsite1, OBLIsite2, OBLIsite3,OBLIsite4, OBLIsite5)
anova(OBLIsite, OBLIsite5) #n.s.



################# END: TRAIT-CLIMATE ####################################














##### ++++++++++++++++++++++++++++++++++++++ #############
################### BEGIN: TRAIT COORDINATION #####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#### Function for fitting MARs

fit.MAR <- function(xvar, yvar, data, method="SMA") {
  if(method=="SMA") meth = 3
  if(method=="MA") meth = 2
  if(method=="OLS") meth =1
  if(length(t(as.vector(data[!(is.na(data[,xvar]) | is.na(data[,yvar])), xvar])))<3){
    return(rep(NA, times=9))
  }
  else{
    if(var(data[,yvar], na.rm=T)==0 | var(data[,xvar],na.rm=T)==0){
      return(rep(NA, times=9))
    }
    else{
      tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), data))
      intercept <- tmp.mod$regression.results$Intercept[meth]
      slope <- tmp.mod$regression.results$Slope[meth]
      slope.lci <- tmp.mod$confidence.intervals[meth,4]
      slope.uci <- tmp.mod$confidence.intervals[meth,5]
      rho <- tmp.mod$r
      r.sq <- tmp.mod$rsquare
      n <- tmp.mod$n
      varX <- var(data[,xvar], na.rm=T)
      varY <- var(data[,yvar], na.rm=T)
      results <- c(intercept, slope,slope.lci, slope.uci, rho, r.sq, n, varX, varY)
      return(results)
    }
  }
}



# calculate significance of a pearson correlation coefficient
rho.sig <- function(rho, n){
  t <- rho/sqrt((1-rho^2)/(n-2))
  pt(-abs(t), n-2)*2
}






############# + WD vs LDMC ###########

###_______________ between Site analysis _______________________________________________


ptm <- proc.time()
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(traits.s$Species)), ncol=18))
colnames(spp.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")

for(i in 1:length(unique(traits.s$Species))){
  species <- levels(traits.s$Species)[i]
  print(species)
  dataz <- traits.s[which(traits.s$Species==species),]
  res <- fit.MAR(xvar='WD',yvar="LDMC",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  spp.results[i,18]<- species
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = niters)  
  #   spp.results[i, 11:16] <- nullbounds
  #   spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  # }
}
proc.time()-ptm
spp.results$sig <- rho.sig(n = spp.results$n, rho = spp.results$Rho)




# Seems to be working!!!

###_______________ Plots w/in Sites level analysis _______________________________________________
t0 <- proc.time()
plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.p$Sitetag)), ncol=18))
colnames(plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.p$Sitetag))){
  Sitetag <- levels(traits.p$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.p[which(traits.p$Sitetag==Sitetag),]
  plot.results[i,1] <- Sitetag
  plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='WD',yvar="LDMC",data=dataz)
  plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   plot.results[i, 11:16] <- nullbounds
  #   plot.results[i,17] <- test.sig(x=plot.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
plot.results$sig <- rho.sig(n = plot.results$n, rho = plot.results$Rho)


###_______________ Trees w/in Sites level analysis _______________________________________________
t0 <- proc.time()
tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Sitetag)), ncol=18))
colnames(tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Sitetag))){
  Sitetag <- levels(traits.t$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.t[which(traits.t$Sitetag==Sitetag),]
  tree.results[i,1] <- Sitetag
  tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='WD',yvar="LDMC",data=dataz)
  tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.results$sig <- rho.sig(n = tree.results$n, rho = tree.results$Rho)

###_______________ Trees w/in Plot level analysis _______________________________________________
t0 <- proc.time()
tree.in.plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Plottag)), ncol=18))
colnames(tree.in.plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Plottag))){
  Plottag <- levels(traits.t$Plottag)[i]
  print(Plottag)
  dataz <- traits.t[which(traits.t$Plottag==Plottag),]
  tree.in.plot.results[i,1] <- Plottag
  tree.in.plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='WD',yvar="LDMC",data=dataz)
  tree.in.plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.in.plot.results$sig <- rho.sig(n = tree.in.plot.results$n, rho = tree.in.plot.results$Rho)



###_______________ Branches w/in Plot level analysis _______________________________________________
t0 <- proc.time()
branch.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Plottag)), ncol=18))
colnames(branch.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Plottag))){
  Plottag <- levels(traits.b$Plottag)[i]
  print(Plottag)
  dataz <- traits.b[which(traits.b$Plottag==Plottag),]
  branch.results[i,1] <- Plottag
  branch.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='WD',yvar="LDMC",data=dataz)
  branch.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.results$sig <- rho.sig(n = branch.results$n, rho = branch.results$Rho)


###_______________ branches w/in Tree level analysis _______________________________________________
t0 <- proc.time()
branch.in.tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Treetag)), ncol=18))
colnames(branch.in.tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Treetag))){
  Treetag <- levels(traits.b$Treetag)[i]
  print(Treetag)
  dataz <- traits.b[which(traits.b$Treetag==Treetag),]
  branch.in.tree.results[i,1] <- Treetag
  branch.in.tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='WD',yvar="LDMC",data=dataz)
  branch.in.tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.in.tree.results$sig <- rho.sig(n = branch.in.tree.results$n, rho = branch.in.tree.results$Rho)


spp.results$Level <- "Btw-site"
plot.results$Level <- "Btw-plot.in.site"
tree.results$Level <- "Btw-tree.in.site"
branch.results$Level<- "Btw-branch.in.plot"
branch.in.tree.results$Level <- "Btw-branch.in.tree"
tree.in.plot.results$Level <- "Btw-tree.in.plot"


WD.LDMC.spp.results <- spp.results
WD.LDMC.plot.results <- plot.results
WD.LDMC.tree.results <- tree.results
WD.LDMC.branch.results<- branch.results
WD.LDMC.branch.in.tree.results <- branch.in.tree.results
WD.LDMC.tree.in.plot.results <- tree.in.plot.results

WD.LDMC.results <- rbind(WD.LDMC.spp.results, WD.LDMC.plot.results, WD.LDMC.tree.results, WD.LDMC.branch.results)





############# + WD vs LMA ###########

###_______________ between Site analysis _______________________________________________


ptm <- proc.time()
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(traits.s$Species)), ncol=18))
colnames(spp.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.s$Species))){
  species <- levels(traits.s$Species)[i]
  print(species)
  dataz <- traits.s[which(traits.s$Species==species),]
  res <- fit.MAR(xvar='WD',yvar="LMA",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  spp.results[i,18] <- species
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = niters)  
  #   spp.results[i, 11:16] <- nullbounds
  #   spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  # }
}
proc.time()-ptm
spp.results$sig <- rho.sig(n = spp.results$n, rho = spp.results$Rho)




# Seems to be working!!!

###_______________ Plots w/in Sites level analysis _______________________________________________
t0 <- proc.time()
plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.p$Sitetag)), ncol=18))
colnames(plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.p$Sitetag))){
  Sitetag <- levels(traits.p$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.p[which(traits.p$Sitetag==Sitetag),]
  plot.results[i,1] <- Sitetag
  plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='WD',yvar="LMA",data=dataz)
  plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   plot.results[i, 11:16] <- nullbounds
  #   plot.results[i,17] <- test.sig(x=plot.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
plot.results$sig <- rho.sig(n = plot.results$n, rho = plot.results$Rho)


###_______________ Trees w/in Sites level analysis _______________________________________________
t0 <- proc.time()
tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Sitetag)), ncol=18))
colnames(tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Sitetag))){
  Sitetag <- levels(traits.t$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.t[which(traits.t$Sitetag==Sitetag),]
  tree.results[i,1] <- Sitetag
  tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='WD',yvar="LMA",data=dataz)
  tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.results$sig <- rho.sig(n = tree.results$n, rho = tree.results$Rho)

###_______________ Trees w/in Plot level analysis _______________________________________________
t0 <- proc.time()
tree.in.plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Plottag)), ncol=18))
colnames(tree.in.plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Plottag))){
  Plottag <- levels(traits.t$Plottag)[i]
  print(Plottag)
  dataz <- traits.t[which(traits.t$Plottag==Plottag),]
  tree.in.plot.results[i,1] <- Plottag
  tree.in.plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='WD',yvar="LMA",data=dataz)
  tree.in.plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.in.plot.results$sig <- rho.sig(n = tree.in.plot.results$n, rho = tree.in.plot.results$Rho)



###_______________ Branches w/in Plot level analysis _______________________________________________
t0 <- proc.time()
branch.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Plottag)), ncol=18))
colnames(branch.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Plottag))){
  Plottag <- levels(traits.b$Plottag)[i]
  print(Plottag)
  dataz <- traits.b[which(traits.b$Plottag==Plottag),]
  branch.results[i,1] <- Plottag
  branch.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='WD',yvar="LMA",data=dataz)
  branch.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.results$sig <- rho.sig(n = branch.results$n, rho = branch.results$Rho)


###_______________ branches w/in Tree level analysis _______________________________________________
t0 <- proc.time()
branch.in.tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Treetag)), ncol=18))
colnames(branch.in.tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varWD","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Treetag))){
  Treetag <- levels(traits.b$Treetag)[i]
  print(Treetag)
  dataz <- traits.b[which(traits.b$Treetag==Treetag),]
  branch.in.tree.results[i,1] <- Treetag
  branch.in.tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='WD',yvar="LMA",data=dataz)
  branch.in.tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.in.tree.results$sig <- rho.sig(n = branch.in.tree.results$n, rho = branch.in.tree.results$Rho)


spp.results$Level <- "Btw-site"
plot.results$Level <- "Btw-plot.in.site"
tree.results$Level <- "Btw-tree.in.site"
branch.results$Level<- "Btw-branch.in.plot"
branch.in.tree.results$Level <- "Btw-branch.in.tree"
tree.in.plot.results$Level <- "Btw-tree.in.plot"


WD.LMA.spp.results <- spp.results
WD.LMA.plot.results <- plot.results
WD.LMA.tree.results <- tree.results
WD.LMA.branch.results<- branch.results
WD.LMA.branch.in.tree.results <- branch.in.tree.results
WD.LMA.tree.in.plot.results <- tree.in.plot.results

WD.LMA.results <- rbind(WD.LMA.spp.results, WD.LMA.plot.results, WD.LMA.tree.results, WD.LMA.branch.results)




############# + LDMC vs LMA ###########

###_______________ between Site analysis _______________________________________________


ptm <- proc.time()
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(traits.s$Species)), ncol=17))
colnames(spp.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLDMC","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(traits.s$Species))){
  species <- levels(traits.s$Species)[i]
  print(species)
  dataz <- traits.s[which(traits.s$Species==species),]
  res <- fit.MAR(xvar='LDMC',yvar="LMA",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = niters)  
  #   spp.results[i, 11:16] <- nullbounds
  #   spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  # }
}
proc.time()-ptm
spp.results$Species <- spp.results$ID
spp.results$sig <- rho.sig(n = spp.results$n, rho = spp.results$Rho)




# Seems to be working!!!

###_______________ Plots w/in Sites level analysis _______________________________________________
t0 <- proc.time()
plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.p$Sitetag)), ncol=18))
colnames(plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLDMC","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.p$Sitetag))){
  Sitetag <- levels(traits.p$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.p[which(traits.p$Sitetag==Sitetag),]
  plot.results[i,1] <- Sitetag
  plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='LDMC',yvar="LMA",data=dataz)
  plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   plot.results[i, 11:16] <- nullbounds
  #   plot.results[i,17] <- test.sig(x=plot.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
plot.results$sig <- rho.sig(n = plot.results$n, rho = plot.results$Rho)


###_______________ Trees w/in Sites level analysis _______________________________________________
t0 <- proc.time()
tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Sitetag)), ncol=18))
colnames(tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLDMC","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Sitetag))){
  Sitetag <- levels(traits.t$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.t[which(traits.t$Sitetag==Sitetag),]
  tree.results[i,1] <- Sitetag
  tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='LDMC',yvar="LMA",data=dataz)
  tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.results$sig <- rho.sig(n = tree.results$n, rho = tree.results$Rho)

###_______________ Trees w/in Plot level analysis _______________________________________________
t0 <- proc.time()
tree.in.plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Plottag)), ncol=18))
colnames(tree.in.plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLDMC","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Plottag))){
  Plottag <- levels(traits.t$Plottag)[i]
  print(Plottag)
  dataz <- traits.t[which(traits.t$Plottag==Plottag),]
  tree.in.plot.results[i,1] <- Plottag
  tree.in.plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='LDMC',yvar="LMA",data=dataz)
  tree.in.plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.in.plot.results$sig <- rho.sig(n = tree.in.plot.results$n, rho = tree.in.plot.results$Rho)



###_______________ Branches w/in Plot level analysis _______________________________________________
t0 <- proc.time()
branch.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Plottag)), ncol=18))
colnames(branch.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLDMC","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Plottag))){
  Plottag <- levels(traits.b$Plottag)[i]
  print(Plottag)
  dataz <- traits.b[which(traits.b$Plottag==Plottag),]
  branch.results[i,1] <- Plottag
  branch.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='LDMC',yvar="LMA",data=dataz)
  branch.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.results$sig <- rho.sig(n = branch.results$n, rho = branch.results$Rho)


###_______________ branches w/in Tree level analysis _______________________________________________
t0 <- proc.time()
branch.in.tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Treetag)), ncol=18))
colnames(branch.in.tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLDMC","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Treetag))){
  Treetag <- levels(traits.b$Treetag)[i]
  print(Treetag)
  dataz <- traits.b[which(traits.b$Treetag==Treetag),]
  branch.in.tree.results[i,1] <- Treetag
  branch.in.tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='LDMC',yvar="LMA",data=dataz)
  branch.in.tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.in.tree.results$sig <- rho.sig(n = branch.in.tree.results$n, rho = branch.in.tree.results$Rho)


spp.results$Level <- "Btw-site"
plot.results$Level <- "Btw-plot.in.site"
tree.results$Level <- "Btw-tree.in.site"
branch.results$Level<- "Btw-branch.in.plot"
branch.in.tree.results$Level <- "Btw-branch.in.tree"
tree.in.plot.results$Level <- "Btw-tree.in.plot"


LDMC.LMA.spp.results <- spp.results
LDMC.LMA.plot.results <- plot.results
LDMC.LMA.tree.results <- tree.results
LDMC.LMA.branch.results<- branch.results
LDMC.LMA.branch.in.tree.results <- branch.in.tree.results
LDMC.LMA.tree.in.plot.results <- tree.in.plot.results

LDMC.LMA.results <- rbind(LDMC.LMA.spp.results, LDMC.LMA.plot.results, LDMC.LMA.tree.results, LDMC.LMA.branch.results)











#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

############# + log.hub vs LMA ###########

###_______________ between Site analysis _______________________________________________


ptm <- proc.time()
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(traits.s$Species)), ncol=17))
colnames(spp.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(traits.s$Species))){
  species <- levels(traits.s$Species)[i]
  print(species)
  dataz <- traits.s[which(traits.s$Species==species),]
  res <- fit.MAR(xvar='log.hub',yvar="LMA",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = niters)  
  #   spp.results[i, 11:16] <- nullbounds
  #   spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  # }
}
proc.time()-ptm
spp.results$Species <- spp.results$ID
spp.results$sig <- rho.sig(n = spp.results$n, rho = spp.results$Rho)




# Seems to be working!!!

###_______________ Plots w/in Sites level analysis _______________________________________________
t0 <- proc.time()
plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.p$Sitetag)), ncol=18))
colnames(plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.p$Sitetag))){
  Sitetag <- levels(traits.p$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.p[which(traits.p$Sitetag==Sitetag),]
  plot.results[i,1] <- Sitetag
  plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="LMA",data=dataz)
  plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   plot.results[i, 11:16] <- nullbounds
  #   plot.results[i,17] <- test.sig(x=plot.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
plot.results$sig <- rho.sig(n = plot.results$n, rho = plot.results$Rho)


###_______________ Trees w/in Sites level analysis _______________________________________________
t0 <- proc.time()
tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Sitetag)), ncol=18))
colnames(tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Sitetag))){
  Sitetag <- levels(traits.t$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.t[which(traits.t$Sitetag==Sitetag),]
  tree.results[i,1] <- Sitetag
  tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="LMA",data=dataz)
  tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.results$sig <- rho.sig(n = tree.results$n, rho = tree.results$Rho)

###_______________ Trees w/in Plot level analysis _______________________________________________
t0 <- proc.time()
tree.in.plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Plottag)), ncol=18))
colnames(tree.in.plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Plottag))){
  Plottag <- levels(traits.t$Plottag)[i]
  print(Plottag)
  dataz <- traits.t[which(traits.t$Plottag==Plottag),]
  tree.in.plot.results[i,1] <- Plottag
  tree.in.plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="LMA",data=dataz)
  tree.in.plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.in.plot.results$sig <- rho.sig(n = tree.in.plot.results$n, rho = tree.in.plot.results$Rho)



###_______________ Branches w/in Plot level analysis _______________________________________________
t0 <- proc.time()
branch.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Plottag)), ncol=18))
colnames(branch.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Plottag))){
  Plottag <- levels(traits.b$Plottag)[i]
  print(Plottag)
  dataz <- traits.b[which(traits.b$Plottag==Plottag),]
  branch.results[i,1] <- Plottag
  branch.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="LMA",data=dataz)
  branch.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.results$sig <- rho.sig(n = branch.results$n, rho = branch.results$Rho)


###_______________ branches w/in Tree level analysis _______________________________________________
t0 <- proc.time()
branch.in.tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Treetag)), ncol=18))
colnames(branch.in.tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLMA", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Treetag))){
  Treetag <- levels(traits.b$Treetag)[i]
  print(Treetag)
  dataz <- traits.b[which(traits.b$Treetag==Treetag),]
  branch.in.tree.results[i,1] <- Treetag
  branch.in.tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="LMA",data=dataz)
  branch.in.tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.in.tree.results$sig <- rho.sig(n = branch.in.tree.results$n, rho = branch.in.tree.results$Rho)


spp.results$Level <- "Btw-site"
plot.results$Level <- "Btw-plot.in.site"
tree.results$Level <- "Btw-tree.in.site"
branch.results$Level<- "Btw-branch.in.plot"
branch.in.tree.results$Level <- "Btw-branch.in.tree"
tree.in.plot.results$Level <- "Btw-tree.in.plot"


log.hub.LMA.spp.results <- spp.results
log.hub.LMA.plot.results <- plot.results
log.hub.LMA.tree.results <- tree.results
log.hub.LMA.branch.results<- branch.results
log.hub.LMA.branch.in.tree.results <- branch.in.tree.results
log.hub.LMA.tree.in.plot.results <- tree.in.plot.results

log.hub.LMA.results <- rbind(log.hub.LMA.spp.results, log.hub.LMA.plot.results, log.hub.LMA.tree.results, log.hub.LMA.branch.results)

quartz(width=5, height=4)
boxplot(Rho~Level, spp.results, xlim=c(0.5,3),ylim=c(-1,1), at=1, main="LMA vs log.hub", ylab="Rho")
boxplot(Rho~Level, plot.results, at=1.5, add=T)
boxplot(Rho~Level, tree.results, at=1.5, add=T, border="blue")
boxplot(Rho~Level, tree.in.plot.results, at=2, add=T)
boxplot(Rho~Level, branch.results, at=2, add=T, border="blue")
boxplot(Rho~Level, branch.in.tree.results, at=2.5, add=T)
abline(h=0)


quartz(width=5, height=4)
boxplot(sig~Level, spp.results, xlim=c(0.5,3),ylim=c(0,1), at=1, main="LMA vs log.hub", ylab="sig")
boxplot(sig~Level, plot.results, at=1.5, add=T)
boxplot(sig~Level, tree.results, at=1.5, add=T, border="blue")
boxplot(sig~Level, tree.in.plot.results, at=2, add=T)
boxplot(sig~Level, branch.results, at=2, add=T, border="blue")
boxplot(sig~Level, branch.in.tree.results, at=2.5, add=T)
abline(h=0.05)











############# + log.hub vs LDMC ###########

###_______________ between Site analysis _______________________________________________


ptm <- proc.time()
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(traits.s$Species)), ncol=17))
colnames(spp.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(traits.s$Species))){
  species <- levels(traits.s$Species)[i]
  print(species)
  dataz <- traits.s[which(traits.s$Species==species),]
  res <- fit.MAR(xvar='log.hub',yvar="LDMC",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
  #   nullbounds <- fit.null(xvar='log.LDMC', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = niters)  
  #   spp.results[i, 11:16] <- nullbounds
  #   spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  # }
}
proc.time()-ptm
spp.results$Species <- spp.results$ID
spp.results$sig <- rho.sig(n = spp.results$n, rho = spp.results$Rho)




# Seems to be working!!!

###_______________ Plots w/in Sites level analysis _______________________________________________
t0 <- proc.time()
plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.p$Sitetag)), ncol=18))
colnames(plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.p$Sitetag))){
  Sitetag <- levels(traits.p$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.p[which(traits.p$Sitetag==Sitetag),]
  plot.results[i,1] <- Sitetag
  plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="LDMC",data=dataz)
  plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LDMC', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   plot.results[i, 11:16] <- nullbounds
  #   plot.results[i,17] <- test.sig(x=plot.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
plot.results$sig <- rho.sig(n = plot.results$n, rho = plot.results$Rho)


###_______________ Trees w/in Sites level analysis _______________________________________________
t0 <- proc.time()
tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Sitetag)), ncol=18))
colnames(tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Sitetag))){
  Sitetag <- levels(traits.t$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.t[which(traits.t$Sitetag==Sitetag),]
  tree.results[i,1] <- Sitetag
  tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="LDMC",data=dataz)
  tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LDMC', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.results$sig <- rho.sig(n = tree.results$n, rho = tree.results$Rho)

###_______________ Trees w/in Plot level analysis _______________________________________________
t0 <- proc.time()
tree.in.plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Plottag)), ncol=18))
colnames(tree.in.plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Plottag))){
  Plottag <- levels(traits.t$Plottag)[i]
  print(Plottag)
  dataz <- traits.t[which(traits.t$Plottag==Plottag),]
  tree.in.plot.results[i,1] <- Plottag
  tree.in.plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="LDMC",data=dataz)
  tree.in.plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LDMC', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.in.plot.results$sig <- rho.sig(n = tree.in.plot.results$n, rho = tree.in.plot.results$Rho)



###_______________ Branches w/in Plot level analysis _______________________________________________
t0 <- proc.time()
branch.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Plottag)), ncol=18))
colnames(branch.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Plottag))){
  Plottag <- levels(traits.b$Plottag)[i]
  print(Plottag)
  dataz <- traits.b[which(traits.b$Plottag==Plottag),]
  branch.results[i,1] <- Plottag
  branch.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="LDMC",data=dataz)
  branch.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LDMC', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.results$sig <- rho.sig(n = branch.results$n, rho = branch.results$Rho)


###_______________ branches w/in Tree level analysis _______________________________________________
t0 <- proc.time()
branch.in.tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Treetag)), ncol=18))
colnames(branch.in.tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varLDMC", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Treetag))){
  Treetag <- levels(traits.b$Treetag)[i]
  print(Treetag)
  dataz <- traits.b[which(traits.b$Treetag==Treetag),]
  branch.in.tree.results[i,1] <- Treetag
  branch.in.tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="LDMC",data=dataz)
  branch.in.tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.LDMC', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.in.tree.results$sig <- rho.sig(n = branch.in.tree.results$n, rho = branch.in.tree.results$Rho)


spp.results$Level <- "Btw-site"
plot.results$Level <- "Btw-plot.in.site"
tree.results$Level <- "Btw-tree.in.site"
branch.results$Level<- "Btw-branch.in.plot"
branch.in.tree.results$Level <- "Btw-branch.in.tree"
tree.in.plot.results$Level <- "Btw-tree.in.plot"


log.hub.LDMC.spp.results <- spp.results
log.hub.LDMC.plot.results <- plot.results
log.hub.LDMC.tree.results <- tree.results
log.hub.LDMC.branch.results<- branch.results
log.hub.LDMC.branch.in.tree.results <- branch.in.tree.results
log.hub.LDMC.tree.in.plot.results <- tree.in.plot.results

log.hub.LDMC.results <- rbind(log.hub.LDMC.spp.results, log.hub.LDMC.plot.results, log.hub.LDMC.tree.results, log.hub.LDMC.branch.results)













############# + log.hub vs WD ###########

###_______________ between Site analysis _______________________________________________


ptm <- proc.time()
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(traits.s$Species)), ncol=17))
colnames(spp.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varWD", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(traits.s$Species))){
  species <- levels(traits.s$Species)[i]
  print(species)
  dataz <- traits.s[which(traits.s$Species==species),]
  res <- fit.MAR(xvar='log.hub',yvar="WD",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
  #   nullbounds <- fit.null(xvar='log.WD', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = niters)  
  #   spp.results[i, 11:16] <- nullbounds
  #   spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  # }
}
proc.time()-ptm
spp.results$Species <- spp.results$ID
spp.results$sig <- rho.sig(n = spp.results$n, rho = spp.results$Rho)




# Seems to be working!!!

###_______________ Plots w/in Sites level analysis _______________________________________________
t0 <- proc.time()
plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.p$Sitetag)), ncol=18))
colnames(plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varWD", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.p$Sitetag))){
  Sitetag <- levels(traits.p$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.p[which(traits.p$Sitetag==Sitetag),]
  plot.results[i,1] <- Sitetag
  plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="WD",data=dataz)
  plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.WD', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   plot.results[i, 11:16] <- nullbounds
  #   plot.results[i,17] <- test.sig(x=plot.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
plot.results$sig <- rho.sig(n = plot.results$n, rho = plot.results$Rho)


###_______________ Trees w/in Sites level analysis _______________________________________________
t0 <- proc.time()
tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Sitetag)), ncol=18))
colnames(tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varWD", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Sitetag))){
  Sitetag <- levels(traits.t$Sitetag)[i]
  print(Sitetag)
  dataz <- traits.t[which(traits.t$Sitetag==Sitetag),]
  tree.results[i,1] <- Sitetag
  tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="WD",data=dataz)
  tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.WD', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.results$sig <- rho.sig(n = tree.results$n, rho = tree.results$Rho)

###_______________ Trees w/in Plot level analysis _______________________________________________
t0 <- proc.time()
tree.in.plot.results <- data.frame(matrix(NA, nrow=length(unique(traits.t$Plottag)), ncol=18))
colnames(tree.in.plot.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varWD", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.t$Plottag))){
  Plottag <- levels(traits.t$Plottag)[i]
  print(Plottag)
  dataz <- traits.t[which(traits.t$Plottag==Plottag),]
  tree.in.plot.results[i,1] <- Plottag
  tree.in.plot.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="WD",data=dataz)
  tree.in.plot.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.WD', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
tree.in.plot.results$sig <- rho.sig(n = tree.in.plot.results$n, rho = tree.in.plot.results$Rho)



###_______________ Branches w/in Plot level analysis _______________________________________________
t0 <- proc.time()
branch.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Plottag)), ncol=18))
colnames(branch.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varWD", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Plottag))){
  Plottag <- levels(traits.b$Plottag)[i]
  print(Plottag)
  dataz <- traits.b[which(traits.b$Plottag==Plottag),]
  branch.results[i,1] <- Plottag
  branch.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="WD",data=dataz)
  branch.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.WD', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.results$sig <- rho.sig(n = branch.results$n, rho = branch.results$Rho)


###_______________ branches w/in Tree level analysis _______________________________________________
t0 <- proc.time()
branch.in.tree.results <- data.frame(matrix(NA, nrow=length(unique(traits.b$Treetag)), ncol=18))
colnames(branch.in.tree.results) <- c("ID", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varlog.hub","varWD", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Species")
for(i in 1:length(unique(traits.b$Treetag))){
  Treetag <- levels(traits.b$Treetag)[i]
  print(Treetag)
  dataz <- traits.b[which(traits.b$Treetag==Treetag),]
  branch.in.tree.results[i,1] <- Treetag
  branch.in.tree.results[i,18] <- unique(dataz$Species)
  res <- fit.MAR(xvar='log.hub',yvar="WD",data=dataz)
  branch.in.tree.results[i,2:10] <- res 
  # if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
  #   nullbounds <- fit.null(xvar='log.WD', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
  #   tree.results[i, 11:16] <- nullbounds
  #   tree.results[i,17] <- test.sig(x=tree.results$Rho[i], test=nullbounds)
  # }
  
}
proc.time()-t0
branch.in.tree.results$sig <- rho.sig(n = branch.in.tree.results$n, rho = branch.in.tree.results$Rho)


spp.results$Level <- "Btw-site"
plot.results$Level <- "Btw-plot.in.site"
tree.results$Level <- "Btw-tree.in.site"
branch.results$Level<- "Btw-branch.in.plot"
branch.in.tree.results$Level <- "Btw-branch.in.tree"
tree.in.plot.results$Level <- "Btw-tree.in.plot"


log.hub.WD.spp.results <- spp.results
log.hub.WD.plot.results <- plot.results
log.hub.WD.tree.results <- tree.results
log.hub.WD.branch.results<- branch.results
log.hub.WD.branch.in.tree.results <- branch.in.tree.results
log.hub.WD.tree.in.plot.results <- tree.in.plot.results

log.hub.WD.results <- rbind(log.hub.WD.spp.results, log.hub.WD.plot.results, log.hub.WD.tree.results, log.hub.WD.branch.results)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## . TABLE S2: Summary of correlations ###########
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

WD.LMA <- WD.LMA.results %>% group_by(Level) %>% summarise(mean = round(mean(Rho, na.rm=T),2), median = round(median(Rho, na.rm=T),2),lq = round(quantile(Rho,.25, na.rm=T ),2), uq = round(quantile(Rho,.75 , na.rm=T),2))
WD.LDMC <- WD.LDMC.results %>% group_by(Level) %>% summarise(mean = round(mean(Rho, na.rm=T),2), median = round(median(Rho, na.rm=T),2),lq = round(quantile(Rho,.25, na.rm=T ),2), uq = round(quantile(Rho,.75 , na.rm=T),2))
LDMC.LMA <- LDMC.LMA.results %>% group_by(Level) %>% summarise(mean = round(mean(Rho, na.rm=T),2), median = round(median(Rho, na.rm=T),2),lq = round(quantile(Rho,.25, na.rm=T ),2), uq = round(quantile(Rho,.75 , na.rm=T),2))
#log.Al_As.LDMC.results %>% group_by(Level) %>% summarise(mean = round(mean(Rho, na.rm=T),2), median = round(median(Rho, na.rm=T),2),lq = round(quantile(Rho,.25, na.rm=T ),2), uq = round(quantile(Rho,.75 , na.rm=T),2))
#log.Al_As.LMA.results %>% group_by(Level) %>% summarise(mean = round(mean(Rho, na.rm=T),2), median = round(median(Rho, na.rm=T),2),lq = round(quantile(Rho,.25, na.rm=T ),2), uq = round(quantile(Rho,.75 , na.rm=T),2))
#log.Al_As.WD.results %>% group_by(Level) %>% summarise(mean = round(mean(Rho, na.rm=T),2), median = round(median(Rho, na.rm=T),2),lq = round(quantile(Rho,.25, na.rm=T ),2), uq = round(quantile(Rho,.75 , na.rm=T),2))
log.hub.LDMC <- log.hub.LDMC.results %>% group_by(Level) %>% summarise(mean = round(mean(Rho, na.rm=T),2), median = round(median(Rho, na.rm=T),2),lq = round(quantile(Rho,.25, na.rm=T ),2), uq = round(quantile(Rho,.75 , na.rm=T),2))
log.hub.LMA <- log.hub.LMA.results %>% group_by(Level) %>% summarise(mean = round(mean(Rho, na.rm=T),2), median = round(median(Rho, na.rm=T),2),lq = round(quantile(Rho,.25, na.rm=T ),2), uq = round(quantile(Rho,.75 , na.rm=T),2))
log.hub.WD <- log.hub.WD.results %>% group_by(Level) %>% summarise(mean = round(mean(Rho, na.rm=T),2), median = round(median(Rho, na.rm=T),2),lq = round(quantile(Rho,.25, na.rm=T ),2), uq = round(quantile(Rho,.75 , na.rm=T),2))

covarsums <- data.frame(trat_pair=c("LDMC-v-LMA","log.HV-vs-LDMC","log.HV-vs-LMA","log.HV-vs-WD","WD-vs-LDMC","WD-vs-LMA"), "Btw_site_means"=rep(NA, 6), "Btw_plot_means"=rep(NA, 6), "Btw_tree_means_in_site"=rep(NA, 6),"Btw_branch_means_in_plot"=rep(NA, 6))

covarsums[1,2:5] <- c(paste0(LDMC.LMA$mean[3]," (",LDMC.LMA$lq[3]," - ",LDMC.LMA$uq[3],")")
                      , paste0(LDMC.LMA$mean[2]," (",LDMC.LMA$lq[2]," - ",LDMC.LMA$uq[2],")")
                      , paste0(LDMC.LMA$mean[4]," (",LDMC.LMA$lq[4]," - ",LDMC.LMA$uq[4],")")
                      , paste0(LDMC.LMA$mean[1]," (",LDMC.LMA$lq[1]," - ",LDMC.LMA$uq[1],")"))

covarsums[2,2:5] <- c(paste0(log.hub.LDMC$mean[3]," (",log.hub.LDMC$lq[3]," - ",log.hub.LDMC$uq[3],")")
                      , paste0(log.hub.LDMC$mean[2]," (",log.hub.LDMC$lq[2]," - ",log.hub.LDMC$uq[2],")")
                      , paste0(log.hub.LDMC$mean[4]," (",log.hub.LDMC$lq[4]," - ",log.hub.LDMC$uq[4],")")
                      , paste0(log.hub.LDMC$mean[1]," (",log.hub.LDMC$lq[1]," - ",log.hub.LDMC$uq[1],")"))

covarsums[3,2:5] <- c(paste0(log.hub.LMA$mean[3]," (",log.hub.LMA$lq[3]," - ",log.hub.LMA$uq[3],")")
                      , paste0(log.hub.LMA$mean[2]," (",log.hub.LMA$lq[2]," - ",log.hub.LMA$uq[2],")")
                      , paste0(log.hub.LMA$mean[4]," (",log.hub.LMA$lq[4]," - ",log.hub.LMA$uq[4],")")
                      , paste0(log.hub.LMA$mean[1]," (",log.hub.LMA$lq[1]," - ",log.hub.LMA$uq[1],")"))

covarsums[4,2:5] <- c(paste0(log.hub.WD$mean[3]," (",log.hub.WD$lq[3]," - ",log.hub.WD$uq[3],")")
                      , paste0(log.hub.WD$mean[2]," (",log.hub.WD$lq[2]," - ",log.hub.WD$uq[2],")")
                      , paste0(log.hub.WD$mean[4]," (",log.hub.WD$lq[4]," - ",log.hub.WD$uq[4],")")
                      , paste0(log.hub.WD$mean[1]," (",log.hub.WD$lq[1]," - ",log.hub.WD$uq[1],")"))

covarsums[5,2:5] <- c(paste0(WD.LDMC$mean[3]," (",WD.LDMC$lq[3]," - ",WD.LDMC$uq[3],")")
                      , paste0(WD.LDMC$mean[2]," (",WD.LDMC$lq[2]," - ",WD.LDMC$uq[2],")")
                      , paste0(WD.LDMC$mean[4]," (",WD.LDMC$lq[4]," - ",WD.LDMC$uq[4],")")
                      , paste0(WD.LDMC$mean[1]," (",WD.LDMC$lq[1]," - ",WD.LDMC$uq[1],")"))

covarsums[6,2:5] <- c(paste0(WD.LMA$mean[3]," (",WD.LMA$lq[3]," - ",WD.LMA$uq[3],")")
                      , paste0(WD.LMA$mean[2]," (",WD.LMA$lq[2]," - ",WD.LMA$uq[2],")")
                      , paste0(WD.LMA$mean[4]," (",WD.LMA$lq[4]," - ",WD.LMA$uq[4],")")
                      , paste0(WD.LMA$mean[1]," (",WD.LMA$lq[1]," - ",WD.LMA$uq[1],")"))

write.csv(covarsums, file=paste0("./",results_dirname,"/TableS4_CorrelationScales.csv"))





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######### . FIG 2new: Trait covariation 6 panels (WD, LMA, LDMC, log.hub) ###################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# where to put the panel labels
letlin <- 0
letadj <- 0

palette(pallight)
# quartz(width=5, height=4.5)
# par(mfrow=c(2,3), mar=c(3.5,3.5,1,1), mgp=c(2,1,0), cex.lab=1.5, oma=c(0,0,2,0))


#quartz(width=7.4, height=4.1)
jpeg(file=paste0("./",results_dirname,"/Fig2_TraitCorrelations.jpg"), width=7.4, height=4.1, res=600, units="in")
mymat <- matrix(c(1,2,3,7
                  ,4,5,6,7), nrow=2, byrow=T)
layout(mat = mymat)#,widths = c(.9,.9,.9,1.3))
par( mar=c(3.6,3.7,1,.5), mgp=c(2,.7,0), cex.lab=1.4, oma=c(0,.2,2,0))

plot(WD~LDMC, traits.t, pch=pchs, cex=.9, col=Species, ylab = expression(paste("WD ",(g/cm^3))), xlab="LDMC (g/g)")
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"LDMC"], traits.t[which(traits.t$Species==i),"WD"])
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(xvar = "LDMC", yvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(xvar = "LDMC", yvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("a)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 3, line= -1.3, adj=.1)
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=3, line=-2,adj=.05, cex=.7)
# legend(x = 0.1, y=.97, legend =c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl")
#        , text.font=3, pch=c(3, rep(16, times=7)), col=pal[c(1,5,4,3,7,8,2,6)], xpd=NA, ncol=8, bty="n")


plot(WD~LMA, traits.t, pch=pchs, cex=.9, col=Species,ylab=expression(paste("WD ",(g/cm^3))), xlab=expression(paste("LMA ",(g/cm^2))))
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"WD"], traits.t[which(traits.t$Species==i),"LMA"])
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(xvar = "LMA", yvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(xvar = "LMA", yvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("b)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 3, line= -1.3, adj=.1)
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=3, line=-2,adj=.05, cex=.7)


plot(LMA~LDMC, traits.t, pch=pchs, cex=.9, col=Species, ylab=expression(paste("LMA ",(g/cm^2))), xlab="LDMC (g/g)")
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"LDMC"], traits.t[which(traits.t$Species==i),"LMA"])
  #print(paste(i, correlation$p.value))
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(xvar = "LDMC", yvar = "LMA",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(xvar = "LDMC", yvar = "LMA",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("c)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 3, line= -1.3, adj=.1)
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=3, line=-2,adj=.05, cex=.7)




plot(log.hub~LDMC, traits.t, pch=pchs, cex=.9, col=Species, ylab=expression(paste(log[10](HV))), xlab="LDMC (g/g)")
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"LDMC"], traits.t[which(traits.t$Species==i),"log.hub"])
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(yvar = "log.hub", xvar = "LDMC",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(yvar = "log.hub", xvar = "LDMC",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("d)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 1, line= -2, adj=.9)# -1.3
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=1, line=-1.1,adj=.95, cex=.7)



plot(log.hub~LMA, traits.t, pch=pchs, cex=.9, col=Species, ylab=expression(paste(log[10](HV))), xlab=expression(paste("LMA ", (g/cm^2))))
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"LMA"], traits.t[which(traits.t$Species==i),"log.hub"])
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(yvar = "log.hub", xvar = "LMA",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(yvar = "log.hub", xvar = "LMA",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("e)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 1, line= -2, adj=.9)
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=1, line=-1.1,adj=.95, cex=.7)



plot(log.hub~WD, traits.t, pch=pchs, cex=.9, col=Species, ylab=expression(paste(log[10](HV))), xlab=expression(paste("WD ", (g/cm^3))))
n_sig <- 0
cors <- rep(NA, times=8)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  correlation <- cor.test(traits.t[which(traits.t$Species==i),"log.hub"], traits.t[which(traits.t$Species==i),"WD"])
  print(paste(i, correlation$p.value))
  cors[j]<- correlation$estimate
  if(correlation$p.value <= 0.05){
    plot.MAR(yvar = "log.hub", xvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2)
    n_sig <- n_sig + 1
  }
  else{
    plot.MAR(yvar = "log.hub", xvar = "WD",data= traits.t[which(traits.t$Species==i),], linecol = pal[j], lwd=2, lty=2)
  }
}
mtext("f)",side=3, line=letlin, adj=letadj)
mtext(text = paste(n_sig,"8", sep="/"), side = 1, line= -2, adj=.9)
mtext(text = paste0("(",round(min(cors),2),"-",round(max(cors),2),")"), side=1, line=-1.1,adj=.95, cex=.7)




## Final panel density plot

par(mar=c(5,3.2,8,1),mgp=c(2,1,0))
plot(density(WD.LMA.spp.results$Rho, na.rm=T, adjust=2),lwd=3, xlim=c(-1,1), main="",xlab='Correlation')
lines(density(WD.LMA.branch.results$Rho, na.rm=T, adjust=2),lwd=1, lty=5)
#lines(density(LMA.WD.branch.in.tree.results$Rho, na.rm=T, adjust=2),lwd=1, lty=3)
lines(density(WD.LMA.tree.results$Rho, na.rm=T, adjust=2),lwd=2, lty=3)
#lines(density(LMA.WD.tree.in.plot.results$Rho, na.rm=T, adjust=2),lwd=1.5, lty=3)
#lines(density(trait.cors$WDvLMA, na.rm=T, adjust=2),lwd=3, lty=3)

lines(density(LDMC.LMA.spp.results$Rho, na.rm=T, adjust=2),lwd=3,col="blue")
lines(density(LDMC.LMA.branch.results$Rho, na.rm=T, adjust=2),lwd=1, lty=5,col="blue")
#lines(density(LDMC.LMA.branch.in.tree.results$Rho, na.rm=T, adjust=2),lwd=1, lty=3,col="blue")
lines(density(LDMC.LMA.tree.results$Rho, na.rm=T, adjust=2),lwd=2, lty=3,col="blue")
#lines(density(LDMC.LMA.tree.in.plot.results$Rho, na.rm=T, adjust=2),lwd=1.5, lty=3,col="blue")
#lines(density(trait.cors$LMAvLDMC, na.rm=T, adjust=2),lwd=3, lty=3,col="blue")
legend('topleft', legend = c("LMA-WD", "LMA-LDMC"), lwd=3, col=c("black","blue"), lty=1, bty="n")
legend('top',inset = -.25,xpd=NA, horiz = F, legend=c("branch","individ","site"), lty=c(5,3,1), lwd=c(1,2,3), bty="n")
mtext("g)",side=3, line=letlin, adj=letadj)
legend(x=-0.25, y=3.5,xjust=.5, legend =c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl")
       , text.font=3, pch=c(3, rep(16, times=7)), col=pal[c(1,5,4,3,7,8,2,6)], xpd=NA, ncol=3)
dev.off()
#quartz.save(file=paste0("./",results_dirname,"/Fig2_TraitCorrelations.pdf"),type = "pdf" )



# END FIG 2new +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### save trait-trait correlations for each species:
trait.cors1 <- data.frame(Species=levels(traits.t$Species), WDvLMA=NA, WDvLDMC=NA, WDvlog.hub=NA,LMAvLDMC=NA, LMAvlog.hub=NA, LDMCvlog.hub=NA)
for (j in 1:length(levels(traits.t$Species))){
  i <- levels(traits.t$Species)[j] 
  tmp <- traits.t[which(traits.t$Species==i),c("WD","LMA","LDMC","log.hub")]
  correlation <- cor(x=tmp, use = 'pairwise')[lower.tri(cor(x=tmp[,c("WD","LMA","LDMC","log.hub")], use = 'pairwise'))]
  trait.cors1[j,-1]<- round(correlation,2)
}






################# END: TRAIT COORDINATION ####################################







##### ++++++++++++++++++++++++++++++++++++++ #############
################### BEGIN: TRAIT PCA #####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++




#####  . branch level PCA on all samples ########

trait.pca <- prcomp(na.omit(traits.b[,c("WD","LMA","LDMC","log.hub")]), scale=T)

# Summary:
#                         PC1    PC2    PC3     PC4
# Proportion of Variance 0.5279 0.2630 0.1363 0.07291
# # Cumulative Proportion  0.5279 0.7908 0.9271 1.00000
#               PC1        PC2        PC3         PC4
# WD      -0.5335178  0.2932623 -0.7022725 -0.36901119
# LMA     -0.6078261 -0.2169686 -0.0298988  0.76326809
# LDMC    -0.4770488 -0.5779293  0.3988094 -0.52855772
# log.hub -0.3439922  0.7300134  0.5889571 -0.04335093


write.csv(round(trait.pca$rotation[,1:2],2),file = paste0(results_dirname,"/Fig3a_PCAloadings.csv"))

### multiply site mean trait values by the PCA rotation to get PC1 and PC2 scores
SitePCs <-  as.matrix(scale(traits.s[,c("WD","LMA","LDMC","log.hub")])) %*% as.matrix(trait.pca$rotation)

traits.s <- data.frame(traits.s, SitePCs)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
####################### . FIG 3: PC1 vs MDc ##############################
quartz(width=4, height=4)
#jpeg(file=paste0("./",results_dirname,"/Fig3_TraitCorrelations.jpg"), width=7.4, height=4.1, res=600, units="in")
par(mar=c(3.5,3.5,1,1), mgp=c(2.5,1,0))
palette(pal)
plot(PC1~MDc, traits.s, pch=ifelse(Species=="ACAC",3,16), col=Species, ylab="site PC1 score", xlab="Moisture Deficit (PET-PPT, mm)")
for (j in 1:length(levels(traits.s$Species))){
  i <- levels(traits.s$Species)[j] 
  (pvalue <- cor.test(x=traits.s[which(traits.s$Species==i),"MDc"], y=traits.s[which(traits.s$Species==i),"PC1"])$p.value)
  plot.MAR(xvar = "MDc", yvar = "PC1",data= traits.s[which(traits.s$Species==i),], linecol = pallight[j], lwd=3)#, lty= ifelse(pvalue<=0.05, 1,2))
  
}
mtext("b)", side=3, line=.1, adj=.05)

legend("bottomleft", legend =c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl")
       , text.font=3, pch=c(3,rep(16, times=7)), lty=1,lwd=1.5, col=pal[c(1,5,4,3,7,8,2,6)], xpd=NA, ncol=2, bty="n", cex=.9)

quartz.save(file=paste0("./",results_dirname,"/Fig3b_PCAClimate_MDchelsa.pdf"),type = "pdf" )






#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##################### . FIG S5: PC2 vs MDc #######################

quartz(width=4, height=4)
par(mar=c(3.5,3.5,1,1), mgp=c(2.5,1,0))
palette(pal)
plot(PC2~MDc, traits.s, pch=ifelse(Species=="ACAC",3,16), col=Species, ylab="site PC2 score", xlab="Moisture Deficit (PET-PPT, mm)")
for (j in 1:length(levels(traits.s$Species))){
  i <- levels(traits.s$Species)[j] 
  (pvalue <- cor.test(x=traits.s[which(traits.s$Species==i),"MDc"], y=traits.s[which(traits.s$Species==i),"PC2"])$p.value)
  plot.MAR(xvar = "MDc", yvar = "PC2",data= traits.s[which(traits.s$Species==i),], linecol = pal[j], lwd=3, lty= ifelse(pvalue<=0.05, 1,2))
}
#mtext("b)", side=3, line=.1, adj=.05)

legend("topleft", legend =c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl")
       , text.font=3, pch=c(3,rep(16, times=7)), lty=1,lwd=1.5, col=pal[c(1,5,4,3,7,8,2,6)], xpd=NA, ncol=2, bty="n", cex=.9)

quartz.save(file=paste0("./",results_dirname,"/FigS6_PC2Climate_MDchelsa.pdf"),type = "pdf" )


##### Test significance to relationship with PCs at site level
testPPT <- lmer(PC1~scale(PPT) + (0 + scale(PPT)|Species), traits.s)
testMD <- lmer(PC1~scale(MD) + (0 + scale(MD)|Species), traits.s)
testPET <- lmer(PC1~scale(PET) + (0 + scale(PET)|Species), traits.s)

testPPTc <- lmer(PC1~scale(PPTc) + (0 + scale(PPTc)|Species), traits.s)
testMDc <- lmer(PC1~scale(MDc) + (0 + scale(MDc)|Species), traits.s)
testPETc <- lmer(PC1~scale(PETc) + (0 + scale(PETc)|Species), traits.s)















##### ++++++++++++++++++++++++++++++++++++++ #############
################### BEGIN: VARIANCE DECOMPOSITION #####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++




#  Approach: use mixed effects models for each species and trait with only random effects (or only branch diameter as fixed effect for WD). Then put the variance components all in one dataframe for additional analysis and visualization

#Note 1: theoretically it shouldn't matter whether these are nested random effects or not, as long as all the lower-level effects are globally unique IDs (according to Zuur et al. 2009).
#In practice, this seems to be qualitatively true but not precisely quantitatively true. It doesn't change any conclusions to define them as nested, but convergence becomes a MAJOR issue for some traits/species. So here I have opted to stick with independant random effects for each organizational level for computation's sake.

#Note 2: after an update of the {lme4} function, some of these models started throwing convergence errors. These models still appear to be fitting appropriately, and have not been tweaked to remove the error. Also 'boundary (singular) fit' warnings indicate that a variance parameter was estimated to be 0.



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########### . Variance Decomposition for all species, all traits #######
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++ first for WD ++++++++++
testACACvd <- lmer(WD~ Diameter + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC")
testESALvd <- lmer(WD~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL")
testEMARGvd <- lmer(WD~ Diameter + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG")
testCOCAvd <- lmer(WD~1 + Diameter + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="COCA")

testAMYGvd <- lmer(WD ~ Diameter + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="AMYG")
# testAMYGvd <- lmer(WD ~ Diameter + (1|Plottag/Sitetag/Treetag), traits.b, subset=Species=="AMYG")
testOBLIvd <- lmer(WD ~ Diameter + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OBLI")
testOVATvd <- lmer(WD ~ Diameter + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OVAT")
testVIMIvd <- lmer(WD ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="VIMI")


## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
Acac <- data.frame(VarCorr(testACACvd))
Esal <- data.frame(VarCorr(testESALvd))
Emarg <- data.frame(VarCorr(testEMARGvd))
Coca <- data.frame(VarCorr(testCOCAvd))
Ovat <- data.frame(VarCorr(testOVATvd))
Vimi <- data.frame(VarCorr(testVIMIvd))
Amyg <- data.frame(VarCorr(testAMYGvd))
Obli <- data.frame(VarCorr(testOBLIvd))

# make a data frame, rows are initially: between tree (within plot), between plot (within site), between site, and within tree
variancesWD <- data.frame(Acac[,4],Esal[,4], Emarg[,4],Coca[,4], Ovat[,4], Vimi[,4], Amyg[,4], Obli[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","wiTree"))
colnames(variancesWD) <- c("ACAC", "ESAL", "EMARG", "COCA", 'OVAT', "VIMI", "AMYG", "OBLI")
# and reorder variance componentsto make them nice plots
variancesWD <- variancesWD[c(3,2,1,4),]
climsens <- c("PET", "MD", "ppt", "pet","PPT+PET", "PPT", "PPT", "pet")







#++++ Next for LMA ++++++++++
# note: raw LMAs and logged SLAs look pretty similar in a model with all spp. So I think I'll go raw LMA for the time being...
# The general inferences from raw LMA and logged SLA are pretty similar,
# with the rank orders being unchanged except for AMYG and OBLI.
testACACvd <- lmer(LMA~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC")
testESALvd <- lmer(LMA~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL")
testEMARGvd <- lmer(LMA~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG")
testCOCAvd <- lmer(LMA~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="COCA")

testAMYGvd <- lmer(LMA ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="AMYG")
testOBLIvd <- lmer(LMA ~ 1 +  (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OBLI")
testOVATvd <- lmer(LMA ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OVAT")
testVIMIvd <- lmer(LMA ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="VIMI")


## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
Acac <- data.frame(VarCorr(testACACvd))
Esal <- data.frame(VarCorr(testESALvd))
Emarg <- data.frame(VarCorr(testEMARGvd))
Coca <- data.frame(VarCorr(testCOCAvd))
Ovat <- data.frame(VarCorr(testOVATvd))
Vimi <- data.frame(VarCorr(testVIMIvd))
Amyg <- data.frame(VarCorr(testAMYGvd))
Obli <- data.frame(VarCorr(testOBLIvd))

# make a data frame, rows are initially: between tree (within plot), between plot (within site), between site, and within tree
variancesLMA <- data.frame(Acac[,4],Esal[,4], Emarg[,4],Coca[,4], Ovat[,4], Vimi[,4], Amyg[,4], Obli[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","wiTree"))
colnames(variancesLMA) <- c("ACAC", "ESAL", "EMARG", "COCA", 'OVAT', "VIMI", "AMYG", "OBLI")
# and reorder variance componentsto make them nice plots
variancesLMA <- variancesLMA[c(3,2,1,4),]







#++++ Then LDMC ++++++++++
testACACvd <- lmer(LDMC~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC")
testESALvd <- lmer(LDMC~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL")
testEMARGvd <- lmer(LDMC~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG")
testCOCAvd <- lmer(LDMC~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="COCA")

testAMYGvd <- lmer(LDMC ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="AMYG")
testOBLIvd <- lmer(LDMC ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OBLI")
testOVATvd <- lmer(LDMC ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OVAT")
testVIMIvd <- lmer(LDMC ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="VIMI")


## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
Acac <- data.frame(VarCorr(testACACvd))
Esal <- data.frame(VarCorr(testESALvd))
Emarg <- data.frame(VarCorr(testEMARGvd))
Coca <- data.frame(VarCorr(testCOCAvd))
Ovat <- data.frame(VarCorr(testOVATvd))
Vimi <- data.frame(VarCorr(testVIMIvd))
Amyg <- data.frame(VarCorr(testAMYGvd))
Obli <- data.frame(VarCorr(testOBLIvd))

# make a data frame, rows are initially: between tree (within plot), between plot (within site), between site, and within tree
variancesLDMC <- data.frame(Acac[,4],Esal[,4], Emarg[,4],Coca[,4], Ovat[,4], Vimi[,4], Amyg[,4], Obli[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","wiTree"))
colnames(variancesLDMC) <- c("ACAC", "ESAL", "EMARG", "COCA", 'OVAT', "VIMI", "AMYG", "OBLI")
# and reorder variance componentsto make them nice plots
variancesLDMC <- variancesLDMC[c(3,2,1,4),]




#++++ Next for Al_As ++++++++++
testACACvd <- lmer(log.Al_As~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC")
testESALvd <- lmer(log.Al_As~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL")
testEMARGvd <- lmer(log.Al_As~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG")
testCOCAvd <- lmer(log.Al_As~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="COCA")

testAMYGvd <- lmer(log.Al_As ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="AMYG")
testOBLIvd <- lmer(log.Al_As ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OBLI")
testOVATvd <- lmer(log.Al_As ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OVAT")
testVIMIvd <- lmer(log.Al_As ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="VIMI")


## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
Acac <- data.frame(VarCorr(testACACvd))
Esal <- data.frame(VarCorr(testESALvd))
Emarg <- data.frame(VarCorr(testEMARGvd))
Coca <- data.frame(VarCorr(testCOCAvd))
Ovat <- data.frame(VarCorr(testOVATvd))
Vimi <- data.frame(VarCorr(testVIMIvd))
Amyg <- data.frame(VarCorr(testAMYGvd))
Obli <- data.frame(VarCorr(testOBLIvd))

# make a data frame, rows are initially: between tree (within plot), between plot (within site), between site, and within tree
variancesAl_As <- data.frame(Acac[,4],Esal[,4], Emarg[,4],Coca[,4], Ovat[,4], Vimi[,4], Amyg[,4], Obli[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","wiTree"))
colnames(variancesAl_As) <- c("ACAC", "ESAL", "EMARG", "COCA", 'OVAT', "VIMI", "AMYG", "OBLI")
# and reorder variance componentsto make them nice plots
variancesAl_As <- variancesAl_As[c(3,2,1,4),]





#++++ Next for hub ++++++++++
# note: actually mathmatically identical to variance decomp of log.Al_As because log.Al_As = -1* log.hub
testACACvd <- lmer(log.hub~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ACAC")
testESALvd <- lmer(log.hub~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="ESAL")
testEMARGvd <- lmer(log.hub~ 1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="EMARG")
testCOCAvd <- lmer(log.hub~1 + (1|Sitetag) + (1|Plottag) + (1|Treetag), data = traits.b, subset=Species=="COCA")

testAMYGvd <- lmer(log.hub ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="AMYG")
testOBLIvd <- lmer(log.hub ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OBLI")
testOVATvd <- lmer(log.hub ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="OVAT")
testVIMIvd <- lmer(log.hub ~ 1 + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b, subset=Species=="VIMI")


## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
Acac <- data.frame(VarCorr(testACACvd))
Esal <- data.frame(VarCorr(testESALvd))
Emarg <- data.frame(VarCorr(testEMARGvd))
Coca <- data.frame(VarCorr(testCOCAvd))
Ovat <- data.frame(VarCorr(testOVATvd))
Vimi <- data.frame(VarCorr(testVIMIvd))
Amyg <- data.frame(VarCorr(testAMYGvd))
Obli <- data.frame(VarCorr(testOBLIvd))

# make a data frame, rows are initially: between tree (within plot), between plot (within site), between site, and within tree
varianceshub <- data.frame(Acac[,4],Esal[,4], Emarg[,4],Coca[,4], Ovat[,4], Vimi[,4], Amyg[,4], Obli[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","wiTree"))
colnames(varianceshub) <- c("ACAC", "ESAL", "EMARG", "COCA", 'OVAT', "VIMI", "AMYG", "OBLI")
# and reorder variance componentsto make them nice plots
varianceshub <- varianceshub[c(3,2,1,4),]




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######### . Full variance decomp of all species ##############

fullWD <- lmer(WD ~ Diameter + (1|Species) + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b[-which(traits.b$Species=="ACAC"),])
fullLDMC <- lmer(LDMC ~ 1 + (1|Species) + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b[-which(traits.b$Species=="ACAC"),])
fullLMA <- lmer(LMA ~ 1 + (1|Species) + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b[-which(traits.b$Species=="ACAC"),])
fullHV <- lmer(log.hub ~ 1 + (1|Species) + (1|Plottag) + (1|Sitetag) + (1|Treetag), traits.b[-which(traits.b$Species=="ACAC"),])
WDvar <- data.frame(VarCorr(fullWD))
LDMCvar <- data.frame(VarCorr(fullLDMC))
LMAvar <- data.frame(VarCorr(fullLMA))
HVvar <- data.frame(VarCorr(fullHV))

variancesall <- data.frame(WDvar[,4],LMAvar[,4], LDMCvar[,4],HVvar[,4], row.names = c("wiPlot", "wiSite","bSite/Clim","bSpecies","wiTree"))
colnames(variancesall) <- c("WD","LMA","LDMC","HV")
# and reorder variance componentsto make them nice plots
variancesall <- variancesall[c(4,3,2,1,5),]
variancesall.scaled <- apply(variancesall, MARGIN=2, FUN = function(x){x/sum(x)})

##### . scaling all variances #############
scaledvariancesWD <- apply(variancesWD, MARGIN=2, FUN= function(x){x/sum(x)})
scaledvariancesLMA <- apply(variancesLMA, MARGIN=2, FUN= function(x){x/sum(x)})
scaledvariancesLDMC <- apply(variancesLDMC, MARGIN=2, FUN= function(x){x/sum(x)})
scaledvarianceshub <- apply(varianceshub, MARGIN=2, FUN= function(x){x/sum(x)})
totalvariancesWD <- colSums(variancesWD)
totalvariancesLMA <- colSums(variancesLMA)
totalvariancesLDMC <- colSums(variancesLDMC)
totalvarianceshub <- colSums(varianceshub)


### .  combining variance estimates into one dataframe:  #################
# # first make them a wide-form dataframe

varswide <- data.frame(rbind(rbind(variancesWD, colSums(variancesWD)), rbind(variancesLMA, colSums(variancesLMA)),
                             rbind(variancesLDMC, colSums(variancesLDMC)), rbind(variancesAl_As, colSums(variancesAl_As)), rbind(varianceshub, colSums(varianceshub))),"trait"=rep(c("WD","LMA","LDMC","Al_As","hub"), each=5))
varswide$level <- rep(c("bSite","wiSite","wiPlot","wiTree","Total"), times=5)

# then melt them into long form
varslong <-melt(data = varswide, id.vars = c("trait","level") , variable.name="Species")
colnames(varslong)[which(colnames(varslong)=="value")] <- "variance"
# now add the species aridity niche
quants$Species <- toupper(quants$spp) # make species codes upper case to match with varslong

# add in species climate niche info from distribution records
# driest range edge
varslong$CMD.90 <- quants$quantCMDcALA.90[match(varslong$Species,quants$Species)] 
varslong$PPT.10 <- quants$quantPPTcALA.10[match(varslong$Species,quants$Species)]
varslong$PET.90 <- quants$quantCMDcALA.90[match(varslong$Species,quants$Species)]
# median
varslong$CMD.50 <- quants$quantCMDcALA.50[match(varslong$Species,quants$Species)]
varslong$PET.50 <- quants$quantPETcALA.50[match(varslong$Species,quants$Species)]
varslong$PPT.50 <- quants$quantPPTcALA.50[match(varslong$Species,quants$Species)]
varslong$level <- factor(varslong$level)
varslong$Species.trait <- with(varslong, paste(Species, trait, sep="."))





## Some might argue that coefficients of variation are more appropriate measures than raw variances for analyzing the change in variance components with aridity. So this code calculates CVs for plotting SI figures.



trait.means <- traits.b %>% group_by(Species) %>% summarise(WD = mean(WD, na.rm=T), LDMC=mean(LDMC, na.rm=T), LMA=mean(LMA, na.rm=T), Al_As=log(mean(Al_As, na.rm=T), base=10), hub=-1*log(mean(hub, na.rm=T), base=10)) #%>% mutate(log.Al_As = log(Al_As,base=10)) %>% select(-Al_As)
# NOTE: Al_As is actually log10(Al_As) but had to keem name for matching with varslong
trait.means.long <- trait.means %>% gather(trait, mean, -Species)
trait.means.long$Species.trait <- paste(trait.means.long$Species, trait.means.long$trait, sep=".")





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## . FIG4: Barplot of variance Decomposition ###########
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## set color choices
#pal.vardecomp <- brewer.pal(n=9, "Set1")
pal.vardecomp <- c(paste0(brewer.pal(n=3, "Set1")[1],"66"), brewer.pal(n=9, "Blues")[c(8,5,2)])
palette(pal.vardecomp)
# set the colors that wil denote within-species, within-genus, within family and across CWMs
colchoices <- c(1,2,4,3,6)


#quartz(width=4.3,height=6.4) 
jpeg(file=paste0("./",results_dirname,"/Fig4_VarainceDecomp_scaled.jpg"), width=4.3, height=6.4, units="in", res=600)
par(mfrow=c(4,1), mar=c(1,3.6,1,3.6), mgp=c(2.3,1,0), oma=c(3.6,0,3,0), cex.lab=1.4, cex.axis=1.1)
p<-barplot(height=as.matrix(scaledvariancesWD)
           , beside=F, names.arg=rep(NA, times=8)
           , col = pal.vardecomp
           , legend.text= c("Btw Sites", "Within Site", "Within Plot", "Within Tree")
           , args.legend=list(bty="n", x=mean(c(1.1,16)), y=1.5, xpd=NA, cex=1.3,xjust=0.5, ncol=4)
           , ylab="% WD Var", las=3
           , xlim=c(1.1,16.2), width=1, space=1)#,  yaxt="n" # ylim=c(0,0.008),
axis(4,at=c(0,0.002,.004,0.006,.008)/max(totalvariancesWD), labels = c("0","0.002","","0.006",""), xpd=NA,col="#333333", col.axis="#333333")
mtext("Tot WD Var", side = 4, line =2.3)
barplot(height=as.matrix(t(totalvariancesWD/max(totalvariancesWD))), add=T
        , space=c(4,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=8))
#text(x=p,y=colSums(scaledvariancesWD)+.001, labels=climsens, srt=90, xpd=T, font=c(2,2,1,1,2,2,2,1))
mtext(text = "a)", side=3, adj=-0.12, line=.2)

p<-barplot(height=as.matrix(scaledvariancesLMA)
           , beside=F, names.arg=rep(NA, times=8)
           , col = pal.vardecomp
           #, legend.text= c("Btw Sites", "Within Site", "Within Plot", "Within Tree"), args.legend=list(bty="n")
           , ylab="% LMA Var", las=3
           , xlim=c(1.1,16.2), width=1, space=1)#, yaxt="n")# , ylim=c(0,2.5e-5)
barplot(height=as.matrix(t(totalvariancesLMA/max(totalvariancesLMA))), add=T
        , space=c(4,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=8))
mtext(text = "b)", side=3, adj=-0.12, line=.2)
mtext("Tot LMA Var", side = 4, line =2.3)
axis(4, at=c(0,.5e-5,1e-5,1.5e-5,2e-5)/max(totalvariancesLMA), labels=c("0","","1e-5","","2e-5"),col="#333333", col.axis="#333333")
p<-barplot(height=as.matrix(scaledvariancesLDMC)
           , beside=F, names.arg=rep(NA, times=8)
           , col = pal.vardecomp
           #, legend.text= c("Btw Sites", "Within Site", "Within Plot", "Within Tree"), args.legend=list(bty="n")
           , ylab=" % LDMC Var", las=3
           , xlim=c(1.1,16.2), width=1, space=1)#, yaxt="n")
barplot(height=as.matrix(t(totalvariancesLDMC/max(totalvariancesLDMC))), add=T
        , space=c(4,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=8))
mtext(text = "c)", side=3, adj=-0.12, line=.2)
mtext("Tot LDMC Var", side = 4, line =2.3)
axis(4,at=c(0,0.001,.002,0.003,.004)/max(totalvariancesLDMC), labels = c("0","","2e-3","","4e-3"),col="#333333", col.axis="#333333" )
# p<-barplot(height=as.matrix(scaledvariancesAl_As)
#            , beside=F, names.arg=rep(NA, times=8)
#            , col = pal.vardecomp
#            #, legend.text= c("Btw Sites", "Within Site", "Within Plot", "Within Tree"), args.legend=list(bty="n")
#            , ylab=expression(paste("Var in ",log[10](A[L]:A[S]))), las=3)
p <- barplot(height=as.matrix(scaledvarianceshub)
             , beside=F, names.arg=rep(NA, times=8)
             , col=pal.vardecomp
             , ylab=expression(paste("% ",log[10](HV)," Var"))
             , xlim=c(1.1,16.2), width=1, space=1)
barplot(height=as.matrix(t(totalvarianceshub/max(totalvarianceshub))), add=T
        , space=c(4,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=8))
mtext(expression(paste("Tot ",log[10](HV)," Var")), 4, line=2.3)
# note: the var decomp of log.Al_As and log.hub are identical, because they are perfectly negatively correlated
mtext(text = "d)", side=3, adj=-0.12, line=.2)
axis(4,at=seq(from=0,to=0.06, by=0.02)/max(totalvarianceshub), labels=c("0","2e-3",NA,"4e-3"),col="#333333", col.axis="#333333")

axis(1,at=p,labels= c("A. acu","E. sal","E. mar","C. cal","E. ova","E. vim","E. amy","E. obl"),font=3, cex.axis = 1.4, tick=F,las=2)
dev.off()
#quartz.save(file=paste0("./",results_dirname,"/Fig4_VarainceDecomp_scaled.pdf"),type = "pdf" )



############### FIGURE S6: Scaled barplots for full dataset #################3
quartz(width=5.3,height=3.5) 
par( mar=c(2,3.6,1,6), mgp=c(2.3,1,0), oma=c(3.6,0,0,0), cex.lab=1.4, cex.axis=1.1)
p<-barplot(height=as.matrix(variancesall.scaled)
           , beside=F, names.arg=c("WD","LMA","LDMC","log10(HV)")
           , col = c(pal[1],pal.vardecomp), 
           , legend.text= c("Btw Species","Btw Sites", "Within Site", "Within Plot", "Within Tree")
           , args.legend=list(bty="n", x=5.8, y=1, xpd=NA, cex=1,xjust=0.5, ncol=1)
           , ylab="% Total Variance", las=3
           )#,  yaxt="n" # ylim=c(0,0.008),
quartz.save(file=paste0("./",results_dirname,"/FigS7_VarainceDecomp_scaled_fulldataset.pdf"),type = "pdf" )





#_______________________________________________________________________________
########### . FIGURE 5: Variance Components f(clim) #####################################
#_______________________________________________________________________________

#quartz(width=7.5, height=4)
jpeg(file=paste0("./",results_dirname,"/Fig5_VaraincebyClimate_MD50.jpg"), height=4, width=7.5, units="in", res=600)
par(mfrow=c(2,4),mar=c(0,3,0,0), oma=c(4,2,2.5,1), mgp=c(3,1,0), cex=.9)
labline = -.9
labadj = 0.04
palette(pal)
## first make a row of the total variances
tmp <- varslong[which(varslong$trait=="WD" & varslong$level=="Total" & varslong$Species!= "ACAC"),]
plot(variance~CMD.50, tmp, pch=16
     , xaxt="n", xlim=c(-250,1180)
     , yaxt="n", ylim=c(min(tmp$variance) - 0.04*min(tmp$variance), max(tmp$variance) + .11*max(tmp$variance)))
axis(2,at=c(0.0015,0.002,0.0025,0.003), labels=c(0.0015,"",0.0025,""))
mtext("Total Variance", side = 2, line=2.5)
mtext("WD", side=3, line=1, font=2)
mtext("a)",side=3, line=labline, adj=labadj)


tmp <- varslong[which(varslong$trait=="LMA" & varslong$level=="Total" & varslong$Species!= "ACAC"),]

plot(variance~CMD.50, tmp, pch=16
     , xaxt="n", xlim=c(-250,1180)
     , yaxt="n", ylim=c(min(tmp$variance) - 0.04*min(tmp$variance), max(tmp$variance) + .11*max(tmp$variance)))
axis(2, at=c(0.9e-5,1e-5,1.1e-5,1.2e-5,1.3e-5), labels=c(NA,"1e-5",NA,"1.2e-5",NA))#,0.00001,"",0.000012))
mtext("LMA", side=3, line=1, font=2)
mtext("b)",side=3, line=labline, adj=labadj)

tmp <- varslong[which(varslong$trait=="LDMC" & varslong$level=="Total" & varslong$Species!= "ACAC"),]
plot(variance~CMD.50, tmp, pch=16
     , xaxt="n", xlim=c(-250,1180)
     , yaxt="n", ylim=c(min(tmp$variance) - 0.04*min(tmp$variance), max(tmp$variance) + .11*max(tmp$variance)))
# plot(variance~CMD.50, varslong[which(varslong$trait=="LDMC" & varslong$level=="Total" & varslong$Species!= "ACAC"),]
#      , pch=16, xaxt="n", xlim=c(-250,1180), ylim=c(0.001,0.0033))
mtext("LDMC", side=3, line=1, font=2)
#abline(lm(variance~CMD.50, varslong[which(varslong$trait=="LDMC" & varslong$level=="Total" & varslong$Species!= "ACAC"),]), lwd=2, lty=2, col="grey")
#mtext(text = paste("p=",round(summary(lm(variance~CMD.50, varslong[which(varslong$trait=="LDMC" & varslong$level=="Total" & varslong$Species!= "ACAC"),]))$coefficients[2,4], digits = 2))
#      , side=3, line = 0, adj = .9)
#legend("topright",legend=c("Total","btw Sites","w/in Site","w/in Plot","w/in Tree"), pch=16, col=c("black",pal[c(1,3,2,4)]),cex=.7)
mtext("c)",side=3, line=labline, adj=labadj)

# plot(variance~CMD.50, varslong[which(varslong$trait=="Al_As" & varslong$level=="Total" & varslong$Species!= "ACAC"),]
#      , pch=16, xaxt="n", xlim=c(-250,1180), yaxt="n")
# mtext(text = expression(bold(log[10](A[L]:A[S]))), side=3, line=1, font=2)
tmp <- varslong[which(varslong$trait=="hub" & varslong$level=="Total" & varslong$Species!= "ACAC"),]
plot(variance~CMD.50, tmp, pch=16
     , xaxt="n", xlim=c(-250,1180)
     , yaxt="n", ylim=c(min(tmp$variance) - 0.04*min(tmp$variance), max(tmp$variance) + .11*max(tmp$variance)))
mtext(text = expression(bold(log[10](HV))), side=3, line=1, font=2)

axis(2,at=c(0.03,0.04,0.05,0.06), labels=c(NA,0.04, NA,0.06))
legend("topright",legend=c("Total","btw Sites","w/in Site","w/in Plot","w/in Tree"), pch=16, col=c("black",pal[c(1,3,2,4)]),cex=.7)
mtext(text = paste("p=",round(summary(lm(variance~CMD.50, varslong[which(varslong$trait=="hub" & varslong$level=="Total" & varslong$Species!= "ACAC"),]))$coefficients[2,4], digits = 2))
      , side=3, line = 0, adj = .9)
abline(lm(variance~CMD.50, varslong[which(varslong$trait=="hub" & varslong$level=="Total" & varslong$Species!= "ACAC"),]), lwd=2, lty=2, col="grey")
mtext("d)",side=3, line=labline, adj=labadj)


### now a row of plots breaking out each variance component
tmp <- varslong[which(varslong$trait=="WD" & varslong$level!="Total" & varslong$Species!= "ACAC"),]
plot(variance*100~CMD.50, tmp, pch=16
     ,  col=factor(level), xlim=c(-250,1180), yaxt="n", ylim=c(min(tmp$variance) - 0.04*min(tmp$variance), max(tmp$variance) + .13*max(tmp$variance))*100)
axis(2, at=c(0.025, 0.05, 0.075, 0.1), labels=c(NA,0.05,NA,0.1))
mtext("Var Components\n(*100)", side = 2, line=2.5)
mtext("e)",side=3, line=labline, adj=labadj)

tmp <- varslong[which(varslong$trait=="LMA" & varslong$level!="Total" & varslong$Species!= "ACAC"),]
plot(variance*100~CMD.50, tmp, pch=16
     ,  col=factor(level), xlim=c(-250,1180), yaxt="n", ylim=c(min(tmp$variance) - 0.04*min(tmp$variance), max(tmp$variance) + .13*max(tmp$variance))*100)
mtext("f)",side=3, line=labline, adj=labadj)

tmp <- varslong[which(varslong$trait=="LDMC" & varslong$level!="Total" & varslong$Species!= "ACAC"),]
plot(variance*100~CMD.50, tmp, pch=16
     ,  col=factor(level), xlim=c(-250,1180), yaxt="n", ylim=c(min(tmp$variance) - 0.04*min(tmp$variance), max(tmp$variance) + .13*max(tmp$variance))*100)
mtext("g)",side=3, line=labline, adj=labadj)
abline(lm(variance*100~CMD.50, varslong[which(varslong$trait=="LDMC" & varslong$level=="bSite" & varslong$Species!= "ACAC"),])
       , col=pal[1], lwd=2, xlim=c(-250,1180), lty=2)
mtext("Species median MD (PET-PPT, mm)", side=1,line = 3, outer=T)

# plot(variance*100~CMD.50, varslong[which(varslong$trait=="Al_As" & varslong$level!="Total" & varslong$Species!= "ACAC"),]
#      , pch=16,  col=factor(level), xlim=c(-250,1180))
# abline(lm(variance*100~CMD.50, varslong[which(varslong$trait=="Al_As" & varslong$level=="bSite" & varslong$Species!= "ACAC"),]), col=pal[1], lwd=2)
tmp <- varslong[which(varslong$trait=="hub" & varslong$level!="Total" & varslong$Species!= "ACAC"),]
plot(variance*100~CMD.50, tmp, pch=16
     ,  col=factor(level), xlim=c(-250,1180), yaxt="n", ylim=c(min(tmp$variance) - 0.04*min(tmp$variance), max(tmp$variance) + .13*max(tmp$variance))*100)
abline(lm(variance*100~CMD.50, varslong[which(varslong$trait=="hub" & varslong$level=="bSite" & varslong$Species!= "ACAC"),]), col=pal[1], lwd=2)

mtext("h)",side=3, line=labline, adj=labadj)
#quartz.save(file=paste0("./",results_dirname,"/Fig5_VaraincebyClimate_MD50.pdf"),type = "pdf" )
dev.off()



####### . Statistical analysis of variance component changes with climate ########

climvartest <- function(trait, varvar = "variance"){
  results <- list()
  tmp <- varslong[which(varslong$trait==trait & varslong$level!= "Total" & varslong$Species!="ACAC"),]
  (results$CMD90 <- summary(lm(get(varvar)~scale(CMD.90)*level, tmp)))
  (results$PET90 <- summary(lm(get(varvar)~scale(PET.90)*level, tmp)))
  (results$PPT10 <- summary(lm(get(varvar)~scale(PPT.10)*level, tmp)))
  (results$CMD50 <- summary(lm(get(varvar)~scale(CMD.50)*level, tmp)))
  (results$PET50 <- summary(lm(get(varvar)~scale(PET.50)*level, tmp)))
  (results$PPT50 <- summary(lm(get(varvar)~scale(PPT.50)*level, tmp)))
  return(results)
}

testWD <- climvartest("WD")
## No significant relationships with anything  
testLMA <- climvartest("LMA")
## sig with PPT.50, marginally sig with PPT.10
testLDMC <- climvartest("LDMC")
## sig with CMD.90, PET.90 and PPT.10, marginal with CMD.50, PPT.50 (.11 w/ PET.50)
testhub <- climvartest("hub")
##  sig with PET.90,CMD.90, sig with CMD.50, PET.50, marg with PPT.50
#lm(get(varvar)~scale(CMD.90)*level, varslong[which(varslong$trait=="WD" & varslong$level!="Total" & varslong$Species!= "ACAC"),])
#summary(testWD) # no significant changes in get(varvar) components with CMD

climvartest.tot <- function(trait, varvar="variance"){
  results <- list()
  tmp <- varslong[which(varslong$trait==trait & varslong$level== "Total" & varslong$Species!="ACAC"),]
  (results$CMD90 <- summary(lm(get(varvar)~scale(CMD.90), tmp)))
  (results$PET90 <- summary(lm(get(varvar)~scale(PET.90), tmp)))
  (results$PPT10 <- summary(lm(get(varvar)~scale(PPT.10), tmp)))
  (results$CMD50 <- summary(lm(get(varvar)~scale(CMD.50), tmp)))
  (results$PET50 <- summary(lm(get(varvar)~scale(PET.50), tmp)))
  (results$PPT50 <- summary(lm(get(varvar)~scale(PPT.50), tmp)))
  return(results)
}


testWD <- climvartest.tot("WD")
## No significant relationships with anything  
testLMA <- climvartest.tot("LMA")
##  nothing
testLDMC <- climvartest.tot("LDMC")
## sig PPT.10
testhub <- climvartest.tot("hub")
##  marg with PET.50, <.15 with CMD.90, PET.90, CMD.50





####### . Calculating CVs #######################

trait.means <- traits.b %>% group_by(Species) %>% summarise(WD = mean(WD, na.rm=T), LDMC=mean(LDMC, na.rm=T), LMA=mean(LMA, na.rm=T), Al_As=log(mean(Al_As, na.rm=T), base=10), hub=-1*log(mean(hub, na.rm=T), base=10)) #%>% mutate(log.Al_As = log(Al_As,base=10)) %>% select(-Al_As)
# NOTE: Al_As is actually log10(Al_As) but had to keem name for matching with varslong
trait.means.long <- trait.means %>% gather(trait, mean, -Species)
trait.means.long$Species.trait <- paste(trait.means.long$Species, trait.means.long$trait, sep=".")

varslong$meanvals <- trait.means.long$mean[match(varslong$Species.trait, trait.means.long$Species.trait)]
varslong$CV <- with(varslong, variance^0.5 / meanvals)



testWD <- climvartest("WD", varvar="CV")
## No significant relationships with anything  
testLMA <- climvartest("LMA", varvar="CV")
## sig with PPT.50,PPT.10, marginally sig with PET.90, CMD.90
testLDMC <- climvartest("LDMC", varvar="CV")
## sig with PPT.50, CMD.50, PET.50 and PPT.10, PET.90, CMD 90
testhub <- climvartest("hub", varvar="CV")
##  sig with PPT.50, PET.50, CMD.50 and PPT.10, PET.90, CMD.90



testWD <- climvartest.tot("WD", varvar="CV")
## No significant relationships with anything  
testLMA <- climvartest.tot("LMA", varvar="CV")
##  sig PPT.50, ms MD.50, ms PPT.10, ms PET.90, MD.90
testLDMC <- climvartest.tot("LDMC", varvar="CV")
## sig PPT.10, ms PET.90, ms CMD.90
testhub <- climvartest.tot("hub", varvar="CV")
##  nothing...



#_______________________________________________________________________________
########### . FIGURE S7: CV Components f(clim) #####################################
#_______________________________________________________________________________

quartz(width=7.5, height=4)
pdf(file=paste0("./",results_dirname,"/FigS7_CV-vs-clim.pdf"), width=7.5, height=4)
par(mfrow=c(2,4),mar=c(0,3,0,0), oma=c(4,1,2.5,1), mgp=c(3,1,0), cex=.9)

palette(pal)
## first make a row of the total CVs
plot(CV~CMD.50, varslong[which(varslong$trait=="WD" & varslong$level=="Total" & varslong$Species!= "ACAC"),]
     , pch=16, xaxt="n", xlim=c(-250,1180))
mtext("Total CV", side = 2, line=2.5)
mtext("WD", side=3, line=1, font=2)

plot(CV~CMD.50, varslong[which(varslong$trait=="LMA" & varslong$level=="Total" & varslong$Species!= "ACAC"),]
     , pch=16, xaxt="n", xlim=c(-250,1180))
mtext("LMA", side=3, line=1, font=2)
abline(lm(CV~CMD.50, varslong[which(varslong$trait=="LMA" & varslong$level=="Total" & varslong$Species!= "ACAC"),]), lwd=2, lty=1)
mtext(text = paste("p=",round(summary(lm(CV~CMD.50, varslong[which(varslong$trait=="LMA" & varslong$level=="Total" & varslong$Species!= "ACAC"),]))$coefficients[2,4], digits = 2))
      , side=3, line = -1, adj = .8)


plot(CV~CMD.50, varslong[which(varslong$trait=="LDMC" & varslong$level=="Total" & varslong$Species!= "ACAC"),]
     , pch=16, xaxt="n", xlim=c(-250,1180))
mtext("LDMC", side=3, line=1, font=2)
abline(lm(CV~CMD.50, varslong[which(varslong$trait=="LDMC" & varslong$level=="Total" & varslong$Species!= "ACAC"),]), lwd=2, lty=2, col="grey")
mtext(text = paste("p=",round(summary(lm(CV~CMD.50, varslong[which(varslong$trait=="LDMC" & varslong$level=="Total" & varslong$Species!= "ACAC"),]))$coefficients[2,4], digits = 2))
      , side=3, line = -1, adj = .8)

# plot(CV~CMD.50, varslong[which(varslong$trait=="Al_As" & varslong$level=="Total" & varslong$Species!= "ACAC"),]
#      , pch=16, xaxt="n", xlim=c(-250,1180))
# mtext(text = expression(bold(log[10](A[L]:A[S]))), side=3, line=1, font=2)
plot(CV~CMD.50, varslong[which(varslong$trait=="hub" & varslong$level=="Total" & varslong$Species!= "ACAC"),]
     , pch=16, xaxt="n", xlim=c(-250,1180))
mtext(text = expression(bold(log[10](HV))), side=3, line=1, font=2)

legend("topright",legend=c("Total","btw Sites","w/in Site","w/in Plot","w/in Tree"), pch=16, col=c("black",pal[c(1,3,2,4)]),cex=.7)


### now a row of plots breaking out each CV component
plot(CV~CMD.50, varslong[which(varslong$trait=="WD" & varslong$level!="Total" & varslong$Species!= "ACAC"),], pch=16,  col=factor(level))
mtext("CV of components", side = 2, line=2.5)

plot(CV~CMD.50, varslong[which(varslong$trait=="LMA" & varslong$level!="Total" & varslong$Species!= "ACAC"),], pch=16,  col=factor(level))

plot(CV~CMD.50, varslong[which(varslong$trait=="LDMC" & varslong$level!="Total" & varslong$Species!= "ACAC"),], pch=16,  col=factor(level))
abline(lm(CV~CMD.50, varslong[which(varslong$trait=="LDMC" & varslong$level=="bSite" & varslong$Species!= "ACAC"),]), col=pal[1], lwd=2)
mtext("Species median MD (PET-PPT, mm)", side=1,line = 3, adj=1,)

# plot(CV~CMD.50, varslong[which(varslong$trait=="Al_As" & varslong$level!="Total" & varslong$Species!= "ACAC"),], pch=16,  col=factor(level))
# abline(lm(CV~CMD.50, varslong[which(varslong$trait=="Al_As" & varslong$level=="bSite" & varslong$Species!= "ACAC"),]), col=pal[1], lwd=2)
plot(CV~CMD.50, varslong[which(varslong$trait=="hub" & varslong$level!="Total" & varslong$Species!= "ACAC"),], pch=16,  col=factor(level))
abline(lm(CV~CMD.50, varslong[which(varslong$trait=="hub" & varslong$level=="bSite" & varslong$Species!= "ACAC"),]), col=pal[1], lwd=2)

dev.off()









###### . FIG 6 example: increasing variance with PPT ############


res.type="pearson"
palvars <- paste0(brewer.pal(n=9, "Set1"), "AA")
colchoices <- c(1,2,4,3,6)

# climvar<- "MDc"
# trait <- "LDMC"
# spp <- "OBLI"
# 
# model <-  lmer(get(trait)~ get(paste(climvar,"scaled", sep="_")) +  (1|Plottag) + (1|Treetag), data = traits.b, subset=Species==spp, REML=F)
# 
# dataz <- traits.b[which(traits.b$Species==spp & traits.b[,trait]>0),]
# ref.group <- ranef(model)
# ref.var.group1 <- tapply(residuals(model, type="pearson", level=1),
#                          dataz[,names(ref.group[1])], var)
# ref.var.group2 <- tapply(residuals(model, type="pearson", level=1),
#                          dataz[,names(ref.group[2])], var)
# 
# mod.resids <- residuals(model, type=res.type)
# #scatter.smooth(abs(mod.resids)~dataz[,paste0("p",climvar)],xlab=climvar)
# 
# 
# 
# ref.group[[2]]$MDc <- traits.p$MDc[match(rownames(ref.group[[2]]),traits.p$Plottag)]
# ref.group[[2]]$PPTc <- traits.p$PPTc[match(rownames(ref.group[[2]]),traits.p$Plottag)]
# ref.group[[2]]$PETc <- traits.p$PETc[match(rownames(ref.group[[2]]),traits.p$Plottag)]
# ref.group[[2]]$var <- ref.var.group2[match(rownames(ref.group[[2]]), names(ref.var.group2))]
# 
# ref.group[[1]]$MDc <- traits.t$MDc[match(rownames(ref.group[[1]]),traits.t$Treetag)]
# ref.group[[1]]$PETc <- traits.t$PETc[match(rownames(ref.group[[1]]),traits.t$Treetag)]
# ref.group[[1]]$PPTc <- traits.t$PPTc[match(rownames(ref.group[[1]]),traits.t$Treetag)]
# 
# ref.group[[1]]$var <- ref.var.group1[match(rownames(ref.group[[1]]), names(ref.var.group1))]
# 
# dataz$TreeEff <- ref.group[[1]]$`(Intercept)`[match(dataz$Treetag,rownames(ref.group[[1]]))]
# #scatter.smooth(abs(ref.group[[2]]$`(Intercept)`)~ref.group[[2]][,paste0("p",climvar)],xlab=climvar, main="Plot Ran Effs")
# #plot(sqrt(abs(ref.group[[1]]$`(Intercept)`))~ref.group[[1]][,paste0("p",climvar)], dataz,xlab=climvar, main="Tree Ran Effs")
# #abline(lm(sqrt(abs(ref.group[[1]]$`(Intercept)`))~ref.group[[1]][,paste0("p",climvar)], dataz))
# 

palvars <- paste0(brewer.pal(n=9, "Set1"), "AA")
palvarsdark <- brewer.pal(n=9, "Set1")
colchoices <- c(1,2,4,3)

# 
# quartz(width=7, height=4)
# par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(0,0,2,5), mgp=c(2,.7,0))

# plot(sqrt(ref.group[[1]]$var)~ref.group[[1]][,climvar], ylab="SD of LDMC (g/g)", col="grey", xlab="Moisture Deficit (mm)")
# points(sqrt(ref.group[[2]]$var)~ref.group[[2]][,climvar], pch=16, cex=1.5)
# 
# # plot(sdLMA~pPPT, traits.t[which(traits.t$Species=="OVAT"),], ylab="SD of LMA (g/cm2)", col="grey", xlab="Mean annual PPT (mm)")
# # points(sdLMA~pPPT, traits.p[which(traits.p$Species=="OVAT"),], pch=16, cex=1.5)
# # summary(lm(sdLMA~pPPT, traits.p[which(traits.p$Species=="OVAT"),]))
# 
# abline(lm(sqrt(ref.group[[2]]$var)~ref.group[[2]][,climvar]), lwd=2)
# abline(lm(sqrt(ref.group[[1]]$var)~ref.group[[1]][,climvar]), lwd=2, col="grey")
# legend("topright",legend=c("w/in Tree", "w/in Plot"), pch=c(1,16), col=c("grey","black"), lty=1, bty="n")
# mtext(3,adj=0, text="a)    E. obliqua", line=.4)

quartz(width=4, height=4)
jpeg(file=paste0("./",results_dirname,"/Fig6_VarPatterns.jpg"), width=4, height=4, units="in", res=600)
par(mar=c(4,4,1,1), oma=c(0,0,2,5), mgp=c(2,.7,0))
b <- barplot(as.matrix(respats[c(4,3,2,1),c("WD","LMA","LDMC","log.hub")]), beside = F, bg=palvars[colchoices], ylab="# Species", names.arg=rep("", times=4), las=2,density = c(30,30,0,30),angle=c(45,90,0,135), col=palvarsdark[colchoices] )
legend(x=4.8, y=4, legend=c("increase","no trend","non-aridity","decrease"), title ="w/in Spp\nvariance trend\nw/ aridity", fill = palvarsdark[rev(colchoices)], xpd=NA, bty="n", angle=c(135,0,90,45), density=c(30,0,30,30))
axis(1,at = b,labels = c("WD","LMA","LDMC","log(HV)"), las=2, tick = F, mgp=c(2,.2,0))
dev.off()

#quartz.save(file=paste0("./",results_dirname,"/Fig6_VarPatterns.pdf"),type = "pdf" )









######### .. Table S2: Twig characteristics ###############################
ts2 <- traits.b %>% group_by(Species) %>% summarise(Diameter = mean(Stem_Di_corrected, na.rm=T), Diameter_sd = sd(Stem_Di_corrected, na.rm=T),
                                                    Length = mean(Stem_Length, na.rm=T), Length_sd = sd(Stem_Length, na.rm=T),       
                                                    N_Leaves = mean(nleaves, na.rm=T), N_Leaves_sd = sd(nleaves, na.rm=T))

ts2.clean <- data.frame(Species = c("Acacia acuminata","E. amygdalina","Corymbia calophylla","E. marginata","E. salmonophloia","E. obliqua","E. ovata","E. viminalis")
                        , Diameter = paste0(round(ts2$Diameter,2)," (",round(ts2$Diameter_sd,2),")")
                        , Length = paste0(round(ts2$Length,2)," (",round(ts2$Length_sd,2),")")
                        , N_Leaves = paste0(round(ts2$N_Leaves,2)," (",round(ts2$N_Leaves_sd,2),")"))

write.csv(ts2.clean, paste0(results_dirname,"/Table_S2_twig_characteristics.csv"))
          
          