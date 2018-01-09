#!/bin/bash/R

############ Analysis for 5 species dataset ############

# Comments so far:
# Histograms show BD load WITHOUT filtering
    # Also, NAs are left as "NAs"; aka, nothing is artificially turned to zero.

# In the table showing BD load of individuals across time, it uses the filtered BD load
# BD load is filtered by the following rules:
    # If any individual read is < 10, then it is set to '0'.
    # If 2 out of 3 samples are '0' and the third sample is <50, then everyting is set to '0'
        # Hopefully, this threshold keeps ones that are truly 'low' load, but discards reads that are noise.


########## OPTPARSE ############










#################################


########### INPUT FOR WORKING ############
otutable <- "/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/MF_and_OTU_edited/otu_table_text.txt"
MFFP <- "/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/MF_and_OTU_edited/MF_withalpha.txt"
dmFP <- "/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/beta_div/bray_curtis_dm.txt"

############## LOAD AND FILTER DATA #####################
library(MASS)
library(vegan)
set.seed(10234)
setwd("/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/")
system("mkdir 5SpeciesAnalysis")
setwd("5SpeciesAnalysis")


# Load mappingfile and dm and otutable
MF <- read.delim(paste0(MFFP), header = TRUE, row.names = 1)
dm_BC <- read.delim(paste0(dmFP), header = TRUE, row.names = 1)
    # Below line is commented out bc it takes a long time to load and I am trying to save time
# otutable <- read.delim(paste0(otutable), header=TRUE, row.names = 1, skip = 1)

##### ADJUSTING DATA: Make same order, make 'time' variable ######
MF.sorted <- MF

# Make a 'time' variable
time <- as.numeric(gsub(".*_","",MF.sorted$PRE_POST_BD_NUM))
MF.sorted <- cbind(MF.sorted, time)
# make a 'time2' variable for tracking across all time
prepos <- MF.sorted$TREATMENT_GROUP
timetotal <- time
for ( i in 1:length(prepos) ) {
    if ( prepos[i] == "Pos" ) {
        timetotal[i] <- as.numeric(timetotal[i]) + 5
    } else if ( prepos[i] == "Tre" ) {
        timetotal[i] <- as.numeric(timetotal[i]) + 16
    }
}
MF.sorted <- cbind(MF.sorted, timetotal)

# Change "NA"s for BD run to '0'
MF.sorted$Bd_Average_Run_3[is.na(MF.sorted$Bd_Average_Run_3)] <- 0
MF.sorted$Bd_Run_1[is.na(MF.sorted$Bd_Run_1)] <- 0
MF.sorted$Bd_Run_2[is.na(MF.sorted$Bd_Run_2)] <- 0


########### Plot BD load in histogram ############
# See if some of them should be zero
Bd1 <- MF.sorted$Bd_Run_1
Bd2 <- MF.sorted$Bd_Run_2
Bd3 <- MF.sorted$Bd_Average_Run_3
allBd <- cbind(Bd1,Bd2,Bd3) # Combine for later
rownames(allBd) <- rownames(MF.sorted)

Bd1[Bd1==0] <- NA
Bd2[Bd2==0] <- NA
Bd3[Bd3==0] <- NA

pdf(file="histogram_bdload.pdf", height = 10, width =5)
par(mfrow = c(3,1))
hist(log10(Bd3+1)
     , xlim = c(0,10)
     , breaks = seq(0,10,by = 0.1))

hist(log10(Bd1+1)
     , xlim = c(0,10)
     , breaks = seq(0,10,by = 0.1))

hist(log10(Bd2+1)
     , xlim = c(0,10)
     , breaks = seq(0,10,by=0.1))
dev.off()

# If 2 BD samples are 0 and the third is less than 50, then get rid of it.
# Also, if anything is less than 10, make it zero anyway.
for ( r in 1:nrow(allBd) ) {
    for ( c in 1:ncol(allBd) ) {
        if ( (allBd[r,c] < 10) & (!is.na(allBd[r,c])) ) {
            allBd[r,c] <- 0
        }
    }
    if ( (sum(allBd[r,] == 0) == 2) & (max(allBd[r,],na.rm=TRUE) < 50) ) {
        allBd[r,] <- c(0,0,0)
    }

}

######## Check to see how many toads/frogs were NEVER infected ########
infectionRate <- aggregate(allBd, by = list(CON = MF.sorted$CONT_NOT_DOSED_BD_UNKNOWNIVE, ID = MF.sorted$ANONYMIZED_NAME, TREAT = MF.sorted$TREATMENT_GROUP), FUN = max)
write.table(infectionRate, file = "infectionRate.txt", sep = "\t")

####### Plot infection per amphibian over time #############

# Get list of species
speciesList <- levels(MF.sorted$SPEC)
speciesList <- speciesList[speciesList != "None"]

# Make a time list
timeList <- levels(factor(MF.sorted$timetotal))

# makes list of species
# In each species, you track individuals (n=6) over time (n=11). value is max BD load for that individual at that time.
bdInfectTime <- list()
MF.sorted[,c("maxInfect","meanInfect")] <- matrix(ncol = 2, nrow = nrow(MF.sorted))
for ( sp in speciesList ) {
    # Make matrix for  individuals in first experiment; first 6 are infected ones
    tempMat <- matrix(ncol = 16, nrow = 10)
    colnames(tempMat) <- seq(1,16,by = 1)
    rownames(tempMat) <- seq(1,10, by = 1)
    tempMatAve <- matrix(ncol = 16, nrow = 10)
    colnames(tempMatAve) <- seq(1,16,by = 1)
    rownames(tempMatAve) <- seq(1,10, by = 1)
    for ( s in 1:10 ) {
        for ( t in 1:16 ) {
            bdload <- which((MF.sorted$timetotal == t) & (MF.sorted$ANONYMIZED_NAME == paste0(sp,"_",s)))
            if ( length(bdload) > 0 ) {
                tempMat[paste0(s),paste0(t)] <- max(allBd[bdload,], na.rm = TRUE)
                tempMatAve[paste0(s),paste0(t)] <- mean(allBd[bdload,], na.rm = TRUE)
                MF.sorted[bdload, c("maxInfect","meanInfect")] <- c(tempMat[paste0(s),paste0(t)],tempMatAve[paste0(s),paste0(t)])
            }
        }
    }
    # Now, make the same table but for "Tre" treatment--- include later if you want. 
    # Commenting out for now bc it's only Bubo
    # tempMat2 <- matrix(ncol = (27-16), nrow = 10)
    # colnames(tempMat2) <- seq(17,27,by=1)
    # rownames(tempMat2) <- seq(1,10,by=1)
    # for ( s in 1:10 ) {
    #     for ( t in 17:27 ) {
    #         bdload <- which((MF.sorted$timetotal == t) & (MF.sorted$ANONYMIZED_NAME == paste0(sp,"_",s)))
    #         MF.sorted[bdload,]
    #         if ( length(bdload) > 0 ) {
    #             tempMat2[paste0(s),paste0(t)] <- max(allBd[bdload,], na.rm = TRUE)
    #         }
    #     }
    # }
    bdInfectTime[[sp]] <- list()
    bdInfectTime[[sp]][["Max"]] <- tempMat
    bdInfectTime[[sp]][["Mean"]] <- tempMatAve
    # bdInfectTime[[sp]][["Exp2"]] <- tempMat2
}

capture.output(bdInfectTime, file = "bdInfecTime.txt")
# Plot each individual taod to see how their infection changes over time

# Each individual is going to be a color of the rainbow
rainbowCol <- c("darkred","red","orange","yellow","green","darkgreen","blue","darkblue","purple","grey")
# Get the max log load of all samples so that we can compare across species
maxLoadLog <- log(max(unlist(bdInfectTime), na.rm=TRUE) +1)

# Make a 6-panel plot for Max and Mean infection over time
pdf("MaxInfectionOverTime_all.pdf", width = 5, height = 8)
par(mfrow=c(3,2))
for ( sp in speciesList ) {
    tableTemp <- bdInfectTime[[sp]][["Max"]]
    minT <- 1
    maxT <- ncol(tableTemp)
    plot(NULL, ylim=c(0,maxLoadLog), xlim=c(minT,maxT),xlab="Time", ylab="BD Load (ln)", xaxt="n", main = paste0(sp))
    axis(side=1, at = seq(minT,maxT), labels = seq(minT,maxT), cex.axis=0.8)
    for ( r in 1:nrow(tableTemp) ) {
        points(log(tableTemp[r,]+1), type = "l", col=rainbowCol[r])
    }
    abline(v = 6, lty = 2, col = "black")
}
# Last panel is the legend
plot(0,0,pch="",axes=FALSE,xlab="",ylab="")
legend("left"
       , legend= c(sapply(seq(1,10),function(x) paste0("Individual ",x)), "(7-10 = Controls)")
       , lty = 1
       , col = c(rainbowCol, NA)
       , bty = "n")
dev.off()

pdf("MeanInfectionOverTime_all.pdf", width = 5, height = 8)
par(mfrow=c(3,2))
for ( sp in speciesList ) {
    tableTemp <- bdInfectTime[[sp]][["Mean"]]
    minT <- 1
    maxT <- ncol(tableTemp)
    plot(NULL, ylim=c(0,maxLoadLog), xlim=c(minT,maxT),xlab="Time", ylab="BD Load (ln)", xaxt="n", main = paste0(sp))
    axis(side=1, at = seq(minT,maxT), labels = seq(minT,maxT), cex.axis=0.8)
    for ( r in 1:nrow(tableTemp) ) {
        points(log(tableTemp[r,]+1), type = "l", col=rainbowCol[r])
    }
    abline(v = 6, lty = 2, col = "black")
}
# Last panel is the legend
plot(0,0,pch="",axes=FALSE,xlab="",ylab="")
legend("left"
       , legend= c(sapply(seq(1,10),function(x) paste0("Individual ",x)), "(7-10 = Controls)")
       , lty = 1
       , col = c(rainbowCol, NA)
       , bty = "n")
dev.off()

########## PLOTTING CONTROLS #########
# First thing to do is to get rid of the reinfected samples
# Get all sample names that are "Treatment"
names_treat <- rownames(MF.sorted)[MF.sorted[,"TREATMENT_GROUP"] == "Tre"]

# Then, filter OTU table, MF, and dm by this
MF.filt <- MF.sorted[-match(names_treat,rownames(MF.sorted)),]
dm_BC.filt <- dm_BC[-which(rownames(dm_BC)%in%names_treat),-which(colnames(dm_BC)%in%names_treat)]
otutable.filt <- otutable[,-which(colnames(otutable)%in%names_treat)]

# Make sure dm and MF.filt is in same order
MF.filt <- MF.filt[match(rownames(dm_BC.filt),rownames(MF.filt) ),]

# Make NMDS
dm_BC_NMDS_wcontrols <- isoMDS(d = dist(dm_BC.filt), k = 2)
# Find max dimensions so I can make all plots the same
maxNMDS <- max(dm_BC_NMDS_wcontrols$points)
minNMDS <- min(dm_BC_NMDS_wcontrols$points)

# Get the species and control name and get rid of all numbers and underscores
MF.filt[,"SampleType"] <- gsub('[0-9]|_','',MF.filt$SPECIES_NAME)
# It also turns out there is one called 'Conrol' instead of "Control" so fix that too
MF.filt[,"SampleType"] <- gsub("Conrol","Control",MF.filt[,"SampleType"])
# Also, "cricket" vs "crickets"
MF.filt[,"SampleType"] <- gsub("Crickets","Cricket",MF.filt[,"SampleType"])
# Also, typo in sterile water
MF.filt[,"SampleType"] <- gsub("Sterlile","Sterile",MF.filt[,"SampleType"])

# Check to make sure it's still same order
rownames(MF.filt) == rownames(dm_BC_NMDS_wcontrols$points)
# Make colors for all of these species and control names
MF.filt[,"SampleType"] <- factor(MF.filt[,"SampleType"]
              , levels = c("Bufo Boreas"
                           ,"Bufo Marinus"
                           ,"Osteopilus septentrionalis"
                           ,"Rana pipiens"
                           ,"Rana Catesbeiana"
                           ,"ControlCricket"
                           ,"ControlHoltfreter"
                           ,"ControlMold"
                           ,"ControlSterile Water")
              )
colSpecies <- cbind(colors = c("red"
                      , "blue"
                      ,"yellow"
                      ,"purple"
                      ,"green"
                      ,"burlywood4"
                      ,"black"
                      ,"darkgreen"
                      ,"lightblue")
                    , type = levels(factor(MF.filt[,"SampleType"]
                                  , levels = c("Bufo Boreas"
                                               ,"Bufo Marinus"
                                               ,"Osteopilus septentrionalis"
                                               ,"Rana pipiens"
                                               ,"Rana Catesbeiana"
                                               ,"ControlCricket"
                                               ,"ControlHoltfreter"
                                               ,"ControlMold"
                                               ,"ControlSterile Water")
                    )))
pdf("Plot_Controls.pdf",7,5)
par(fig=c(0,0.7,0,1))
plot(dm_BC_NMDS_wcontrols$points
     , col = colSpecies[,"colors"][MF.filt[,"SampleType"]]
     , pch = c(21,21,21,21,21,19,19,19,19)[MF.filt[,"SampleType"]])
par(fig=c(0.68,1,0,1), mar=c(0,0,0,0),new=TRUE)
plot(0,0,xlab="",ylab="",pch="",axes=FALSE)
legend("left"
       , bty = 'n'
       , legend = levels(MF.filt[,"SampleType"])
       , col = colSpecies[,"colors"]
       , pch = c(21,21,21,21,21,19,19,19,19)
)
dev.off()


########## PLOT TREATMENTS ###########
# First, filter out controls
MF.filt.nocon <- MF.filt[-grep("(c|C)ontrol",rownames(MF.filt)),]
dm_BC.filt.nocon <- dm_BC.filt[rownames(MF.filt.nocon),rownames(MF.filt.nocon)]
dm_BC_NMDS <- isoMDS(dist(dm_BC.filt.nocon), k = 2)
minNMDS_nocon <- min(dm_BC_NMDS$points)
maxNMDS_nocon <- max(dm_BC_NMDS$points)
stress_nocon <- dm_BC_NMDS$stress

# List of names "before" 
preNames <- rownames(MF.filt.nocon)[grep("Pre",MF.filt.nocon$TREATMENT_GROUP)]

# find mean dispersion for each species in the PRE group
MF_pre <- MF.filt.nocon[MF.filt.nocon$TREATMENT_GROUP == "Pre",]
dm_BC_pre <- dm_BC.filt.nocon[match(rownames(MF_pre), rownames(dm_BC.filt.nocon)),match(rownames(MF_pre), rownames(dm_BC.filt.nocon))]
betadisper_pre <- betadisper(dist(dm_BC_pre), group = MF_pre$SPEC)
distances_pre <- betadisper_pre$distances

# Gamma distribution
gammaNLL <- function(pars, data) {
    alpha <- pars[1]
    theta <- pars[2]
    return(-sum(dgamma(data,shape=alpha, scale = theta,log=TRUE)))
}

preDispersion <- list()
for ( sp in speciesList ) {
    # Make list for species for each individual.
    preDispersion[[paste0(sp)]] <- list()
    # Get distances for each species: This is the data.
    sp_temp <- distances_pre[grep(sp,names(distances_pre))]
    u <- mean(log(sp_temp))
    v <- sd(log(sp_temp))
    # Find fit for gamma distribution-- even though I KNOW it's not normal
    gamma_optim <- optim(par=c(3,0.5), fn = gammaNLL, data = sp_temp)
    # Get parameters
    fit_alpha <- gamma_optim$par[1]
    fit_theta <- gamma_optim$par[2]
    preDispersion[[paste0(sp)]][["par"]] <- c(fit_alpha,fit_theta)
}



############# Pre/Post dispersion: are they different? ###############
pre.post.Disp <- list()
for ( sp in speciesList ) {
    pre.post.Disp[[sp]] <- list()
    
    # Keep only Pre/post and the species we want to look at
    toKeep <- MF.filt.nocon$SPEC == sp
    # Apply to all data
    MF.temp <- MF.filt.nocon[toKeep,]
    dm_BC.temp <- dm_BC.filt.nocon[toKeep,toKeep]
    dm_BC_NMDS.temp <- dm_BC_NMDS$points[toKeep,]
    
    # Now, split into infected or not infected
    infect <- MF.temp$BD100KZSP_3122011 == "y"
    noninfect <- MF.temp$BD100KZSP_3122011 == "n"
    # Apply to all data
    MF.temp.infect <- MF.temp[infect,] 
    MF.temp.noninfect <- MF.temp[noninfect,]
    dm_BC.temp.infect <- dm_BC.temp[infect,infect]
    dm_BC.temp.noninfect <- dm_BC.temp[noninfect,noninfect]
    dm_BC_NMDS.temp.infect <- dm_BC_NMDS.temp[infect,]
    dm_BC_NMDS.temp.noninfect <- dm_BC_NMDS.temp[noninfect,]
    
    # Now, split into pre- and post- infection
    infect.pre <- MF.temp.infect$TREATMENT_GROUP == "Pre"
    infect.post <- MF.temp.infect$TREATMENT_GROUP == "Pos"
    noninfect.pre <- MF.temp.noninfect$TREATMENT_GROUP == "Pre"
    noninfect.post <- MF.temp.noninfect$TREATMENT_GROUP == "Pos"
    # Apply to data
    MF.temp.infect.pre <- MF.temp.infect[infect.pre,]
    MF.temp.noninfect.pre <- MF.temp.noninfect[noninfect.pre,]
    MF.temp.infect.post <- MF.temp.infect[infect.post,]
    MF.temp.noninfect.post <- MF.temp.noninfect[noninfect.post,]
    
    dm_BC.temp.infect.pre <- dm_BC.temp.infect[infect.pre,infect.pre]
    dm_BC.temp.noninfect.pre <- dm_BC.temp.noninfect[noninfect.pre,noninfect.pre]
    dm_BC.temp.infect.post <- dm_BC.temp.infect[infect.post,infect.post]
    dm_BC.temp.noninfect.post <- dm_BC.temp.noninfect[noninfect.post,noninfect.post]
    
    # Now, run betadisper on all combos
    # Infected Pre
    infect.pre.betadisper <- betadisper(d = dist(dm_BC.temp.infect.pre), group = as.factor(MF.temp.infect.pre$time) )
    infect.pre.lm <- lm(infect.pre.betadisper$distances  ~ as.numeric(as.character(MF.temp.infect.pre$time)))
    pre.post.Disp[[sp]][["Infected:Pre"]] <- anova(infect.pre.lm)
    # Infected Post
    infect.post.betadisper <- betadisper(d = dist(dm_BC.temp.infect.post), group = as.factor(MF.temp.infect.post$time) )
    infect.post.lm <- lm(infect.post.betadisper$distances  ~ as.numeric(as.character(MF.temp.infect.post$time)))
    pre.post.Disp[[sp]][["Infected:Post"]] <- anova(infect.post.lm)
    # Noninfected Pre
    noninfect.pre.betadisper <- betadisper(d = dist(dm_BC.temp.noninfect.pre), group = as.factor(MF.temp.noninfect.pre$time) )
    noninfect.pre.lm <- lm(noninfect.pre.betadisper$distances  ~ as.numeric(as.character(MF.temp.noninfect.pre$time)))
    pre.post.Disp[[sp]][["Noninfected:Pre"]] <- anova(noninfect.pre.lm)
    # Noninfected Post
    noninfect.post.betadisper <- betadisper(d = dist(dm_BC.temp.noninfect.post), group = as.factor(MF.temp.noninfect.post$time) )
    noninfect.post.lm <- lm(noninfect.post.betadisper$distances  ~ as.numeric(as.character(MF.temp.noninfect.post$time)))
    pre.post.Disp[[sp]][["Noninfected:Post"]] <- anova(noninfect.post.lm)
    
    # Compare pre vs post: infected
    infect.betadisper <- betadisper(d = dist(dm_BC.temp.infect), group = factor(MF.temp.infect$TREATMENT_GROUP))
    pre.post.Disp[[sp]][["Infected:prevspos"]] <- anova(infect.betadisper)
    # Compare pre vs post: noninfected
    noninfect.betadisper <- betadisper(d = dist(dm_BC.temp.noninfect), group = factor(MF.temp.noninfect$TREATMENT_GROUP))
    pre.post.Disp[[sp]][["Noninfected:prevspos"]] <- anova(noninfect.betadisper)
    
    # Now, save all this information!
    pre.post.Disp[[sp]][["distanceMat"]] <- list()
    pre.post.Disp[[sp]][["distanceMat"]][["all"]] <- dm_BC.temp
    pre.post.Disp[[sp]][["distanceMat"]][["infect"]] <- dm_BC.temp.infect
    pre.post.Disp[[sp]][["distanceMat"]][["noninfect"]] <- dm_BC.temp.noninfect
    pre.post.Disp[[sp]][["distanceMat"]][["infectpre"]] <- dm_BC.temp.infect.pre
    pre.post.Disp[[sp]][["distanceMat"]][["infectpost"]] <- dm_BC.temp.infect.post
    pre.post.Disp[[sp]][["distanceMat"]][["noninfectpre"]] <- dm_BC.temp.noninfect.pre
    pre.post.Disp[[sp]][["distanceMat"]][["noninfectpost"]] <- dm_BC.temp.noninfect.post
    
    pre.post.Disp[[sp]][["NMDS"]] <- list()
    pre.post.Disp[[sp]][["NMDS"]][["all"]] <- dm_BC_NMDS.temp
    pre.post.Disp[[sp]][["NMDS"]][["infect"]] <- dm_BC_NMDS.temp.infect
    pre.post.Disp[[sp]][["NMDS"]][["noninfect"]] <- dm_BC_NMDS.temp.noninfect
    
    pre.post.Disp[[sp]][["MF"]] <- list()
    pre.post.Disp[[sp]][["MF"]][["all"]] <- MF.temp
    pre.post.Disp[[sp]][["MF"]][["infect"]] <- MF.temp.infect
    pre.post.Disp[[sp]][["MF"]][["noninfect"]] <- MF.temp.noninfect
    pre.post.Disp[[sp]][["MF"]][["infectpre"]] <- MF.temp.infect.pre
    pre.post.Disp[[sp]][["MF"]][["infectpost"]] <- MF.temp.infect.post
    pre.post.Disp[[sp]][["MF"]][["noninfectpre"]] <- MF.temp.noninfect.pre
    pre.post.Disp[[sp]][["MF"]][["noninfectpost"]] <- MF.temp.noninfect.post
    
    pre.post.Disp[[sp]][["betadisp"]] <- list()
    pre.post.Disp[[sp]][["betadisp"]][["infectpre"]] <- infect.pre.betadisper
    pre.post.Disp[[sp]][["betadisp"]][["infectpost"]] <- infect.post.betadisper
    pre.post.Disp[[sp]][["betadisp"]][["noninfectpre"]] <- noninfect.pre.betadisper
    pre.post.Disp[[sp]][["betadisp"]][["noninfectpost"]] <- noninfect.post.betadisper
}

### Make plots of all species, infected or not, pre and post. Highlight ones with BD found on them.
for ( sp in speciesList ) {
    pdf(file = paste0(sp,"_betaplots.pdf"), height = 10, width = 5)
    par(fig = c(0,1,0.8,1), xpd = TRUE)
    plot(0,0, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "", sub = paste0("Stress: ",signif(stress_nocon/100,2)))
    legend("center"
           , legend = c("Pre infection","Post infection", "BD detected", "no BD detected")
           , pch = 21
           , pt.bg = c("blue", "red", "grey","grey")
           , col = c("white","white","black","white")
           , ncol = 2
    )
    par(fig = c(0,1,0.4,0.8), new = TRUE)
    plot(pre.post.Disp[[sp]][["NMDS"]][["infect"]]
         , bg = c("red","blue")[factor(pre.post.Disp[[sp]][["MF"]][["infect"]]$TREATMENT_GROUP)]
         , col = c("black", "white")[factor(pre.post.Disp[[sp]][["MF"]][["infect"]]$maxInfect > 0, levels = c("TRUE","FALSE"))]
         , pch = 21
         # , cex = as.numeric(as.character((pre.post.Disp[[sp]][["MF"]][["infect"]]$time)))^(1/2)
         , xlim = c(minNMDS_nocon, maxNMDS_nocon)
         , ylim = c(minNMDS_nocon, maxNMDS_nocon)
         , xlab = "NMDS1"
         , ylab = "NMDS2"
         , main = "Infected: Pre and Post"
         , sub = paste0("p=",signif(pre.post.Disp[[sp]][["Infected:prevspos"]]$`Pr(>F)`[1],2))
    )
    par(fig = c(0,1,0,0.4), new = TRUE)
    plot(pre.post.Disp[[sp]][["NMDS"]][["noninfect"]]
         , bg = c("red","blue")[factor(pre.post.Disp[[sp]][["MF"]][["noninfect"]]$TREATMENT_GROUP)]
         , col = c("black", "white")[factor(pre.post.Disp[[sp]][["MF"]][["noninfect"]]$maxInfect > 0, levels = c("TRUE","FALSE"))]
         , pch = 21
         # , cex = as.numeric(as.character((pre.post.Disp[[sp]][["MF"]][["noninfect"]]$time)))^(1/2)
         , xlim = c(minNMDS_nocon, maxNMDS_nocon)
         , ylim = c(minNMDS_nocon, maxNMDS_nocon)
         , xlab = "NMDS1"
         , ylab = "NMDS2"
         , main = "Not infected: Pre and Post"
         , sub = paste0("p=",signif(pre.post.Disp[[sp]][["Noninfected:prevspos"]]$`Pr(>F)`[1],2))
         
    )
    dev.off()
}


######## AVE LOAD VS AVE DIST ##########
#### Calculate how much individuals travel in pre-inoculation populations for EACH species- NMDS space ######
correlationData <- list()
for ( sp in speciesList ) {
    correlationData[[sp]] <- list()
    travelDist <- matrix(ncol = 6, nrow = length(timeList))
    NMDS.temp <- pre.post.Disp[[sp]][["NMDS"]][["all"]]
    for ( n in 1:6 ) {
        t <- 1
        individual <- pre.post.Disp[[sp]][["MF"]][["infect"]]$ANONYMIZED_NAME == paste0(sp,"_",n)
        first <- c()
        while ( t < length(timeList) ) {
            t <- (t + 1)
            newfirst <- NMDS.temp[which(individual & (pre.post.Disp[[sp]][["MF"]][["infect"]]$timetotal == (t-1))),]
            newsecond <- NMDS.temp[which(individual & (pre.post.Disp[[sp]][["MF"]][["infect"]]$timetotal == t)),]
            if ( length(newfirst) > 0 ) {
                first <- newfirst
            }
            if ( length(first) == 0 ) {
                next
            }
            while (( length(newsecond) == 0 ) & (t < length(timeList))) {
                t <- t + 1
                newsecond <- NMDS.temp[which(individual & (pre.post.Disp[[sp]][["MF"]][["infect"]]$timetotal == (t))),]
            }
            second <- newsecond
            if ( length(second) == 0 ) {
                next
            }
            riserun <- second-first
            distance <- sqrt(riserun[1]^2 + riserun[2]^2)
            
            travelDist[t,n] <- distance
        }
    }
    correlationData[[sp]][["travelDist"]] <- travelDist
    
    # Now that we have the distance between each timepoint, we should plot the distance travelled by infection rate
    correlation <- matrix(ncol = 2, nrow = 6)
    colnames(correlation) <- c("distance","load")
    for ( n in 1:6 ) {
        individual <- pre.post.Disp[[sp]][["MF"]][["infect"]]$ANONYMIZED_NAME == paste0(sp,"_",n)
        postInfectIndiv <- pre.post.Disp[[sp]][["MF"]][["infect"]]$TREATMENT_GROUP == "Pos"
        loadAve <- mean(pre.post.Disp[[sp]][["MF"]][["infect"]]$meanInfect[which(individual & postInfectIndiv)])
        distAve <- mean(travelDist[,n],na.rm = TRUE)
        correlation[n,] <- c(distAve, log(loadAve + 1))
    }
    correlationData[[sp]][["correlation"]] <- correlation
}


# Colors for species, taken from presentation
colSpecies <- c("red", "blue","yellow","purple","green")

# Find min and max of all correlation data
minmaxLoad <- matrix(ncol = 2, nrow = 5)
minmaxDist <- matrix(ncol = 2, nrow = 5)
allCorrelation <- matrix(ncol = 2)
for (s in 1: length(speciesList)) {
    minmaxLoad[s,] <- range(correlationData[[speciesList[s]]][["correlation"]][,2])
    minmaxDist[s,] <- range(correlationData[[speciesList[s]]][["correlation"]][,1])
    allCorrelation <-  rbind(allCorrelation, correlationData[[speciesList[s]]][["correlation"]])
}
minPlotLoad <- min(minmaxLoad[,1], na.rm = TRUE)
maxPlotLoad <- max(minmaxLoad[,2], na.rm = TRUE)
minPlotDist <- min(minmaxDist[,1], na.rm = TRUE)
maxPlotDist <- max(minmaxDist[,2], na.rm = TRUE)

# Statistical test for regression of allCorrelation
correlation.lm <- lm(allCorrelation[,2] ~ allCorrelation[,1])
summary(correlation.lm)
correlation.lm.anova <- anova(correlation.lm)

pdf(file="DistancevsLoad_allspecies_NMDSspace.pdf")
plot(NA,NA, type = "n", xlim = c(minPlotDist, maxPlotDist), ylim = c(minPlotLoad, maxPlotLoad)
     , xlab = "Average Distance Travelled per individual in NMDS space"
     , ylab = "Average Load (log)"
     , sub = paste0("Regression: p = ", signif(correlation.lm.anova$`Pr(>F)`[1],2))
)
for (s in 1:length(speciesList)) {
    points(correlationData[[speciesList[s]]][["correlation"]]
           , pch = 19
           , col = colSpecies[s])
}
dev.off()



#### Calculate how much individuals travel in pre-inoculation populations for EACH species- BC space ######
for ( sp in speciesList ) {
    # Make enpty list
    correlationData[[sp]] <- list()
    # Make matrix to record how much they've travelled
    travelDist <- matrix(ncol = 6, nrow = length(6:11))
    # Get the dm with post-infected individuals
    betadisp.temp <- pre.post.Disp[[sp]][["betadisp"]][["infectpost"]]$distances
    
    for ( n in 1:6 ) { # For each individual taod
        individual <- pre.post.Disp[[sp]][["MF"]][["infectpost"]]$ANONYMIZED_NAME == paste0(sp,"_",n)
        for ( t in 1:6 ) {
            timet <- pre.post.Disp[[sp]][["MF"]][["infectpost"]]$time == t+5
            dist.temp <- betadisp.temp[which(individual & timet)]
            if (length(dist.temp) > 0) {
                travelDist[t,n] <- dist.temp
            }
        }
    }
    correlationData[[sp]][["travelDistBC"]] <- travelDist
    
    # Now that we have the distance between each timepoint, we should plot the distance travelled by infection rate
    correlation <- matrix(ncol = 2, nrow = 6)
    colnames(correlation) <- c("distance","load")
    for ( n in 1:6 ) {
        individual <- pre.post.Disp[[sp]][["MF"]][["infectpost"]]$ANONYMIZED_NAME == paste0(sp,"_",n)
        loadAve <- mean(pre.post.Disp[[sp]][["MF"]][["infectpost"]]$meanInfect[individual])
        distAve <- mean(travelDist[,n],na.rm = TRUE)
        correlation[n,] <- c(distAve, log(loadAve + 1))
    }
    correlationData[[sp]][["correlationBC"]] <- correlation
}


# Colors for species, taken from presentation
colSpecies <- c("red", "blue","yellow","purple","green")

# Find min and max of all correlation data
minmaxLoad <- matrix(ncol = 2, nrow = 5)
minmaxDist <- matrix(ncol = 2, nrow = 5)
allCorrelation <- matrix(ncol = 2)
for (s in 1: length(speciesList)) {
    minmaxLoad[s,] <- range(correlationData[[speciesList[s]]][["correlationBC"]][,2])
    minmaxDist[s,] <- range(correlationData[[speciesList[s]]][["correlationBC"]][,1])
    allCorrelation <-  rbind(allCorrelation, correlationData[[speciesList[s]]][["correlationBC"]])
}
minPlotLoad <- min(minmaxLoad[,1], na.rm = TRUE)
maxPlotLoad <- max(minmaxLoad[,2], na.rm = TRUE)
minPlotDist <- min(minmaxDist[,1], na.rm = TRUE)
maxPlotDist <- max(minmaxDist[,2], na.rm = TRUE)
# Statistical test for regression of allCorrelation
correlation.lm <- lm(allCorrelation[,2] ~ allCorrelation[,1])
summary(correlation.lm)
correlation.lm.anova <- anova(correlation.lm)

pdf(file="DistancevsLoad_allspecies_BCspace.pdf")
plot(NA,NA, type = "n", xlim = c(minPlotDist, maxPlotDist), ylim = c(minPlotLoad, maxPlotLoad)
     , xlab = "Average Distance Travelled per individual in Bray-curtis space"
     , ylab = "Average Load (log)"
     , sub = paste0("Regression: p = ", signif(correlation.lm.anova$`Pr(>F)`[1],2))
     
)
for (s in 1:length(speciesList)) {
    points(correlationData[[speciesList[s]]][["correlationBC"]]
           , pch = 19
           , col = colSpecies[s])
}
dev.off()

######## IMMEDIATE LOAD vs IMMEDIATE DISTANCE ##########