#!/bin/bash/R

############ Analysis for 5 species dataset ############

########## OPTPARSE ############










#################################


########### INPUT FOR WORKING ############
otutable <- ""
MFFP <- "/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/MF_and_OTU_edited/MF_withalpha.txt"
dmFP <- "/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/beta_div/bray_curtis_dm.txt"

############## LOAD DATA #####################
library(MASS)
library(vegan)
set.seed(10234)
setwd("/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/TESTING\ PLOTS")


MF <- read.delim(paste0(MFFP), header = TRUE, row.names = 1)

dm_BC <- read.delim(paste0(dmFP), header = TRUE, row.names = 1)

# Make all the same order; also, inadvertently filters MF to only include dm_BC values
MF.sorted <- MF[match(rownames(dm_BC), rownames(MF)),]

# Make NMDS
dm_BC_NMDS <- isoMDS(d = dist(dm_BC), k = 2)
# Find max dimensions so I can make all plots the same
maxNMDS <- max(dm_BC_NMDS$points)
minNMDS <- min(dm_BC_NMDS$points)

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

# Okay, for each species...let's plot this out.
speciesList <- levels(MF.sorted$SPEC)
speciesList <- speciesList[speciesList != "None"]

# Make a time list
timeList <- levels(factor(MF.sorted$timetotal))



########### Plot BD load in histogram ############
# See if some of them should be zero
Bd1 <- MF.sorted$Bd_Run_1
Bd2 <- MF.sorted$Bd_Run_2
Bd3 <- MF.sorted$Bd_Average_Run_3
allBd <- cbind(Bd1,Bd2,Bd3) # Combine for later

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


########## Fix BD loads if inconsistent ##############
# If there are 2 reads that are '0' and the third one is less than 50, set it all to zero
for ( r in 1:nrow(allBd) ) {
    if ( (sum(allBd[r,] == 0) == 2) & (max(allBd[r,],na.rm=TRUE) < 50) ) {
        print(allBd[r,])
        allBd[r,] <- c(0,0,0)
    }
}

######## Check to see how many toads/frogs were NEVER infected ########
infectionRate <- aggregate(allBd, by = list(CON = MF.sorted$CONT_NOT_DOSED_BD_UNKNOWNIVE, ID = MF.sorted$ANONYMIZED_NAME, TREAT = MF.sorted$TREATMENT_GROUP), FUN = max)


####### Plot infection per amphibian over time #############

# This makes a list of each species and each individual over time to see what the max BD laod was

bdInfectTime <- list()
for ( sp in speciesList ) {
    # Set up a matrix where columns time points sampled and rows is the frog #
    tempMat <- matrix(ncol = 11, nrow = 6)
    colnames(tempMat) <- seq(1,11,by = 1)
    rownames(tempMat) <- seq(1,6, by = 1)
    for ( s in 1:6 ) {
        for ( t in 1:11 ) {
            bdload <- which((MF.sorted$timetotal == t) & (MF.sorted$ANONYMIZED_NAME == paste0(sp,"_",s)))
            if ( length(bdload) > 0 ) {
                tempMat[s,t] <- max(allBd[bdload,], na.rm = TRUE)
            }
        }
    }
    bdInfectTime[[sp]] <- tempMat
}
bdInfectTime

############# Pre/Post dispersion: are they different? ###############
pre.post.Disp <- list()
for ( sp in speciesList ) {
    pre.post.Disp[[sp]] <- list()
    
    # Keep only Pre/post and the species we want to look at
    toKeep <- MF.sorted$TREATMENT_GROUP != "Tre" & MF.sorted$SPEC == sp
    # Apply to all data
    MF.temp <- MF.sorted[toKeep,]
    dm_BC.temp <- dm_BC[toKeep,toKeep]
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
}

### Make plots of all species, infected or not, pre and post. Highlight ones with BD found on them.
for ( sp in speciesList ) {
    pdf(file = paste0(sp,"_betaplots.pdf"), height = 10, width = 5)
    par(fig = c(0,1,0.8,1), xpd = TRUE)
    plot(0,0, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "")
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
         , col = c("black", "white")[factor(pre.post.Disp[[sp]][["MF"]][["infect"]]$Bd_Average_Run_3 > 0, levels = c("TRUE","FALSE"))]
         , pch = 21
         # , cex = as.numeric(as.character((pre.post.Disp[[sp]][["MF"]][["infect"]]$time)))^(1/2)
         , xlim = c(minNMDS, maxNMDS)
         , ylim = c(minNMDS, maxNMDS)
         , xlab = "NMDS1"
         , ylab = "NMDS2"
         , main = "Infected: Pre and Post"
         , sub = paste0("p=",signif(pre.post.Disp[[sp]][["Infected:prevspos"]]$`Pr(>F)`[1],2))
    )
    par(fig = c(0,1,0,0.4), new = TRUE)
    plot(pre.post.Disp[[sp]][["NMDS"]][["noninfect"]]
         , bg = c("red","blue")[factor(pre.post.Disp[[sp]][["MF"]][["noninfect"]]$TREATMENT_GROUP)]
         , col = c("black", "white")[factor(pre.post.Disp[[sp]][["MF"]][["noninfect"]]$Bd_Average_Run_3 > 0, levels = c("TRUE","FALSE"))]
         , pch = 21
         # , cex = as.numeric(as.character((pre.post.Disp[[sp]][["MF"]][["noninfect"]]$time)))^(1/2)
         , xlim = c(minNMDS, maxNMDS)
         , ylim = c(minNMDS, maxNMDS)
         , xlab = "NMDS1"
         , ylab = "NMDS2"
         , main = "Not infected: Pre and Post"
         , sub = paste0("p=",signif(pre.post.Disp[[sp]][["Noninfected:prevspos"]]$`Pr(>F)`[1],2))
         
    )
    dev.off()
}


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
        loadAve <- mean(pre.post.Disp[[sp]][["MF"]][["infect"]]$Bd_Average_Run_3[which(individual & postInfectIndiv)])
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
        loadAve <- mean(pre.post.Disp[[sp]][["MF"]][["infect"]]$Bd_Average_Run_3[individual])
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
for (s in 1: length(speciesList)) {
    minmaxLoad[s,] <- range(correlationData[[speciesList[s]]][["correlation"]][,2])
    minmaxDist[s,] <- range(correlationData[[speciesList[s]]][["correlation"]][,1])
}
minPlotLoad <- min(minmaxLoad[,1], na.rm = TRUE)
maxPlotLoad <- max(minmaxLoad[,2], na.rm = TRUE)
minPlotDist <- min(minmaxDist[,1], na.rm = TRUE)
maxPlotDist <- max(minmaxDist[,2], na.rm = TRUE)

pdf(file="DistancevsLoad_allspecies_NMDSspace.pdf")
plot(NA,NA, type = "n", xlim = c(minPlotDist, maxPlotDist), ylim = c(minPlotLoad, maxPlotLoad)
     , xlab = "Average Distance Travelled per individual in NMDS space"
     , ylab = "Average Load (log)"
)
for (s in 1:length(speciesList)) {
    points(correlationData[[speciesList[s]]][["correlation"]]
           , pch = 19
           , col = colSpecies[s])
}
dev.off()

#### Test: Is there a relationship between how much you move in PCOA space and how likely you are to be infected?
