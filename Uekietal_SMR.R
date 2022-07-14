## title: 'SMR: regularization path plots'
## author: 'H.Hino'
## coded date:2021/04/27
library(MASS)
library(dplyr)
library(ggplot2)
library(caret)
library(glmnet)
library(car)
library(corrplot)
library(fields)
library(plotmo)
library(doParallel)

## specify number of cores available
registerDoParallel(8)

## create directory to store the results
dir.name <- paste0("./result_poly_path_scirep/")
if (!file.exists(dir.name)) {
    dir.create(dir.name)
}

obsData <- read.csv("./Uekietal_SMR_input.csv")

## remove 'NA'
obsData <- na.omit(obsData)

## assign class labels (Tectonic settings)
obsData$Cluster_ID[which(obsData$Cluster_ID == 1)] <- "CA"
obsData$Cluster_ID[which(obsData$Cluster_ID == 2)] <- "IOA"
obsData$Cluster_ID[which(obsData$Cluster_ID == 3)] <- "IA"
obsData$Cluster_ID[which(obsData$Cluster_ID == 4)] <- "BAB"
obsData$Cluster_ID[which(obsData$Cluster_ID == 5)] <- "CFB"
obsData$Cluster_ID[which(obsData$Cluster_ID == 6)] <- "MOR"
obsData$Cluster_ID[which(obsData$Cluster_ID == 7)] <- "OP"
obsData$Cluster_ID[which(obsData$Cluster_ID == 8)] <- "OI"

obsData$Cluster_ID <- as.factor(obsData$Cluster_ID)

cNum <- length(levels(obsData$Cluster_ID))  ## number of classes

## order of classes in accordance with Petrelli and Perugini (2016)
ppOrder <- c("CA", "IOA", "IA", "CFB", "OP", "OI", "BAB", "MOR")

## consider 2nd order combination
cands <- c("SiO2", "TiO2", "Al2O3", "Fe2O3", "MgO", "CaO", "Na2O", "K2O", 
    "La", "Ce", "Nd", "Sm", "Gd", "Yb", "Lu", "Ba", "Hf", "Nb", "Rb", 
    "Sr", "Ta", "Th", "Y", "Zr")
addVar <- NULL
## multiplication
for (i in 1:(length(cands) - 1)) {
    for (j in (i + 1):length(cands)) {
        tmp1 <- obsData %>% dplyr::select(cands[i])
        tmp2 <- obsData %>% dplyr::select(cands[j])
        addVar <- cbind(addVar, as.numeric(unlist(tmp1 * tmp2)))
        colnames(addVar)[dim(addVar)[2]] <- paste(cands[i], cands[j], 
            sep = "*")
    }
}
## division
for (i in 1:(length(cands))) {
    for (j in seq(1:length(cands))[-i]) {
        tmp1 <- obsData %>% dplyr::select(cands[i])
        tmp2 <- obsData %>% dplyr::select(cands[j])
        addVar <- cbind(addVar, as.numeric(unlist(tmp2/tmp1)))
        colnames(addVar)[dim(addVar)[2]] <- paste(cands[j], cands[i], 
            sep = "/")
    }
}

obsData <- cbind(obsData, addVar)

featN <- dim(obsData)[2] - 2
obsData <- dplyr::select(obsData, -Number)

X <- dplyr::select(obsData, -Cluster_ID)
PPvars <- colnames(X)

## remove feature with variance zero, which is not informative
if (sum(is.na(apply(obsData[, -1], 2, var))) >= 1) {
    obsData <- obsData[, -which(is.na(apply(obsData[, -1], 2, var))) + 
        1]
}

## split dataset for cross validation
set.seed(123)  ## set the random number seed to asure reproducibility
cvIDs <- createFolds(obsData$Cluster_ID, k = 10)

### box-cox transformation followed by scaling
cvErrsSVMBC <- NULL
BCdata <- obsData
BCdata[, -1] <- scale(BCdata[, -1])
for (i in 1:dim(X)[2]) {
    tmp <- powerTransform(BCdata[, i + 1] - min(BCdata[, i + 1]) + 0.01)
    
    BCdata[, i + 1] <- bcPower(BCdata[, i + 1] - min(BCdata[, i + 1]) + 
        0.1, tmp$roundlam)
}
save.image(file = paste0(dir.name,"PolyEnv.RData"))


################### sparse multinomial regression
Fullmodel <- glmnet(x = as.matrix(BCdata[, -1]), y = BCdata[, 1], fam = "mult")
save(Fullmodel, file = paste0(dir.name,"SMRFullmodel.RData"))


cat("First k-Coeffs. \n")
for (k in c(5, 10, 15)) {
    print(k)
    for (j in 1:length(names(Fullmodel$beta))) {
        cf <- Fullmodel$beta[[j]]
        print(names(Fullmodel$beta)[j])
        numO <- num <- 0
        for (i in 1:length(Fullmodel$lambda)) {
            if (sum(abs(cf[, i])) > 0) {
                numO <- num
                num <- sum(cf[, i] != 0)
                
                myID <- which(abs(cf[, i]) != 0)
                if (num >= k) {
                  print(cf[myID[order(-abs(cf[myID, i]))], i])
                  break
                }
            }
        }
    }
    cat("\n\n")
}


for (i in 1:length(Fullmodel$classname)){
    pdf(paste0(dir.name,"Path_",Fullmodel$classname[i],".pdf"))
    plot_glmnet(Fullmodel, label = 5, nresponse = i, lwd = 2)
    dev.off()
}

