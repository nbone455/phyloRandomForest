# packages
library(ggplot2)
library(phytools)
library(geiger)
#

set.seed(7)

data <- list()
al <- list()
bl <- list()
cl <- list()
dl <- list()
el <- list()
fl <- list()
n <- 200 

tree <- geiger::sim.bdtree(n=n)
tree <- force.ultrametric(tree)


# phylo-correlated trait (P)
matrix <- matrix(byrow = T, nrow = 2, ncol = 2, data = c(NA,0.8,0.8,NA))
Q <- matrix * 0.1
qq <- Q 
qq[is.na(qq)] <- 0
diag(qq) <- rowSums(qq)*-1
datP <- sim.char(tree, par=qq, model = "discrete", root = 1)[,,1]

xx <- setNames(datP, (paint$phy$tip.label))
s <- make.simmap(paint$phy, xx, model = "ER")
cols <- setNames(c("darkgoldenrod1","green4"), c("1","2"))
plotSimmap(s, colors = cols )

# correleatd with clades maybe?

# lets say if state = 1, then there is a positive weight on Y, and a negtive weight for =2
datP_new <- case_when(
  datP == 1 ~ 0.45,
  datP == 2 ~ -0.45
  )




alpha <- 0.95

for (i in 1:n) {
  # six uncorrelated traits
al[[i]] <- sample(runif(100,max = 10, min = -5), 1)
bl[[i]] <- sample(rnorm(100,mean = 15, sd = 1.2), 1)
cl[[i]] <- sample(rnorm(100,mean= 10, sd = 3.5),1)
dl[[i]] <- sample(rnorm(100,mean= 3, sd = 1),1)
el[[i]] <- sample(runif(100, max = 10, min = 0),1)
fl[[i]] <- sample(rnorm(100, mean = 8.5, sd = 0.5),1)
  
  # this is making the classifer variable, Y 
  temp  <- (alpha*al[[i]] + alpha*bl[[i]] + alpha*cl[[i]] + alpha*dl[[i]] + alpha*el[[i]] + alpha*fl[[i]] +  alpha * datP_new[[i]])
  data[[i]] <- temp
}

#mergeCols <- c("unlist.data.","unlist.al.","unlist.b.","unlist.c.")
# uncorrelated traits 'A' - 'F'
a <- data.frame(unlist(al))
b <- data.frame(unlist(bl))
c <- data.frame(unlist(cl))
d <- data.frame(unlist(dl))
e <- data.frame(unlist(el))
f <- data.frame(unlist(fl))

# phylocorrelated trait 'P' 
datP <- data.frame(datP)
# response variable 'Y'
y <- data.frame(unlist(data))
full <- cbind(y,a,b,c,d,e,f,datP)
colnames(full) <- c("Y","A","B","C","D","E","F","P")
full_new <- full

#plot(full_new$Y)
# td object
td_full <- make.treedata(tree, full_new)
# figure out clades
# https://github.com/uyedaj/treeplyr/wiki
library(treeplyr)
paint <- paint_clades(td_full,nclades = 4,name = 'clade',interactive = TRUE)
paint$dat$Y <- factor(paint$dat$Y)

# need a trait that has a correlation with the value for 'clade' 
phylo_corr <- vector()

phylo_corr <- case_when(
  paint$dat$clade == 1 ~ 0,
  paint$dat$clade == 2 ~ 0.25,
  paint$dat$clade == 3 ~ 0.50,
  paint$dat$clade == 4 ~ 0.75,
  paint$dat$clade == 5 ~ 1,
)

paint$dat$phy_corr <- phylo_corr



Y_discrete <- vector()
Y_discrete <- case_when(
  Y_new  < 50 ~ 0,
 Y_new  > 50 ~ 1,
)

paint$dat$Y <- factor(Y_discrete)




# three tests 
## without the trait 
noCorrTrait_dat <- paint$dat
noCorrTrait_dat$P <- NULL
noCorrTrait_dat$phy_corr <- NULL
## without clade
noClade_dat <- paint$dat
noClade_dat$clade <- NULL

noRando_dat <- paint$dat 
noRando_dat$B <- NULL


one.error_data <- replicate(n = 4, expr = list())
zero.error_data <- replicate(n = 4, expr = list())
oob.error_data <- replicate(n = 4, expr = list())
error_noClade <- list()
error_noTrait <- list()
error_full <- list()
error_random <- list()

for (i in 1:500)  {
  library(randomForest)
  
    # Random forest what if a random uncorrelated trait was missing 
    inde <- sample(1:nrow(full_new), size = floor(nrow(full_new)/2))
    training_noR <- noRando_dat[inde,]
    training_noR <- na.omit(training_noR)
    rf_classifer_noR <- randomForest(Y ~., data = training_noR, ntree = 6000, mtry = 4, importance = TRUE)
    one.error_data[[4]][[i]] <- rf_classifer_noR$confusion[,3][2]
    zero.error_data[[4]][[i]] <- rf_classifer_noR$confusion[,3][1]
    oob.error_data[[4]][[i]] <- mean(rf_classifer_noR$err.rate[,1])
    
    # full RF model 
    pred_noR = predict(rf_classifer_noR, newdata = noRando_dat)
    cm_noR = table(noRando_dat$Y, pred_noR)
    #confusion matrix (this is the FP + FN / N)
    error_random[i] = (cm_noR[1,2] + cm_noR[2,1]) / (sum(cm_noR))
      
    # RANDOM FOREST (everything)
    inde <- sample(1:nrow(full_new), size = floor(nrow(full_new)/2))
    training <- paint$dat[inde,]
    training <- na.omit(training)
    # random forests threshold 
    library(randomForest)
    rf_classifer_full <- randomForest(Y ~., data = training, ntree = 10000, mtry = 4, importance = TRUE)
    #varImpPlot(rf_classifer_full, main = "everything")
    
    one.error_data[[1]][[i]] <- rf_classifer_full$confusion[,3][2]
    zero.error_data[[1]][[i]] <- rf_classifer_full$confusion[,3][1]
    oob.error_data[[1]][[i]] <- mean(rf_classifer_full$err.rate[,1])
   
    # full RF model 
    pred_full = predict(rf_classifer_full, newdata = paint$dat)
    cm_full = table(paint$dat$Y, pred_full)
    #confusion matrix (this is the FP + FN / N)
    error_full[i] = (cm_full[1,2] + cm_full[2,1]) / (sum(cm_full))
    
    
     # WITHOUT HAVING THE TRAIT 
    
    inde <- sample(1:nrow(full_new), size = floor(nrow(full_new)/2))
    training_noT <- noCorrTrait_dat[inde,]
    training_noT <- na.omit(training_noT)
    # random forests threshold 
    library(randomForest)
    rf_classifer_noT <- randomForest(Y ~., data = training_noT, ntree = 6000, mtry = 4, importance = TRUE)
    #varImpPlot(rf_classifer_noT, main = "no phylo Trait")
    
    one.error_data[[2]][[i]] <- rf_classifer_noT$confusion[,3][2]
    zero.error_data[[2]][[i]] <- rf_classifer_noT$confusion[,3][1]
    oob.error_data[[2]][[i]] <- mean(rf_classifer_noT$err.rate[,1])
    
    pred_noTrait = predict(rf_classifer_noT, newdata = noCorrTrait_dat)
    cm_noTrait = table(noCorrTrait_dat$Y, pred_noTrait)
    error_noTrait[i] = (cm_noTrait[1,2] + cm_noTrait[2,1]) / (sum(cm_noTrait))
    
    # WITHOUT CLADE
    inde <- sample(1:nrow(full_new), size = floor(nrow(full_new)/2))
    training_noC <- noClade_dat[inde,]
    training_noC <- na.omit(training_noC)
    # random forests threshold 
    library(randomForest)
    rf_classifer_noC <- randomForest(Y ~., data = training_noC, ntree = 6000, mtry = 4, importance = TRUE)
    #varImpPlot(rf_classifer_noC, main = "no phylo")
    
    pred_noClade = predict(rf_classifer_noC, newdata = noClade_dat)
    cm_noClade = table(noClade_dat$Y, pred_noClade)
    error_noClade[i] = (cm_noClade[1,2] + cm_noClade[2,1]) / (sum(cm_noClade))
    
    one.error_data[[3]][[i]] <- rf_classifer_noC$confusion[,3][2]
    zero.error_data[[3]][[i]] <- rf_classifer_noC$confusion[,3][1]
    oob.error_data[[3]][[i]] <- mean(rf_classifer_noC$err.rate[,1])
    
    }
    
    
#summarize each of these
no_rando <- (unlist(oob.error_data[[4]]))
no_clade <- (unlist(oob.error_data[[3]]))
no_trait <- (unlist(oob.error_data[[2]]))
all_err <- (unlist(oob.error_data[[1]]))

boxplot(no_rando, no_clade, no_trait, all_err, names = c("No Random Trait", "No Clade Info", "No Tree-Correlated", "All"),  main = "OOB ERROR", ylab = 
          "OOB error rate")


no_clade_test <- (unlist(error_noClade))
no_trait_test <- (unlist(error_noTrait))
no_rando_test <- (unlist(error_random))
full_test     <- (unlist(error_full))

boxplot(no_rando_test, no_clade_test, no_trait_test, full_test, names = c("No Random Trait", "No Clade Info", "No Tree-Correlated", "All"),  main = "Overall ERROR", ylab = 
          "Overall error rate")
