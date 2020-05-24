# On-imputation-methods
# Estimating a regression analysis from a data matrix with missing data (MCAR); using the MICE package
# Estimating one true model and with Complete Case Analysis, Mean Imputation, Stochastic Regression Imputation, Random Forest Imputation, Predictive Mean Matching

###library 
library(mice)
library(ggplot2)
library(tidyverse)
library(xtable)

#################################
## in this section we control all important variables/parameters of the simulation study 
## sample size
iN = 1000

## no. of independent variables
ik = 1

## no. of repititions in the Monte Carlo
iM = 100

## no. of regression models we will estimate
## will set the size of arrays to store results 
rm = 6 

## standard deviation of the errors
sigma =  1

## the proportion of missing values in the experiement
pp = 0.1

## set min and max for IVs
xMax <- 20
xMin <- 0

#############################
## here we set all the parameters of the model

## parameters 
beta = rep(1, ik)

#############################
## here we take samples of draws from some distributions
## and control the seed
set.seed(2020)

#error = rnorm(iN, mean=0, sd=1)
m_error = matrix(rnorm(iN*iM, sd=sigma), iN, iM)

#matrix of independent variables simulated 
#simulated from uniform distribution
mX = matrix(runif(iN*ik, min=xMin, max=xMax), iN, ik) 

pos_missing = matrix(sample(c(TRUE, FALSE),
                            iN*(ik+1), replace=T, 
                            prob=c(pp, 1-pp)), iN, ik+1)

#### from here we will not do any sampling. 

#############################
## Monte Carlo and the experiment

## we need some container to hold the output
## we initialize the results matrix to be all zeros!
### matrix to store coefficent estimates
betas = array(0, dim=c(ik, rm, iM)) 

### matrix to store coefficent standard errors
SE = array(0, dim=c(ik, rm, iM)) 
###matrix to store 95% CIs 
CI = array(0, dim=c(rm, 2*ik, iM)) 

## the body of the experiment
for(iter in 1:iM){
  
  ##generate Y from beta and the independent variable
  vy = mX %*% beta + m_error[,iter]
  
  ## here we make the missing observations!!!
  ## make a big matrix of the data
  ##containing the dependent variables and the independent ones
  m_data = cbind(vy, mX)
  ##make a second big matrix of the data 
  ##containing the dependent and independent variables
  ##for use without creating any missing observations
  full_data <- cbind(vy,mX)
  
  ## now we make the misssing obs.
  m_data[pos_missing] = NA
  model_data <- as.data.frame(m_data)
  
  
  ### Giving the variables more intuitive names
  ### Y and X1, X2.... X:iK
  colnames(m_data) <- c("Y", paste0("X", 1:ik))
  colnames(full_data) <- c("Y", paste0("X", 1:ik))
  
  ### set some data frames for each dataset 
  #for the "full model" analysis
  original_data <- as.data.frame(full_data) 
  #for the complete case analysis
  CCA_data <- as.data.frame(m_data) 
  
  ### m = imputations no. imputations in the MI
  ### imputing mean of each variable with missing data
  imp_mean <- mice(m_data, method = "mean", m = 3, seed = 1)

  ### stochastic regression imputation
  imp_sto <- mice(m_data, method = "norm.nob", m = 3, seed = 1)
  
  ### random forest imputation
  imp_rf <- mice(m_data, method = "rf", m = 3, seed = 1) 

  ### Predictive Mean Matching
  imp_PMM <- mice(m_data, m = 3, method = "pmm",seed = 1) 
 
  
  ### run regression models with the dependent variable Y and iK 
  model_fulldata <- lm(data=original_data, 
                       formula(paste0("Y ~ ",
                                str_c(paste0("X",1:ik), 
                                     collapse=" + "))))
  model_CCA <- lm(data=CCA_data, 
                  formula(paste0("Y ~ ",
                                 str_c(paste0("X",1:ik), 
                                       collapse=" + "))))
  ####
  mean <- with(data = imp_mean,
               expr= lm( formula(paste0("Y ~ ",
                                        str_c(paste0("X",1:ik), 
                                              collapse=" + ")))))
  stochastic <- with(data = imp_sto,
                     expr = lm(formula(paste0("Y ~ ",
                                       str_c(paste0("X",1:ik), 
                                            collapse=" + ")))))
  randomforest <- with(data = imp_rf,
                       expr = lm(formula(paste0("Y ~ ",
                                        str_c(paste0("X",1:ik), 
                                              collapse=" + ")))))
  PMM <- with(data = imp_PMM,
              expr = lm(formula(paste0("Y ~ ",
                                str_c(paste0("X",1:ik),
                                      collapse=" + ")))))
  
  ##pooling results from the eastimated regression models of the imputed data
  model_mean <- summary(pool(mean), conf.int = TRUE)
  model_stochastic <- summary(pool(stochastic), conf.int = TRUE)
  model_rf <- summary(pool(randomforest), conf.int = TRUE)
  model_PMM <- summary(pool(PMM), conf.int = TRUE)
  
  ##If running the code with iK > 1 , change the code to also save these in the result containers
  ##Storing Betas for each  model monte carlo
  betas[,1,iter] = coef(model_fulldata)[-1]
  betas[,2,iter] = coef(model_CCA)[-1]
  betas[,3,iter] = model_mean[2,2]
  betas[,4,iter] = model_stochastic[2,2]
  betas[,5,iter] = model_rf[2,2]
  betas[,6,iter] = model_PMM[2,2]

  ##Storing standards errors 
  SE[,1,iter] = coef(summary(model_fulldata))[2, "Std. Error"]
  SE[,2,iter] = coef(summary(model_CCA))[2, "Std. Error"]
  SE[,3,iter] = model_mean[2,3]
  SE[,4,iter] = model_stochastic[2,3]
  SE[,5,iter] = model_rf[2,3]
  SE[,6,iter] = model_PMM[2,3]
  model_mean
  ##Storing Confidence intervals for betas
  CI[1,,iter] = confint(model_fulldata)[-1,] 
  CI[2,,iter] = confint(model_CCA)[-1,]
  CI[3,,iter] = c(model_mean[2,7], model_mean[2,8])
  CI[4,,iter] = c(model_stochastic[2,7], model_stochastic[2,8])
  CI[5,,iter] = c(model_rf[2,7], model_rf[2,8])
  CI[6,,iter] = c(model_PMM[2,7], model_PMM[2,8])
}

##From here we will analyse and use code to present the results 
##Extracting the first imputed dataset of each MI to plot
mean_data <- complete(imp_mean)
stochastic_data <- complete(imp_sto)
rf_data <- complete(imp_rf)
PMM_data <- complete(imp_PMM)

##Raw Bias, Percentage Bias and RMSE for Betas
RB <- rowMeans(betas[,,] - beta)
PB <- 100 * abs((rowMeans(betas[,,]) - beta) / beta)
RMSE <- sqrt(rowMeans((betas)[,,] - beta)^2)

##Coverage Rate and Average Width for Betas
CR <- rowMeans(CI[,1,] < beta & beta < CI[,2,])
AW <- rowMeans(CI[,2,] - CI[,1,])

##Extracting pooled R squared from each MI case
R_1 <- (pool.r.squared(mean))
R_2 <- (pool.r.squared(stochastic))
R_3 <- (pool.r.squared(randomforest))
R_4 <- (pool.r.squared(PMM))

### Pooling betas and standard errors from the Monte Carlo
SE_pooled <- rowMeans(SE[,,])
betas_pooled <- rowMeans(betas[,,])

##Combining R squared in one array
CD <- array(cbind(summary(model_fulldata)$r.squared,
                  summary(model_CCA)$r.squared, R_1[1,1],
                  R_2[1,1],
                  R_3[1,1],
                  R_4[1,1]))  

##All in one dataframe to put in table 
StatsTable <-  as.data.frame(cbind(betas_pooled, SE_pooled, RB, PB, RMSE, CR, AW, CD))

#setting column  and rownames for the table
#use as numeric to 
colnames(StatsTable) <- c("Estimates","Std.Error", "RB", "PB", "RMSE", "CR", "AW", "R.squared") 
rownames(StatsTable) <- c("Full","CCA","Mean", "Sto", "RF", "PMM")

## creating latex file to present the results in a table 
print(xtable(StatsTable, digits=3, type = "latex",), file = "file2.tex")


######################################################
## Here we make some examples how to plot the results
######################################################

#Plots of confidence intervals of beta in the Monte Carlo
### set plot window 
par(mfrow=c(2,3))
#### Full Model
ID <- which(!(CI[1,1,] <= 1 & 1 <= CI[1,2,]))
#initializing the plot
plot(0, 
     xlim = c(.96, 1.04), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(beta), 
     main = "Full")
# colors vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"

#reference for true beta
abline(v = 1, lty = 2)

for(j in 1:iM){
  lines(c(CI[1,1,j], CI[1,2,j]), 
        c(j,j), 
        col = colors[j], 
        lwd = 2)
}

#### Comple Caste
ID <- which(!(CI[2,1,] <= 1 & 1 <= CI[2,2,]))
#initializing the plot
plot(0, 
     xlim = c(.93, 1.07), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(beta), 
     main = "CCA")
# colors vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"

# reference for true beta
abline(v = 1, lty = 2)

for(j in 1:iM){
  lines(c(CI[2,1,j], CI[2,2,j]), 
        c(j,j), 
        col = colors[j], 
        lwd = 2)
}

#### Mean Model
ID <- which(!(CI[3,1,] <= 1 & 1 <= CI[3,2,]))
#initializing the plot
plot(0, 
     xlim = c(.8, 1.05), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(beta), 
     main = "Mean")
# colors vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"
# reference for beta = 1
# draw reference line at mu=5
abline(v = 1, lty = 2)

for(j in 1:iM){
  lines(c(CI[3,1,j], CI[3,2,j]), 
        c(j,j), 
        col = colors[j], 
        lwd = 2)
}


#### Stochastic Regression Imp 
ID <- which(!(CI[4,1,] <= 1 & 1 <= CI[4,2,]))
#initializing the plot
plot(0, 
     xlim = c(.95, 1.05), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(beta), 
     main = "Stochastic")
# colors vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"
# reference for beta = 1
# draw reference line at mu=5
abline(v = 1, lty = 2)

for(j in 1:iM){
  lines(c(CI[4,1,j], CI[4,2,j]), 
        c(j,j), 
        col = colors[j], 
        lwd = 2)
}

#### Determenistic Regression Imp 
ID <- which(!(CI[5,1,] <= 1 & 1 <= CI[5,2,]))
#initializing the plot
plot(0, 
     xlim = c(.95, 1.05), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(beta), 
     main = "Random Forest")
# colors vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"
# reference for beta = 1
# draw reference line at mu=5
abline(v = 1, lty = 2)

for(j in 1:iM){
  lines(c(CI[5,1,j], CI[5,2,j]), 
        c(j,j), 
        col = colors[j], 
        lwd = 2)
}

#### PMM 
ID <- which(!(CI[6,1,] <= 1 & 1 <= CI[6,2,]))
#initializing the plot
plot(0, 
     xlim = c(.95, 1.05), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(beta), 
     main = "PMM")
# colors vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"
# reference for beta = 1
# draw reference line at mu=5
abline(v = 1, lty = 2)

for(j in 1:iM){
  lines(c(CI[6,1,j], CI[6,2,j]), 
        c(j,j), 
        col = colors[j], 
        lwd = 2)
}

#########################################
#Plots of imputed values vs observed values, first imputed dataset in the Multiple Imputation

## set plot window
par(mfrow=c(2,2))
# mean imputation, observed vs imputed values

plot(original_data$X1[pos_missing[,2]],  #observed values
     xlim = c(0, iN*pp), ylim = c(xMin-5, xMax+5),
     main = "Mean",
     xlab = "X", ylab = "Y")
points((mean_data$X1[pos_missing[,2]]), #imputed values
       col = "red")
legend("topleft",   # Legend
       c("Observed Values", "Imputed Values"),
       pch = c(1, 1, NA),
       lty = c(NA, NA, 1),
       cex = .5,
       col = c("black", "red"))

#stochastic regression, observed vs imputed values
plot(original_data$X1[pos_missing[,2]], #observed values
     xlim = c(0, iN*pp), ylim = c(xMin-5, xMax+5),
     main = "Stochastic",
     xlab = "X", ylab = "Y")

points((stochastic_data$X1[pos_missing[,2]]), #imputed values
       col = "red")
legend("topleft",  # Legend
       c("Observed Values", "Imputed Values"),
       pch = c(1, 1, NA),
       lty = c(NA, NA, 1),
       cex = .5,
       col = c("black", "red"))

#Random Forest, observed vs imputed values 
plot(original_data$X1[pos_missing[,2]], #observed values
     xlim = c(0, iN*pp), ylim = c(xMin-5, xMax+5),
     main = "Random Forest",
     xlab = "X", ylab = "Y")

points((rf_data$X1[pos_missing[,2]]), #imputed values
       col = "red")
legend("topleft",   # Legend
       c("Observed Values", "Imputed Values"),
       pch = c(1, 1, NA),
       lty = c(NA, NA, 1),
       cex = .5,
       col = c("black", "red"))

#PMM, observed vs imputed values
plot(original_data$X1[pos_missing[,2]], #observed values
     xlim = c(0, iN*pp), ylim = c(xMin-5, xMax+5),
     main = "PMM",
     xlab = "X", ylab = "Y")
points((PMM_data$X1[pos_missing[,2]]),#imputed values
       col = "red")
legend("topleft",  # Legend
       c("Observed Values", "Imputed Values"),
       pch = c(1, 1, NA),
       lty = c(NA, NA, 1),
       cex = .5,
       col = c("black", "red"))
