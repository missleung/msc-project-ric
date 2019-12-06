

#When we do simulations, biomarker values are assumed to be non-negative for the definition to work!!

########################
####### RIC_exp function
########################

library(DescTools)
library(data.table)

#RIC_exp is a function used to plot the RIC curve on simulated dataset

RIC_exp <- function(proportion, biomarker, benefit){ 
  
  #Grabbing biomarker threshold such that it matches the proportion, in this case, our biomarker is y_no_trt
  k <- quantile(biomarker, probs = 1 - proportion)

  #Creating treatment vectors such that patients are assigned to treatment if biomarker is greater than threshold
  new_trt <- sapply(k, function (x) as.numeric(biomarker >= x))
  
  #relative benefit
  b1 <- apply(new_trt, 2, function(x) sum((benefit)*x)/
               sum(benefit))
  
  #return relative benefit
  b1
}

################################################
###### Simulation for Multiple Linear Regression
################################################

#Simulation on one sample data for biomarker function data with n subjects and two risk factors.

n <- 50
int <- 0
b1 <- 1
b2 <- 0.5
b3 <- -0.25
b4 <- -0.5
b5 <- -0.25

#INITIATE a vector of aucs
auc_vec_lm_1<-c()
auc_vec_lm_2<-c()
auc_vec_lm_3<-c()
auc_vec_lm_4<-c()

for (i in 1: 100){

set.seed(i)
x1 <- runif(n, 0, 10)
x2 <- rnorm(n, 10, sd=6)
e <- rnorm(n,sd=5)

#Treatment was randomly assigned 
trt <- rbinom(n, size=1, prob=0.5)

#if biomarker and benefit is calculated  with errors
y_err <- int + b1*x1 + b2*x2 + b3*trt + I(b4*trt*x1) + I(b4*trt*x2) + e

#Make a dataframe out of the whole data
mat_err = data.frame(x1=x1, x2=x2, trt=trt, y=y_err)

##Prediction model
lm_err <- lm(y ~ x1 + x2 + trt + I(trt*x1) + I(trt*x2), data=mat_err)

###NOTE: Marker values were calculated as the predicted number of exacerbations without treatment by setting, 
##for all individuals, the treatment variable to 0. (Sadatsafavi, In Review)

#Creating data set that has all trt variable as 0
no_trt_data <- as.data.frame(cbind(mat_err[,!names(mat_err) %in% "trt"], trt = 0))

#Creating data set that has all trt variable as 1
all_trt_data <- as.data.frame(cbind(mat_err[,!names(mat_err) %in% "trt"], trt = 1))

#Biomarker is the predicted y values when treatment is 0
y_hat_no_trt <- predict(lm_err, newdata = no_trt_data)

#Predict outcome when all trt variable is 1
y_hat_all_trt <- predict(lm_err, newdata = all_trt_data)

#Predict Benefit
benefit_predict <- y_hat_no_trt - y_hat_all_trt

#True y when treatment is 0 (true biomarker)
y_true_no_trt <- int + b1*x1 + b2*x2 
y_true_all_trt <- int + b1*x1 + b2*x2 + b3 + I(b4*x1) + I(b4*x2)
benefit_true <- y_true_no_trt - y_true_all_trt

#Case 1: both biomarker and benefits are predicted

#CALCULATING AUC
ric_x <- seq(0, 1, by = 0.01)
ric_y_1 <- RIC_exp(proportion=ric_x, biomarker = y_hat_no_trt, benefit = benefit_predict)
auc_vec_lm_1 <- append(auc_vec_lm_1, AUC(ric_x,ric_y_1))

#Case 2: predicted biomarker and true benefit

#CALCULATING AUC
ric_y_2 <- RIC_exp(proportion=ric_x, biomarker = y_hat_no_trt, benefit = benefit_true)
auc_vec_lm_2 <- append(auc_vec_lm_2, AUC(ric_x,ric_y_2))

#Case 3: true biomarker and predicted benefit

#CALCULATING AUC
ric_y_3 <- RIC_exp(proportion=ric_x, biomarker = y_true_no_trt, benefit = benefit_predict)
auc_vec_lm_3 <- append(auc_vec_lm_3, AUC(ric_x,ric_y_3))

#Case 4: both biomarker and benefits is true

#CALCULATING AUC
ric_y_4 <- RIC_exp(proportion=ric_x, biomarker = y_true_no_trt, benefit = benefit_true)
auc_vec_lm_4 <- append(auc_vec_lm_4, AUC(ric_x,ric_y_4))

cat("Done:", i,"\n")

}

#Summary of AUC values on simulated data
cat("Summary statistics on AUC values for Case 1 data \n",
    "Mean AUC: ", mean(auc_vec_lm_1), "\n",
    "Median AUC: ", median(auc_vec_lm_1), "\n",
    "SD on AUC: ",sd(auc_vec_lm_1), "\n",
    "50% range: ", quantile(auc_vec_lm_1,0.25), " to ",quantile(auc_vec_lm_1,0.75),"\n",
    "Range on AUC", min(auc_vec_lm_1), " to ", max(auc_vec_lm_1))

jpeg("Rplot_LM_1.jpeg")
plot(density(auc_vec_lm_1), main = "Density Plot of AUCs in Case 1", xlab="AUC values", xlim=c(0.4,0.9), ylim=c(0,25))
dev.off()

cat("Summary statistics on AUC values for Case 2 data \n",
    "Mean AUC: ", mean(auc_vec_lm_2), "\n",
    "Median AUC: ", median(auc_vec_lm_2), "\n",
    "SD on AUC: ",sd(auc_vec_lm_2),
    "50% range: ", quantile(auc_vec_lm_2,0.25), " to ",quantile(auc_vec_lm_2,0.75),"\n",
    "Range on AUC", min(auc_vec_lm_2), " to ", max(auc_vec_lm_2))

jpeg("Rplot_LM_2.jpeg")
plot(density(auc_vec_lm_2), main = "Density Plot of AUCs in Case 2", xlab="AUC values", xlim=c(0.4,0.9), ylim=c(0,25))
dev.off()

cat("Summary statistics on AUC values for Case 3 data \n",
    "Mean AUC: ", mean(auc_vec_lm_3), "\n",
    "Median AUC: ", median(auc_vec_lm_3), "\n",
    "SD on AUC: ",sd(auc_vec_lm_3),
    "50% range: ", quantile(auc_vec_lm_3,0.25), " to ",quantile(auc_vec_lm_3,0.75),"\n",
    "Range on AUC", min(auc_vec_lm_3), " to ", max(auc_vec_lm_3))

jpeg("Rplot_LM_3.jpeg")
plot(density(auc_vec_lm_3), main = "Density Plot of AUCs in Case 3", xlab="AUC values", xlim=c(0.4,0.9), ylim=c(0,25))
dev.off()

cat("Summary statistics on AUC values for Case 4 data \n",
    "Mean AUC: ", mean(auc_vec_lm_4), "\n",
    "Median AUC: ", median(auc_vec_lm_4), "\n",
    "SD on AUC: ",sd(auc_vec_lm_4),
    "50% range: ", quantile(auc_vec_lm_4,0.25), " to ",quantile(auc_vec_lm_4,0.75),"\n",
    "Range on AUC", min(auc_vec_lm_4), " to ", max(auc_vec_lm_4))

jpeg("Rplot_LM_4.jpeg")
plot(density(auc_vec_lm_4), main = "Density Plot of AUCs in Case 4", xlab="AUC values", xlim=c(0.4,0.9), ylim=c(0,25))
dev.off()
################################################
###### Simulation for Poisson Regression
################################################
#Simulation on one sample data for biomarker function data with n subjects and two risk factors.

n <- 20
int <- 0.2
b1 <- 0.1
b2 <- 0.25
b3 <- -0.25
b4 <- -0.05
b5 <- -0.125


#INITIATE a vector of aucs
auc_vec_glm_1<-c()
auc_vec_glm_2<-c()
auc_vec_glm_3<-c()
auc_vec_glm_4<-c()

for (i in 1: 100){
  
  set.seed(i)
  x1 <- runif(n, 0, 10)
  x2 <- rnorm(n, 10, 6)

  #Treatment was randomly assigned
  trt <- rbinom(n, size=1, prob=0.5)

  #Biomarker and benefit are with errors
  y_err <- rpois(n=n, lambda = exp(int + b1*x1 + b2*x2 + b3*trt + I(b4*trt*x1) + I(b4*trt*x2)))
  
  #Make a dataframe out of the whole data
  mat_err = data.frame(x1=x1, x2=x2, trt=trt, y=y_err)
  
  #Prediction model
  glm_err <- glm(y ~ -1 + x1 + x2 + I(trt*x1) , family = "poisson", data=mat_err)
  
  #Creating data set that has all trt variable as 0
  no_trt_data <- as.data.frame(cbind(mat_err[,!names(mat_err) %in% "trt"], trt = 0))
  
  #Creating data set that has all trt variable as 1
  all_trt_data <- as.data.frame(cbind(mat_err[,!names(mat_err) %in% "trt"], trt = 1))
  
  #Biomarker is the predicted y values when treatment is 0
  y_hat_no_trt <- predict(glm_err, newdata = no_trt_data, type = "response")
  
  #Predict outcome when all trt variable is 1
  y_hat_all_trt <- predict(glm_err, newdata = all_trt_data, type = "response")
  
  #Predict benefit 
  benefit_predict <- y_hat_no_trt - y_hat_all_trt
  
  #True benefit  
  y_true_no_trt <- exp(b1*x1 + b2*x2)
  y_true_all_trt <- exp(b1*x1 + b2*x2 + b3 + b4*x1)
  benefit_true <- y_true_no_trt - y_true_all_trt 
  
  #Case 1: both biomarker and benefits are predicted
  
  #CALCULATING AUC
  ric_x <- seq(0, 1, by = 0.01)
  ric_y_1 <- RIC_exp(proportion=ric_x, biomarker = y_hat_no_trt, benefit = benefit_predict)
  auc_vec_glm_1 <- append(auc_vec_glm_1, AUC(ric_x,ric_y_1))
  
  #Case 2: predicted biomarker and true benefit
  
  ric_y_2 <- RIC_exp(proportion=ric_x, biomarker = y_hat_no_trt, benefit = benefit_true)
  auc_vec_glm_2 <- append(auc_vec_glm_2, AUC(ric_x,ric_y_2))
  
  #Case 3: true biomarker and predicted benefit
  
  #CALCULATING AUC
  ric_y_3 <- RIC_exp(proportion=ric_x, biomarker = y_true_no_trt, benefit = benefit_predict)
  auc_vec_glm_3 <- append(auc_vec_glm_3, AUC(ric_x,ric_y_3))
  
  cat("Done:", i,"\n")
    
  #Case 4: both biomarker and benefits is true
  
  #CALCULATING AUC
  ric_y_4 <- RIC_exp(proportion=ric_x, biomarker = y_true_no_trt, benefit = benefit_true)
  auc_vec_glm_4 <- append(auc_vec_glm_4, AUC(ric_x,ric_y_4))
  
}

par(mfrow=c(2,2)) 

#Summary of AUC values on simulated data
cat("Summary statistics on AUC values for Case 1 data \n",
    "Mean AUC: ", mean(auc_vec_glm_1), "\n",
    "Median AUC: ", median(auc_vec_glm_1), "\n",
    "SD on AUC: ",sd(auc_vec_glm_1), "\n",
    "50% range: ", quantile(auc_vec_glm_1,0.25), " to ",quantile(auc_vec_glm_1,0.75),"\n",
    "Range on AUC", min(auc_vec_glm_1), " to ", max(auc_vec_glm_1))


jpeg("poisson_case1.jpg")
plot(density(auc_vec_glm_1), main = "Density Plot of AUCs in Case 1", xlab="AUC values", xlim=c(0.5,1), ylim=c(0,10))
dev.off()

cat("Summary statistics on AUC values for Case 2 data \n",
    "Mean AUC: ", mean(auc_vec_glm_2), "\n",
    "Median AUC: ", median(auc_vec_glm_2), "\n",
    "SD on AUC: ",sd(auc_vec_glm_2), "\n",
    "50% range: ", quantile(auc_vec_glm_2,0.25), " to ",quantile(auc_vec_glm_2,0.75),"\n",
    "Range on AUC", min(auc_vec_glm_2), " to ", max(auc_vec_glm_2))

jpeg("poisson_case2.jpg")
plot(density(auc_vec_glm_2), main = "Density Plot of AUCs in Case 2", xlab="AUC values", xlim=c(0.5,1), ylim=c(0,10))
dev.off()

cat("Summary statistics on AUC values for Case 3 data \n",
    "Mean AUC: ", mean(auc_vec_glm_3), "\n",
    "Median AUC: ", median(auc_vec_glm_3), "\n",
    "SD on AUC: ",sd(auc_vec_glm_3), "\n",
    "50% range: ", quantile(auc_vec_glm_3,0.25), " to ",quantile(auc_vec_glm_3,0.75),"\n",
    "Range on AUC", min(auc_vec_glm_3), " to ", max(auc_vec_glm_3))

jpeg("poisson_case3.jpg")
plot(density(auc_vec_glm_3), main = "Density Plot of AUCs in Case 3", xlab="AUC values",  xlim=c(0.5,1), ylim=c(0,10))
dev.off()

cat("Summary statistics on AUC values for Case 4 data \n",
    "Mean AUC: ", mean(auc_vec_glm_4), "\n",
    "Median AUC: ", median(auc_vec_glm_4), "\n",
    "SD on AUC: ",sd(auc_vec_glm_4), "\n",
    "50% range: ", quantile(auc_vec_glm_4,0.25), " to ",quantile(auc_vec_glm_4,0.75),"\n",
    "Range on AUC", min(auc_vec_glm_4), " to ", max(auc_vec_glm_4))

jpeg("poisson_case4.jpg")
plot(density(auc_vec_glm_4), main = "Density Plot of AUCs in Case 4", xlab="AUC values", xlim=c(0.5,1), ylim=c(0,10))
dev.off()


