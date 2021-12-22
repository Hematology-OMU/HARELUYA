################################################################################################################ 
## Describe Calibration plots and ROC curves in RSF and DeepSurv models
## This code requires the CSV file which include the predictive probabilities of 1-year OS after allo-HSCT in the test cohort using Deepsurv with python
################################################################################################################ 

library(gtools)
library(timeROC)
library(randomForestSRC)
library(tidyverse)
library(survival)

createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}
FILEPATH = './data/analysed_df.csv'
Analysed_df_imp = read.csv(FILEPATH, fileEncoding = 'cp932', stringsAsFactors=TRUE)
Analysed_df_imp <- select(Analysed_df_imp, -c(.Sex.Mismatch3,dx_to_sct_day_trump))
Analysed_df_imp <- Analysed_df_imp %>% mutate(.Event = .OS)
Analysed_df_imp <- Analysed_df_imp %>% mutate(.MonthEvent = .MonthOS)

Analysed_df_imp$.HCT.CI <- as.factor(Analysed_df_imp$.HCT.CI)
Analysed_df_imp$.PS24 <- as.factor(Analysed_df_imp$.PS24)
Analysed_df_imp$.Analysed_Disease <- as.factor(Analysed_df_imp$.Analysed_Disease)
Analysed_df_imp$.Disease.satatus.all <- as.factor(Analysed_df_imp$.Disease.satatus.all)
Analysed_df_imp$Tx_pattern <- as.factor(Analysed_df_imp$Tx_pattern)
Analysed_df_imp$R_CMVAntiB <- as.factor(Analysed_df_imp$R_CMVAntiB)


karafuto <- 1808

train_df <- Analysed_df_imp %>% filter(.SCT.Year <= 2016)
test_df <- Analysed_df_imp %>% filter(.SCT.Year >= 2017)
Outcome_Month <- 12

# Set hyper-parameters of RSF
mtry = 2
nodesize = 80
nodedepth = 59
ntree = 500

# Develop predictive RSF model
modelRFSRC <- rfsrc(Surv(.MonthEvent, .Event) ~ 
                      .Age
                    +.HCT.CI
                    +.PS24
                    +.Analysed_Disease
                    +.Disease.satatus.all
                    +R_CMVAntiB
                    +Tx_pattern,
                    data = train_df, seed = -karafuto, ntree = ntree, 
                    mtry =mtry, nodesize = nodesize, nodedepth = nodedepth)

# Calculate the predictive probabilities of 1-year OS after allo-HSCT in the test cohort using RSF
predict_outcome_rfsrc <- function(x,fu_month){
  y <- which( predict(modelRFSRC, x)$time.interest <= fu_month)
  return (predict(modelRFSRC, x)$survival[,max(y)])
}
rfsrc_predict = createEmptyDf(nrow(test_df), 1, colnames = c( "study_ID" ) )
for( i in 1:nrow(test_df) ){
  rfsrc_predict[ i, 1 ] = as.character(test_df[i,"study_ID"])
}
temp_p <- predict_outcome_rfsrc(test_df,Outcome_Month)
rfsrc_predict <- cbind(rfsrc_predict,temp_p)
colnames(rfsrc_predict)[2] <- c("rfsrc_predict_OS")
test_df <- merge(test_df, rfsrc_predict, by="study_ID", all=T)
rm(rfsrc_predict)

# calculate timeROC for the predictive probability of 1-year OS after allo-HSCT in test cohort using RSF
test_rfsrc <- test_df %>% select(c(.Event,.MonthEvent,rfsrc_predict_OS))

rfsrc_tROC <-timeROC(T=test_rfsrc$.MonthEvent,
                         delta=test_rfsrc$.Event,marker=(1-test_rfsrc$rfsrc_predict_OS),
                         cause=1,weighting="marginal",
                         times=c(Outcome_Month),
                         iid=TRUE)
print(confint(rfsrc_tROC))

# Cut the predictive probabilities of 1-year OS after allo-HSCT in test cohort using RSF into 10 bins
devide_n <- 10
quantiles_risk <- quantcut(test_rfsrc$rfsrc_predict_OS, q=devide_n) 
test_rfsrc$quantiles_risk <- quantiles_risk 
for (i in 1:devide_n){
  assign(paste0("KM_estimates_",i),survfit(Surv(.MonthEvent, .Event)~1, data=test_rfsrc[test_rfsrc$quantiles_risk==levels(test_rfsrc$quantiles_risk)[i],]))
}
mean_survival_risk_rsf <-tapply(test_rfsrc$rfsrc_predict_OS, test_rfsrc$quantiles_risk, mean) 

estimate_KMvalue <- function(x,fu_month){
  y <- which(x$time <= fu_month)
  return (x$surv[max(y)])
}

# Kaplan - Meier estimates (observed survival probability) for each bin
KM_est_rsf <- c(estimate_KMvalue(KM_estimates_1,Outcome_Month), 
            estimate_KMvalue(KM_estimates_2,Outcome_Month), 
            estimate_KMvalue(KM_estimates_3,Outcome_Month), 
            estimate_KMvalue(KM_estimates_4,Outcome_Month), 
            estimate_KMvalue(KM_estimates_5,Outcome_Month), 
            estimate_KMvalue(KM_estimates_6,Outcome_Month), 
            estimate_KMvalue(KM_estimates_7,Outcome_Month), 
            estimate_KMvalue(KM_estimates_8,Outcome_Month), 
            estimate_KMvalue(KM_estimates_9,Outcome_Month),
            estimate_KMvalue(KM_estimates_10,Outcome_Month))

KM_est_event <- 1 - KM_est_rsf ## observed event risk
mean_event_risk <- 1 - mean_survival_risk_rsf ## predicted event risk

rsf_x <- data.frame(
  predicted_survival = mean_survival_risk_rsf,
  KM_survival  = KM_est_rsf
)

# Read csv which include the predictive probabilities of 1-year OS after allo-HSCT in the test cohort using Deepsurv with python
FILEPATH = './data/deep_s_12prediction.csv'
DeepS_12predct = read.csv(FILEPATH, fileEncoding = 'cp932')
DeepS_12predct$.HCT.CI <- as.factor(DeepS_12predct$.HCT.CI)
DeepS_12predct$.PS24 <- as.factor(DeepS_12predct$.PS24)
DeepS_12predct$.Analysed_Disease <- as.factor(DeepS_12predct$.Analysed_Disease)
DeepS_12predct$.Disease.satatus.all <- as.factor(DeepS_12predct$.Disease.satatus.all)
DeepS_12predct$Tx_pattern <- as.factor(DeepS_12predct$Tx_pattern)
DeepS_12predct$R_CMVAntiB <- as.factor(DeepS_12predct$R_CMVAntiB)

test_DeepS <- DeepS_12predct %>% select(c(.OS,.MonthOS,deepS_p_predict))

# calculate timeROC for the predictive probability of 1-year OS after allo-HSCT in test cohort using DeepSurv
DeepS_tROC <-timeROC(T=test_DeepS$.MonthOS,
                     delta=test_DeepS$.OS,marker=(1-test_DeepS$deepS_p_predict),
                     cause=1,weighting="marginal",
                     times=c(Outcome_Month),
                     iid=TRUE)
print(confint(DeepS_tROC))

# Cut the predictive probabilities of 1-year OS after allo-HSCT in test cohort using DeepSurv into 10 bins
devide_n <- 10
quantiles_risk <- quantcut(test_DeepS$deepS_p_predict, q=devide_n) 
test_DeepS$quantiles_risk <- quantiles_risk 
for (i in 1:devide_n){
  assign(paste0("KM_estimates_",i),survfit(Surv(.MonthOS, .OS)~1, data=test_DeepS[test_DeepS$quantiles_risk==levels(test_DeepS$quantiles_risk)[i],]))
}
mean_survival_risk_deepS <-tapply(test_DeepS$deepS_p_predict, test_DeepS$quantiles_risk, mean) 
estimate_KMvalue <- function(x,fu_month){
  y <- which(x$time <= fu_month)
  return (x$surv[max(y)])
}

# Kaplan - Meier estimates (observed survival probability) for each bin
KM_est_deepS <- c(estimate_KMvalue(KM_estimates_1,Outcome_Month), 
            estimate_KMvalue(KM_estimates_2,Outcome_Month), 
            estimate_KMvalue(KM_estimates_3,Outcome_Month), 
            estimate_KMvalue(KM_estimates_4,Outcome_Month), 
            estimate_KMvalue(KM_estimates_5,Outcome_Month), 
            estimate_KMvalue(KM_estimates_6,Outcome_Month), 
            estimate_KMvalue(KM_estimates_7,Outcome_Month), 
            estimate_KMvalue(KM_estimates_8,Outcome_Month), 
            estimate_KMvalue(KM_estimates_9,Outcome_Month),
            estimate_KMvalue(KM_estimates_10,Outcome_Month))

KM_est_event <- 1- KM_est_deepS ## observed event risk
mean_event_risk <- 1- mean_survival_risk_deepS ## predicted event risk

deepS_x <- data.frame(
  predicted_survival = mean_survival_risk_deepS,
  KM_survival  = KM_est_deepS
)


# Describe ROC curves of RSF and DeepSurv
plot(rfsrc_tROC,time=Outcome_Month,col="goldenrod",lwd=2,title=FALSE)
plot(DeepS_tROC,time=Outcome_Month,col="navy",add=TRUE,lwd=2,lty=2)
legend("bottomright", legend = c("DeepSurv","RSF"), col = c("navy","goldenrod"), lty = c(2,1),
       text.col = "black",
       bg = "white")

# Describe calibration plot of RSF and DeepSurv
x <-  data.frame(
  predicted_surv_rsf = mean_survival_risk_rsf,
  KM_surv_rsf  = KM_est_rsf,
  predicted_surv_deepS = mean_survival_risk_deepS,
  KM_surv_deepS  = KM_est_deepS
)

h <- ggplot(x, aes())
h <- h + theme_light()
h <- h + geom_line(aes(x = predicted_surv_rsf, y = KM_surv_rsf, colour='RSF'), size=1.0)
h <- h + geom_line(aes(x = predicted_surv_deepS, y = KM_surv_deepS, colour='DeepSurv'),linetype="longdash", size=1.0)
h <- h + geom_point(aes(x = predicted_surv_rsf, y = KM_surv_rsf, colour='RSF'),size = 1.5)
h <- h + geom_point(aes(x = predicted_surv_deepS, y = KM_surv_deepS, colour='DeepSurv'),size = 1.5)
h <- h + xlim(0.2,1.0)+ylim(0.2,1.0)
h <- h + geom_abline(intercept = 0, slope = 1, 
                     colour = "blue",
                     size = 0.5,
                     alpha = 0.5,
                     linetype = 2) 
h <- h + labs(x="Predicted probability of Overall Survival", y="Observed probability of Overall Survival")
h <- h + theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
h <- h + theme(legend.title=element_blank())
h <- h + scale_color_manual(values = c("navy","goldenrod"))
plot(h)
