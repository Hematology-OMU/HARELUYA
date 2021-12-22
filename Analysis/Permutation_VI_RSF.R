################################################################################################################ 
## Calculate Permutation Variable Importance for the predictive probability of 1-year OS after allo-HSCT using RSF
################################################################################################################ 

library(survival)
library(timeROC)
library(randomForestSRC)
library(tidyverse)
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
set.seed(karafuto)
Outcome_Month <-12
train_df <- Analysed_df_imp %>% filter(.SCT.Year <= 2016)
test_df <- Analysed_df_imp %>% filter(.SCT.Year >= 2017)

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


# Function calculating permutation variable importance based on Arguments (modelRFSRC: RSF model, test_df: test deta set)
calc_permVI_RSF <- function(modelRFSRC,test_df){
  
  # At first, calculate baseline AUC value before permutation
  predict_EFS_rfsrc <- function(x,fu_month){
    y <- which( predict(modelRFSRC, x)$time.interest <= fu_month)
    return (predict(modelRFSRC, x)$survival[,max(y)])
  }
  rfsrc_predict = createEmptyDf(nrow(test_df), 1, colnames = c( "study_ID" ) )
  for( i in 1:nrow(test_df) ){
    rfsrc_predict[ i, 1 ] = as.character(test_df[i,"study_ID"])
  }
  temp_p <- predict_EFS_rfsrc(test_df,Outcome_Month)
  rfsrc_predict <- cbind(rfsrc_predict,temp_p)
  colnames(rfsrc_predict)[2] <- c("rfsrc_predict_OS")
  
  test_rfsrc <- test_df %>% select(c(study_ID, .Event, .MonthEvent))
  test_rfsrc <- merge(test_rfsrc, rfsrc_predict, by="study_ID", all=T)
  rm(rfsrc_predict)
  
  # Calculate AUC for the probability of 1 year OS after allo-HSCT in test cohort as a baseline AUC value
  rfsrc_tROC <-timeROC(T=test_rfsrc$.MonthEvent,
                       delta=test_rfsrc$.Event,marker=(1-test_rfsrc$rfsrc_predict_OS),
                       cause=1,weighting="marginal",
                       times=c(Outcome_Month),
                       iid=FALSE)
  test_AUC <- rfsrc_tROC$AUC[[2]]
  
  # Function calculating AUC after permutation of test_column
  calc_permAUC <- function(modelRFSRC,test_df,test_column){
    
    test_df_perm <- test_df
    temp <- test_df_perm[,paste(test_column)] 
    test_df_perm[,paste(test_column)] <- sample(temp,length(temp),replace = FALSE)
    
    perm_predict = createEmptyDf(nrow(test_df_perm), 1, colnames = c("study_ID"))
    for( i in 1:nrow(test_df_perm) ){
      perm_predict[ i, 1 ] = as.character(test_df_perm[i,"study_ID"])
    }
    temp_p <- predict_EFS_rfsrc(test_df_perm,Outcome_Month)
    perm_predict <- cbind(perm_predict,temp_p)
    colnames(perm_predict)[2] <- c("perm_predict_OS")
    
    perm_test_rfsrc <- test_df_perm %>% select(c(study_ID, .Event, .MonthEvent))
    perm_test_rfsrc <- merge(perm_test_rfsrc, perm_predict, by="study_ID", all=T)
    
    perm_tROC <-timeROC(T=perm_test_rfsrc$.MonthEvent,
                         delta=perm_test_rfsrc$.Event,marker=(1-perm_test_rfsrc$perm_predict_OS),
                         cause=1,weighting="marginal",
                         times=c(Outcome_Month),
                         iid=FALSE)
    perm_AUC <- perm_tROC$AUC[[2]]
    return(perm_AUC)
  }
  
  Age_permAUC <- calc_permAUC(modelRFSRC,test_df, ".Age")
  HCTCI_permAUC <- calc_permAUC(modelRFSRC,test_df, ".HCT.CI")
  PS_permAUC <- calc_permAUC(modelRFSRC,test_df, ".PS24")
  Disease_permAUC <- calc_permAUC(modelRFSRC,test_df, ".Analysed_Disease")
  Dis_Status_permAUC <- calc_permAUC(modelRFSRC,test_df, ".Disease.satatus.all")
  Tx_pattern_permAUC <- calc_permAUC(modelRFSRC,test_df, "Tx_pattern")
  RCMV_permAUC <- calc_permAUC(modelRFSRC,test_df, "R_CMVAntiB")
  
  # Difference between baseline AUC and permutation AUC values 
  permVI_df <- rbind(c("Age", (test_AUC - Age_permAUC)),
                     c("HCT-CI", (test_AUC - HCTCI_permAUC)),
                     c("PS", (test_AUC - PS_permAUC)),
                     c("Disease", (test_AUC - Disease_permAUC)),
                     c("Disease status",(test_AUC - Dis_Status_permAUC)),
                     c("Transplant procedure",(test_AUC - Tx_pattern_permAUC)),
                     c("CMV serostatus",(test_AUC - RCMV_permAUC)))
  colnames(permVI_df ) <- c( "Var_name", "Var_Importance")
  return(permVI_df)
}

permVI_df <- calc_permVI_RSF(modelRFSRC,test_df)
permVI_df <- as.data.frame(permVI_df)
permVI_df$Var_Importance <- as.numeric(permVI_df$Var_Importance)
permVI_df <- permVI_df %>% arrange(-Var_Importance)

g = ggplot(permVI_df)
g = g + aes(x = reorder(Var_name, Var_Importance), y = Var_Importance)
g = g + geom_bar(stat="identity")
g = g + theme_bw() + theme(legend.position = "none")
g = g + xlab("Variable Name") + ylab("Variable Importance")
g = g + coord_flip(ylim=c(-0.01,0.125))
g = g + scale_fill_manual(values=c("#80c0ff"))
plot(g)
