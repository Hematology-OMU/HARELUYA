################################################################################################################ 
## Calculate Permutation Variable Importance for the predictive probability of 1-year OS after allo-HSCT using DeepSurv
## This code requires the CSV files before and after variable permutation which include the predictive probabilities of 1-year OS after allo-HSCT in the test cohort using DeepSurv with python
################################################################################################################ 

library(survival)
library(timeROC)
library(tidyverse)

Outcome_Month <-12

DS_base_pred_df = read.csv('./data/deepS_test_base_prediction.csv', fileEncoding = 'cp932')
DS_Age_pred_df = read.csv('./data/deep_srv_perm_Age_prediction.csv', fileEncoding = 'cp932')
DS_HCTCI_pred_df = read.csv('./data/deep_srv_perm_HCTCI_prediction.csv', fileEncoding = 'cp932')
DS_PS_pred_df = read.csv('./data/deep_srv_perm_PS_prediction.csv', fileEncoding = 'cp932')
DS_Dis_pred_df = read.csv('./data/deep_srv_perm_Dis_prediction.csv', fileEncoding = 'cp932')
DS_DisSt_pred_df = read.csv('./data/deep_srv_perm_DisSt_prediction.csv', fileEncoding = 'cp932')
DS_TxP_pred_df = read.csv('./data/deep_srv_perm_TxP_prediction.csv', fileEncoding = 'cp932')
DS_RCMV_pred_df = read.csv('./data/deep_srv_perm_RCMV_prediction.csv', fileEncoding = 'cp932')

# Function calculating AUC based on Arguments (df, Outcome_Month) using timeROC 
calc_AUC <- function(df,Outcome_Month){
  tROC <-timeROC(T=df$.MonthOS,
                delta=df$.OS,marker=(1-df$predicted_1y_OS),
                cause=1,weighting="marginal",
                times=c(Outcome_Month),
                iid=FALSE)
  AUC <- tROC$AUC[[2]]
  return(AUC)
}

DS_base_AUC <- calc_AUC(DS_base_pred_df,Outcome_Month)
DS_perm_Age_AUC <- calc_AUC(DS_Age_pred_df,Outcome_Month)
DS_perm_HCTCI_AUC <- calc_AUC(DS_HCTCI_pred_df,Outcome_Month)
DS_perm_PS_AUC <- calc_AUC(DS_PS_pred_df,Outcome_Month)
DS_perm_Dis_AUC <- calc_AUC(DS_Dis_pred_df,Outcome_Month)
DS_perm_DisSt_AUC <- calc_AUC(DS_DisSt_pred_df,Outcome_Month)
DS_perm_TxP_AUC <- calc_AUC(DS_TxP_pred_df,Outcome_Month)
DS_perm_RCMV_AUC <- calc_AUC(DS_RCMV_pred_df,Outcome_Month)

permVI_df <- rbind(c("Age", (DS_base_AUC - DS_perm_Age_AUC)),
                   c("HCT-CI", (DS_base_AUC - DS_perm_HCTCI_AUC)),
                   c("PS", (DS_base_AUC - DS_perm_PS_AUC)),
                   c("Disease", (DS_base_AUC - DS_perm_Dis_AUC)),
                   c("Disease status",(DS_base_AUC - DS_perm_DisSt_AUC)),
                   c("Transplant procedure",(DS_base_AUC - DS_perm_TxP_AUC)),
                   c("CMV serostatus",(DS_base_AUC - DS_perm_RCMV_AUC)))
colnames(permVI_df ) <- c( "Var_name", "Var_Importance")
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

