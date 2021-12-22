################################################################################################################ 
## Comparison of the OS between recommendation and non-recommendation groups using RSF
################################################################################################################ 

library(survival)
library(randomForestSRC)
library(tidyverse)
library(survminer)
library(car)

karafuto <- 1808
set.seed(karafuto)

FILEPATH = './data/analysed_df.csv'
Analysed_df_imp = read.csv(FILEPATH, fileEncoding = 'cp932', stringsAsFactors=TRUE)

Analysed_df_imp <- Analysed_df_imp %>% rename(.Event = .OS)
Analysed_df_imp <- Analysed_df_imp %>% rename(.MonthEvent = .MonthOS)
Analysed_df_imp$.HCT.CI <- as.factor(Analysed_df_imp$.HCT.CI)
Analysed_df_imp$.PS24 <- as.factor(Analysed_df_imp$.PS24)
Analysed_df_imp$.Analysed_Disease <- as.factor(Analysed_df_imp$.Analysed_Disease)
Analysed_df_imp$.Disease.satatus.all <- as.factor(Analysed_df_imp$.Disease.satatus.all)
Analysed_df_imp$Tx_pattern <- as.factor(Analysed_df_imp$Tx_pattern)
Analysed_df_imp$.Sex.Mismatch3 <- as.factor(Analysed_df_imp$.Sex.Mismatch3)
Analysed_df_imp$R_CMVAntiB <- as.factor(Analysed_df_imp$R_CMVAntiB)
Outcome_Month <-12

# Set hyper-parameters of RSF
mtry = 2
nodesize = 80
nodedepth = 59
ntree = 500

train_df <- Analysed_df_imp %>% filter(.SCT.Year <= 2016)
test_df <- Analysed_df_imp %>% filter(.SCT.Year >= 2017)

train_df %>% group_by(Tx_pattern) %>% summarize(count =n())
test_df %>% group_by(Tx_pattern) %>% summarize(count =n())

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
                    mtry = mtry, nodesize = nodesize, nodedepth = nodedepth)

# Replace from Tx_pattern to Tx_pattern_actual
select_colmn_test <- c("study_ID",".Age",".HCT.CI",".PS24",".Analysed_Disease",".Disease.satatus.all",
                       "Tx_pattern","R_CMVAntiB")
test_df <- test_df %>% select(select_colmn_test)
test_df <- rename(test_df, Tx_pattern_actual = Tx_pattern)

# Create a function to return predicted outcome value at fu_month for x (dataframe)
predict_outcome_rsf <- function(x,fu_month){
  y <- which( predict(modelRFSRC, x)$time.interest <= fu_month)
  return (predict(modelRFSRC, x)$survival[,max(y)])
}
createEmptyDf = function( nrow, ncol, colnames = c() ){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}
rfsrc_p_month = createEmptyDf(nrow(test_df), 1, colnames = c("study_ID"))

# Calculate 10 predictive probabilities of 1-year OS after allo-HSCT for each Tx_pattern in test cohort
for( i in 1:nrow(test_df) ){
  rfsrc_p_month[ i, 1 ] = as.character(test_df[i,"study_ID"])
}
for( j in 1:10){
  test_df_temp <- test_df %>% mutate(Tx_pattern = factor(j,levels=c(1,2,3,4,5,6,7,8,9,10)))
  
  start.time<-proc.time()
  rfsrc_p_month <- cbind(rfsrc_p_month, predict_outcome_rsf(test_df_temp,Outcome_Month))
  end.time<-proc.time()
  print(paste("Tx_pattern_no: ",j, "process_time: ", round((end.time-start.time)[[3]]/60,2),"min"))
  
}
colnames(rfsrc_p_month)[2:11] <- c("rfsrc_p_1","rfsrc_p_2","rfsrc_p_3","rfsrc_p_4","rfsrc_p_5",
                                   "rfsrc_p_6","rfsrc_p_7","rfsrc_p_8","rfsrc_p_9","rfsrc_p_10")

# Calculate SD of 10 predictive probabilities of 1-year OS after allo-HSCT
rfsrc_p_month <- rfsrc_p_month %>% mutate(p_sd = apply(rfsrc_p_month %>% select(-"study_ID"),1,sd))

# Merge 10 predictive probabilities of 1-year OS after allo-HSCT in each patient to test_df
test_df <- merge(test_df, rfsrc_p_month, by="study_ID", all=T)

# When the actual donor is unrelated (Tx_pattern_actual = 2,7), predictive probabilities of 1-year OS for related donors (Tx_pattern_actual = 1,6) are put into -1 not to recommend them as the optimal donor.
test_df <- test_df %>% mutate( rfsrc_p_1  = if_else( Tx_pattern_actual == 2 | Tx_pattern_actual == 7 , -1, rfsrc_p_1) )
test_df <- test_df %>% mutate( rfsrc_p_6  = if_else( Tx_pattern_actual == 2 | Tx_pattern_actual == 7 , -1, rfsrc_p_1) )

# When the actual donor is CB or Haplo (Tx_pattern_actual = 3,4,5,8,9,10), predictive probabilities of 1-year OS for related and unrelated donor (Tx_pattern_actual = 1,2,6,7) are put into -1 not to recommend them as the the optimal donor.
test_df <- test_df %>% mutate( rfsrc_p_1  = if_else( Tx_pattern_actual == 3 | Tx_pattern_actual == 4 |
                                                       Tx_pattern_actual == 5 | Tx_pattern_actual == 8 |
                                                       Tx_pattern_actual == 9 | Tx_pattern_actual == 10, -1, rfsrc_p_1) )
test_df <- test_df %>% mutate( rfsrc_p_2  = if_else( Tx_pattern_actual == 3 | Tx_pattern_actual == 4 |
                                                       Tx_pattern_actual == 5 | Tx_pattern_actual == 8 |
                                                       Tx_pattern_actual == 9 | Tx_pattern_actual == 10, -1, rfsrc_p_2) )
test_df <- test_df %>% mutate( rfsrc_p_6  = if_else( Tx_pattern_actual == 3 | Tx_pattern_actual == 4 |
                                                       Tx_pattern_actual == 5 | Tx_pattern_actual == 8 |
                                                       Tx_pattern_actual == 9 | Tx_pattern_actual == 10, -1, rfsrc_p_6) )
test_df <- test_df %>% mutate( rfsrc_p_7  = if_else( Tx_pattern_actual == 3 | Tx_pattern_actual == 4 |
                                                       Tx_pattern_actual == 5 | Tx_pattern_actual == 8 |
                                                       Tx_pattern_actual == 9 | Tx_pattern_actual == 10, -1, rfsrc_p_7) )

# Tx_pattern showing the highest predictive probability of 1-year OS after allo-HSCT among 10 predictive probabilities is put into Tx_pattern_best column
test_df <- cbind(test_df, Tx_pattern_best = apply(test_df, 1, function(x){which.max(x[9:18])}))
test_df$Tx_pattern_best <- factor(test_df$Tx_pattern_best, levels = c(1,2,3,4,5,6,7,8,9,10))
# Discriminate recommendation and non-recommendation groups based on match or mismatch between Tx_pattern_actual and Tx_pattern_best
test_df <- test_df %>% mutate( recomm_Tx  = if_else( Tx_pattern_actual == Tx_pattern_best, 1, 0) )

# recomm_Tx_caliper takes a specific caliper range into consideration in recommendation group within caliper range 
caliper <- 0.25
test_df <- cbind(test_df, p_max = apply(test_df[9:18], 1, max))
test_df <- cbind(test_df, p_min = apply(test_df,1,function(x){return(as.numeric(x["p_max"])-as.numeric(x["p_sd"])*caliper)}))
test_df <- test_df %>% mutate(recomm_Tx_caliper = case_when( Tx_pattern_actual == 1 ~ ifelse(test_df$rfsrc_p_1 <= p_max & test_df$rfsrc_p_1 >= p_min,1,0),
                                                   Tx_pattern_actual == 2 ~ ifelse(test_df$rfsrc_p_2 <= p_max & test_df$rfsrc_p_2 >= p_min,1,0),
                                                   Tx_pattern_actual == 3 ~ ifelse(test_df$rfsrc_p_3 <= p_max & test_df$rfsrc_p_3 >= p_min,1,0),
                                                   Tx_pattern_actual == 4 ~ ifelse(test_df$rfsrc_p_4 <= p_max & test_df$rfsrc_p_4 >= p_min,1,0),
                                                   Tx_pattern_actual == 5 ~ ifelse(test_df$rfsrc_p_5 <= p_max & test_df$rfsrc_p_5 >= p_min,1,0),
                                                   Tx_pattern_actual == 6 ~ ifelse(test_df$rfsrc_p_6 <= p_max & test_df$rfsrc_p_6 >= p_min,1,0),
                                                   Tx_pattern_actual == 7 ~ ifelse(test_df$rfsrc_p_7 <= p_max & test_df$rfsrc_p_7 >= p_min,1,0),
                                                   Tx_pattern_actual == 8 ~ ifelse(test_df$rfsrc_p_8 <= p_max & test_df$rfsrc_p_8 >= p_min,1,0),
                                                   Tx_pattern_actual == 9 ~ ifelse(test_df$rfsrc_p_9 <= p_max & test_df$rfsrc_p_9 >= p_min,1,0),
                                                   Tx_pattern_actual == 10 ~ ifelse(test_df$rfsrc_p_10 <= p_max & test_df$rfsrc_p_10 >= p_min,1,0)))
test_df <- merge(test_df %>% select(study_ID,recomm_Tx, recomm_Tx_caliper, Tx_pattern_actual, Tx_pattern_best),Analysed_df_imp %>% filter(.SCT.Year >= 2017), by = 'study_ID', all.x = T )

test_df$recomm_Tx <- as.factor(test_df$recomm_Tx)
test_df$recomm_Tx_caliper <- as.factor(test_df$recomm_Tx_caliper)

# Univariate analysis for OS between recommendation and non-recommendation groups in all cases using log-lank and cox PH model
survdiff(Surv(.MonthEvent,.Event) ~ recomm_Tx, data=test_df)
COX_comp_recomm <- coxph(Surv(.MonthEvent, .Event) ~ recomm_Tx, data = test_df)
summary(COX_comp_recomm)
# Plot KM OS curve
survfit_uni <-survfit(Surv(.MonthEvent,.Event)~recomm_Tx, data=test_df)
j_univ <- ggsurvplot(survfit_uni, data=test_df, size=2,
                       legend.title = element_blank(),
                       legend=c(.8,.2), legend.labs = c("Non-recommendation", "Recommendation"),
                       xlab="Months after Transplantation", ylab="Probability of Overall Survival",censor=F, 
                       ggtheme = theme_light(base_size = 13, base_line_size = 1.5, base_rect_size = 0.8),
                       xlim = c(0,24),break.time.by = 6,
                       font.x = c(12, "bold", "black"),
                       font.y = c(12, "bold", "black"),
                       palette ="jco",
                       surv.scale = "percent",
                       ########## risk table #########,
                       risk.table = TRUE,
                       risk.table.title= "Number at risk",
                       risk.table.subtitle = "and remember about censoring.",
                       risk.table.height = 0.12,
                       tables.theme = theme_void(),
                       tables.y.text = FALSE,
                       fontsize = 4
)
j_univ

# Univariate analysis for OS between recommendation group within caliper range and non-recommendation group without caliper range in all cases using log-lank and cox PH model
survdiff(Surv(.MonthEvent,.Event) ~ recomm_Tx_caliper, data=test_df)
COX_comp_recomm_cal <- coxph(Surv(.MonthEvent, .Event) ~ recomm_Tx_caliper, data = test_df)
summary(COX_comp_recomm_cal)
# Plot KM OS curve
survfit_uni_cal <-survfit(Surv(.MonthEvent,.Event)~recomm_Tx_caliper, data=test_df)
j_univ_cal <- ggsurvplot(survfit_uni_cal, data=test_df, size=2,
                     legend.title = element_blank(),
                     legend=c(.8,.2), legend.labs = c("Non-recommendation without caliper range", "Recommendation within caliper range"),
                     xlab="Months after Transplantation", ylab="Probability of Overall Survival",censor=F, 
                     ggtheme = theme_light(base_size = 13, base_line_size = 1.5, base_rect_size = 0.8),
                     xlim = c(0,24),break.time.by = 6,
                     font.x = c(12, "bold", "black"),
                     font.y = c(12, "bold", "black"),
                     palette ="jco",
                     surv.scale = "percent",
                     ########## risk table #########,
                     risk.table = TRUE,
                     risk.table.title= "Number at risk",
                     risk.table.subtitle = "and remember about censoring.",
                     risk.table.height = 0.12,
                     tables.theme = theme_void(),
                     tables.y.text = FALSE,
                     fontsize = 4
)
j_univ_cal


# Multivariate analysis between recommendation and non-recommendation groups in all cases using cox PH model
COX_comp_recomm_multi<- coxph(Surv(.MonthEvent, .Event) ~ 
                               .Age
                             +.HCT.CI
                             +.PS24
                             +.Analysed_Disease
                             +.Disease.satatus.all
                             +.Sex.Mismatch3
                             +R_CMVAntiB
                             +dx_to_sct_day_trump 
                             +Tx_pattern_actual
                             +.SCT.Year
                             +recomm_Tx, data = test_df)
summary(COX_comp_recomm_multi)

# Plot cox PH survival curves stratified by recommendation and non-recommendation and adjusted by multivariate factors in all cases
COX_comp_match_multi_strata<- coxph(Surv(.MonthEvent, .Event) ~ 
                                      .Age
                                    +.HCT.CI
                                    +.PS24
                                    +.Analysed_Disease
                                    +.Disease.satatus.all
                                    +.Sex.Mismatch3
                                    +R_CMVAntiB
                                    +dx_to_sct_day_trump 
                                    +Tx_pattern_actual
                                    +.SCT.Year
                                    +strata(recomm_Tx), data = test_df)

survfit_strata <-survfit(COX_comp_match_multi_strata)
j_strata <- ggsurvplot(survfit_strata, data=test_df, size=2,
                       legend.title = element_blank(),
                       legend=c(.8,.2), legend.labs = c("Non-recommendation", "Recommendation"),
                       xlab="Months after Transplantation", ylab="Probability of Overall Survival",censor=F, 
                       ggtheme = theme_light(base_size = 13, base_line_size = 1.5, base_rect_size = 0.8),
                       xlim = c(0,24),break.time.by = 6,
                       font.x = c(12, "bold", "black"),
                       font.y = c(12, "bold", "black"),
                       palette ="jco",
                       surv.scale = "percent",
                       ########## risk table #########,
                       risk.table = TRUE,
                       risk.table.title= "Number at risk",
                       risk.table.subtitle = "and remember about censoring.",
                       risk.table.height = 0.12,
                       tables.theme = theme_void(),
                       tables.y.text = FALSE,
                       fontsize = 4
)
j_strata


# Multivariate analysis between recommendation group within caliper range and non-recommendation group without caliper range in all cases using cox PH model
COX_comp_recomm_multi_caliper<- coxph(Surv(.MonthEvent, .Event) ~
                                .Age
                              +.HCT.CI
                              +.PS24
                              +.Analysed_Disease
                              +.Disease.satatus.all
                              +.Sex.Mismatch3
                              +R_CMVAntiB
                              +dx_to_sct_day_trump
                              +Tx_pattern_actual
                              +.SCT.Year
                              +recomm_Tx_caliper, data = test_df)
summary(COX_comp_recomm_multi_caliper)

# Plot cox PH survival curves stratified by recommendation group within caliper range and non-recommendation group without caliper range and adjusted by multivariate factors in all cases
COX_comp_match_multi_strata_caliper<- coxph(Surv(.MonthEvent, .Event) ~
                                      .Age
                                    +.HCT.CI
                                    +.PS24
                                    +.Analysed_Disease
                                    +.Disease.satatus.all
                                    +.Sex.Mismatch3
                                    +R_CMVAntiB
                                    +dx_to_sct_day_trump
                                    +Tx_pattern_actual
                                    +.SCT.Year
                                    +strata(recomm_Tx_caliper), data = test_df)

survfit_strata_caliper <-survfit(COX_comp_match_multi_strata_caliper)
j_strata_caliper <- ggsurvplot(survfit_strata_caliper, data=test_df, size=2,
                       legend.title = element_blank(),
                       legend=c(.8,.2), legend.labs = c("Non-recommendation without caliper range", "Recommendation within caliper range"),
                       xlab="Months after Transplantation", ylab="Probability of Overall Survival",censor=F,
                       ggtheme = theme_light(base_size = 13, base_line_size = 1.5, base_rect_size = 0.8),
                       xlim = c(0,24),break.time.by = 6,
                       font.x = c(12, "bold", "black"),
                       font.y = c(12, "bold", "black"),
                       palette ="jco",
                       surv.scale = "percent",
                       ########## risk table #########,
                       risk.table = TRUE,
                       risk.table.title= "Number at risk",
                       risk.table.subtitle = "and remember about censoring.",
                       risk.table.height = 0.12,
                       tables.theme = theme_void(),
                       tables.y.text = FALSE,
                       fontsize = 4
)
j_strata_caliper


# Describe Chord diagram of the relationship between actual and RSF-recommended allo-HSCT procedures in all cases
library(circlize)
library(RColorBrewer)
ColColor <- colorRampPalette(c("#2ca9e1", "#e7609e", "#f7c114", "#543f32")) #行の色
RowColor <- colorRampPalette(c("#ec6d51", "#3f312b", "#674196", "#82ae46")) #列の色

test_df <- test_df %>% mutate(Tx_pattern_best = case_when(
  Tx_pattern_best == 1 ~ "R-MRD",
  Tx_pattern_best == 2 ~ "R-MUD",
  Tx_pattern_best == 3 ~ "R-UCB",
  Tx_pattern_best == 4 ~ "R-Hap-CY",
  Tx_pattern_best == 5 ~ "R-Hap-nCY",
  Tx_pattern_best == 6 ~ "M-MRD",
  Tx_pattern_best == 7 ~ "M-MUD",
  Tx_pattern_best == 8 ~ "M-UCB",
  Tx_pattern_best == 9 ~ "M-Hap-CY",
  Tx_pattern_best == 10 ~ "M-Hap-nCY"
))
test_df$Tx_pattern_best <- factor(test_df$Tx_pattern_best, 
                                  levels=c("R-MRD","R-MUD","R-UCB","R-Hap-CY","R-Hap-nCY","M-MRD","M-MUD","M-UCB","M-Hap-CY","M-Hap-nCY"))

test_df <- test_df %>% mutate(Tx_pattern_actual = case_when(
  Tx_pattern_actual == 1 ~ "R-MRD",
  Tx_pattern_actual == 2 ~ "R-MUD",
  Tx_pattern_actual == 3 ~ "R-UCB",
  Tx_pattern_actual == 4 ~ "R-Hap-CY",
  Tx_pattern_actual == 5 ~ "R-Hap-nCY",
  Tx_pattern_actual == 6 ~ "M-MRD",
  Tx_pattern_actual == 7 ~ "M-MUD",
  Tx_pattern_actual == 8 ~ "M-UCB",
  Tx_pattern_actual == 9 ~ "M-Hap-CY",
  Tx_pattern_actual == 10 ~ "M-Hap-nCY"
))
test_df$Tx_pattern_actual <- factor(test_df$Tx_pattern_actual, 
                                  levels=c("R-MRD","R-MUD","R-UCB","R-Hap-CY","R-Hap-nCY","M-MRD","M-MUD","M-UCB","M-Hap-CY","M-Hap-nCY"))

real_best_tx_matrix <- with(test_df, table(Tx_pattern_actual, Tx_pattern_best)) 

par(cex = 0.9,mar=c(0,0,0,0))
chordDiagram(real_best_tx_matrix,
             directional = TRUE,direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             transparency = 0.5,
             row.col = RowColor(nrow(real_best_tx_matrix)),
             grid.col = RowColor(nrow(real_best_tx_matrix))
)

# Describe chord diagram of the relationship between actual and recommended donors (Std and Alt)
test_df <- test_df %>% mutate(donor_best_std_alt = case_when(
  Tx_pattern_best == "R-MRD" ~ "Std",
  Tx_pattern_best == "R-MUD" ~ "Std",
  Tx_pattern_best == "R-UCB" ~ "Alt",
  Tx_pattern_best == "R-Hap-CY" ~ "Alt",
  Tx_pattern_best == "R-Hap-nCY" ~ "Alt",
  Tx_pattern_best == "M-MRD" ~ "Std",
  Tx_pattern_best == "M-MUD" ~ "Std",
  Tx_pattern_best == "M-UCB" ~ "Alt",
  Tx_pattern_best == "M-Hap-CY" ~ "Alt",
  Tx_pattern_best == "M-Hap-nCY" ~ "Alt"
))
test_df$donor_best_std_alt <- factor(test_df$donor_best_std_alt, levels=c("Std","Alt"))

test_df <- test_df %>% mutate(donor_real_std_alt = case_when(
  Tx_pattern_actual == "R-MRD" ~ "Std",
  Tx_pattern_actual == "R-MUD" ~ "Std",
  Tx_pattern_actual == "R-UCB" ~ "Alt",
  Tx_pattern_actual == "R-Hap-CY" ~ "Alt",
  Tx_pattern_actual == "R-Hap-nCY" ~ "Alt",
  Tx_pattern_actual == "M-MRD" ~ "Std",
  Tx_pattern_actual == "M-MUD" ~ "Std",
  Tx_pattern_actual == "M-UCB" ~ "Alt",
  Tx_pattern_actual == "M-Hap-CY" ~ "Alt",
  Tx_pattern_actual == "M-Hap-nCY" ~ "Alt"
))
test_df$donor_real_std_alt <- factor(test_df$donor_real_std_alt, levels=c("Std","Alt"))

std_alt_matrix <- with(test_df, table(donor_real_std_alt, donor_best_std_alt)) 

par(cex = 0.9,mar=c(0,0,0,0))
chordDiagram(std_alt_matrix,
             directional = TRUE,direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             transparency = 0.5, 
             row.col = RowColor(nrow(std_alt_matrix)),
             grid.col = RowColor(nrow(std_alt_matrix)) 
)

# Describe chord diagram of the relationship between actual and recommended intensities (RIC and MAC)
test_df <- test_df %>% mutate(best_RIC_MAC = case_when(
  Tx_pattern_best == "R-MRD" ~ "RIC",
  Tx_pattern_best == "R-MUD" ~ "RIC",
  Tx_pattern_best == "R-UCB" ~ "RIC",
  Tx_pattern_best == "R-Hap-CY" ~ "RIC",
  Tx_pattern_best == "R-Hap-nCY" ~ "RIC",
  Tx_pattern_best == "M-MRD" ~ "MAC",
  Tx_pattern_best == "M-MUD" ~ "MAC",
  Tx_pattern_best == "M-UCB" ~ "MAC",
  Tx_pattern_best == "M-Hap-CY" ~ "MAC",
  Tx_pattern_best == "M-Hap-nCY" ~ "MAC"
))
test_df$best_RIC_MAC <- factor(test_df$best_RIC_MAC, levels=c("RIC","MAC"))

test_df <- test_df %>% mutate(real_RIC_MAC = case_when(
  Tx_pattern_actual == "R-MRD" ~ "RIC",
  Tx_pattern_actual == "R-MUD" ~ "RIC",
  Tx_pattern_actual == "R-UCB" ~ "RIC",
  Tx_pattern_actual == "R-Hap-CY" ~ "RIC",
  Tx_pattern_actual == "R-Hap-nCY" ~ "RIC",
  Tx_pattern_actual == "M-MRD" ~ "MAC",
  Tx_pattern_actual == "M-MUD" ~ "MAC",
  Tx_pattern_actual == "M-UCB" ~ "MAC",
  Tx_pattern_actual == "M-Hap-CY" ~ "MAC",
  Tx_pattern_actual == "M-Hap-nCY" ~ "MAC"
))
test_df$real_RIC_MAC <- factor(test_df$real_RIC_MAC, levels=c("RIC","MAC"))

RIC_MAC_matrix <- with(test_df, table(real_RIC_MAC, best_RIC_MAC)) 

par(cex = 0.9,mar=c(0,0,0,0))
chordDiagram(RIC_MAC_matrix,
             directional = TRUE,direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             transparency = 0.5,
             row.col = RowColor(nrow(RIC_MAC_matrix)),
             grid.col = RowColor(nrow(RIC_MAC_matrix))
)

# Ignore conditioning intensity category
test_df$Tx_pattern_actual <- as.character(test_df$Tx_pattern_actual)
test_df$Tx_pattern_best <- as.character(test_df$Tx_pattern_best)
test_df$Tx_pattern_actual  <- str_remove(test_df$Tx_pattern_actual, "R-")
test_df$Tx_pattern_actual  <- str_remove(test_df$Tx_pattern_actual, "M-")
test_df$Tx_pattern_best  <- str_remove(test_df$Tx_pattern_best, "R-")
test_df$Tx_pattern_best  <- str_remove(test_df$Tx_pattern_best, "M-")
test_df$Tx_pattern_actual <- factor(test_df$Tx_pattern_actual,
                                    levels=c("MRD","MUD","UCB","Hap-CY","Hap-nCY"))
test_df$Tx_pattern_best <- factor(test_df$Tx_pattern_best,
                                 levels=c("MRD","MUD","UCB","Hap-CY","Hap-nCY"))

test_df <- test_df %>% mutate( donor_matched  = if_else( Tx_pattern_actual == Tx_pattern_best, 1, 0) )

# Describe chord diagram of the relationship between actual and RSF-recommended allo-HSCT donors
donor_matrix <- with(test_df, table(Tx_pattern_actual, Tx_pattern_best)) 

par(cex = 0.9,mar=c(0,0,0,0))
chordDiagram(donor_matrix,
             directional = TRUE,direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             transparency = 0.5,
             row.col = RowColor(nrow(donor_matrix)),
             grid.col = RowColor(nrow(donor_matrix))
)

# Sub-group analysis in patients transplanted from alternative donor (CB and Haplo)
# Comparison of OS between the recommendation and non-recommendation groups in the patient subgroup transplanted from alternative donors
# filter cases transplanted from alternative donors
test_df_real_alt <- filter(test_df, donor_real_std_alt == "Alt")
test_df_real_alt$Tx_pattern_actual <- factor(test_df_real_alt$Tx_pattern_actual, 
                                                levels=c("UCB","Hap-CY","Hap-nCY"))


# Univariate analysis for OS between recommendation and non-recommendation groups (CB, Haplo-CY, Haplo-nCY) in the patients subgroup transplanted from alternative donor using log-lank test and cox PH model
survdiff(Surv(.MonthEvent,.Event) ~ donor_matched, data=test_df_real_alt)
COX_comp_match_real_alt_univ <- coxph(Surv(.MonthEvent, .Event) ~ donor_matched, data = test_df_real_alt)
summary(COX_comp_match_real_alt_univ)
# Plot KM OS curve
survfit_comp_match_real_alt_uni <-survfit(Surv(.MonthEvent,.Event)~donor_matched, data=test_df_real_alt)
j_comp_match_real_alt_univ <- ggsurvplot(survfit_comp_match_real_alt_uni, data=test_df_real_alt, size=2,
                              legend.title = element_blank(),
                              legend=c(.8,.2), legend.labs = c("Non-recommendation", "Recommendation"),
                              xlab="Months after Transplantation", ylab="Probability of Overall Survival",censor=F, 
                              ggtheme = theme_light(base_size = 13, base_line_size = 1.5, base_rect_size = 0.8),
                              xlim = c(0,24),break.time.by = 6,
                              font.x = c(12, "bold", "black"),
                              font.y = c(12, "bold", "black"),
                              palette ="jco",
                              surv.scale = "percent",
                              ########## risk table #########,
                              risk.table = TRUE,
                              risk.table.title= "Number at risk",
                              risk.table.subtitle = "and remember about censoring.",
                              risk.table.height = 0.12,
                              tables.theme = theme_void(),
                              tables.y.text = FALSE,
                              fontsize = 4
)
j_comp_match_real_alt_univ

# Multivariate analysis between recommendation and non-recommendation groups in the patients subgroup transplanted from alternative donor
COX_comp_match_real_alt_multi<- coxph(Surv(.MonthEvent, .Event) ~ .Age
                                      +.HCT.CI
                                      +.PS24
                                      +.Analysed_Disease
                                      +.Disease.satatus.all
                                      +.Sex.Mismatch3
                                      +R_CMVAntiB
                                      +dx_to_sct_day_trump
                                      +real_RIC_MAC
                                      +Tx_pattern_actual
                                      +.SCT.Year
                                      +donor_matched, data = test_df_real_alt)
summary(COX_comp_match_real_alt_multi)

# Plot cox PH survival curves stratified by recommendation and non-recommendation and adjusted by multivariate factors in the patients subgroup transplanted from alternative donor
COX_comp_match_real_alt_multi_strata<- coxph(Surv(.MonthEvent, .Event) ~  
                                               .Age
                                             +.HCT.CI
                                             +.PS24
                                             +.Analysed_Disease
                                             +.Disease.satatus.all
                                             +.Sex.Mismatch3
                                             +R_CMVAntiB
                                             +dx_to_sct_day_trump 
                                             +real_RIC_MAC
                                             +Tx_pattern_actual
                                             +.SCT.Year
                                              +strata(donor_matched), data = test_df_real_alt)

survfit_real_alt_matched_strata <-survfit(COX_comp_match_real_alt_multi_strata)
j_real_alt_matched_strata <- ggsurvplot(survfit_real_alt_matched_strata, data=test_df_real_alt, size=2,
                                       legend.title = element_blank(),
                                       legend=c(.8,.2), legend.labs = c("Non-recommendation", "Recommendation"),
                                       xlab="Months after Transplantation", ylab="Probability of Overall Survival",censor=F, 
                                       ggtheme = theme_light(base_size = 13, base_line_size = 1.5, base_rect_size = 0.8),
                                       xlim = c(0,24),break.time.by = 6,
                                       font.x = c(12, "bold", "black"),
                                       font.y = c(12, "bold", "black"),
                                       palette ="jco",
                                       surv.scale = "percent",
                                       ########## risk table #########,
                                       risk.table = TRUE,
                                       risk.table.title= "Number at risk",
                                       risk.table.subtitle = "and remember about censoring.",
                                       risk.table.height = 0.12,
                                       tables.theme = theme_void(),
                                       tables.y.text = FALSE,
                                       fontsize = 4
)
j_real_alt_matched_strata

