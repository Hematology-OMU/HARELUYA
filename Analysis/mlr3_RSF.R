################################################################################################################ 
## Random search for optimal hyper-parameters (myry, nodesize, and nodedepth) in RSF using mlr3 packeage
################################################################################################################ 

library(tidyverse)

FILEPATH = './data/analysed_df.csv'
Analysed_df_imp = read.csv(FILEPATH, fileEncoding = 'cp932', stringsAsFactors=TRUE)
Analysed_df_imp <- select(Analysed_df_imp, -c(study_ID,.Sex.Mismatch3,dx_to_sct_day_trump))
Analysed_df_imp <- Analysed_df_imp %>% rename(.Event = .OS)
Analysed_df_imp <- Analysed_df_imp %>% rename(.MonthEvent = .MonthOS)
Analysed_df_imp$.HCT.CI <- as.factor(Analysed_df_imp$.HCT.CI)
Analysed_df_imp$.PS24 <- as.factor(Analysed_df_imp$.PS24)
Analysed_df_imp$.Analysed_Disease <- as.factor(Analysed_df_imp$.Analysed_Disease)
Analysed_df_imp$.Disease.satatus.all <- as.factor(Analysed_df_imp$.Disease.satatus.all)
Analysed_df_imp$Tx_pattern <- as.factor(Analysed_df_imp$Tx_pattern)
Analysed_df_imp$R_CMVAntiB <- as.factor(Analysed_df_imp$R_CMVAntiB)
Outcome_Month <-12

karafuto <- 1808
set.seed(karafuto)

train_df <- Analysed_df_imp %>% filter(.SCT.Year <= 2016)
test_df <- Analysed_df_imp %>% filter(.SCT.Year >= 2017)
train_df <- train_df %>% select(-.SCT.Year)

library(survivalmodels)
library(mlr3)
library(mlr3proba)
library(paradox)
library(mlr3pipelines)
library(mlr3extralearners)
library(mlr3tuning)

TrumpML <- TaskSurv$new("TrumpML", train_df, time = ".MonthEvent", event = ".Event")

search_space <- ps(
  ## p_int for integer valued parameters
  mtry = p_int(lower = 1, upper = 7),
  nodesize = p_int(lower = 10, upper = 80),
  nodedepth = p_int(lower = 10, upper = 100)
)

learner <- lrn("surv.rfsrc",ntree = 500)

instance <- TuningInstanceSingleCrit$new(
  task = TrumpML,
  learner = learner,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  search_space = search_space,
  terminator = trm("evals", n_evals = 100),
  store_benchmark_result = FALSE,
  store_models = FALSE
)
tt = tnr("random_search")

# perform random search for optimal hyper-parameters
tt$optimize(instance)

# Show the result
instance$result_learner_param_vals
instance$result_y
instance$archive
instance$archive$benchmark_result
instance$archive$benchmark_result$score(msr("surv.cindex"))
