#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn_pandas import DataFrameMapper
import torch
import torchtuples as tt
from pycox.datasets import metabric
from pycox.models import CoxPH
from pycox.evaluation import EvalSurv
import os

#%%
karafuto = 1
epochs = 300
batch_size = 128
lr = 0.001
dropout = 0.534917552
batch_norm = True
n_nodes = 64
frac=0.2

# %%
np.random.seed(karafuto)
_ = torch.manual_seed(karafuto)
dataset_fp = os.path.join(os.path.dirname(__file__),'data\\analysed_df.csv')
analysed_df = pd.read_csv(dataset_fp, encoding = 'cp932')
analysed_df = analysed_df[['.SCT.Year','.Age','.HCT.CI','.PS24','.Analysed_Disease','.Disease.satatus.all',
                            'Tx_pattern','R_CMVAntiB','.Sex.Mismatch3','dx_to_sct_day_trump','.MonthOS','.OS']]
analysed_df.head()

# %%
analysed_df_dum = pd.get_dummies(analysed_df, columns=['.HCT.CI','.PS24','.Analysed_Disease','.Disease.satatus.all','Tx_pattern','R_CMVAntiB'])
# %%
df_train = analysed_df_dum[analysed_df_dum['.SCT.Year'] <= 2016]
df_val = df_train.sample(frac=frac, random_state = karafuto)
df_train = df_train.drop(df_val.index)
df_test = analysed_df_dum[analysed_df_dum['.SCT.Year'] >= 2017]

print(len(df_train))
print(len(df_val))
print(len(df_test))

# %%
cols_standardize = ['.Age']
cols_leave = ['.Analysed_Disease_ALL', '.Analysed_Disease_AML',
       '.Analysed_Disease_ATL', '.Analysed_Disease_MDS',
       '.Analysed_Disease_ML', '.Analysed_Disease_MPN', 
       '.Disease.satatus.all_CR', '.Disease.satatus.all_nonCR', 'Tx_pattern_1',
       'Tx_pattern_2', 'Tx_pattern_3', 'Tx_pattern_4', 'Tx_pattern_5',
       'Tx_pattern_6', 'Tx_pattern_7', 'Tx_pattern_8', 'Tx_pattern_9',
       'Tx_pattern_10', '.HCT.CI_0', '.HCT.CI_1', '.HCT.CI_2', '.HCT.CI_3', '.PS24_0',
       '.PS24_1', '.PS24_2', '.PS24_3', '.PS24_4','R_CMVAntiB_no', 'R_CMVAntiB_yes']

standardize = [([col], StandardScaler()) for col in cols_standardize]
leave = [(col, None) for col in cols_leave]
x_mapper = DataFrameMapper(standardize + leave)

# %%
x_train = x_mapper.fit_transform(df_train).astype('float32')
x_val = x_mapper.transform(df_val).astype('float32')
x_test = x_mapper.transform(df_test).astype('float32')

# %%
get_target = lambda df: (df['.MonthOS'].values, df['.OS'].values)
y_train = get_target(df_train)
y_val = get_target(df_val)
durations_test, events_test = get_target(df_test)
val = x_val, y_val
# %%
in_features = x_train.shape[1]
num_nodes = [n_nodes,n_nodes,n_nodes,n_nodes,n_nodes,n_nodes]
out_features = 1
batch_norm = batch_norm
dropout = dropout
output_bias = False
net = tt.practical.MLPVanilla(in_features, num_nodes, out_features, batch_norm,
                              dropout, activation=torch.nn.ReLU, output_bias=output_bias)
# %%
model = CoxPH(net, tt.optim.Adam())

# %%
batch_size = batch_size
model.optimizer.set_lr(lr)
epochs = epochs
callbacks = [tt.callbacks.EarlyStopping()]
verbose = True

#%%time
log = model.fit(x_train, y_train, batch_size, epochs, callbacks, verbose,
                val_data=val, val_batch_size=batch_size)
_ = log.plot()
model.partial_log_likelihood(*val).mean()
_ = model.compute_baseline_hazards()

# %%
# Calculate the predictive probabilities of 1-year OS after allo-HSCT in test cohort using DeepSurv and output results to CSV file
surv_performacne = model.predict_surv_df(x_test)
temp_surv_12_performance = surv_performacne.loc[[12]].transpose()
temp_surv_12_performance = temp_surv_12_performance.reset_index()
temp_surv_12_performance = temp_surv_12_performance[12]
temp_surv_12_performance = temp_surv_12_performance.rename('deepS_p_predict')

ev = EvalSurv(surv_performacne, durations_test, events_test, censor_surv='km')
ev.concordance_td()
# %%
testdf_performance = analysed_df[analysed_df['.SCT.Year'] >= 2017]
testdf_performance = testdf_performance.reset_index()
testdf_performance = pd.concat([testdf_performance, temp_surv_12_performance], axis=1)

# %%
testdf_performance.to_csv(os.path.join(os.path.dirname(__file__),'data\\deep_s_12prediction.csv'),index = False)

# %%
# Calculate the predictive probabilities of 1-year OS after allo-HSCT for ten allo-HSCT procedures in each test cohort patient
df_12_month_p = pd.DataFrame(index = range(0, len(df_test)-1 ), columns=[])
for j in range(1, 11):

       df_test_temp = analysed_df[analysed_df['.SCT.Year'] >= 2017]
       df_test_temp = df_test_temp.drop(columns=['Tx_pattern','.MonthOS','.OS'])
       df_test_temp['Tx_pattern_1'] = 0
       df_test_temp['Tx_pattern_2'] = 0
       df_test_temp['Tx_pattern_3'] = 0
       df_test_temp['Tx_pattern_4'] = 0
       df_test_temp['Tx_pattern_5'] = 0
       df_test_temp['Tx_pattern_6'] = 0
       df_test_temp['Tx_pattern_7'] = 0
       df_test_temp['Tx_pattern_8'] = 0
       df_test_temp['Tx_pattern_9'] = 0
       df_test_temp['Tx_pattern_10'] = 0

       text = 'Tx_pattern_{}'
       print(text.format(j))
       df_test_temp[text.format(j)] = 1

       df_test_temp = pd.get_dummies(df_test_temp, columns=['.HCT.CI','.PS24','.Analysed_Disease', '.Disease.satatus.all','R_CMVAntiB'])
       x_test_np = x_mapper.transform(df_test_temp).astype('float32')

       # Input the predictive probabilities of 1-year OS for 'Tx_pattern_j' into 'df_12_month_p' in test cohort
       surv = model.predict_surv_df(x_test_np)
       temp_surv_12 = surv.loc[[12]].transpose()
       temp_surv_12 = temp_surv_12.reset_index()
       temp_surv_12 = temp_surv_12[12]

       text = 'deepS_p_{}'
       print(text.format(j))
       temp_surv_12 = temp_surv_12.rename(text.format(j))
       df_12_month_p = pd.concat([df_12_month_p, temp_surv_12], axis=1)

# %%
df = analysed_df[analysed_df['.SCT.Year'] >= 2017]
df = df.reset_index()
df = pd.concat([df, df_12_month_p], axis=1)
df.to_csv(os.path.join(os.path.dirname(__file__),'data\\deep_srv_Tx10pattern.csv'),index = False)
