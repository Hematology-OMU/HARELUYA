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
                            'Tx_pattern','R_CMVAntiB','.MonthOS','.OS']]                        
analysed_df = analysed_df.dropna()
df_train = analysed_df[analysed_df['.SCT.Year'] <= 2016]
df_val = df_train.sample(frac=frac, random_state = karafuto)
df_train = df_train.drop(df_val.index)
df_test = analysed_df[analysed_df['.SCT.Year'] >= 2017]

# %%
cols_standardize = ['.Age']
cols_leave = ['.Analysed_Disease_ALL', '.Analysed_Disease_AML',
       '.Analysed_Disease_ATL', '.Analysed_Disease_MDS',
       '.Analysed_Disease_ML', '.Analysed_Disease_MPN', 
       '.Disease.satatus.all_CR', '.Disease.satatus.all_nonCR', 'Tx_pattern_1',
       'Tx_pattern_2', 'Tx_pattern_3', 'Tx_pattern_4', 'Tx_pattern_5',
       'Tx_pattern_6', 'Tx_pattern_7', 'Tx_pattern_8', 'Tx_pattern_9',
       'Tx_pattern_10', '.HCT.CI_0', '.HCT.CI_1', '.HCT.CI_2', '.HCT.CI_3', '.PS24_0',
       '.PS24_1', '.PS24_2', '.PS24_3', '.PS24_4', 'R_CMVAntiB_no', 'R_CMVAntiB_yes']

standardize = [([col], StandardScaler()) for col in cols_standardize]
leave = [(col, None) for col in cols_leave]
x_mapper = DataFrameMapper(standardize + leave)
get_target = lambda df: (df['.MonthOS'].values, df['.OS'].values)

# %%
df_train = pd.get_dummies(df_train, columns=['.HCT.CI','.PS24','.Analysed_Disease', '.Disease.satatus.all','Tx_pattern','R_CMVAntiB'])
df_val = pd.get_dummies(df_val, columns=['.HCT.CI','.PS24','.Analysed_Disease', '.Disease.satatus.all','Tx_pattern','R_CMVAntiB'])

x_df_train = x_mapper.fit_transform(df_train).astype('float32')
y_df_train = get_target(df_train)
x_df_val = x_mapper.transform(df_val).astype('float32')
y_df_val = get_target(df_val)
val = x_df_val, y_df_val

# %%
# Develop DeepSurv model in train & val data
in_features = x_df_train.shape[1]
num_nodes = [n_nodes,n_nodes,n_nodes,n_nodes,n_nodes,n_nodes]
out_features = 1
batch_norm = batch_norm
dropout = dropout
output_bias = False

net = tt.practical.MLPVanilla(in_features, num_nodes, out_features, batch_norm,
                              dropout, activation=torch.nn.ReLU, output_bias=output_bias)
model = CoxPH(net, tt.optim.Adam())
batch_size = batch_size
model.optimizer.set_lr(lr)
epochs = epochs
callbacks = [tt.callbacks.EarlyStopping()]
verbose = True

#%%time
log = model.fit(x_df_train, y_df_train, batch_size, epochs, callbacks, verbose,
                val_data=val, val_batch_size=batch_size)
_ = model.compute_baseline_hazards()

# %%
# Create 6 kinds of permutated test data by shuffling variable
temp = df_test['.Age'].sample(frac = 1).reset_index(drop = True)
df_test_perm_Age = pd.concat([df_test.drop(".Age", axis=1).reset_index(drop = True), temp], axis = 1)

temp = df_test['.HCT.CI'].sample(frac = 1).reset_index(drop = True)
df_test_perm_HCTCI = pd.concat([df_test.drop(".HCT.CI", axis=1).reset_index(drop = True), temp], axis = 1)

temp = df_test['.PS24'].sample(frac = 1).reset_index(drop = True)
df_test_perm_PS = pd.concat([df_test.drop(".PS24", axis=1).reset_index(drop = True), temp], axis = 1)

temp = df_test['.Analysed_Disease'].sample(frac = 1).reset_index(drop = True)
df_test_perm_Dis = pd.concat([df_test.drop(".Analysed_Disease", axis=1).reset_index(drop = True), temp], axis = 1)

temp = df_test['.Disease.satatus.all'].sample(frac = 1).reset_index(drop = True)
df_test_perm_Dis_St = pd.concat([df_test.drop(".Disease.satatus.all", axis=1).reset_index(drop = True), temp], axis = 1)

temp = df_test['Tx_pattern'].sample(frac = 1).reset_index(drop = True)
df_test_perm_Txpattern = pd.concat([df_test.drop("Tx_pattern", axis=1).reset_index(drop = True), temp], axis = 1)

temp = df_test['R_CMVAntiB'].sample(frac = 1).reset_index(drop = True)
df_test_perm_RCMV = pd.concat([df_test.drop("R_CMVAntiB", axis=1).reset_index(drop = True), temp], axis = 1)

# %%
# Predict 1-year OS for each test patient
def create_df_1yOS_prediction_df(model,df):
    df_test = df
    df_dum = pd.get_dummies(df, columns=['.HCT.CI','.PS24','.Analysed_Disease', '.Disease.satatus.all','Tx_pattern','R_CMVAntiB'])
    x_df_dum = x_mapper.transform(df_dum).astype('float32')
    surv = model.predict_surv_df(x_df_dum)
    temp_surv_12 = surv.loc[[12]].transpose()
    temp_surv_12 = temp_surv_12.reset_index()
    temp_surv_12 = temp_surv_12[12]
    temp_surv_12 = temp_surv_12.rename('predicted_1y_OS')
    df_test = df_test.reset_index()
    df_test = df_test.drop("index", axis=1)
    df_test = pd.concat([df_test, temp_surv_12], axis=1)
    return df_test

# %%
df_base_test = create_df_1yOS_prediction_df(model, df_test)
df_test_perm_Age = create_df_1yOS_prediction_df(model, df_test_perm_Age)
df_test_perm_HCTCI = create_df_1yOS_prediction_df(model, df_test_perm_HCTCI)
df_test_perm_PS = create_df_1yOS_prediction_df(model, df_test_perm_PS)
df_test_perm_Dis = create_df_1yOS_prediction_df(model, df_test_perm_Dis)
df_test_perm_Dis_St = create_df_1yOS_prediction_df(model, df_test_perm_Dis_St)
df_test_perm_Txpattern = create_df_1yOS_prediction_df(model, df_test_perm_Txpattern)
df_test_perm_RCMV = create_df_1yOS_prediction_df(model, df_test_perm_RCMV)
#%%
df_base_test.to_csv(os.path.join(os.path.dirname(__file__),'data\\deepS_test_base_prediction.csv'),index = False)
df_test_perm_Age.to_csv(os.path.join(os.path.dirname(__file__),'data\\deep_srv_perm_Age_prediction.csv'),index = False)
df_test_perm_HCTCI.to_csv(os.path.join(os.path.dirname(__file__),'data\\deep_srv_perm_HCTCI_prediction.csv'),index = False)
df_test_perm_PS.to_csv(os.path.join(os.path.dirname(__file__),'data\\deep_srv_perm_PS_prediction.csv'),index = False)
df_test_perm_Dis.to_csv(os.path.join(os.path.dirname(__file__),'data\\deep_srv_perm_Dis_prediction.csv'),index = False)
df_test_perm_Dis_St.to_csv(os.path.join(os.path.dirname(__file__),'data\\deep_srv_perm_DisSt_prediction.csv'),index = False)
df_test_perm_Txpattern.to_csv(os.path.join(os.path.dirname(__file__),'data\\deep_srv_perm_TxP_prediction.csv'),index = False)
df_test_perm_RCMV.to_csv(os.path.join(os.path.dirname(__file__),'data\\deep_srv_perm_RCMV_prediction.csv'),index = False)