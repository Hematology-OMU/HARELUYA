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
import random
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.model_selection import KFold
import os
karafuto = 1808

#%%
dataset_fp = os.path.join(os.path.dirname(__file__),'data\\analysed_df.csv')
analysed_df = pd.read_csv(dataset_fp, encoding = 'cp932')
analysed_df = analysed_df[['.SCT.Year','.Age','.HCT.CI','.PS24','.Analysed_Disease','.Disease.satatus.all',
                            'Tx_pattern','R_CMVAntiB','.MonthOS','.OS']]
analysed_df = analysed_df.dropna()
analysed_df.head()
analysed_df = pd.get_dummies(analysed_df, columns=['.HCT.CI','.PS24','.Analysed_Disease', '.Disease.satatus.all','Tx_pattern','R_CMVAntiB'])

#%%
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
df_CV = analysed_df[analysed_df['.SCT.Year'] <= 2016]

#%%
random.seed(karafuto)

#%%
def CrossValidation():
    n_seed = random.randint(1,1000)
    nodes =  [8, 16, 32, 64, 128]
    n_nodes = random.choice(nodes)
    n_layers = random.randint(1,6)
    bins = [16, 32, 64, 128, 256, 512]
    batch_size = random.choice(bins)
    dropout = random.random()
    batch_norm = True
    lr = 0.001 # defalt value in pytorch packege
    frac=0.2
    epochs = 300
    print('(batch_size)', batch_size)
    print('(n_nodes)', n_nodes)
    print('(n_layers)', n_layers)
    print('-----')
    if n_layers == 1:
        num_nodes = [n_nodes]
    if n_layers == 2:
        num_nodes = [n_nodes,n_nodes]
    if n_layers == 3:
        num_nodes = [n_nodes,n_nodes,n_nodes]
    if n_layers == 4:
        num_nodes = [n_nodes,n_nodes,n_nodes,n_nodes]
    if n_layers == 5:
        num_nodes = [n_nodes,n_nodes,n_nodes,n_nodes,n_nodes]
    if n_layers == 6:
        num_nodes = [n_nodes,n_nodes,n_nodes,n_nodes,n_nodes,n_nodes]

    np.random.seed(n_seed)
    _ = torch.manual_seed(n_seed)
    kf = KFold(n_splits=5, shuffle=True, random_state=karafuto)
    c_index = np.array([])
    for train, test in kf.split(df_CV):
        train_df = df_CV.iloc[train]
        val_df = train_df.sample(frac=frac, random_state = karafuto)
        train_df = train_df.drop(val_df.index)
        test_df = df_CV.iloc[test]

        x_train_df = x_mapper.fit_transform(train_df).astype('float32')
        y_train_df = get_target(train_df)
        x_val_df = x_mapper.transform(val_df).astype('float32')
        y_val_df = get_target(val_df)
        val = x_val_df, y_val_df
        x_test_df = x_mapper.transform(test_df).astype('float32')
        durations_test_df, events_test_df = get_target(test_df)

        in_features = x_train_df.shape[1]
        out_features = 1
        net = tt.practical.MLPVanilla(in_features, num_nodes, out_features, batch_norm,
                                dropout, activation=torch.nn.ReLU, output_bias=False)
        model = CoxPH(net, tt.optim.Adam())
        model.optimizer.set_lr(lr)
        callbacks = [tt.callbacks.EarlyStopping()]
        verbose = False

        log = model.fit(x_train_df, y_train_df, batch_size, epochs, callbacks, verbose,
                    val_data=val, val_batch_size=batch_size)

        _ = model.compute_baseline_hazards()

        surv_performacne = model.predict_surv_df(x_test_df)
        ev = EvalSurv(surv_performacne, durations_test_df, events_test_df, censor_surv='km')
        ev.concordance_td()

        c_index = np.append(c_index, ev.concordance_td())

    temp_result = pd.DataFrame({"n_seed" : n_seed,
                            "n_nodes" : n_nodes,
                            "n_layers": n_layers,
                            "batch_size": batch_size,
                            "dropout": dropout,
                            "c-index": np.mean(c_index)},index=[1])
    return temp_result

#%%
RS_result = pd.DataFrame(index = [], columns=["n_seed","n_nodes","n_layers","batch_size","dropout","c-index"])
for j in range(0, 100):
    print('(cycle number)', j)
    temp = CrossValidation()
    RS_result = pd.concat([RS_result, temp], axis=0, ignore_index=True)
    print(temp)

    if j%20 == 19:
        RS_result.to_csv(os.path.join(os.path.dirname(__file__),'data\\DS_RandomSearch_test.csv'),index = False)