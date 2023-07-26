#!/usr/bin/python3
import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
import json
import scipy.optimize as opt
import scipy.stats as stats
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml import plot_utils
from matplotlib.backends.backend_pdf import PdfPages
from apply_model import apply_ML_model
from HypertritonRestframeBoost import get_df

pt_min = 4
pt_max = 5

BDT_cut = -1.7346938775510203

data_file = '../Data/DataTable_18_pass3.root'
data_tree = 'DataTable'
MC_file = f'../Data/SignalTable_20g7.root'
bkg_file = f'../Data/DataTable_18LS_pass3.root'
model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'
output_file = f"../Output/ML_Hypertriton_output_{pt_min}<pt<{pt_max}_all.pdf"


selected_data_hndl, dataH =  apply_ML_model(model_file, data_file, tree_name=data_tree, output_file=None, BDT_cut=BDT_cut, pt_cut=(pt_min, pt_max))
   

df = selected_data_hndl.get_data_frame()
print('Number of selected candidates:', len(df))
print(df.columns)

plot_utils.plot_distr(
        [selected_data_hndl], 
        column='m', 
        bins=100, 
        labels=['selected data'], 
        colors=['orange'], 
        density=True,
        fill=True, 
        histtype='step', 
        alpha=0.5,
        )
plt.show()
exit()
df_cos_theta = get_df(selected_data_hndl)
print(df_cos_theta[:10])


