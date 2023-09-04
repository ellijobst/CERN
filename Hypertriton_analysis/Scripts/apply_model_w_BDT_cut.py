#!/usr/bin/python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import scipy.optimize as opt
import scipy.stats as stats
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml import plot_utils
from apply_model import apply_ML_model, save_output_as_pdf
from HypertritonRestframeBoost import get_df



# select pt cut
print("Pt_min:")
pt_min = input()
print("pt_max:")
pt_max = input()
print("MC?[y/n]")
MC_yn = input()


print(f'{pt_min}<pt<{pt_max}')

# get best BDT cut
with open(f"../Output/SignEff/signxeff_{pt_min}_pt_{pt_max}.json", 'r') as f:
        signxeff = json.load(f)

BDT_cuts = np.linspace(-15,10,100)
Best_BDT = BDT_cuts[signxeff.index(max(signxeff))]

print('Best BDT cut:', Best_BDT)

# load data
#MCs
if MC_yn == "y":
       MC = "MC"
       data_file = '../Data/SignalTable_B_20g7.root'
       data_tree = 'SignalTable'

# data
elif MC_yn == "n":
        MC = ""
        data_file = '../Data/DataTable_18_pass3.root'
        data_tree = "DataTable"

#model
model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'



# apply model to data (cut on ct is included in this function)
selected_data_hndl, dataH =  apply_ML_model(model_file, data_file, tree_name=data_tree, output_file=None, BDT_cut=Best_BDT, pt_cut=(pt_min, pt_max))


# get BDT cut efficiency
# NOTE: this is just for MC files!!!

before=dataH.get_n_cand()
after=selected_data_hndl.get_n_cand()

print("BDT CUT EFFICIENCY:", after/before)


# plot signal (for crosscheck)
'''plot_utils.plot_distr(
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
'''

# Alternative method: Convert data to a dictionary with numpy arrays
# data = {key: df[key].values for key in df.keys()}
# rdf = ROOT.RDF.FromNumpy(data)


# calculate cos theta*
m_He = 2.80839 #GeV
m_Pi = 0.139570 #GeV

# apply hypertriton restframe boost to get cos(theta*) columns
df_with_cos_theta = get_df(selected_data_hndl)
print(df_with_cos_theta[:10])

# save as root file
selected_data_hndl.set_data_frame(df_with_cos_theta)
selected_data_hndl.write_df_to_root_files(f"SelectedDataFrame{MC}_{pt_min}_pt_{pt_max}", "df")
print("df saved to .root file.")



