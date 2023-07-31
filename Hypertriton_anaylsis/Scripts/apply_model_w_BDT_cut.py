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
import ROOT

def signal_and_bkg_fit(m, loc_sig, scale_sig, A_sig, loc_bkg, scale_bkg, A_bkg):
    return  A_sig*stats.norm.pdf(m, loc=loc_sig, scale=scale_sig)+A_bkg*stats.expon.pdf(m, loc=loc_bkg, scale=scale_bkg)

# select pt cut
print("Pt_min:")
pt_min = input()
print("pt_max:")
pt_max = input()

# pt_min = 3
# pt_max = 4
print(f'{pt_min}<pt<{pt_max}')

# get best BDT cut
with open(f"../Output/SignEff/significance_{pt_min}_pt_{pt_max}.json", 'r') as f:
        significance = json.load(f)
with open(f"../Output/SignEff/efficiency_{pt_min}_pt_{pt_max}.json", 'r') as f:
        efficiency = json.load(f)

signxeff = [x*y for x,y in zip(significance, efficiency)]
BDT_cuts = np.linspace(-15,10,100)
Best_BDT = BDT_cuts[signxeff.index(max(signxeff))]
print('Best BDT cut:', Best_BDT)

# load data
data_file = '../Data/DataTable_18_pass3.root'
data_tree = 'DataTable'
model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'

# apply model to data (cut on ct is included in this function)
selected_data_hndl, dataH =  apply_ML_model(model_file, data_file, tree_name=data_tree, output_file=None, BDT_cut=Best_BDT, pt_cut=(pt_min, pt_max))
   
# turn into df
df = selected_data_hndl.get_data_frame()
print('Number of selected candidates:', len(df))
print(df.columns)

# plot signal
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


# Convert data to a dictionary with numpy arrays
data = {key: df[key].values for key in df.keys()}
print(data)

rdf = ROOT.RDF.FromNumpy(data)
rdf.Display().Print()

# save as .root file
rdf.Snapshot('df', 'SelectedData_{pt_min}_pt_{pt_max}.root')



# plotting with scipy
'''

plt.figure()
bins_m = np.linspace(np.amin(df['m']), np.amax(df['m']), 101)
count_m, bins_m = np.histogram(df['m'], bins=bins_m)
centered_bins = (bins_m[:-1] + bins_m[1:]) / 2

m_Hyp = 2.991
p0=[m_Hyp, 3*0.0017, 50, 3, 1, 5]
# p0=StartParamsNorm(df)
optimizedParameters, pcov = opt.curve_fit(signal_and_bkg_fit, centered_bins, count_m, p0=p0) #TODO:das funktioniert nicht !

plt.figure()
plot_utils.plot_distr(selected_data_hndl, column=['m'], bins=100)
plt.plot(np.linspace(bins_m[0], bins_m[-1], 100), signal_and_bkg_fit(np.linspace(bins_m[0], bins_m[-1], 100), *optimizedParameters), label='Fit')
save_output_as_pdf(f'../Output/Distribution_m_{pt_min}_pt_{pt_max}.pdf')
plt.show()
plt.close()


exit()
'''

# calculate cos theta*
m_He = 2.80839 #GeV
m_Pi = 0.139570 #GeV

df_with_cos_theta = get_df(selected_data_hndl)
print(df_with_cos_theta[:10])

# plot cos theta distribution of Hypertriton Candidates
bins = np.linspace(-1, 1, 6)
count_data, bins = np.histogram(df_with_cos_theta['cos_theta_beam'], bins=bins)

plt.stairs(count_data, bins, label=f'Selected Candidates \n {pt_min}<'+r'$p_T$'+f'<{pt_max}', color='orangered')
plt.ylabel('Counts')
plt.xlabel(r'$\cos{\theta *}_{beam}$')
plt.legend()
save_output_as_pdf(f'../Output/He3Distribution_CosThetaBeam_{pt_min}_pt_{pt_max}.pdf')  
print('output_saved')
plt.show()

