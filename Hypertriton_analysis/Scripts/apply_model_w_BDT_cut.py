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
# import ROOT

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
start = np.amin(df["m"])
stop = np.amax(df["m"])



#TODO: figure out how to plot bkg/signal
'''
# plot selected candidates in hisotgram
c = ROOT.TCanvas()
rdf = ROOT.RDataFrame("df", "SelectedDataFrame.root")
inv_mass = rdf.Histo1D(("InvariantMassHistogram", "3<pt<6; m[GeV]", 100, start, stop),"m")
inv_mass.Draw()
c.Draw()

# fit a exponential background an a gaussian signal
# fit_func = ROOT.TF1("fit_func", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Exp(x, [4], [5])")
fit_func = ROOT.TF1("fit_func", "gaus(0) + [3]*expo(4)")
fit_func.SetParameter(1, 2.993)
fit_func.SetParameter(2, 0.005)
fit_func.SetNpx(1000)
inv_mass.Fit(fit_func)
pars = fit_func.GetParameters()
c.Draw()

bkg_func = ROOT.TF1("bkg_func", "[0]*expo(1)")
bkg_func.SetLineColor(ROOT.kGreen)
bkg_func.SetNpx(1000)
bkg_func.Draw("Same")
bkg_func.SetParameters(pars[0], pars[1], pars[2])
c.Draw()

signal_func = ROOT.TF1("signal_func", "gaus(0)")
signal_func.SetLineColor(ROOT.kBlue)
signal_func.SetNpx(1000)
signal_func.Draw("Same")
signal_func.SetParameters(pars[3], pars[4], pars[5])
c.Draw()

# legend = ROOT.TLegend(0.45, 0.65, 0.73, 0.85)
# legend.SetTextFont(72)
# legend.SetTextSize(0.04)
# legend.AddEntry(inv_mass, "Data", "l")
# legend.AddEntry(bkg_func, "Background fit", "l")
# legend.AddEntry(signal_func, "Signal fit", "l")
# legend.AddEntry(fit_func, "Global Fit", "l")
# legend.Draw("Same")
# c.Draw()

c.Print(f"SelectedDataDistribution_{pt_min}_{pt_max}.pdf", "Title: m")
print("Canvas saved as SelectedDataDistribution.pdf")
# print('Number of selected candidates:', len(df))
# print(df.columns)
'''
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

#Alternative method: ROOFit

# x = ROOT.RooRealVar("x", "invariantMass", 2.7, 3.1)
# mean = ROOT.RooReaLVar("mean", "meanOfGaussian", 2.992, 2.9, 3.0)
# sigma = ROOT.RooRealVar("sigma", "StdOfGaussian", 0.005, 0.001, 0.01)
# A = ROOT.RooRealVar("")


#cut on m?
# rdf.Filter()#TODO: FILTER EINFÜGEN!


# calculate cos theta*
m_He = 2.80839 #GeV
m_Pi = 0.139570 #GeV

df_with_cos_theta = get_df(selected_data_hndl)
print(df_with_cos_theta[:10])

# plot cos theta distribution of Hypertriton Candidates
bins = np.linspace(-1, 1, 6)
count_data, bins = np.histogram(df_with_cos_theta['cos_theta_beam'], bins=bins)

selected_data_hndl.set_data_frame(df_with_cos_theta)
selected_data_hndl.write_df_to_root_files("SelectedDataFrame", "df")
print("df saved to .root file.")

exit()
plt.stairs(count_data, bins, label=f'Selected Candidates \n {pt_min}<'+r'$p_T$'+f'<{pt_max}', color='orangered')
plt.ylabel('Counts')
plt.xlabel(r'$\cos{\theta *}_{beam}$')
plt.legend()
save_output_as_pdf(f'../Output/He3Distribution_CosThetaBeam_{pt_min}_pt_{pt_max}.pdf')  
print('output_saved')
plt.show()

