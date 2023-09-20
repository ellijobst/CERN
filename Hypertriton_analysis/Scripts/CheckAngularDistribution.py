import ROOT
import json
import numpy as np

def Procedure(cmin, nbins, inplane, model):
    EP = "IP" if inplane else "OP"

    
    # fit the difference between matter/antimatter
    with open(Inputdir+f"NhypCorrected3.0pt6.4centr{cmin}_{EP}_M.json", 'r') as f:
        nhyp_m = json.load(f) 
    with open(Inputdir+f"NhypCorrectedErr3.0pt6.4centr{cmin}_{EP}_M.json", 'r') as f:
        err = json.load(f)
    with open(Inputdir+f"NhypCorrected3.0pt6.4centr{cmin}_{EP}_AM.json", 'r') as f:
        nhyp_am = json.load(f)
    with open(Inputdir+f"NhypCorrectedErr3.0pt6.4centr{cmin}_{EP}_AM.json", 'r') as f:
        err_am = json.load(f)
    
    
    diff = [x-y for x,y in zip(nhyp_m, nhyp_am)]
    diff_err = [np.sqrt(x**2+y**2) for x,y in zip(err, err_am)]

    bins = np.linspace(-1,1,nbins+1)
    centered_bins = (bins[:-1] + bins[1:]) / 2
    
    c = ROOT.TCanvas()
    h = ROOT.TH1D("hist", f"M-AM, Inplane:{inplane}, Cmin:{cmin}, Model:{model}; cos(#theta *); Counts/bin", nbins, -1, 1)
    h.FillN(nbins, np.array(centered_bins), np.array(diff) )
    for i in range(1,nbins+1,1):
        h.SetBinError(i, diff_err[i-1])
    
        #define fit function
    fit_func = ROOT.TF1("fit_func", model)
    if model=="pol2":
        fit_func.FixParameter(1,0)


    # fit, here we use a chi2 fit, 
    h.Fit(fit_func)

    h.GetYaxis().SetTitle("Counts per bin")


    h.Draw()
    h.GetFunction("fit_func").SetLineColor(ROOT.kRed)

    legend = ROOT.TLegend(.12, .85, .3, .935)
    legend.SetTextSize(0.025)
    legend.SetBorderSize(1)
    legend.AddEntry(h, "matter-antimatter", "l")
    legend.AddEntry(fit_func, "linear Fit", "l").SetLineColor(ROOT.kRed)

    legend.Draw("same")
    pars = fit_func.GetParameters()
    errors = fit_func.GetParErrors()
    chisq = fit_func.GetChisquare()
    
    pave2 = ROOT.TPaveText(0.32, 0.8, 0.55, 0.935, "NDC")
    pave2.AddText("#Chi^{2}"+f" = {str(chisq)[:10]}")
    # degrees of freedom is number of bins minus number of fit parameters
    pave2.AddText(f"Probabilty = {str(ROOT.TMath.Prob(chisq, nbins-1))[:10]}")
    pave2.AddText(f"slope: {str(pars[int(model[-1])])[:5]}#pm{str(errors[int(model[-1])])[:5]}")
    pave2.SetBorderSize(1)
    pave2.SetFillColor(ROOT.kWhite)
    pave2.SetTextFont(42)

    pave2.Draw()

    c.Draw()
    c.SaveAs(Outputdir + "Polarization.pdf")

if __name__ == "__main__":
    # Input files directory
    Inputdir = "../Output/TESTTEST/"
    Outputdir = "../Output/TESTTEST/"

    # centrality cuts
    cmin_list = [10]#[0,5,10]
    
    # models
    model_list=["pol0", "pol1", "pol2"]

    # event plane cut
    inplane_list =[True, False]

    #number of cos(theta*) bins
    nbins = 5

    c = ROOT.TCanvas()
    pave = ROOT.TPaveText(0.2, 0.2, 0.8, 0.8, "NDC")
    pave.AddText("Output for:")
    pave.AddText("Centrality >0,5,10")
    pave.AddText("In plane and Out of plane")
    pave.AddText("pol0, pol1 and pol2")
    pave.Draw()
    c.Draw()
    c.SaveAs("Test/Polarization.pdf(")

    nbins = 5

    for cmin in cmin_list:
        print(f"Centrality > {cmin}")
        for inplane in inplane_list:
            print(f"Inplane:{inplane}")
            for model in model_list:
                print(f"Model:{model}")
                Procedure(cmin, nbins, inplane, model)
                
    c = ROOT.TCanvas()
    c.Draw()
    c.SaveAs(Outputdir + "AngularDistribution.pdf)")