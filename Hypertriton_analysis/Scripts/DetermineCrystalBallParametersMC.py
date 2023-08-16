import ROOT
import numpy as np
import json

#-------------SETTINGS------------------
# matter?
matter = "false"

#output file
Output= "./Output/RooFitCB3pt6_AM"

#pt cuts
pt_min = 3
pt_max = 6


print("------------SETTINGS---------------")
print(f" matter = {matter}")
print(f"{pt_min}<pt<{pt_max}")
print("------------------------------------")


#--------------FUNCTIONS--------------

def CreatePlots(i):
    # create canvas
    c1 = ROOT.TCanvas()
    mframe = m.frame(Title="m histogram + Fit")
    
   
    bins = np.linspace(-1,1,8)

    # create histogram
    inv_mass_df = rdf.Filter(f"{bins[i]} < cos_theta_beam and cos_theta_beam < {bins[i+1]} and Matter == {matter} ")
    inv_mass = inv_mass_df.Histo1D(("InvariantMassHistogram", f"{pt_min}"+"<p_{T}<"+f"{pt_max}, {str(bins[i])[:5]}"+"< cos(#theta_{beam})<"+f"{str(bins[i+1])[:5]}; m[GeV]", 80, 2.96, 3.04),"m")
    inv_mass_copy = inv_mass.Clone()
    
    # get RooDataHist
    # NOTE: use the .GetPtr() since inv_mass is just a ptr to a TH1D and not the TH1D itself!!
    dh = ROOT.RooDataHist("dh", "dataset with m", m, inv_mass.GetPtr())
    
    # draw histogram using again 80 bins
    # NOTE: plotting options: https://root.cern/doc/master/group__Plotting.html 
    dh.plotOn(mframe, MarkerSize=0, Binning=80)
    
    # fit pdf and plot fit 
    fit_func.fitTo(dh)
    fit_func.plotOn(mframe, ROOT.RooFit.Precision(1e-5),  ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineColor(ROOT.kRed)) 
    
    # plot signal and background
    fit_func.plotOn(mframe, ROOT.RooFit.Components("Expo"), ROOT.RooFit.LineColor(ROOT.kGreen+1) , LineWidth=2)
    fit_func.plotOn(mframe, ROOT.RooFit.Components("gauss"), ROOT.RooFit.LineColor(ROOT.kBlue+1) , LineWidth=2)
    
    legend = ROOT.TLegend(.12, .8, .49, .935)
    legend.SetBorderSize(1)
    legend.AddEntry(dh, "Data", "le")
    legend.AddEntry(fit_func, "Fit Crystal Ball Shape", "l").SetLineColor(ROOT.kRed)

    legend.SetTextSize(0.025)
    
    # Add Fit Parameters to canvas 
    pave = ROOT.TPaveText(0.5, 0.8, 0.77, 0.935, "NDC")
    pave.SetBorderSize(1)
    pave.SetFillColor(ROOT.kWhite)
    pave.SetTextFont(42)

    t0 = pave.AddText(f"#mu = {str(m0.getVal())[:10]} #pm {str(m0.getError())[:10]}")
    t1 = pave.AddText(f"#sigma = {str(sigma.getVal())[:10]} #pm {str(sigma.getError())[:10]}")
    t2 = pave.AddText("#alpha_{L} "+f"= {str(alphal.getVal())[:10]} #pm {str(alphal.getError())[:10]}")
    t3 = pave.AddText("#alpha_{R}"+f" = {str(alphar.getVal())[:10]} #pm {str(alphar.getError())[:10]}")
    t4 = pave.AddText("n_{L}"+f"=  {str(nl.getVal())[:10]} #pm {str(nl.getError())[:10]}")
    t5= pave.AddText("n_{R}"+f"=  {str(nr.getVal())[:10]} #pm {str(nr.getError())[:10]}")
 

    pave2 = ROOT.TPaveText(0.78, 0.75, 0.95, 0.935, "NDC")
    pave2.AddText("#Chi^{2}"+f" = {str(mframe.chiSquare())[:10]}")
    pave2.SetBorderSize(1)
    pave2.SetFillColor(ROOT.kWhite)
    pave2.SetTextFont(42)
    
    
    mframe.Draw()
    mframe.GetYaxis().SetRangeUser(0,20000)
    legend.Draw("same")
       
    pave.Draw("same")
    pave2.Draw("same")
    

    c1.Draw()

    c1.SaveAs(f"{Output}.pdf")
    
#     return Nhyp, error
    return nl.getVal(), nr.getVal(), alphal.getVal(), alphar.getVal(), nl.getError(), nr.getError(), alphal.getError(), alphar.getError()

#--------------CODE-----------------------
if __name__ == "__main__":
    
    
    # define variables
    # NOTE: we dont need amplitudes just a fraction since this is pdfs
    m = ROOT.RooRealVar("m", "m [GeV]", 2.96, 3.04)
    m0 = ROOT.RooRealVar("m0", "mean of CB", 2.992,  2.990, 2.994)
    sigma = ROOT.RooRealVar("sigma", "sigma of CB", 0.001,  0.0001, 0.0015)

    
    alphal = ROOT.RooRealVar("alphal", "alpha L", 2.05,1, 2.1)
    nl = ROOT.RooRealVar("nl", "n L", 4.7,  3.9, 7.2)
    alphar = ROOT.RooRealVar("alphar", "alpha right", 2.08, 1, 2.2)
    nr = ROOT.RooRealVar("nr", "n right", 3.7, 2.9 , 6.5)


    # define pdf model
    fit_func = ROOT.RooCrystalBall("fit_func", "Crystal Ball PDF", m, m0, sigma, alphal, nl, alphar, nr)
    # gaus = ROOT.RooGaussian("gaus", "gaus PDF", m, m0, sigma)


    # import data
    rdf = ROOT.RDataFrame("df", f"./SelectedDataFrames/SelectedDataFrameMC_{pt_min}_pt_{pt_max}.root")

    # i goes from 0 to 6, since there are 7 bins
    bins = np.linspace(-1,1,8)
        

    # plots cos(theta*) distribution-----------------------------------------------
    c2 = ROOT.TCanvas()
    cos_theta = rdf.Filter(f"Matter == {matter}").Histo1D(("CosThetaHistogram", "cos(theta*)_wrt_beam", 100, -1, 1),"cos_theta_beam")
    cos_theta.Draw()
    c2.Draw()
    c2.SaveAs(f"{Output}.pdf(")


    # create plots with signal fit-------------------------------------------------
    results = [CreatePlots(i) for i in range(7)]
    
    print(results)
    nl = [results[i][0] for i in range(len(results))]
    nr = [results[i][1] for i in range(len(results))]
    alphal = [results[i][2] for i in range(len(results))]
    alphar = [results[i][3] for i in range(len(results))]
    
    nl_err = [results[i][4] for i in range(len(results))]
    nr_err = [results[i][5] for i in range(len(results))]
    alphal_err = [results[i][6] for i in range(len(results))]
    alphar_err = [results[i][7] for i in range(len(results))]
    

    if matter == "true":
        matter3 = "M"
    elif matter == "false":
        matter3 = "AM"
    with open(f"FixedCBParams{pt_min}pt{pt_max}_{matter3}.json", 'w') as f:
        json.dump(results, f, indent=2) 

    print("results saved.")

    
        
    c4 = ROOT.TCanvas()
    centered_bins = (bins[:-1] + bins[1:]) / 2


    gr0 = ROOT.TGraphErrors(7, centered_bins, np.array(nl), np.array([1/7]*7), np.array(nl_err))
    gr1 = ROOT.TGraphErrors(7, centered_bins, np.array(nr), np.array([1/7]*7), np.array(nr_err))
    gr1.SetLineColor(ROOT.kRed)
#     gr1.SetLineWidth(1)
    gr0.GetYaxis().SetRangeUser(5.5,8)
    
    gr0.SetTitle("Parameters from CB MC Fit")
    gr0.GetXaxis().SetTitle("cos(#theta_{beam})")
    
    ROOT.gStyle.SetTitleFontSize(0.045)

    gr0.Draw()
    gr1.Draw("Same")
    
    legend4 = ROOT.TLegend(.7, .7, .85, .85)
    legend4.SetBorderSize(0)
    legend4.AddEntry(gr0, "n_{L}", "l")
    legend4.AddEntry(gr1, "n_{R}", "l").SetLineColor(ROOT.kRed)

    legend4.SetTextSize(0.025)
    legend4.Draw()
    
    c4.Draw()
    c4.SaveAs(f"{Output}.pdf")
    
    
    c5 = ROOT.TCanvas()


    gr0 = ROOT.TGraphErrors(7, centered_bins, np.array(alphal), np.array([1/7]*7), np.array(alphal_err))
    gr1 = ROOT.TGraphErrors(7, centered_bins, np.array(alphar), np.array([1/7]*7), np.array(alphar_err))
    gr1.SetLineColor(ROOT.kRed)
    
    gr0.SetTitle("Parameters from CB MC Fit")
    gr0.GetXaxis().SetTitle("cos(#theta_{beam})")
    
    ROOT.gStyle.SetTitleFontSize(0.045)

    gr0.Draw()
    gr1.Draw("Same")
    
    legend5 = ROOT.TLegend(.7, .7, .85, .85)
    legend5.SetBorderSize(0)
    legend5.AddEntry(gr0, "#alpha_{L}", "l")
    legend5.AddEntry(gr1, "#alpha_{R}", "l").SetLineColor(ROOT.kRed)

    legend5.SetTextSize(0.025)
    legend5.Draw()
    
    c5.SaveAs(f"{Output}.pdf)")
    