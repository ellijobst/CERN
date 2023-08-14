import ROOT
import numpy as np
import json

# matter?
matter = "true"

#number of bins
nbins = 7

#centrality cut
cmin = 0

#pt cut
pt_min = 3
pt_max = 6


#output file: Method, PDF Form, pT range,centrality cut matter?
Output= f"./Output/RooFitCBExpo{pt_min}pt{pt_max}centr{cmin}_M"


print("------------SETTINGS---------------")
print(f"  Matter:             {matter}")
print(f"  Selected pT range:  {pt_min}<pt<{pt_max}")
print(f"  Centrality Cut:     centrality > {cmin}")
print(f"  Number of bins:     {nbins}")
print(f"  Event Plane Cut: -- ")
print("------------------------------------")

#---------------------------------------------------------------------------------------------
#parameters from Fitting MCs with RooCrystalBall

if matter == "true":
        matter3 = "M"
elif matter == "false":
        matter3 = "AM"
with open(f"./CBFixedMCParameters/FixedCBParams{pt_min}pt{pt_max}_{matter3}.json", 'r') as f:
        results = json.load(f)
        
nl_list = [results[i][0] for i in range(len(results))]
nr_list = [results[i][1] for i in range(len(results))]
alphal_list = [results[i][2] for i in range(len(results))]
alphar_list = [results[i][3] for i in range(len(results))]

nl_err = [results[i][4] for i in range(len(results))]
nr_err = [results[i][5] for i in range(len(results))]
alphal_err = [results[i][6] for i in range(len(results))]
alphar_err = [results[i][7] for i in range(len(results))]


#------------------------------------------------------------------------------------
def CreatePlots(i):
    # fix parameters based on MCS
    alphal = ROOT.RooRealVar("alphal", "alpha L", alphal_list[i], alphal_list[i], alphal_list[i])
    nl = ROOT.RooRealVar("nl", "n L", nl_list[i], nl_list[i], nl_list[i])
    alphar = ROOT.RooRealVar("alphar", "alpha right", alphar_list[i], alphar_list[i], alphar_list[i])
    nr = ROOT.RooRealVar("nr", "n right", nr_list[i], nr_list[i], nr_list[i])
    
    CB = ROOT.RooCrystalBall("CB", "Crystal Ball PDF", m, m0, sigma, alphal, nl, alphar, nr)
    expo = ROOT.RooExponential("Expo", "expo pdf", m, lam )
    fit_func = ROOT.RooAddPdf("fit_func", "Fit_function_pdf", ROOT.RooArgList(CB, expo), ROOT.RooArgList(frac))
    
    # create canvas
    c1 = ROOT.TCanvas()
    mframe = m.frame(Title="m histogram + Fit")
    
   
    bins = np.linspace(-1,1,nbins+1)

    # create histogram
    inv_mass_df = rdf.Filter(f"{bins[i]} < cos_theta_beam and cos_theta_beam < {bins[i+1]} and Matter == {matter}  and centrality > {cmin} ")
    inv_mass = inv_mass_df.Histo1D(("InvariantMassHistogram", f"{pt_min}"+"<p_{T}<"+f"{pt_max}, {str(bins[i])[:5]}"+"< cos(#theta_{beam})<"+f"{str(bins[i+1])[:5]}; m[GeV]", 80, 2.96, 3.04),"m")
    inv_mass_copy = inv_mass.Clone()
    
    # get RooDataHist
    # NOTE: use the .GetPtr() since inv_mass is just a ptr to a TH1D and not the TH1D itself!!
    dh = ROOT.RooDataHist("dh", "dataset with m", m, inv_mass.GetPtr())
    
    # draw histogram using again 80 bins
    # NOTE: plotting options: https://root.cern/doc/master/group__Plotting.html 
    # important: in order to be able to access the object later on, it needs to be explicitly named!!
    dh.plotOn(mframe, ROOT.RooFit.Name("dh"), ROOT.RooFit.MarkerSize(0))#, Binning=80)
    
    # fit pdf and plot fit 
    fit_func.fitTo(dh)
    fit_func.plotOn(mframe, ROOT.RooFit.Precision(1e-5),  ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("fit_func")) 
    
    # plot signal and background
    fit_func.plotOn(mframe, ROOT.RooFit.Components("Expo"), ROOT.RooFit.LineColor(ROOT.kGreen+1) , LineWidth=2)
    fit_func.plotOn(mframe, ROOT.RooFit.Components("CB"), ROOT.RooFit.LineColor(ROOT.kBlue+1) , LineWidth=2)
    
    legend = ROOT.TLegend(.12, .8, .49, .935)
    legend.SetBorderSize(1)
    legend.AddEntry(dh, "Data", "le")

    # fit_func.SetLineColor(ROOT.kRed)
    legend.AddEntry(fit_func, "Fit = Signal + Background", "l").SetLineColor(ROOT.kRed)
    legend.AddEntry(expo, "Background = (1-a)*exp(x #lambda))", "l").SetLineColor(ROOT.kGreen+1)
    legend.AddEntry(CB, "Signal = a*CB(x, #mu, #sigma, n_{L}, #alpha_{L}, n_{R}, #alpha_{R})", "l").SetLineColor(ROOT.kBlue)

    legend.SetTextSize(0.025)
    
    # Add Fit Parameters to canvas 
    pave = ROOT.TPaveText(0.5, 0.75, 0.77, 0.935, "NDC")
    pave.SetBorderSize(1)
    pave.SetFillColor(ROOT.kWhite)
    pave.SetTextFont(42)

    t0 = pave.AddText(f"#mu = {str(m0.getVal())[:10]} #pm {str(m0.getError())[:10]}")
    t1 = pave.AddText(f"#sigma = {str(sigma.getVal())[:10]} #pm {str(sigma.getError())[:10]}")
#     t2 = pave.AddText("#alpha_{L} "+f"= {str(alphal.getVal())[:10]} #pm {str(alphal.getError())[:10]}")
#     t3 = pave.AddText("#alpha_{R}"+f" = {str(alphar.getVal())[:10]} #pm {str(alphar.getError())[:10]}")
#     t4 = pave.AddText("n_{L}"+f"=  {str(nl.getVal())[:10]} #pm {str(nl.getError())[:10]}")
#     t5= pave.AddText("n_{R}"+f"=  {str(nr.getVal())[:10]} #pm {str(nr.getError())[:10]}")
    t6 = pave.AddText(f"a = {str(frac.getVal())[:10]} #pm {str(frac.getError())[:10]}")
    t7 = pave.AddText(f"#lambda = {str(lam.getVal())[:10]} #pm {str(lam.getError())[:10]}")
    
    chisq = mframe.chiSquare("fit_func", "dh", 4)
        
    pave2 = ROOT.TPaveText(0.78, 0.8, 0.95, 0.935, "NDC")
    pave2.AddText("#Chi^{2}/ndf"+f" = {str(chisq)[:10]}")
    # degrees of freedom is number of bins minus number of fit parameters
    pave2.AddText(f"Probabilty = {str(ROOT.TMath.Prob(chisq*(80-4), 80-4))[:10]}")
    pave2.SetBorderSize(1)
    pave2.SetFillColor(ROOT.kWhite)
    pave2.SetTextFont(42)
    
    Nhyp = inv_mass.GetEntries()*frac.getVal()
    error = inv_mass.GetEntries()*frac.getError()
    
    mframe.Draw()
#     mframe.GetYaxis().SetRangeUser(0,50)
    legend.Draw("same")
       
    pave.Draw("same")
    pave2.Draw("same")
    
    c1.Draw()
    c1.SaveAs(f"{Output}.pdf")
    
    return Nhyp, error



if __name__ == "__main__":
    
    
    #----------------------------------------------------------------------------------------
    # define variables
    # NOTE: we dont need amplitudes just a fraction since this is pdfs
    m = ROOT.RooRealVar("m", "m [GeV]", 2.96, 3.04)
    m0 = ROOT.RooRealVar("m0", "m0", 2.992,  2.990, 2.994) #mean of CB
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.001,  0.0001, 0.002) #standard deviation of CB
    frac = ROOT.RooRealVar("frac", "fraction", 0.5, 0, 1) #ratio between Signal and Bkg
    lam = ROOT.RooRealVar("Lambda", "slope of expo", -7, -50, 0) #parameter for the exponential
    
    
    # define pdf model
    # it looks like this: a*CrystalBall+(1-a)*exponential, where a is the frac parameter
    # see: https://root.cern.ch/doc/master/classRooAddPdf.html 

    
    # import data
    rdf = ROOT.RDataFrame("df", f"./SelectedDataFrames/SelectedDataFrame_{pt_min}_pt_{pt_max}.root")

    # i goes from 0 to 6, since there are 7 bins
    bins = np.linspace(-1,1,nbins+1)
        

    # plots cos(theta*) distribution-----------------------------------------------
    c2 = ROOT.TCanvas()
    cos_theta = rdf.Filter(f"Matter == {matter}").Histo1D(("CosThetaHistogram", "cos(theta*)_wrt_beam", 100, -1, 1),"cos_theta_beam")
    cos_theta.Draw()
    c2.Draw()
    c2.SaveAs(f"{Output}.pdf(")


    # create plots with signal fit-------------------------------------------------
    results = [CreatePlots(i) for i in range(nbins)]
    number_of_Hypertritons = [results[i][0] for i in range(nbins)]
    errorbars = [results[i][1] for i in range(nbins)]

    # plot number of Hypertritons (raw) -------------------------------------------
    c_Eff = ROOT.TCanvas()
    centered_bins = (bins[:-1] + bins[1:]) / 2
    
    
    h_raw = ROOT.TH1D("hist", f"Number of Hypertritons(raw) for {pt_min}"+"<p_{T}<"+f"{pt_max}; cos(#theta*) wrt beam", 7, -1, 1)
    h_raw.FillN(7, np.array(centered_bins), np.array(number_of_Hypertritons) )
    for i in range(1,8,1):
        h_raw.SetBinError(i, errorbars[i-1])

    h_raw.GetYaxis().SetTitle("Counts per bin")
#     h_raw.GetYaxis().SetRangeUser(0,150)

    ROOT.gStyle.SetTitleFontSize(0.045)

    h_raw.Draw()

    c_Eff.Draw()
    c_Eff.SaveAs(f"{Output}.pdf")

    # draw pt distr ---------------------------------------------------------------

    
    matter2 = 0.0
    if matter == "true":
        matter2 = 1.0
    

    # load the files that have been created with HypertritonRestframeBoost.py
    # since they are too big, we actually use the same files but just read in some of the columns: pt, Matter, Centrality, cos(theta*) wrt beam
    rdf_gen = ROOT.RDF.FromCSV("generatedHypertritons_cut.csv")
    rdf_reco = ROOT.RDF.FromCSV("reconstructedHypertritons_cut.csv")
    print(rdf_gen.GetColumnNames())

    pt_distr = ROOT.TCanvas()
    # NOTE: here the and must be used otherwise it will be interpreted wrong!!
    # NOTE: in the original files it is "matter" for gen, "Matter" for reco. In the *_cut.csv files it is "Matter" for both!
    rdf_gen_cut = rdf_gen.Filter(f"{pt_min} < pt and pt < {pt_max} and Matter == {matter2} and Centrality > {cmin}")
    rdf_reco_cut = rdf_reco.Filter(f"{pt_min} < pt and pt < {pt_max} and Matter == {matter2} and Centrality > {cmin}")

    pthistgen = rdf_gen_cut.Histo1D(("ptHistogram", "p_{T} Distribution for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max}"+"; p_{T}[GeV]", 30, pt_min, pt_max),"pt")
    pthistreco = rdf_reco_cut.Histo1D(("ptHistogram", "p_{T} Distribution for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max}"+"; p_{T}[GeV]", 30, pt_min, pt_max),"pt")

    pthistgen.SetFillColor(ROOT.kBlue-10)
    pthistreco.SetFillColor(ROOT.kRed-10)
    pthistgen.Draw("HIST")
    pthistreco.Draw("Same HIST")
    ROOT.gStyle.SetOptStat(0)

    legend3 = ROOT.TLegend(.7, .7, .85, .85)
    legend3.SetBorderSize(0)
    legend3.AddEntry(pthistgen.GetPtr(), "generated", "f")
    legend3.AddEntry(pthistreco.GetPtr(), "reconstruced", "f")

    legend3.SetTextSize(0.025)
    legend3.Draw()

    pt_distr.RedrawAxis()
    pt_distr.Draw()
    pt_distr.SaveAs(f"{Output}.pdf")
    
    #plot efficiency for whole pt range -------------------------------------------
    histgen = rdf_gen_cut.Histo1D(("cosThetaHistogram", "Cos(#theta_{beam}) Distribution for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max}"+"; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")
    histreco = rdf_reco_cut.Histo1D(("cosThetaHistogram", "Cos(#theta_{beam}) Distribution for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max}"+"; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")

    histgen.SetFillColor(ROOT.kBlue-10)
    histreco.SetFillColor(ROOT.kRed-10)
    histgen.Draw("hist")
    histreco.Draw("Same hist")

    histgen.GetYaxis().SetRangeUser(0,250000)

    pt_distr.Draw()
    pt_distr.SaveAs(f"{Output}.pdf")

    # calulate detector efficiency
    d = ROOT.TCanvas()
    efficiency = histreco.Clone()

    efficiency.Divide(histreco.GetPtr(), histgen.GetPtr(), 1, 1, "B")
    efficiency.Draw("E")
    efficiency.SetTitle("Detector Efficiency for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max}")
    efficiency.GetYaxis().SetTitle("Efficiency")
    efficiency.GetYaxis().SetRangeUser(0.2, 0.8)
    d.Draw()
    d.SaveAs(f"{Output}.pdf")

    count = [efficiency.GetBinContent(i) for i in range(1,nbins+1,1)]
    eff_err = [efficiency.GetBinError(i) for i in range(1,nbins+1,1)]
    print(count)


    #plot number of Hypertritons (corrected)---------------------------------------
    yerr = [np.sqrt((errorbars[i]/count[i])**2+((number_of_Hypertritons[i]*eff_err[i])/count[i]**2)**2) for i in range(nbins)]

    c_Eff = ROOT.TCanvas()
    centered_bins = (bins[:-1] + bins[1:]) / 2
    # TODO
    h_cor = ROOT.TH1D("hist", "Number of Hypertritons(corrected) for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max};"+" cos(#theta_{beam})", 7, -1, 1)
    
    h_cor.FillN(7, np.array(centered_bins), np.array(number_of_Hypertritons)/np.array(count) )
    for i in range(1,8,1):
        h_cor.SetBinError(i, yerr[i-1])

    h_cor.GetYaxis().SetTitle("Counts per bin")
#     h_cor.GetYaxis().SetRangeUser(0,300)

    ROOT.gStyle.SetTitleFontSize(0.045)

    h_cor.Draw()

    
    #save number of hypertritons after correction in file
    Nhyp_corrected = np.array(number_of_Hypertritons)/np.array(count)
    Nhyp_corrected = Nhyp_corrected.tolist()
    with open(f"NhypCorrected{pt_min}pt{pt_max}centr{cmin}_{matter3}.json", 'w') as f:
        json.dump(Nhyp_corrected, f, indent=2) 
    with open(f"NhypCorrectedErr{pt_min}pt{pt_max}centr{cmin}_{matter3}.json", 'w') as f:
        json.dump(yerr, f, indent=2) 
        
    ROOT.gStyle.SetTitleFontSize(0.045)

    h_cor.Draw()

    c_Eff.Draw()
    c_Eff.SaveAs(f"{Output}.pdf")


    # plot efficiencies for different pt ranges (with uncertainties) --------------
    c4 = ROOT.TCanvas()

    rdf_gen_cut_34 = rdf_gen.Filter(f"3 < pt and pt < 4")
    rdf_reco_cut_34 = rdf_reco.Filter(f"3 < pt and pt < 4")
    rdf_gen_cut_45 = rdf_gen.Filter(f"4 < pt and pt < 5")
    rdf_reco_cut_45 = rdf_reco.Filter(f"4 < pt and pt < 5")
    rdf_gen_cut_56 = rdf_gen.Filter(f"5 < pt and pt < 6")
    rdf_reco_cut_56 = rdf_reco.Filter(f"5 < pt and pt < 6")
    rdf_gen_cut_36 = rdf_gen.Filter(f"3 < pt and pt < 6")
    rdf_reco_cut_36 = rdf_reco.Filter(f"3 < pt and pt < 6")

    histgen34 = rdf_gen_cut_34.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")
    histreco34 = rdf_reco_cut_34.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")

    efficiency34 = histreco34.Clone()
    efficiency34.Divide(histreco34.GetPtr(), histgen34.GetPtr(), 1, 1, "B")
    efficiency34.SetTitle("Detector Efficiency for different p_{T} ranges.")
    efficiency34.SetLineColor(9)
    efficiency34.Draw("E")
    efficiency34.GetYaxis().SetTitle("Efficiency")

    histgen45 = rdf_gen_cut_45.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")
    histreco45 = rdf_reco_cut_45.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")

    efficiency45 = histreco45.Clone()
    efficiency45.Divide(histreco45.GetPtr(), histgen45.GetPtr(), 1, 1, "B")
    efficiency45.Draw("Same E")
    efficiency45.SetLineColor(8)

    histgen56 = rdf_gen_cut_56.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")
    histreco56 = rdf_reco_cut_56.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")

    efficiency56 = histreco45.Clone()
    efficiency56.Divide(histreco56.GetPtr(), histgen56.GetPtr(), 1, 1, "B")
    efficiency56.Draw("Same E")
    efficiency56.SetLineColor(6)

    histgen36 = rdf_gen_cut_36.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")
    histreco36 = rdf_reco_cut_36.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")

    efficiency36 = histreco36.Clone()
    efficiency36.Divide(histreco36.GetPtr(), histgen36.GetPtr(), 1, 1, "B")
    efficiency36.Draw("Same E")
    efficiency36.SetLineColor(1)

    efficiency34.GetXaxis().SetRangeUser(-1,1)
    efficiency34.GetYaxis().SetRangeUser(0.25,0.7)

    legend3 = ROOT.TLegend(.73, .7, .86, .85)
    legend3.SetBorderSize(1)
    legend3.AddEntry(efficiency36, "3 < p_{T} < 6", "l")
    legend3.AddEntry(efficiency34, "3 < p_{T} < 4", "l")
    legend3.AddEntry(efficiency45, "4 < p_{T} < 5", "l")
    legend3.AddEntry(efficiency56, "5 < p_{T} < 6", "l")
    legend3.SetTextSize(0.025)
    legend3.Draw()

    c4.SetGrid()
    c4.Draw()
    c4.SaveAs(f"{Output}.pdf)")