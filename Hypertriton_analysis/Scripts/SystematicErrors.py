#!/usr/bin/python3
import ROOT
import numpy as np
import json
import logging
import datetime


# TODO: add option to not save output,since it will be a lot of files, remove plotting if possible to save computing time
# TODO: check that everything is done for matter and antimatter

#--------------FUNCTIONS--------------

def CreatePlotsMC(i, fit_func, rdf, m, m0, sigma, alphal, nl, alphar, nr, pt_min, pt_max, matter, Output, nbins):
    '''
    This function is fitting a RooCrystalBall Shape on MC SignalTable and storing the fit parameters in a .json file, and plots in .pdf file.
    Returns: fit parameter 
    '''
    # create canvas
    c1 = ROOT.TCanvas()
    mframe = m.frame(Title="m histogram + Fit")
    
   
    bins = np.linspace(-1,1,nbins+1)

    # create histogram
    # inv_mass_df = rdf.Filter(f"{bins[i]} < CosThetaWrtBeam and CosThetaWrtBeam < {bins[i+1]} and Matter == {matter} ")
    # inv_mass = inv_mass_df.Histo1D(("InvariantMassHistogram", f"{pt_min}"+"<p_{T}<"+f"{pt_max}, {str(bins[i])[:5]}"+"< cos(#theta_{beam})<"+f"{str(bins[i+1])[:5]}; m[GeV]", 80, 2.96, 3.04),"m")
    inv_mass = rdf.Filter(f"{bins[i]} < CosThetaWrtBeam and CosThetaWrtBeam < {bins[i+1]} and Matter == {matter} ")\
        .Histo1D(("InvariantMassHistogram", f"{pt_min}"+"<p_{T}<"+f"{pt_max}, {str(bins[i])[:5]}"+"< cos(#theta_{beam})<"+f"{str(bins[i+1])[:5]}; m[GeV]", 80, 2.96, 3.04),"m")
    
    # inv_mass_copy = inv_mass.Clone()
    
    # get RooDataHist
    dh = ROOT.RooDataHist("dh", "dataset with m", m, inv_mass.GetPtr())
    
    # draw histogram using again 80 bins
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
#     mframe.GetYaxis().SetRangeUser(0,20000)
    legend.Draw("same")
       
    pave.Draw("same")
    pave2.Draw("same")
    

    c1.Draw()
    if i < nbins-1:
        c1.SaveAs(f"{Output}.pdf")
    else:
        c1.SaveAs(f"{Output}.pdf)")
    
    del c1
    del dh
    del inv_mass

    return nl.getVal(), nr.getVal(), alphal.getVal(), alphar.getVal(), nl.getError(), nr.getError(), alphal.getError(), alphar.getError(), sigma.getVal(), sigma.getError()


def DetermineCBParametersFromMC(rdfMC, Output,  pt_min, pt_max, matter, nbins):
    '''
    This Function Defines the Variables and PDF for the Fit and stores the obtained Parameters for each bin in a list.
    Needs to be done for Matter and Antimatter and each BDTEfficiency
    Returns: Fit parameter List
    '''
    
    # define variables
    m = ROOT.RooRealVar("m", "m [GeV]", 2.96, 3.04)
    m0 = ROOT.RooRealVar("m0", "mean of CB", 2.992,  2.990, 2.994)
    sigma = ROOT.RooRealVar("sigma", "sigma of CB", 0.001,  0.0001, 0.0015)

    
    alphal = ROOT.RooRealVar("alphal", "alpha L", 2.05,1, 2.1)
    nl = ROOT.RooRealVar("nl", "n L", 4.7,  3.9, 7.2)
    alphar = ROOT.RooRealVar("alphar", "alpha right", 2.08, 1, 2.2)
    nr = ROOT.RooRealVar("nr", "n right", 3.7, 2.9 , 6.5)

    # define pdf model
    fit_func = ROOT.RooCrystalBall("fit_func", "Crystal Ball PDF", m, m0, sigma, alphal, nl, alphar, nr)


    # import data
    rdf = rdfMC
        

    # plots cos(theta*) distribution-----------------------------------------------
    c2 = ROOT.TCanvas()
    cos_theta = rdf.Filter(f"Matter == {matter}").Histo1D(("CosThetaHistogram", "cos(theta*)_wrt_beam", 100, -1, 1),"CosThetaWrtBeam")
    cos_theta.Draw()
    c2.Draw()
    c2.SaveAs(f"{Output}.pdf(")


    # create plots with signal fit-------------------------------------------------
    results = [CreatePlotsMC(i, fit_func, rdf, m, m0, sigma, alphal, nl, alphar, nr, pt_min, pt_max, matter, Output, nbins) for i in range(nbins)]
    
    '''
    # nl = [results[i][0] for i in range(len(results))]
    # nr = [results[i][1] for i in range(len(results))]
    # alphal = [results[i][2] for i in range(len(results))]
    # alphar = [results[i][3] for i in range(len(results))]
    # sigma = [results[i][8] for i in range(len(results))]
    
    # nl_err = [results[i][4] for i in range(len(results))]
    # nr_err = [results[i][5] for i in range(len(results))]
    # alphal_err = [results[i][6] for i in range(len(results))]
    # alphar_err = [results[i][7] for i in range(len(results))]
    # sigma_err = [results[i][9] for i in range(len(results))]
    '''
    if matter == "true":
        matter3 = "M"
    elif matter == "false":
        matter3 = "AM"
#     with open(f"./CBFixedMCParameters/FixedCBParams{pt_min}pt{pt_max}_{matter3}.json", 'w') as f:
    with open(f"./forSystematicErrors/FixedCBParams{pt_min}pt{pt_max}_{matter3}.json", 'w') as f:
        json.dump(results, f, indent=2) 

    logging.info("MC fit parameters saved.")
    del c2

    return results


def CreatePlots(i, alphal_list, alphar_list, nl_list, nr_list, m, m0, sigma, lam, a, b, frac, nbins, cmin, matter, rdf, pt_min, pt_max, Output, model):
    
    # Fix parameters based on MC fit
    alphal = ROOT.RooRealVar("alphal", "alpha L", alphal_list[i], alphal_list[i], alphal_list[i])
    nl = ROOT.RooRealVar("nl", "n L", nl_list[i], nl_list[i], nl_list[i])
    alphar = ROOT.RooRealVar("alphar", "alpha right", alphar_list[i], alphar_list[i], alphar_list[i])
    nr = ROOT.RooRealVar("nr", "n right", nr_list[i], nr_list[i], nr_list[i])
    
    # Define PDF for Fit
    CB = ROOT.RooCrystalBall("CB", "Crystal Ball PDF", m, m0, sigma, alphal, nl, alphar, nr)
    expo = ROOT.RooExponential("Expo", "expo pdf", m, lam )
    pol1 = ROOT.RooPolynomial("linear", "linear", m, ROOT.RooArgList(b))
    pol2 = ROOT.RooPolynomial("linear", "linear", m, ROOT.RooArgList(b, a))

    if model == "ExpoCB":
        fit_func = ROOT.RooAddPdf("fit_func", "Fit_function_pdf", ROOT.RooArgList(CB, expo), ROOT.RooArgList(frac))
    elif model == "Pol1CB":
        fit_func = ROOT.RooAddPdf("fit_func", "Fit_function_pdf", ROOT.RooArgList(pol1, expo), ROOT.RooArgList(frac))
    elif model == "Pol2CB":
        fit_func = ROOT.RooAddPdf("fit_func", "Fit_function_pdf", ROOT.RooArgList(pol2, expo), ROOT.RooArgList(frac))
    else:
        print("Not a valid model selected. Exiting...")
        exit()

    # create canvas
    c1 = ROOT.TCanvas()
    mframe = m.frame(Title="m histogram + Fit")
    
   
    bins = np.linspace(-1,1,nbins+1)

    # create histogram
    inv_mass_df = rdf.Filter(f"{bins[i]} < CosThetaWrtBeam and CosThetaWrtBeam < {bins[i+1]} and Matter == {matter}  and centrality > {cmin} ")
    inv_mass = inv_mass_df.Histo1D(("InvariantMassHistogram", f"{pt_min}"+"<p_{T}<"+f"{pt_max}, {str(bins[i])[:5]}"+"< cos(#theta_{beam})<"+f"{str(bins[i+1])[:5]}; m[GeV]", 80, 2.96, 3.04),"m")
    inv_mass_copy = inv_mass.Clone()
    
    # get RooDataHist
    dh = ROOT.RooDataHist("dh", "dataset with m", m, inv_mass.GetPtr())
    
    # draw histogram using again 80 bins
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


def GetNHypertritonCandidates(rdfData,  Output, results, pt_min, pt_max, cmin, matter, nbins = 7, fitsigma=True):
    #parameters from Fitting MCs with RooCrystalBall

    if matter == "true":
            matter3 = "M"
    elif matter == "false":
            matter3 = "AM"
    # with open(f"./forSystematicErrors/FixedCBParams{pt_min}pt{pt_max}_{matter3}.json", 'r') as f:
    #         results = json.load(f)
    

    nl_list = [results[i][0] for i in range(len(results))]
    nr_list = [results[i][1] for i in range(len(results))]
    alphal_list = [results[i][2] for i in range(len(results))]
    alphar_list = [results[i][3] for i in range(len(results))]

    nl_err = [results[i][4] for i in range(len(results))]
    nr_err = [results[i][5] for i in range(len(results))]
    alphal_err = [results[i][6] for i in range(len(results))]
    alphar_err = [results[i][7] for i in range(len(results))]

    sigma_list = [results[i][8] for i in range(len(results))]
    sigma_err = [results[i][9] for i in range(len(results))]
    
    #----------------------------------------------------------------------------------------
    # define variables
    # NOTE: we dont need amplitudes just a fraction since this is pdfs
    m = ROOT.RooRealVar("m", "m [GeV]", 2.96, 3.04)
    m0 = ROOT.RooRealVar("m0", "m0", 2.992,  2.990, 2.994) #mean of CB
   
    frac = ROOT.RooRealVar("frac", "fraction", 0.5, 0, 1) #ratio between Signal and Bkg
    lam = ROOT.RooRealVar("Lambda", "slope of expo", -7, -50, 0) #parameter for the exponential
    # TODO: diese parameter ranges müssen angepasst werden
    a = ROOT.RooRealVar("a", "f(x)=ax²+bx+c", -7, -50, 0)
    b = ROOT.RooRealVar("b", "f(x)=ax²+bx+c", -7, -50, 0)
    # NOTE: Constant doesnt have to be added

    if fitsigma == True:
        # sigma is determined by fit of Data
        sigma = ROOT.RooRealVar("sigma", "sigma", 0.001,  0.0001, 0.002) #standard deviation of CB
    else:
        # sigma is determined by MC fit
        sigma = ROOT.RooRealVar("sigma", "sigma", sigma_list[i], sigma_list[i], sigma_list[i])

    
    # import data
    rdf = rdfData

    # i goes from 0 to 6, since there are 7 bins
    bins = np.linspace(-1,1,nbins+1)
        

    # plots cos(theta*) distribution-----------------------------------------------
    c2 = ROOT.TCanvas()
    cos_theta = rdf.Filter(f"Matter == {matter}").Histo1D(("CosThetaHistogram", "cos(theta*)_wrt_beam", 100, -1, 1),"CosThetaWrtBeam")
    cos_theta.Draw()
    c2.Draw()
    c2.SaveAs(f"{Output}.pdf(")


    # create plots with signal fit-------------------------------------------------
    results = [CreatePlots(i, alphal_list, alphar_list, nl_list, nr_list, m, m0, sigma, lam, a, b, frac, nbins, cmin, matter, rdf, pt_min, pt_max, Output) for i in range(nbins)]
    number_of_Hypertritons = [results[i][0] for i in range(nbins)]
    errorbars = [results[i][1] for i in range(nbins)]

    # plot number of Hypertritons (raw) -------------------------------------------
    c_Eff = ROOT.TCanvas()
    centered_bins = (bins[:-1] + bins[1:]) / 2
    
    
    h_raw = ROOT.TH1D("hist raw", f"Number of Hypertritons(raw) for {pt_min}"+"<p_{T}<"+f"{pt_max}; cos(#theta*) wrt beam", 7, -1, 1)
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
    h_cor = ROOT.TH1D("hist corrected", "Number of Hypertritons(corrected) for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max};"+" cos(#theta_{beam})", 7, -1, 1)
    
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
    with open(f"./forSystematicErrors/NhypCorrected{pt_min}pt{pt_max}centr{cmin}_{matter3}.json", 'w') as f:
        json.dump(Nhyp_corrected, f, indent=2) 
    with open(f"./forSystematicErrors/NhypCorrectedErr{pt_min}pt{pt_max}centr{cmin}_{matter3}.json", 'w') as f:
        json.dump(yerr, f, indent=2) 
        
    ROOT.gStyle.SetTitleFontSize(0.045)

    h_cor.Draw()

    c_Eff.Draw()
    c_Eff.SaveAs(f"{Output}.pdf)")

    return Nhyp_corrected, yerr
    

def ProcessForOneBDTEff(rdfData, rdfMC, pt_min, pt_max, BDTEfficiency, nbins, matter):
    #-------------SETTINGS------------------
    # TODO: Add centrality cut, ggf.
    matter2 = "AM"
    if matter == "true":
        matter2 = "M"

    # Specify Output File Paths
    OutputMC= f"./forSystematicErrors/RooFitCB{pt_min}pt{pt_max}_BDT{BDTEfficiency}_{matter2}"
    OutputData= f"./forSystematicErrors/RooFitCBExpo{pt_min}pt{pt_max}centr{cmin}_BDT{BDTEfficiency}_{matter2}"

    #centrality cut
    cmin = 0
    
    print("------------MC SETTINGS-------------")
    print(f"  Matter:             {matter}")
    print(f"  Selected pT range:  {pt_min}<pt<{pt_max}")
    print("------------------------------------")

    
    # Create the Dataframes and apply cuts

    rdfData_new = rdfData.Filter(f"BDTEfficiency < {BDTEfficiency}")#TODO: plus ggf. Centrality cut
    rdfMC_new = rdfMC.Filter(f"BDTEfficiency < {BDTEfficiency}")#TODO: plus ggf. Centrality cut
    
    # Calculate the Parameters for CB from fitting the MCs
    results = DetermineCBParametersFromMC(rdfMC, OutputMC, matter, pt_min=pt_min, pt_max=pt_max)


    print("------------Data SETTINGS-----------")
    print(f"  Matter:             {matter}")
    print(f"  Selected pT range:  {pt_min}<pt<{pt_max}")
    print(f"  Centrality Cut:     centrality > {cmin}")
    print(f"  Number of bins:     {nbins}")
    print(f"  Event Plane Cut: -- ")
    print("------------------------------------")

    # do it for every model (Pol1, Pol2, Expo)+CB(with and without fixed sigma)
    ExpoNhyp, ExpoYerr = GetNHypertritonCandidates(rdfData, OutputData, results, pt_min, pt_max, cmin, matter, nbins, fitsigma=True, model="ExpoCB")
    Pol1Nhyp, Pol1Yerr = GetNHypertritonCandidates(rdfData, OutputData, results, pt_min, pt_max, cmin, matter, nbins, fitsigma=True, model="Pol1CB")
    Pol2Nhyp, Pol2Yerr = GetNHypertritonCandidates(rdfData, OutputData, results, pt_min, pt_max, cmin, matter, nbins, fitsigma=True, model="Pol2CB")  
    ExpoNhypS, ExpoYerrS = GetNHypertritonCandidates(rdfData, OutputData, results, pt_min, pt_max, cmin, matter, nbins, fitsigma=False, model="ExpoCB")
    Pol1NhypS, Pol1YerrS = GetNHypertritonCandidates(rdfData, OutputData, results, pt_min, pt_max, cmin, matter, nbins, fitsigma=False, model="Pol1CB")
    Pol2NhypS, Pol2YerrS = GetNHypertritonCandidates(rdfData, OutputData, results, pt_min, pt_max, cmin, matter, nbins, fitsigma=False, model="Pol2CB")        


    Nhyp = [ExpoNhyp, ExpoNhypS, Pol1Nhyp, Pol1NhypS, Pol2Nhyp, Pol2NhypS]
    Yerr = [ExpoYerr, ExpoYerrS, Pol1Yerr, Pol1YerrS, Pol2Yerr, Pol2YerrS]
    return Nhyp, Yerr


def StoreSlope(ResultsMatter, ResultsAntimatter, nbins, OutputSlope, sample_size):
    np.random.seed(2023)

    #TODO: check if replace True or False! (mit Zurücklegen oder ohne)

    # Select a random Element for each bin 
    Nhyp_m = [np.random.choice(ResultsMatter[f"bin{i}"], size=sample_size, replace=False) for i in range(0,7,1)]
    Nhyp_am = [np.random.choice(ResultsAntimatter[f"bin{i}"], size=sample_size, replace=False) for i in range(0,7,1)]

    Yerr_m = [np.random.choice(ResultsMatter[f"bin{i}yerr"], size=sample_size, replace=False) for i in range(0,7,1)]
    Yerr_am = [np.random.choice(ResultsAntimatter[f"bin{i}yerr"], size=sample_size, replace=False) for i in range(0,7,1)]

    #define fit function
    fit_func0 = ROOT.TF1("fit_func0", "pol0")
    fit_func1 = ROOT.TF1("fit_func1", "pol1")
    fit_func2 = ROOT.TF1("fit_func2", "pol2")
    fit_func2.FixParameter(1,0)

    slopes = []
    slope_errors = []
    # fit difference of Matter and Antimatter
    for i in range(sample_size):

        diff = [x[i]-y[i] for x,y in zip(Nhyp_m, Nhyp_am)]
        diff_err = [np.sqrt(x[i]**2+y[i]**2) for x,y in zip(Yerr_m, Yerr_am)]

        bins = np.linspace(-1,1,nbins+1)
        centered_bins = (bins[:-1] + bins[1:]) / 2



        c = ROOT.TCanvas()

        h = ROOT.TH1D("hist", "Nhyp matter-antimatter", 7, -1, 1)
        h.FillN(7, np.array(centered_bins), np.array(diff) )
        for j in range(1,8,1):
            h.SetBinError(j, diff_err[j-1])





    # fit, here we use a chi2 fit, 
    # since the counts are corrected and therefore the uncertainties are gaussian
    # chi2 method is default, "F" means taking the errors of the data into account
    #TODO: select different Fitting functions, no polarisation:pol0, Spin1/2:pol1, Spin3/2:pol2
        fit_func = fit_func1

        h.Fit(fit_func)
        h.GetXaxis().SetTitle("cos(#theta_{beam}*)")
        h.GetYaxis().SetTitle("Counts per bin")

        h.Draw()

        # fit_func.Draw("same")
        fit_func.SetFillColor(ROOT.kRed-10)

        legend = ROOT.TLegend(.12, .85, .3, .935)
        legend.SetTextSize(0.025)
        legend.SetBorderSize(1)
        legend.AddEntry(h, "matter-antimatter", "l")
        legend.AddEntry(fit_func, "linear Fit", "l").SetLineColor(ROOT.kRed)

        legend.Draw()

    

        pars = fit_func.GetParameters()
        errors = fit_func.GetParErrors()
        chisq = fit_func.GetChisquare()

        pave2 = ROOT.TPaveText(0.32, 0.8, 0.55, 0.935, "NDC")
        pave2.AddText("#Chi^{2}"+f" = {str(chisq)[:10]}")
        # degrees of freedom is number of bins minus number of fit parameters
        pave2.AddText(f"Probabilty = {str(ROOT.TMath.Prob(chisq, 7-1))[:10]}")
        pave2.SetBorderSize(1)
        pave2.SetFillColor(ROOT.kWhite)
        pave2.SetTextFont(42)

        pave2.Draw("Same")

        c.Draw()
        c.SaveAs(f"{OutputSlope}.pdf")

        # TODO:change when checking Pol2, Pol0
        slopes.append(pars[1])
        slope_errors.append(errors[1])

    return slopes, slope_errors


def GetResults(rdfData, rdfMC, matter, pt_min, pt_max, BDTEfficiencies, nbins):

    matter2 = "AM"
    if matter == "true":
        matter2 = "M"


    Results = {}
    Results["bin0"]=[]
    Results["bin1"]=[]
    Results["bin2"]=[]
    Results["bin3"]=[]
    Results["bin4"]=[]
    Results["bin5"]=[]
    Results["bin6"]=[]

    Results["bin0yerr"] = []
    Results["bin1yerr"] = []
    Results["bin2yerr"] = []
    Results["bin3yerr"] = []
    Results["bin4yerr"] = []
    Results["bin5yerr"] = []
    Results["bin6yerr"] = []


    for BDTEff in BDTEfficiencies:
        logging.info(f"BDTEfficiency: {BDTEff}")
        Nhyp, Yerr = ProcessForOneBDTEff(rdfData, rdfMC, pt_min, pt_max, BDTEff, nbins, matter)
        Results["bin0"] += [x[0] for x in Nhyp]
        Results["bin1"] += [x[1] for x in Nhyp]
        Results["bin2"] += [x[2] for x in Nhyp]
        Results["bin3"] += [x[3] for x in Nhyp]
        Results["bin4"] += [x[4] for x in Nhyp]
        Results["bin5"] += [x[5] for x in Nhyp]
        Results["bin6"] += [x[6] for x in Nhyp]

        Results["bin0yerr"] += [x[0] for x in Yerr]
        Results["bin1yerr"] += [x[1] for x in Yerr]
        Results["bin2yerr"] += [x[2] for x in Yerr]
        Results["bin3yerr"] += [x[3] for x in Yerr]
        Results["bin4yerr"] += [x[4] for x in Yerr]
        Results["bin5yerr"] += [x[5] for x in Yerr]
        Results["bin6yerr"] += [x[6] for x in Yerr]

    with open(f"./forSystematicErrors/Results{matter2}.json", 'w') as f:
            json.dump(Results, f, indent=2) 



#--------CODE-----------------
#TODO: mache dasselbe für antimatter nochmal damit wir später die Differenz fitten können
#TODO: change Input
if __name__ == "__main__":
    logging.basicConfig(filename='SystematicErrors-31082023.log', level=logging.DEBUG)
    logging.info(f"{datetime.datetime.now()}")
 

    pt_min = 3.0
    pt_max = 6.0

    BestBDTEff = 0.5

    #number of bins
    nbins = 7

    # Outfput file paths
    OutputSlope = 

    BDTEfficiencies = np.linspace(BestBDTEff-10, BestBDTEff+10, 21)

    rdfData = ROOT.RDataFrame("df", f"forSystematicErrors/SystematicsTestData").Filter(f"{pt_min} < pt and pt < {pt_max}")
    rdfMC = ROOT.RDataFrame("df", f"forSystematicErrors/SystematicsTestMC").Filter(f"{pt_min} < pt and pt < {pt_max}")

    # extract Number of Hypertritons for Matter and Antimatter
    ResultsMatter = GetResults(rdfData, rdfMC, matter="true", pt_min=pt_min, pt_max=pt_max, BDTEfficiencies=BDTEfficiencies)
    ResultsAntimatter = GetResults(rdfData, rdfMC, matter="false", pt_min=pt_min, pt_max=pt_max, BDTEfficiencies=BDTEfficiencies)


    h = ROOT.TH1D("hist", "Nhyp matter-antimatter", 100, -1, 1)
    h.SetCanExtend(ROOT.TH1.kAllAxes)
    # h.FillN(10000, BINS, StoreSlope(ResultsMatter, ResultsAntimatter, nbins, OutputSlope, sample_size=10000))
    # TODO: fill the slopes in a histogram
    slopes, slope_errors = StoreSlope(ResultsMatter, ResultsAntimatter, nbins, OutputSlope, sample_size=10000)
    for x in slopes:
        h.Fill(x)

    b = ROOT.RooRealVar("b", "slope", -100, 100, 40)

    dh = ROOT.RooDataHist("dh", "dh", b, h.GetPtr())