#!/usr/bin/python3
import ROOT
import numpy as np
import json
import logging
import datetime


#--------------FUNCTIONS--------------

# These two function determine the CB parameters from the MC fit
def CreatePlotsMC(i, fit_func, rdf, m, m0, sigma, alphal, nl, alphar, nr, cmin, cmax, matter, Output, nbins, save_Output):
    '''
    This function is fitting a RooCrystalBall Shape on MC SignalTable and storing the fit parameters in a .json file, and plots in .pdf file.
    Returns: fit parameter 
    '''

    # create canvas
    if save_Output == True:
        c1 = ROOT.TCanvas()
        mframe = m.frame(Title="m histogram + Fit")
    
    # pt bins
    # bins = np.linspace(2,6,nbins+1)
    bins = np.array([2,3,4,6])


    # create histogram
    inv_mass = rdf.Filter(f"{bins[i]} < pt and pt < {bins[i+1]} and Matter == {matter} and {cmin} < centrality and centrality < {cmax}")\
        .Histo1D(("InvariantMassHistogram", f"{cmin}"+"< centrality <"+f"{cmax}, {str(bins[i])[:5]}"+"< pt <"+f"{str(bins[i+1])[:5]}; m[GeV]", binnumber, 2.96, 3.04),"m")

    # inv_mass_copy = inv_mass.Clone()
    
    # get RooDataHist
    dh = ROOT.RooDataHist("dh", "dataset with m", m, inv_mass.GetPtr())
    del inv_mass
    # fit pdf and plot fit 
    fit_func.fitTo(dh)

    if save_Output == True:
        # draw histogram and fit
        dh.plotOn(mframe, MarkerSize=0, Binning=binnumber)    
        fit_func.plotOn(mframe, ROOT.RooFit.Precision(1e-5),  ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineColor(ROOT.kRed)) 
        
        # plot signal and background
        # fit_func.plotOn(mframe, ROOT.RooFit.Components("Expo"), ROOT.RooFit.LineColor(ROOT.kGreen+1) , LineWidth=2)
        # fit_func.plotOn(mframe, ROOT.RooFit.Components("gauss"), ROOT.RooFit.LineColor(ROOT.kBlue+1) , LineWidth=2)
        
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


    return nl.getVal(), nr.getVal(), alphal.getVal(), alphar.getVal(), nl.getError(), nr.getError(), alphal.getError(), alphar.getError(), sigma.getVal(), sigma.getError()

def DetermineCBParametersFromMC(rdfMC, Output,  cmin, cmax, matter, nbins, save_Output):
    '''
    This Function Defines the Variables and PDF for the Fit and stores the obtained Parameters for each bin in a list.
    Needs to be done for Matter and Antimatter and each BDTEfficiency
    Returns: Fit parameter List
    '''
    
    # define variables
    m = ROOT.RooRealVar("m", "m [GeV]", 2.96, 3.04)
    m0 = ROOT.RooRealVar("m0", "mean of CB", 2.992,  2.990, 2.994)
    sigma = ROOT.RooRealVar("sigma", "sigma of CB", 0.001,  0.0001, 0.002)

    
    alphal = ROOT.RooRealVar("alphal", "alpha L", 2.05,0.8, 2.1)
    nl = ROOT.RooRealVar("nl", "n L", 4.7,  3, 7.2)
    alphar = ROOT.RooRealVar("alphar", "alpha right", 2.08, 0.8, 2.2)
    nr = ROOT.RooRealVar("nr", "n right", 3.7, 2 , 6.5)

    # define pdf model
    fit_func = ROOT.RooCrystalBall("fit_func", "Crystal Ball PDF", m, m0, sigma, alphal, nl, alphar, nr)


    # import data
    rdf = rdfMC
        

    # plots cos(theta*) distribution-----------------------------------------------
    if save_Output == True:
        c2 = ROOT.TCanvas()
        cos_theta = rdf.Filter(f"Matter == {matter}").Histo1D(("CosThetaHistogram", "pt", 100, 2, 6),"pt")
        cos_theta.Draw()
        c2.Draw()
        c2.SaveAs(f"{Output}.pdf(")


    # create plots with signal fit-------------------------------------------------
    results = [CreatePlotsMC(i, fit_func, rdf, m, m0, sigma, alphal, nl, alphar, nr, cmin, cmax, matter, Output, nbins, save_Output) for i in range(nbins)]
    
    if matter == "true":
        matter3 = "M"
    elif matter == "false":
        matter3 = "AM"

    with open(f"../Output/EPCutOutput/FixedCBParamsCentr{cmin}-{cmax}_{matter3}.json", 'w') as f:
        json.dump(results, f, indent=2) 

    logging.info("MC fit parameters saved.")
    del c2
    del rdf
    del cos_theta

    return results

# these two functions obtain the Number of Hypertriton Candidates from the Data
def CreatePlots(i, alphal_list, alphar_list, nl_list, nr_list, m, m0, sigma_list, lam, lin, quad, frac, nbins, cmin, matter, rdf, pt_min, pt_max, Output, model, fitsigma, save_Output):
    
    # Fix parameters based on MC fit
    alphal = ROOT.RooRealVar("alphal", "alpha L", alphal_list[i], alphal_list[i], alphal_list[i])
    nl = ROOT.RooRealVar("nl", "n L", nl_list[i], nl_list[i], nl_list[i])
    alphar = ROOT.RooRealVar("alphar", "alpha right", alphar_list[i], alphar_list[i], alphar_list[i])
    nr = ROOT.RooRealVar("nr", "n right", nr_list[i], nr_list[i], nr_list[i])


    if fitsigma == True:
        # sigma is determined by fit of Data
        sigma = ROOT.RooRealVar("sigma", "sigma", sigma_list[i],  sigma_list[i]-0.0008, sigma_list[i]+0.0005) #standard deviation of CB
    else:
        # sigma is determined by MC fit
        sigma = ROOT.RooRealVar("sigma", "sigma", sigma_list[i], sigma_list[i], sigma_list[i])
    
    # Define PDF for Fit
    CB = ROOT.RooCrystalBall("CB", "Crystal Ball PDF", m, m0, sigma, alphal, nl, alphar, nr)
    expo = ROOT.RooExponential("Expo", "expo pdf", m, lam )
    pol1 = ROOT.RooPolynomial("Pol1", "linear", m, ROOT.RooArgList(lin))
    pol2 = ROOT.RooPolynomial("Pol2", "quadratic", m, ROOT.RooArgList(lin, quad))

    # bins = np.linspace(2,6,nbins+1)
    bins = np.array([2,3,4,6])

    # create histogram
    inv_mass = rdf.Filter(f"{bins[i]} < pt and pt < {bins[i+1]} and Matter == {matter} and {cmin} < centrality and centrality < {cmax}")\
        .Histo1D(("InvariantMassHistogram", f"{pt_min}"+"<p_{T}<"+f"{pt_max}, {str(bins[i])[:5]}"+"< cos(#theta_{beam})<"+f"{str(bins[i+1])[:5]}; m[GeV]", binnumber, 2.96, 3.04),"m")
    # inv_mass_copy = inv_mass.Clone()
    
    # get RooDataHist
    dh = ROOT.RooDataHist("dh", "dataset with m", m, inv_mass.GetPtr())
    

    mframe = m.frame(Title="m histogram + Fit")

    # Create Canvas
    if save_Output == True:
        c1 = ROOT.TCanvas()

        # Draw histogram 
        dh.plotOn(mframe, ROOT.RooFit.Name("dh"), ROOT.RooFit.MarkerSize(0))#, Binning=80)




    # Fit using different PDFs
    if model == "ExpoCB":

        fit_func = ROOT.RooAddPdf("fit_func", "Fit_function_pdf", ROOT.RooArgList(CB, expo), ROOT.RooArgList(frac))

        fit_func.fitTo(dh)
        if save_Output == True:
            fit_func.plotOn(mframe, ROOT.RooFit.Precision(1e-5),  ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("fit_func")) 
            fit_func.plotOn(mframe, ROOT.RooFit.Components("Expo"), ROOT.RooFit.LineColor(ROOT.kGreen+1) , LineWidth=2)
            fit_func.plotOn(mframe, ROOT.RooFit.Components("CB"), ROOT.RooFit.LineColor(ROOT.kBlue+1) , LineWidth=2)
            
            # Plot Legend
            legend = ROOT.TLegend(.12, .8, .49, .935)
            legend.SetBorderSize(1)
            legend.AddEntry(dh, "Data", "le")
            legend.AddEntry(fit_func, "Fit = Signal + Background", "l").SetLineColor(ROOT.kRed)
            legend.AddEntry(expo, "Background = (1-a)*exp(x #lambda))", "l").SetLineColor(ROOT.kGreen+1)
            legend.AddEntry(CB, "Signal = a*CB(x, #mu, #sigma, n_{L}, #alpha_{L}, n_{R}, #alpha_{R})", "l").SetLineColor(ROOT.kBlue)
            legend.SetTextSize(0.025)

            pave = ROOT.TPaveText(0.5, 0.75, 0.77, 0.935, "NDC")
            pave.SetBorderSize(1)
            pave.SetFillColor(ROOT.kWhite)
            pave.SetTextFont(42)

            t0 = pave.AddText(f"#mu = {str(m0.getVal())[:10]} #pm {str(m0.getError())[:10]}")
            t1 = pave.AddText(f"#sigma = {str(sigma.getVal())[:10]} #pm {str(sigma.getError())[:10]}")
            t6 = pave.AddText(f"a = {str(frac.getVal())[:10]} #pm {str(frac.getError())[:10]}")
            t7 = pave.AddText(f"#lambda = {str(lam.getVal())[:10]} #pm {str(lam.getError())[:10]}")

    elif model == "Pol1CB":
        fit_func = ROOT.RooAddPdf("fit_func", "Fit_function_pdf", ROOT.RooArgList(CB, pol1), ROOT.RooArgList(frac))

        fit_func.fitTo(dh)
        
        if save_Output == True:
            fit_func.plotOn(mframe, ROOT.RooFit.Precision(1e-5),  ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("fit_func")) 
            fit_func.plotOn(mframe, ROOT.RooFit.Components("Pol1"), ROOT.RooFit.LineColor(ROOT.kGreen+1) , LineWidth=2)
            fit_func.plotOn(mframe, ROOT.RooFit.Components("CB"), ROOT.RooFit.LineColor(ROOT.kBlue+1) , LineWidth=2)
            
            # Plot Legend
            legend = ROOT.TLegend(.12, .8, .49, .935)
            legend.SetBorderSize(1)
            legend.AddEntry(dh, "Data", "le")
            legend.AddEntry(fit_func, "Fit = Signal + Background", "l").SetLineColor(ROOT.kRed)
            legend.AddEntry(pol1, "Background = (1-a)*(c_{1}x))", "l").SetLineColor(ROOT.kGreen+1)
            legend.AddEntry(CB, "Signal = a*CB(x, #mu, #sigma, n_{L}, #alpha_{L}, n_{R}, #alpha_{R})", "l").SetLineColor(ROOT.kBlue)
            legend.SetTextSize(0.025)

            pave = ROOT.TPaveText(0.5, 0.75, 0.77, 0.935, "NDC")
            pave.SetBorderSize(1)
            pave.SetFillColor(ROOT.kWhite)
            pave.SetTextFont(42)

            t0 = pave.AddText(f"#mu = {str(m0.getVal())[:10]} #pm {str(m0.getError())[:10]}")
            t1 = pave.AddText(f"#sigma = {str(sigma.getVal())[:10]} #pm {str(sigma.getError())[:10]}")
            t6 = pave.AddText(f"a = {str(frac.getVal())[:10]} #pm {str(frac.getError())[:10]}")
            t7 = pave.AddText("c_{1}"+f" = {str(lin.getVal())[:10]} #pm {str(lin.getError())[:10]}")

    elif model == "Pol2CB":
        fit_func = ROOT.RooAddPdf("fit_func", "Fit_function_pdf", ROOT.RooArgList(CB, pol2), ROOT.RooArgList(frac))

        fit_func.fitTo(dh)

        if save_Output == True:
            fit_func.plotOn(mframe, ROOT.RooFit.Precision(1e-5),  ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("fit_func")) 
            fit_func.plotOn(mframe, ROOT.RooFit.Components("Pol2"), ROOT.RooFit.LineColor(ROOT.kGreen+1) , LineWidth=2)
            fit_func.plotOn(mframe, ROOT.RooFit.Components("CB"), ROOT.RooFit.LineColor(ROOT.kBlue+1) , LineWidth=2)
            
            # Plot Legend
            legend = ROOT.TLegend(.12, .8, .49, .935)
            legend.SetBorderSize(1)
            legend.AddEntry(dh, "Data", "le")
            legend.AddEntry(fit_func, "Fit = Signal + Background", "l").SetLineColor(ROOT.kRed)
            legend.AddEntry(pol2, "Background = (1-a)*(c_{2}x^2+c_{1}x))", "l").SetLineColor(ROOT.kGreen+1)
            legend.AddEntry(CB, "Signal = a*CB(x, #mu, #sigma, n_{L}, #alpha_{L}, n_{R}, #alpha_{R})", "l").SetLineColor(ROOT.kBlue)
            legend.SetTextSize(0.025)

            pave = ROOT.TPaveText(0.5, 0.75, 0.77, 0.935, "NDC")
            pave.SetBorderSize(1)
            pave.SetFillColor(ROOT.kWhite)
            pave.SetTextFont(42)

            t0 = pave.AddText(f"#mu = {str(m0.getVal())[:10]} #pm {str(m0.getError())[:10]}")
            t1 = pave.AddText(f"#sigma = {str(sigma.getVal())[:10]} #pm {str(sigma.getError())[:10]}")
            t6 = pave.AddText(f"a = {str(frac.getVal())[:10]} #pm {str(frac.getError())[:10]}")
            t7 = pave.AddText("c_{2}"+f" = {str(quad.getVal())[:10]} #pm {str(quad.getError())[:10]}")
            t8 = pave.AddText("c_{1}"+f" = {str(lin.getVal())[:10]} #pm {str(lin.getError())[:10]}")

    else:
        print("Not a valid model selected. Exiting...")
        exit()

    
    chisq = mframe.chiSquare("fit_func", "dh", 4)
    if save_Output == True:
        pave2 = ROOT.TPaveText(0.78, 0.8, 0.95, 0.935, "NDC")
        pave2.AddText("#Chi^{2}/ndf"+f" = {str(chisq)[:10]}")
        # degrees of freedom is number of bins minus number of fit parameters
        pave2.AddText(f"Probabilty = {str(ROOT.TMath.Prob(chisq*(binnumber-4), binnumber-4))[:10]}")
        pave2.SetBorderSize(1)
        pave2.SetFillColor(ROOT.kWhite)
        pave2.SetTextFont(42)
    
    Nhyp = inv_mass.GetEntries()*frac.getVal()
    error = inv_mass.GetEntries()*frac.getError()
    
    if save_Output == True:
        mframe.Draw()
    #     mframe.GetYaxis().SetRangeUser(0,50)
        legend.Draw("same")
        
        pave.Draw("same")
        pave2.Draw("same")
        
        c1.Draw()
        c1.SaveAs(f"{Output}.pdf")
    del dh
    del mframe
    del inv_mass
    
    return Nhyp, error

def GetNHypertritonCandidates(rdfData, results, BDTEfficiency, cmin, cmax, matter, nbins, fitsigma=True, model="ExpoCB", save_Output=False):
    #parameters from Fitting MCs with RooCrystalBall
    if matter == "true":
            matter3 = "M"
    elif matter == "false":
            matter3 = "AM"
    fitsigma2 = ""
    if fitsigma == True:
        fitsigma2 = "S"

    BDTEff = "%.2f" %BDTEfficiency
    Output= Outputdir + f"RooFit{model}{fitsigma2}Centr{cmin}-{cmax}_BDT{BDTEff}_{matter3}"#TODO: mache es wie in cpp file wo man das ende in der config bestimmt und dann einfach anhängt
    # with open(f"../Output/EPCutOutput/FixedCBParams{pt_min}pt{pt_max}_{matter3}.json", 'r') as f:
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
    quad = ROOT.RooRealVar("c2", "f(x)=const+lin*x+quad*x²", -7, -50, 0)
    lin = ROOT.RooRealVar("c1", "f(x)=const+lin*x+quad*x²", -7, -50, 50)
    # NOTE: Constant doesnt have to be added

    
    # import data
    rdf = rdfData

    # i goes from 0 to 6, since there are 7 bins
    # bins = np.linspace(2,6,nbins+1)
    bins = np.array([2,3,4,6])
        

    # plots cos(theta*) distribution-----------------------------------------------
    if save_Output == True:
        c2 = ROOT.TCanvas()
        cos_theta = rdf.Filter(f"Matter == {matter}").Histo1D(("CosThetaHistogram", "pt distribution", 100, 2, 6),"pt")
        cos_theta.Draw()
        # s = cos_theta.FindObject("stats")
        # s.SetX1NDC(0.7)
        # s.SetX2NDC(0.9)
        c2.Draw()
        c2.SaveAs(f"{Output}.pdf(")


    # create plots with signal fit-------------------------------------------------
    results = [CreatePlots(i, alphal_list, alphar_list, nl_list, nr_list, m, m0, sigma_list, lam, lin, quad, frac, nbins, cmin, matter, rdf, cmin, cmax, Output, model, fitsigma, save_Output) for i in range(nbins)]
    number_of_Hypertritons = [results[i][0] for i in range(nbins)]
    errorbars = [results[i][1] for i in range(nbins)]

    # plot number of Hypertritons (raw) -------------------------------------------
    if save_Output == True:
        c_Eff = ROOT.TCanvas()
        centered_bins = (bins[:-1] + bins[1:]) / 2
        
        
        h_raw = ROOT.TH1D("hist raw", f"Number of Hypertritons(raw) for {cmin}"+"<centrality<"+f"{cmax}; pt; counts/bin", nbins, 2, 6)
        h_raw.FillN(nbins, np.array(centered_bins), np.array(number_of_Hypertritons) )
        for i in range(1,nbins+1,1):
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

        
    # NOTE: here the and must be used otherwise it will be interpreted wrong!!
    # NOTE: in the original files it is "matter" for gen, "Matter" for reco. In the *_cut.csv files it is "Matter" for both!
    rdf_gen_cut = rdf_gen.Filter(f"{cmin} < Centrality and Centrality < {cmax} and Matter == {matter2} and 2 < pt and pt < 6")
    rdf_reco_cut = rdf_reco.Filter(f"{cmin} < Centrality and Centrality < {cmax} and Matter == {matter2} and 2 < pt and pt < 6")

    pthistgen = rdf_gen_cut.Histo1D(("ptHistogram", "p_{T} Distribution for"+f"{cmin}"+"< Centrality <"+f"{cmax}"+"; p_{T}[GeV]", 30, 2, 6),"pt")
    pthistreco = rdf_reco_cut.Histo1D(("ptHistogram", "p_{T} Distribution for"+f"{cmin}"+"< Centrality <"+f"{cmax}"+"; p_{T}[GeV]", 30, 2, 6),"pt")
        
    if save_Output == True:
        pt_distr = ROOT.TCanvas()
        
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
    histgen = rdf_gen_cut.Histo1D(("cosThetaHistogram", "p_{T} Distribution for"+f"{cmin}"+"< Centrality <"+f"{cmax}"+"; p_{T}", nbins, 2, 6),"pt")
    histreco = rdf_reco_cut.Histo1D(("cosThetaHistogram", "p_{T} Distribution for"+f"{cmin}"+"< Centrality <"+f"{cmax}"+"; p_{T}", nbins, 2, 6),"pt")

    if save_Output == True:
        histgen.SetFillColor(ROOT.kBlue-10)
        histreco.SetFillColor(ROOT.kRed-10)
        histgen.Draw("hist")
        histreco.Draw("Same hist")

        histgen.GetYaxis().SetRangeUser(0,250000)

        pt_distr.Draw()
        pt_distr.SaveAs(f"{Output}.pdf")

    # calulate detector efficiency
    
    efficiency = histreco.Clone()

    efficiency.Divide(histreco.GetPtr(), histgen.GetPtr(), 1, 1, "B")
    if save_Output == True:
        d = ROOT.TCanvas()
        efficiency.Draw("E")
        efficiency.SetTitle("Detector Efficiency for"+f"{cmin}"+"<centrality<"+f"{cmax}")
        efficiency.GetYaxis().SetTitle("Efficiency")
        efficiency.GetYaxis().SetRangeUser(0.2, 0.8)
        d.Draw()
        d.SaveAs(f"{Output}.pdf")

    count = [efficiency.GetBinContent(i) for i in range(1,nbins+1,1)]
    eff_err = [efficiency.GetBinError(i) for i in range(1,nbins+1,1)]
    print(count)


    #plot number of Hypertritons (corrected)---------------------------------------
    yerr = [np.sqrt((errorbars[i]/count[i])**2+((number_of_Hypertritons[i]*eff_err[i])/count[i]**2)**2) for i in range(nbins)]
    if save_Output == True:
        # c_Eff = ROOT.TCanvas()
        centered_bins = (bins[:-1] + bins[1:]) / 2

        h_cor = ROOT.TH1D("hist corrected", "Number of Hypertritons(corrected) for"+f"{cmin}"+"< centrality <"+f"{cmax};"+" p_{T}", nbins, 2, 6)
        h_cor.FillN(nbins, np.array(centered_bins), np.array(number_of_Hypertritons)/np.array(count) )
        for i in range(1,nbins+1,1):
            h_cor.SetBinError(i, yerr[i-1])

        h_cor.GetYaxis().SetTitle("Counts per bin")

        ROOT.gStyle.SetTitleFontSize(0.045)

        h_cor.Draw()
        d.Draw()
        d.SaveAs(f"{Output}.pdf)")
        print("output saved!")

    
    #save number of hypertritons after correction in file
    Nhyp_corrected = np.array(number_of_Hypertritons)/np.array(count)
    Nhyp_corrected = Nhyp_corrected.tolist()
    with open(Outputdir + f"NhypCorrectedCentr{cmin}-{cmax}_{matter3}.json", 'w') as f:
        json.dump(Nhyp_corrected, f, indent=2) 
    with open(Outputdir + f"NhypCorrectedErrCentr{cmin}-{cmax}_{matter3}.json", 'w') as f:
        json.dump(yerr, f, indent=2) 

   
    return Nhyp_corrected, yerr
    
# These functions are used for the rest
def ProcessForOneBDTEff(rdfData, rdfMC, cmin, cmax, BDTEfficiency, nbins, matter, inplane, save_output):
    #-------------SETTINGS------------------
    # TODO: Add centrality cut, ggf.
    matter2 = "AM"
    if matter == "true":
        matter2 = "M"


    # Specify Output File Paths
    OutputMC= Outputdir + f"RooFitCBcentr{cmin}-{cmax}_BDT{BDTEfficiency}_{matter2}"
   
    
    print("------------MC SETTINGS-------------")
    print(f"  Matter:             {matter}")
    print(f"  Selected Centrality range:  {cmin}-{cmax}%")
    print("------------------------------------")

    
    # Create the Dataframes and apply cuts
    print("Filter BDTEfficiency = ", BDTEfficiency)
    rdfData_new = rdfData.Filter(f"BDTEfficiency < {BDTEfficiency}")#TODO: plus ggf. Centrality cut
    rdfMC_new = rdfMC.Filter(f"BDTEfficiency < {BDTEfficiency}")#TODO: plus ggf. Centrality cut
    
    # Calculate the Parameters for CB from fitting the MCs
    results = DetermineCBParametersFromMC(rdfMC_new, OutputMC, cmin, cmax, matter, nbins, save_Output=save_output)


    print("------------Data SETTINGS-----------")
    print(f"  Matter:             {matter}")
    print(f"  Selected Centrality range:  {cmin}-{cmax}%")
    print(f"  Number of bins:     {nbins}")
    print(f"  Event Plane Cut:    inplane: {inplane} ")
    print("------------------------------------")

    # do it for every model (Pol1, Pol2, Expo)+CB(with and without fixed sigma)
    ExpoNhyp, ExpoYerr = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, cmin, cmax, matter, nbins, fitsigma=True, model="ExpoCB", save_Output=save_output)
    # Pol1Nhyp, Pol1Yerr = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, cmin, cmax, matter, nbins, fitsigma=True, model="Pol1CB", save_Output=False)
    # Pol2Nhyp, Pol2Yerr = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, cmin, cmax, matter, nbins, fitsigma=True, model="Pol2CB", save_Output=False)  
    # ExpoNhypS, ExpoYerrS = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, cmin, cmax, matter, nbins, fitsigma=False, model="ExpoCB", save_Output=False)
    # Pol1NhypS, Pol1YerrS = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, cmin, cmax, matter, nbins, fitsigma=False, model="Pol1CB", save_Output=False)
    # Pol2NhypS, Pol2YerrS = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, cmin, cmax, matter, nbins, fitsigma=False, model="Pol2CB", save_Output=False)        


    Nhyp = [ExpoNhyp]#, ExpoNhypS, Pol1Nhyp, Pol1NhypS, Pol2Nhyp, Pol2NhypS]
    Yerr = [ExpoYerr]#, ExpoYerrS, Pol1Yerr, Pol1YerrS, Pol2Yerr, Pol2YerrS]
    return Nhyp, Yerr


def GetResults(rdfData, rdfMC, matter,cmin, cmax, BDTEfficiencies, nbins, inplane, save_output):

    matter2 = "AM"
    if matter == "true":
        matter2 = "M"


    Results = {}
    Results["bin0"]=[]
    Results["bin1"]=[]
    Results["bin2"]=[]
    # Results["bin3"]=[]
    # Results["bin4"]=[]
    # Results["bin5"]=[]
    # Results["bin6"]=[]

    Results["bin0yerr"] = []
    Results["bin1yerr"] = []
    Results["bin2yerr"] = []
    # Results["bin3yerr"] = []
    # Results["bin4yerr"] = []
    # Results["bin5yerr"] = []
    # Results["bin6yerr"] = []


    for BDTEff in BDTEfficiencies:
        logging.info(f"BDTEfficiency: {BDTEff}")
        Nhyp, Yerr = ProcessForOneBDTEff(rdfData, rdfMC, cmin, cmax, BDTEff, nbins, matter, inplane, save_output)
        Results["bin0"] += [x[0] for x in Nhyp]
        Results["bin1"] += [x[1] for x in Nhyp]
        Results["bin2"] += [x[2] for x in Nhyp]
        # Results["bin3"] += [x[3] for x in Nhyp]
        # Results["bin4"] += [x[4] for x in Nhyp]
        # Results["bin5"] += [x[5] for x in Nhyp]
        # Results["bin6"] += [x[6] for x in Nhyp]

        Results["bin0yerr"] += [x[0] for x in Yerr]
        Results["bin1yerr"] += [x[1] for x in Yerr]
        Results["bin2yerr"] += [x[2] for x in Yerr]
        # Results["bin3yerr"] += [x[3] for x in Yerr]
        # Results["bin4yerr"] += [x[4] for x in Yerr]
        # Results["bin5yerr"] += [x[5] for x in Yerr]
        # Results["bin6yerr"] += [x[6] for x in Yerr]

    with open(Outputdir + f"Results{matter2}.json", 'w') as f:
            json.dump(Results, f, indent=2) 
    del Nhyp
    del Yerr
    return Results

def StoreSlope(ResultsMatter, ResultsAntimatter, nbins, OutputSlope, sample_size, cmin):
    np.random.seed(1773)

    #NOTE: Replace has to be True, since the nubmer of points per bin is less than the sample size

    # Select a random Element for each bin 
    Nhyp_m = [np.random.choice(np.array(ResultsMatter[f"bin{i}"]), size=sample_size, replace=True) for i in range(0,nbins,1)]
    Nhyp_am = [np.random.choice(np.array(ResultsAntimatter[f"bin{i}"]), size=sample_size, replace=True) for i in range(0,nbins,1)]

    Yerr_m = [np.random.choice(np.array(ResultsMatter[f"bin{i}yerr"]), size=sample_size, replace=True) for i in range(0,nbins,1)]
    Yerr_am = [np.random.choice(np.array(ResultsAntimatter[f"bin{i}yerr"]), size=sample_size, replace=True) for i in range(0,nbins,1)]

    #define fit function
    fit_func0 = ROOT.TF1("fit_func0", "pol0")
    fit_func1 = ROOT.TF1("fit_func1", "pol1")
    fit_func2 = ROOT.TF1("fit_func2", "pol2")
    fit_func2.FixParameter(1,0)

    slopes = []
    slope_errors = []

    # fit difference of Matter and Antimatter
    h = ROOT.TH1D("hist", "Nhyp matter-antimatter", nbins, 2, 6)

    for i in range(sample_size):
        diff = [x[i]-y[i] for x,y in zip(Nhyp_m, Nhyp_am)]
        diff_err = [np.sqrt(x[i]**2+y[i]**2) for x,y in zip(Yerr_m, Yerr_am)]

        # bins = np.linspace(2, 6, nbins+1)
        bins = np.array([2,3,4,6])
        centered_bins = (bins[:-1] + bins[1:]) / 2

        # c = ROOT.TCanvas()

        h.Reset()
        h.FillN(nbins, np.array(centered_bins), np.array(diff) )
        for j in range(1,nbins+1,1):
            h.SetBinError(j, diff_err[j-1])

        #TODO when running: select different Fitting functions, no polarisation:pol0 (=fit_func0), Spin1/2:pol1(=fit_func1), Spin3/2:pol2(=fit_func2)
        fit_func = fit_func1

        h.Fit(fit_func)
        '''
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
        '''
    

        pars = fit_func.GetParameters()
        errors = fit_func.GetParErrors()
        chisq = fit_func.GetChisquare()
        '''
        pave2 = ROOT.TPaveText(0.32, 0.8, 0.55, 0.935, "NDC")
        pave2.AddText("#Chi^{2}"+f" = {str(chisq)[:10]}")
        # degrees of freedom is number of bins minus number of fit parameters
        pave2.AddText(f"Probabilty = {str(ROOT.TMath.Prob(chisq, nbins-1))[:10]}")
        pave2.SetBorderSize(1)
        pave2.SetFillColor(ROOT.kWhite)
        pave2.SetTextFont(42)

        pave2.Draw("Same")

        c.Draw()
        c.SaveAs(f"../Output/EPCutOutput/TestingNhypSlopeFit.pdf")
        '''
        logging.info(f"Chi2:{chisq}, Prob:{ROOT.TMath.Prob(chisq, nbins-1)}")

        # TODO when running:change when checking Pol2, Pol0
        slopes.append(pars[1])
        slope_errors.append(errors[1])

    del h
    with open(Outputdir+f"SlopesCentr{cmin}-{cmax}.json", 'w') as f:
        json.dump(slopes, f, indent=2) 
    with open(Outputdir+f"SlopesErrCentr{cmin}-{cmax}.json", 'w') as f:
        json.dump(slope_errors, f, indent=2) 

    return slopes, slope_errors


#--------CODE-----------------


if __name__ == "__main__":
    logging.basicConfig(filename='v2.log', level=logging.DEBUG)
    logging.info(f"{datetime.datetime.now()}")
 
    #---------------------SETTINGS------------------
    cmin = 0
    cmax = 10
    # and 30-50%

    # BDT Cut
    BestBDTEff = 0.81

    #number of pt bins
    nbins = 3

    # number of bins in histogram per cos(theta) bin
    binnumber = 50

    #inplane?
    inplane = True
    #TODO:Change OUTPUTDIRECTORY!!

    #save_output?
    save_output = True

    logging.info("------CONFIG-----------")
    logging.info(f"cmin: {cmin}")
    logging.info(f"cmax: {cmax}")
    logging.info(f"BDTCut: {BestBDTEff}")
    logging.info(f"Number of bins: {nbins}")
    logging.info(f"Inplane = {inplane}")
    logging.info("-----------------------")


    #-----------------PATHS FOR OUTPUT-------------------
    # Outfput file paths
    OutputSlope = "../Output/V2/SlopeOutput"

    #Outputdir path
    Outputdir = "../Output/V2/In/"

    #Input dir path
    Inputdir = "../Output/V2/"

    InputDataFileName = "2.0_pt_6.0_SystematicsEPangleData"
    InputMCFileName = "2.0_pt_6.0_SystematicsEPangleMC"


    #---------------------------------------------------------
    # Cut depending on if it is inplane or out of plane
    if inplane:
        EPanglecut = f"(Phi > EPangle-{np.pi/4} and Phi < EPangle + {np.pi/4}) or Phi > EPangle + {3*np.pi/4} or Phi < EPangle - {3*np.pi/4}"
    else:
        EPanglecut = f"not((Phi > EPangle-{np.pi/4} and Phi < EPangle + {np.pi/4}) or Phi > EPangle + {3*np.pi/4} or Phi < EPangle - {3*np.pi/4})"
        # EPanglecut = f"(Phi < EPangle - {np.pi/4} and Phi > EPangle - {3*np.pi/4}) or (Phi < EPangle + {3*np.pi/4} or Phi > EPangle + {np.pi/4})"  
    
    # BDTEfficiencies that are used
    BDTEfficiencies = np.linspace(BestBDTEff-0.0, BestBDTEff+0.0, 1)
    print("BDTEfficiencies:", BDTEfficiencies)
    
    # Load dataframes
    rdfData = ROOT.RDataFrame("df", Inputdir+InputDataFileName).Filter(f"{cmin} < centrality and centrality < {cmax} and 1 < ct and ct < 35 and {EPanglecut}")
    
    rdfMC = ROOT.RDataFrame("df", Inputdir+InputMCFileName).Filter(f"{cmin} < centrality and centrality < {cmax} and 1 < ct and ct < 35 ")#and {EPanglecut}")

    # Extract Number of Hypertritons for Matter and Antimatter
    ResultsMatter = GetResults(
        rdfData, 
        rdfMC, 
        matter="true", 
        cmin=cmin, 
        cmax=cmax, 
        BDTEfficiencies=BDTEfficiencies, 
        nbins=nbins, 
        inplane=inplane, 
        save_output=save_output)
    ResultsAntimatter = GetResults(rdfData, rdfMC, matter="false",  cmin=cmin, cmax=cmax,  BDTEfficiencies=BDTEfficiencies, nbins=nbins, inplane=inplane, save_output=save_output)

    # with open(f"../Output/EPCutOutput/ResultsM.json", 'r') as f:
    #     ResultsMatter = json.load(f)

    # with open(f"../Output/EPCutOutput/ResultsAM.json", 'r') as f:
    #     ResultsAntimatter = json.load(f)
    exit()
    # Create Histogram
    h = ROOT.TH1D("hist", "Nhyp matter-antimatter", 100, -100, 100)
    h.SetCanExtend(ROOT.TH1.kAllAxes)

    # Fill the slopes into histogram
    slopes, slope_errors = StoreSlope(ResultsMatter, ResultsAntimatter, nbins, OutputSlope, sample_size=10000, cmin=cmin)#TODO: change when running

    # with open(Outputdir+f"Slopes{pt_min}pt{pt_max}centr{cmin}.json", 'r') as f:
    #     slopes = json.load(f)
    # with open(Outputdir + f"SlopesErr{pt_min}pt{pt_max}centr{cmin}.json", 'r') as f:
    #     slope_errors = json.load(f)

    del ResultsMatter
    del ResultsAntimatter

    for x in slopes:
        h.Fill(x)

    # Define Variables
    # TODO: Adjust ranges!
    b = ROOT.RooRealVar("b", "slope", -100, 100)
    b0 = ROOT.RooRealVar("b0", "mean of Gaussian", 30, -100,  100)
    sigma = ROOT.RooRealVar("sigma", "sigma of Gaussian", 25, 1,  100,)

    dh = ROOT.RooDataHist("dh", "dh", b, h)

    gauss = ROOT.RooGaussian("gauss", "gaussian PDF", b, b0, sigma)
    gauss.fitTo(dh)
    

    c1 = ROOT.TCanvas()
    mframe = b.frame(Title="Slopes Histogram + Gaussian Fit")  # RooPlot
    dh.plotOn(mframe, ROOT.RooFit.MarkerSize(0), ROOT.RooFit.Binning(40), ROOT.RooFit.Name("dh"))
    gauss.plotOn(mframe, ROOT.RooFit.Name("gauss"))

    pave = ROOT.TPaveText(0.78, 0.75, 0.95, 0.935, "NDC")
    pave.AddText(f"#mu = {str(b0.getVal())[:10]} #pm {str(b0.getError())[:10]}")
    pave.AddText(f"#sigma = {str(sigma.getVal())[:10]} #pm {str(sigma.getError())[:10]}")
    pave.SetBorderSize(1)
    pave.SetFillColor(ROOT.kWhite)
    pave.SetTextFont(42)

    mframe.Draw()
    pave.Draw("same")

    c1.Draw()
    c1.SaveAs(Outputdir+"SlopeHistogram.pdf")