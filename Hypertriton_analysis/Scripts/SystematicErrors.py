#!/usr/bin/python3
import ROOT
import numpy as np
import json
import logging
import datetime


#--------------FUNCTIONS--------------

def CreatePlotsMC(i, fit_func, rdf, m, m0, sigma, alphal, nl, alphar, nr, pt_min, pt_max, matter, Output, nbins, save_Output):
    '''
    This function is fitting a RooCrystalBall Shape on MC SignalTable and storing the fit parameters in a .json file, and plots in .pdf file.
    Returns: fit parameter 
    '''

    # create canvas
    if save_Output == True:
        c1 = ROOT.TCanvas()
        mframe = m.frame(Title="InvMassHistogram + CrystalBallFit")
    
   
    bins = np.linspace(-1,1,nbins+1)

    # create invariant mass histogram
    inv_mass = rdf.Filter(f"{bins[i]} < CosThetaWrtBeam and CosThetaWrtBeam < {bins[i+1]} and Matter == {matter} ")\
        .Histo1D(("InvariantMassHistogram", f"{pt_min}"+"<p_{T}<"+f"{pt_max}, {str(bins[i])[:5]}"+"< cos(#theta_{beam})<"+f"{str(bins[i+1])[:5]}; m[GeV]; Counts/bin", binnumber, 2.96, 3.04),"m")
    
    # get RooDataHist
    dh = ROOT.RooDataHist("dh", "dataset with m", m, inv_mass.GetPtr())
    del inv_mass

    # fit pdf and plot fit 
    fit_func.fitTo(dh)

    if save_Output == True:
        # draw histogram and fit
        dh.plotOn(mframe, MarkerSize=0, Binning=binnumber)    
        fit_func.plotOn(mframe, ROOT.RooFit.Precision(1e-5),  ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineColor(ROOT.kRed)) 
        
        # Add legend
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

def DetermineCBParametersFromMC(rdf, Output,  pt_min, pt_max, matter, nbins, save_Output):
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


    # plots cos(theta*) distribution-----------------------------------------------
    if save_Output == True:
        c2 = ROOT.TCanvas()
        cos_theta = rdf.Filter(f"Matter == {matter}").Histo1D(("CosThetaHistogram", "cos(theta*)_wrt_beam", 100, -1, 1),"CosThetaWrtBeam")
        cos_theta.Draw()
        c2.Draw()
        c2.SaveAs(f"{Output}.pdf(")


    # create plots with signal fit-------------------------------------------------
    results = [CreatePlotsMC(i, fit_func, rdf, m, m0, sigma, alphal, nl, alphar, nr, pt_min, pt_max, matter, Output, nbins, save_Output) for i in range(nbins)]
    
    matter3 = "M" if matter == "true" else "AM"

    with open(f"../Output/EPCutOutput/FixedCBParams{pt_min}pt{pt_max}_{matter3}.json", 'w') as f:
        json.dump(results, f, indent=2) 

    logging.info("MC fit parameters saved.")
    del c2
    del rdf
    del cos_theta

    return results

def CreatePlots(i, alphal_list, alphar_list, nl_list, nr_list, m, m0, sigma_list, lam, lin, quad, frac, nbins, cmin, matter, rdf, pt_min, pt_max, Output, model, fitsigma, save_Output):
    
    # Fix Crystal Ball parameters based on MC fit
    alphal = ROOT.RooRealVar("alphal", "alpha L", alphal_list[i], alphal_list[i], alphal_list[i])
    nl = ROOT.RooRealVar("nl", "n L", nl_list[i], nl_list[i], nl_list[i])
    alphar = ROOT.RooRealVar("alphar", "alpha right", alphar_list[i], alphar_list[i], alphar_list[i])
    nr = ROOT.RooRealVar("nr", "n right", nr_list[i], nr_list[i], nr_list[i])


    if fitsigma:
        # sigma is determined by fit of Data
        sigma = ROOT.RooRealVar("sigma", "sigma", sigma_list[i],  sigma_list[i]-0.0008, sigma_list[i]+0.0005) 
    else:
        # sigma is determined by MC fit
        sigma = ROOT.RooRealVar("sigma", "sigma", sigma_list[i], sigma_list[i], sigma_list[i])
    
    # Define PDF for Fit
    CB = ROOT.RooCrystalBall("CB", "Crystal Ball PDF", m, m0, sigma, alphal, nl, alphar, nr)
    expo = ROOT.RooExponential("Expo", "expo pdf", m, lam )
    pol1 = ROOT.RooPolynomial("Pol1", "linear", m, ROOT.RooArgList(lin))
    pol2 = ROOT.RooPolynomial("Pol2", "quadratic", m, ROOT.RooArgList(lin, quad))

    bins = np.linspace(-1,1,nbins+1)

    # create histogram
    inv_mass = rdf.Filter(f"{bins[i]} < CosThetaWrtBeam and CosThetaWrtBeam < {bins[i+1]} and Matter == {matter} and centrality > {cmin} ")\
        .Histo1D(("InvariantMassHistogram", f"{pt_min}"+"<p_{T}<"+f"{pt_max}, {str(bins[i])[:5]}"+"< cos(#theta_{beam})<"+f"{str(bins[i+1])[:5]}; m[GeV]", binnumber, 2.96, 3.04),"m")
    
    # Get RooDataHist
    dh = ROOT.RooDataHist("dh", "dataset with m", m, inv_mass.GetPtr())
    

    mframe = m.frame(Title=f"{pt_min}"+"<p_{T}<"+f"{pt_max}, {str(bins[i])[:5]} <"+"< cos(#theta_{beam})<"+f"< {str(bins[i+1])[:5]}")

    # Create Canvas
    if save_Output == True:
        c1 = ROOT.TCanvas()

        # Draw histogram 
        dh.plotOn(mframe, ROOT.RooFit.Name("dh"), ROOT.RooFit.MarkerSize(0))



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
        legend.Draw("same")   
        pave.Draw("same")
        pave2.Draw("same")
        
        c1.Draw()
        c1.SaveAs(f"{Output}.pdf")
    del dh
    del mframe
    del inv_mass
    
    return Nhyp, error

def GetNHypertritonCandidates(rdf, results, BDTEfficiency, pt_min, pt_max, cmin, matter, nbins, fitsigma=True, model="ExpoCB", save_Output=False):
    #parameters from Fitting MCs with RooCrystalBall
    matter3 = "M" if matter == "true" else "AM"

    fitsigma2 = "S" if fitsigma == True else ""
    	
    inplane2 = "IP" if inplane == True else "OP"

    BDTEff = "%.2f" %BDTEfficiency
    Output= Outputdir + f"RooFit{model}{fitsigma2}{pt_min}pt{pt_max}centr{cmin}_BDT{BDTEff}_{inplane2}_{matter3}"

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
    m = ROOT.RooRealVar("m", "m [GeV]", 2.96, 3.04)
    m0 = ROOT.RooRealVar("m0", "m0", 2.992,  2.990, 2.994) 
   
    frac = ROOT.RooRealVar("frac", "fraction", 0.5, 0, 1) 
    lam = ROOT.RooRealVar("Lambda", "slope of expo", -7, -50, 0) 
    quad = ROOT.RooRealVar("c2", "f(x)=const+lin*x+quad*x²", -7, -50, 0)
    lin = ROOT.RooRealVar("c1", "f(x)=const+lin*x+quad*x²", -7, -50, 50)

    bins = np.linspace(-1,1,nbins+1)
        

    # plots cos(theta*) distribution-----------------------------------------------
    if save_Output == True:
        c2 = ROOT.TCanvas()
        cos_theta = rdf.Filter(f"Matter == {matter}").Histo1D(("CosThetaHistogram", "cos(theta*)_wrt_beam", 100, -1, 1),"CosThetaWrtBeam")
        cos_theta.Draw()
        c2.Draw()
        c2.SaveAs(f"{Output}.pdf(")


    # create plots with signal fit-------------------------------------------------
    results = [CreatePlots(i, alphal_list, alphar_list, nl_list, nr_list, m, m0, sigma_list, lam, lin, quad, frac, nbins, cmin, matter, rdf, pt_min, pt_max, Output, model, fitsigma, save_Output) for i in range(nbins)]
    number_of_Hypertritons = [results[i][0] for i in range(nbins)]
    errorbars = [results[i][1] for i in range(nbins)]

    # plot number of Hypertritons (raw) -------------------------------------------
    if save_Output == True:
        c_Eff = ROOT.TCanvas()
        centered_bins = (bins[:-1] + bins[1:]) / 2
        
        
        h_raw = ROOT.TH1D("hist raw", f"Number of Hypertritons(raw) for {pt_min}"+"<p_{T}<"+f"{pt_max}; cos(#theta*) wrt beam", nbins, -1, 1)
        h_raw.FillN(nbins, np.array(centered_bins), np.array(number_of_Hypertritons) )
        for i in range(1,nbins+1,1):
            h_raw.SetBinError(i, errorbars[i-1])

        h_raw.GetYaxis().SetTitle("Counts per bin")

        ROOT.gStyle.SetTitleFontSize(0.045)

        h_raw.Draw()

        c_Eff.Draw()
        c_Eff.SaveAs(f"{Output}.pdf")

    # draw pt distr ---------------------------------------------------------------

    
    matter2 = 1.0 if matter == "true" else 0.0
    
   
    # load the files that have been created with HypertritonRestframeBoost.py
    rdf_gen = ROOT.RDF.FromCSV("generatedHypertritons_cut.csv")
    rdf_reco = ROOT.RDF.FromCSV("reconstructedHypertritons_cut.csv")
    print(rdf_gen.GetColumnNames())


    rdf_gen_cut = rdf_gen.Filter(f"{pt_min} < pt and pt < {pt_max} and Matter == {matter2} and Centrality > {cmin}")
    rdf_reco_cut = rdf_reco.Filter(f"{pt_min} < pt and pt < {pt_max} and Matter == {matter2} and Centrality > {cmin}")

    pthistgen = rdf_gen_cut.Histo1D(("ptHistogram", "p_{T} Distribution for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max}"+"; p_{T}[GeV]", 30, pt_min, pt_max),"pt")
    pthistreco = rdf_reco_cut.Histo1D(("ptHistogram", "p_{T} Distribution for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max}"+"; p_{T}[GeV]", 30, pt_min, pt_max),"pt")
        
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
    histgen = rdf_gen_cut.Histo1D(("cosThetaHistogram", "Cos(#theta_{beam}) Distribution for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max}"+"; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")
    histreco = rdf_reco_cut.Histo1D(("cosThetaHistogram", "Cos(#theta_{beam}) Distribution for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max}"+"; cos(#theta_{beam})", nbins, -1, 1),"cos_theta_beam")

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
    if save_Output == True:
        centered_bins = (bins[:-1] + bins[1:]) / 2

        h_cor = ROOT.TH1D("hist corrected", "Number of Hypertritons(corrected) for"+f"{pt_min}"+"<p_{T}<"+f"{pt_max};"+" cos(#theta_{beam})", nbins, -1, 1)
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
    with open(Outputdir + f"NhypCorrected{pt_min}pt{pt_max}centr{cmin}_{inplane2}_{matter3}.json", 'w') as f:
        json.dump(Nhyp_corrected, f, indent=2) 
    with open(Outputdir + f"NhypCorrectedErr{pt_min}pt{pt_max}centr{cmin}_{inplane2}_{matter3}.json", 'w') as f:
        json.dump(yerr, f, indent=2) 

   
    return Nhyp_corrected, yerr
    
def ProcessForOneBDTEff(rdfData, rdfMC, pt_min, pt_max, BDTEfficiency, nbins, matter, cmin, inplane, save_output):
    #-------------SETTINGS------------------
    # TODO: Add centrality cut, ggf.
    matter2 = "M" if matter == "true" else "AM"
    inplane2 = "IP" if inplane == True else "OP"


    # Specify Output File Paths
    OutputMC= Outputdir + f"RooFitCB{pt_min}pt{pt_max}_BDT{BDTEfficiency}_{matter2}" #for MCs we do not differentiate between matter and antimatter
    
    print("Filter BDTEfficiency = ", BDTEfficiency)

    # Apply BDTEfficiency Cut
    rdfData_new = rdfData.Filter(f"BDTEfficiency < {BDTEfficiency}")
    rdfMC_new = rdfMC.Filter(f"BDTEfficiency < {BDTEfficiency}")
    
    # Calculate the Parameters for CB from fitting the MCs
    results = DetermineCBParametersFromMC(rdfMC_new, OutputMC, pt_min, pt_max, matter, nbins, save_Output=save_output)


    print("-------------SETTINGS---------------")
    print(f"  Matter:             {matter}")
    print(f"  Selected pT range:  {pt_min}<pt<{pt_max}")
    print(f"  Centrality Cut:     centrality > {cmin}")
    print(f"  Number of bins:     {nbins}")
    print(f"  Event Plane Cut:    inplane: {inplane} ")
    print("------------------------------------")

    # do it for every model (Pol1, Pol2, Expo)+CB(with and without fixed sigma)
    if CalculateSystematicErrors:
        ExpoNhyp, ExpoYerr = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, pt_min, pt_max, cmin, matter, nbins, fitsigma=True, model="ExpoCB", save_Output=False)
        Pol1Nhyp, Pol1Yerr = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, pt_min, pt_max, cmin, matter, nbins, fitsigma=True, model="Pol1CB", save_Output=False)
        Pol2Nhyp, Pol2Yerr = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, pt_min, pt_max, cmin, matter, nbins, fitsigma=True, model="Pol2CB", save_Output=False)  
        ExpoNhypS, ExpoYerrS = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, pt_min, pt_max, cmin, matter, nbins, fitsigma=False, model="ExpoCB", save_Output=False)
        Pol1NhypS, Pol1YerrS = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, pt_min, pt_max, cmin, matter, nbins, fitsigma=False, model="Pol1CB", save_Output=False)
        Pol2NhypS, Pol2YerrS = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, pt_min, pt_max, cmin, matter, nbins, fitsigma=False, model="Pol2CB", save_Output=False)        

        Nhyp = [ExpoNhyp, ExpoNhypS, Pol1Nhyp, Pol1NhypS, Pol2Nhyp, Pol2NhypS]
        Yerr = [ExpoYerr, ExpoYerrS, Pol1Yerr, Pol1YerrS, Pol2Yerr, Pol2YerrS]
    else:
        ExpoNhyp, ExpoYerr = GetNHypertritonCandidates(rdfData_new, results, BDTEfficiency, pt_min, pt_max, cmin, matter, nbins, fitsigma=True, model="ExpoCB", save_Output=save_output)
        Nhyp = [ExpoNhyp]
        Yerr = [ExpoYerr]

    return Nhyp, Yerr

def GetResults(rdfData, rdfMC, matter, pt_min, pt_max, BDTEfficiencies, nbins, cmin, inplane, save_output):

    matter2 = "M" if matter == "true" else "AM"
    inplane2 = "IP" if inplane == True else "OP"

    # Create and Fill Dictionary
    Results = {}
    for i in range(nbins):
        Results[f"bin{i}"]=[]
        Results[f"bin{i}yerr"] = []

    for BDTEff in BDTEfficiencies:
        logging.info(f"BDTEfficiency: {BDTEff}")
        Nhyp, Yerr = ProcessForOneBDTEff(rdfData, rdfMC, pt_min, pt_max, BDTEff, nbins, matter, cmin, inplane, save_output)
        for i in range(nbins):
            Results[f"bin{i}"] += [x[i] for x in Nhyp]
            Results[f"bin{i}yerr"] += [x[i] for x in Yerr]

    # save to file
    with open(Outputdir + f"Results{pt_min}pt{pt_max}_BDT{BestBDTEff}_{inplane2}_{matter2}.json", 'w') as f:
            json.dump(Results, f, indent=2) 

    del Nhyp
    del Yerr
    return Results

def StoreSlope(ResultsMatter, ResultsAntimatter, nbins, OutputSlope, sample_size, cmin):
    np.random.seed(1773)

    inplane2 = "IP" if inplane == True else "OP"

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
    h = ROOT.TH1D("hist", "Nhyp matter-antimatter", nbins, -1, 1)

    for i in range(sample_size):
        diff = [x[i]-y[i] for x,y in zip(Nhyp_m, Nhyp_am)]
        diff_err = [np.sqrt(x[i]**2+y[i]**2) for x,y in zip(Yerr_m, Yerr_am)]

        bins = np.linspace(-1,1,nbins+1)
        centered_bins = (bins[:-1] + bins[1:]) / 2

        h.Reset()
        h.FillN(nbins, np.array(centered_bins), np.array(diff) )
        for j in range(1,nbins+1,1):
            h.SetBinError(j, diff_err[j-1])

        #TODO when running: select different Fitting functions, no polarisation:pol0 (=fit_func0), Spin1/2:pol1(=fit_func1), Spin3/2:pol2(=fit_func2)
        fit_func = fit_func1 if AngularDistribution == "Pol1" else fit_func2

        h.Fit(fit_func)

        pars = fit_func.GetParameters()
        errors = fit_func.GetParErrors()
        chisq = fit_func.GetChisquare()
            
        logging.info(f"Chi2:{chisq}, Prob:{ROOT.TMath.Prob(chisq, nbins-1)}")

        slopes.append(pars[1])
        slope_errors.append(errors[1])

    del h
    with open(Outputdir+f"Slopes{pt_min}pt{pt_max}centr{cmin}_{inplane2}.json", 'w') as f:
        json.dump(slopes, f, indent=2) 
    with open(Outputdir+f"SlopesErr{pt_min}pt{pt_max}centr{cmin}_{inplane2}.json", 'w') as f:
        json.dump(slope_errors, f, indent=2) 

    return slopes, slope_errors

def CreateSystematicErrorsHisto(slopes):
    # Create Histogram
    h = ROOT.TH1D("hist", "Nhyp matter-antimatter", 100, -100, 100)
    h.SetCanExtend(ROOT.TH1.kAllAxes)

    inplane2 = "IP" if inplane == True else "OP"

    # Fill Histogram with slopes
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
    c1.SaveAs(Outputdir+f"SlopeHistogram_centr{cmin}_{inplane2}.pdf")


#--------CODE-----------------


if __name__ == "__main__":
    logging.basicConfig(filename="SystematicErrors-"+datetime.datetime.today().strftime('%Y-%m-%d')+".log", level=logging.DEBUG)
    logging.info(f"{datetime.datetime.now()}")
 
    #---------------------SETTINGS------------------
    pt_min = 3.0
    pt_max = 6.4


    # BDT Cut
    BestBDTEff = 0.81

    #number of cos(theta) bins
    nbins = 5

    # number of bins in histogram per cos(theta) bin
    binnumber = 50

    #centrality cut
    cmin = 10 

    #inplane?
    inplane = True

    #Calculate Systematic Errors?
    CalculateSystematicErrors = False

    # trying to fit with?
    AngularDistribution = "Pol1"

    #save_output?
    save_output = True

    logging.info("------CONFIG-----------")
    logging.info(f"pt_min: {pt_min}")
    logging.info(f"pt_max: {pt_max}")
    logging.info(f"BDTCut: {BestBDTEff}")
    logging.info(f"Number of bins: {nbins}")
    logging.info(f"Centrality > {cmin}")
    logging.info(f"Inplane = {inplane}")
    logging.info("-----------------------")


    #-----------------PATHS FOR OUTPUT-------------------
    # Outfput file paths
    OutputSlope = "../Output/SlopeOutput"

    #Outputdir path
    Outputdir = "../Output/TESTTEST/"

    #Input dir path
    Inputdir = "../Output/"

    InputDataFileName = "3.0_pt_6.4_SystematicsEPangleData_3"
    InputMCFileName = "3.0_pt_6.4_SystematicsEPangleMC_3"


    #---------------------------------------------------------
    # Cut depending on if it is inplane or out of plane
    if inplane:
        EPanglecut = f"(Phi > EPangle-{np.pi/4} and Phi < EPangle + {np.pi/4}) or Phi > EPangle + {3*np.pi/4} or Phi < EPangle - {3*np.pi/4}"
    else:
        EPanglecut = f"not((Phi > EPangle-{np.pi/4} and Phi < EPangle + {np.pi/4}) or Phi > EPangle + {3*np.pi/4} or Phi < EPangle - {3*np.pi/4})"
    
    # BDTEfficiencies that are used
    if CalculateSystematicErrors == True:
        BDTEfficiencies = np.linspace(BestBDTEff-0.1, BestBDTEff+0.1, 21)
    else:
        BDTEfficiencies = np.array([BestBDTEff])

    print("BDTEfficiencies:", BDTEfficiencies)
    
    # Load dataframes
    rdfData = ROOT.RDataFrame("df", Inputdir+InputDataFileName).Filter(f"{pt_min} < pt and pt < {pt_max} and 1 < ct and ct < 35 and {EPanglecut}")
    
    rdfMC = ROOT.RDataFrame("df", Inputdir+InputMCFileName).Filter(f"{pt_min} < pt and pt < {pt_max} and 1 < ct and ct < 35 ")#and {EPanglecut}")

    # Extract Number of Hypertritons for Matter and Antimatter
    ResultsMatter = GetResults(
        rdfData, 
        rdfMC, 
        matter="true", 
        pt_min=pt_min, 
        pt_max=pt_max, 
        BDTEfficiencies=BDTEfficiencies, 
        nbins=nbins, 
        cmin=cmin, 
        inplane=inplane, 
        save_output=save_output
        )
    ResultsAntimatter = GetResults(
        rdfData, 
        rdfMC, 
        matter="false", 
        pt_min=pt_min, 
        pt_max=pt_max, 
        BDTEfficiencies=BDTEfficiencies, 
        nbins=nbins, 
        cmin=cmin, 
        inplane=inplane, 
        save_output=save_output
        )

    # with open(f"../Output/EPCutOutput/ResultsM.json", 'r') as f:
    #     ResultsMatter = json.load(f)

    # with open(f"../Output/EPCutOutput/ResultsAM.json", 'r') as f:
    #     ResultsAntimatter = json.load(f)

    if CalculateSystematicErrors == False:
        print("Number of Hypertritons succesfully calculated. Exiting...")
        exit()


    # Calculate the slopes for every case
    slopes, slope_errors = StoreSlope(ResultsMatter, ResultsAntimatter, nbins, OutputSlope, sample_size=10000, cmin=cmin)

    del ResultsMatter
    del ResultsAntimatter

    # with open(Outputdir+f"Slopes{pt_min}pt{pt_max}centr{cmin}.json", 'r') as f:
    #     slopes = json.load(f)
    # with open(Outputdir + f"SlopesErr{pt_min}pt{pt_max}centr{cmin}.json", 'r') as f:
    #     slope_errors = json.load(f)

    CreateSystematicErrorsHisto(slopes)

    