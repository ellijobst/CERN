import ROOT
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json

#matter? options: "true", "false"  all lowercase!!
matter = "true"


def CreatePlots(i, output="Output.pdf"):
    c = ROOT.TCanvas()
    inv_mass = rdf.Filter(f"{bins[i]} < cos_theta_beam and cos_theta_beam < {bins[i+1]} and Matter == {matter} ")\
    .Histo1D(("InvariantMassHistogram", "3<p_{T}<6,  "+f"{str(bins[i])[:5]}"+"< cos(#theta_{beam})<"+f"{str(bins[i+1])[:5]}; m[GeV]", 100, 2.96, 3.04),"m")
    
    inv_mass.SetFillColor(ROOT.kBlue-10)
    inv_mass_copy = inv_mass.Clone()

    inv_mass.Draw()

    fit_func = ROOT.TF1("fit_func", "gausn(0) + [3]*expo(4)")
    fit_func.SetParLimits(1, 2.98, 2.998)
    fit_func.SetParLimits(2, 0.0001, 0.002)
    fit_func.SetNpx(1000)
    
    inv_mass.Fit(fit_func)
    pars = fit_func.GetParameters()
    errors = fit_func.GetParErrors()
    
    c.SetLogy()

    # add Background and Signal
    bkg_func = ROOT.TF1("bkg_func", "[0]*expo(1)", 2.96, 3.04)
    bkg_func.SetLineColor(ROOT.kGreen)
    bkg_func.SetNpx(1000)
    bkg_func.Draw("Same")
    bkg_func.SetParameters(pars[3], pars[4], pars[5])

    signal_func = ROOT.TF1("signal_func", "gausn(0)", 2.96, 3.04)
    signal_func.SetLineColor(ROOT.kBlue)
    signal_func.SetNpx(1000)
    signal_func.Draw("Same")
    signal_func.SetParameters(pars[0], pars[1], pars[2])


    # add legend
    legend = ROOT.TLegend(.12, .8, .49, .935)
    legend.SetBorderSize(1)
    legend.AddEntry(inv_mass_copy, "Data", "f")
    legend.AddEntry(fit_func, "Fit = A gausn(x, #mu, #sigma) + B expo(x, y, #lambda)", "l")
    legend.AddEntry(signal_func, "Signal = A gausn(x, #mu, #sigma)", "l")
    legend.AddEntry(bkg_func, "Background = B expo(x, y, #lambda))", "l")
    legend.SetTextSize(0.025)
    legend.Draw()
    
    # get output parameters
    pars = fit_func.GetParameters()
    parerrors = fit_func.GetParErrors()
    parameters = [pars[i] for i in range(6)]
    errors = [parerrors[i] for i in range(6)]
    output = [f"{parameters[i]} +- {errors[i]}" for i in range(6)]
    
    number_of_Hypertritons.append(pars[0]/(3.04-2.96)*100)
    errorbars.append(parerrors[0]/(3.04-2.96)*100)

    # Add Fit Parameters to canvas
    pave = ROOT.TPaveText(0.5, 0.7, 0.77, 0.935, "NDC")
    pave.SetBorderSize(1)
    pave.SetFillColor(ROOT.kWhite)
    pave.SetTextFont(42)

    t0 = pave.AddText(f"A = {str(parameters[0])[:10]} #pm {str(errors[0])[:10]}")
    t1 = pave.AddText(f"#mu = {str(parameters[1])[:10]} #pm {str(errors[1])[:10]}")
    t2 = pave.AddText(f"#sigma = {str(parameters[2])[:10]} #pm {str(errors[2])[:10]}")
    t3 = pave.AddText(f"B = {str(parameters[3])[:10]} #pm {str(errors[3])[:10]}")
    t4 = pave.AddText(f"y = {str(parameters[4])[:10]} #pm {str(errors[4])[:10]}")
    t5 = pave.AddText(f"#lambda = {str(parameters[5])[:10]} #pm {str(errors[5])[:10]}")

    pave.Draw()

    c.Draw()
    c.SaveAs("Output.pdf")

if __name__ == "__main__":
    # load data
    rdf = ROOT.RDataFrame("df", "SelectedDataFrame_3_pt_6.root")

    number_of_Hypertritons = []
    errorbars = []

    bins = np.linspace(-1,1,8)
        

    # plots cos(theta*) distribution-----------------------------------------------
    c1 = ROOT.TCanvas()
    cos_theta = rdf.Filter(f"Matter == {matter}").Histo1D(("CosThetaHistogram", "cos(theta*)_wrt_beam", 100, -1, 1),"cos_theta_beam")
    cos_theta.Draw()
    c1.Draw()
    c1.SaveAs("Output.pdf(")


    # create plots with signal fit-------------------------------------------------
    for i in range(7):
        CreatePlots(i, "Output.pdf")

    # plot number of Hypertritons (raw) -------------------------------------------
    
    c_Eff = ROOT.TCanvas()
    centered_bins = (bins[:-1] + bins[1:]) / 2


    gr = ROOT.TGraphErrors(7, centered_bins, np.array(number_of_Hypertritons), np.array([1/7]*7), np.array(errorbars))
    gr.SetFillColor(ROOT.kOrange+6)
    gr.SetTitle("Number of Hypertritons(raw) for 3<p_{T}<6")#, "cos(#theta_{beam})", "Counts")
    gr.GetXaxis().SetTitle("cos(#theta_{beam})")
    gr.GetYaxis().SetTitle("Counts")
    gr.GetYaxis().SetRangeUser(0,150)
    # gr.SetMarkerColor(4)
    # gr.SetMarkerStyle(21)

    ROOT.gStyle.SetBarWidth(1.17)
    ROOT.gStyle.SetTitleFontSize(0.045)

    gr.Draw("ABZ")

    c_Eff.Draw()
    c_Eff.SaveAs("Output.pdf")

    # draw pt distr ---------------------------------------------------------------

    pt_min = 3
    pt_max = 6
    
    matter2 = 0.0
    if matter == "true":
        matter2 = 1.0
    

    # load the files that have been created with HypertritonRestframeBoost.py
    rdf_gen = ROOT.RDF.FromCSV("generatedHypertritons_cut.csv")
    rdf_reco = ROOT.RDF.FromCSV("reconstructedHypertritons_cut.csv")
    print(rdf_reco.GetColumnNames())

    pt_distr = ROOT.TCanvas()
    # NOTE: here the and must be used otherwise it will be interpreted wrong!!
    # NOTE: in the original files it is "matter" for gen, "Matter" for reco. In the *_cut.csv files it is "Matter" for both!
    rdf_gen_cut = rdf_gen.Filter(f"{pt_min} < pt and pt < {pt_max} and Matter == {matter2}")
    rdf_reco_cut = rdf_reco.Filter(f"{pt_min} < pt and pt < {pt_max} and Matter == {matter2}")

    pthistgen = rdf_gen_cut.Histo1D(("ptHistogram", "p_{T} Distribution for 3<p_{T}<6; p_{T}[GeV]", 30, 3, 6),"pt")
    pthistreco = rdf_reco_cut.Histo1D(("ptHistogram", "p_{T} Distribution for 3<p_{T}<6; p_{T}[GeV]", 30, 3, 6),"pt")

   
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
    pt_distr.SaveAs("Output.pdf")
    
    #plot efficiency for whole pt range -------------------------------------------
    histgen = rdf_gen_cut.Histo1D(("cosThetaHistogram", "Cos(#theta_{beam}) Distribution for 3<p_{T}<6; cos(#theta_{beam})", 7, -1, 1),"cos_theta_beam")
    histreco = rdf_reco_cut.Histo1D(("cosThetaHistogram", "Cos(#theta_{beam}) Distribution for 3<p_{T}<6; cos(#theta_{beam})", 7, -1, 1),"cos_theta_beam")
    histgen.SetFillColor(ROOT.kBlue-10)
    histreco.SetFillColor(ROOT.kRed-10)
    histgen.Draw("hist")
    histreco.Draw("Same hist")
    histreco.SetDefaultSumw2()
    histgen.GetYaxis().SetRangeUser(0,250000)

    pt_distr.Draw()
    pt_distr.SaveAs("Output.pdf")

    d = ROOT.TCanvas()
    efficiency = histreco.Clone()
    efficiency.Sumw2()
    efficiency.Divide(histreco.GetPtr(), histgen.GetPtr(), 1, 1, "B")
    efficiency.Draw("E")
    efficiency.SetTitle("Detector Efficiency for 3<p_{T}<6")
    efficiency.GetYaxis().SetTitle("Efficiency")
    efficiency.GetYaxis().SetRangeUser(0.2, 0.8)
    d.Draw()
    d.SaveAs("Output.pdf")

    count = [efficiency.GetBinContent(i) for i in range(1,8,1)]
    eff_err = [efficiency.GetBinError(i) for i in range(1,8,1)]
    print(count)


    #plot number of Hypertritons (corrected)---------------------------------------
    yerr = [np.sqrt((errorbars[i]/count[i])**2+((number_of_Hypertritons[i]*eff_err[i])/count[i]**2)**2) for i in range(7)]

    c_Eff = ROOT.TCanvas()
    centered_bins = (bins[:-1] + bins[1:]) / 2


    gr = ROOT.TGraphErrors(7, centered_bins, np.array(number_of_Hypertritons)/np.array(count),\
                        np.array([1/7]*7), np.array(yerr))
    gr.SetFillColor(ROOT.kOrange+6)
    gr.SetTitle("Number of Hypertritons(corrected) for 3<p_{T}<6")#, "cos(#theta_{beam})", "Counts")
    gr.GetXaxis().SetTitle("cos(#theta_{beam})")
    gr.GetYaxis().SetTitle("Counts")
    gr.GetYaxis().SetRangeUser(0,300)

    ROOT.gStyle.SetBarWidth(1.17)
    ROOT.gStyle.SetTitleFontSize(0.045)

    gr.Draw("ABZ")

    c_Eff.Draw()
    c_Eff.SaveAs("Output.pdf")

    # plot efficiencies for different pt ranges------------------------------------

    '''
    with open("./DetectorEfficiency_3_pt_4.json", 'r') as f:
        detector_efficiency34 = np.array(json.load(f))
    with open("./DetectorEfficiency_4_pt_5.json", 'r') as f:
        detector_efficiency45 = np.array(json.load(f))
    with open("./DetectorEfficiency_5_pt_6.json", 'r') as f:
        detector_efficiency56 = np.array(json.load(f))
    
    c_Eff_2 = ROOT.TCanvas()

    # plot efficiency for the whole range
    
    gr_3_pt_6 = ROOT.TGraphErrors(7, centered_bins, detector_efficiency, np.array([1/7]*7), np.array([0.00001]*7))
    gr_3_pt_6.GetXaxis().SetTitle("cos(#theta_{beam})")
    gr_3_pt_6.GetYaxis().SetTitle("Efficency")
    gr_3_pt_6.SetTitle("Detector Efficiencies for different p_{T} bins")
    gr_3_pt_6.GetXaxis().SetRangeUser(-1,1)
    gr_3_pt_6.GetYaxis().SetRangeUser(0.2,0.8)

    gr_3_pt_6.SetMarkerStyle(8)
    gr_3_pt_6.SetMarkerSize(0.5)
    gr_3_pt_6.Draw("APZ")

    gr_3_pt_4 = ROOT.TGraphErrors(7, centered_bins, detector_efficiency34, np.array([1/7]*7), np.array([0.00001]*7))
    gr_3_pt_4.SetMarkerStyle(8)
    gr_3_pt_4.SetMarkerSize(0.5)
    gr_3_pt_4.SetMarkerColor(ROOT.kBlue-4)
    gr_3_pt_4.SetLineColor(ROOT.kBlue-4)
    gr_3_pt_4.Draw("PZ")

    gr_4_pt_5 = ROOT.TGraphErrors(7, centered_bins, detector_efficiency45, np.array([1/7]*7), np.array([0.00001]*7))
    gr_4_pt_5.SetMarkerStyle(8)
    gr_4_pt_5.SetMarkerSize(0.5)
    gr_4_pt_5.SetMarkerColor(ROOT.kGreen-2)
    gr_4_pt_5.SetLineColor(ROOT.kGreen-2)
    gr_4_pt_5.Draw("PZ")

    gr_5_pt_6 = ROOT.TGraphErrors(7, centered_bins, detector_efficiency56, np.array([1/7]*7), np.array([0.00001]*7))
    gr_5_pt_6.SetMarkerStyle(8)
    gr_5_pt_6.SetMarkerSize(0.5)
    gr_5_pt_6.SetMarkerColor(ROOT.kMagenta-4)
    gr_5_pt_6.SetLineColor(ROOT.kMagenta-4)
    gr_5_pt_6.Draw("PZ")

    # add legend
    legend2 = ROOT.TLegend(.73, .7, .86, .85)
    legend2.SetBorderSize(1)
    legend2.AddEntry(gr_3_pt_6, "3 < p_{T} < 6", "l")
    legend2.AddEntry(gr_3_pt_4, "3 < p_{T} < 4", "l")
    legend2.AddEntry(gr_4_pt_5, "4 < p_{T} < 5", "l")
    legend2.AddEntry(gr_5_pt_6, "5 < p_{T} < 6", "l")
    legend2.SetTextSize(0.025)
    legend2.Draw()

    c_Eff_2.SetGrid()

    c_Eff_2.Draw()
    c_Eff_2.SaveAs("Output.pdf)")
    '''
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


    histgen34 = rdf_gen_cut_34.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", 7, -1, 1),"cos_theta_beam")
    histreco34 = rdf_reco_cut_34.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", 7, -1, 1),"cos_theta_beam")
    histreco34.SetDefaultSumw2()

    efficiency34 = histreco34.Clone()
    efficiency34.Sumw2()
    efficiency34.Divide(histreco34.GetPtr(), histgen34.GetPtr(), 1, 1, "B")
    efficiency34.SetTitle("Detector Efficiency for different p_{T} ranges.")
    efficiency34.SetLineColor(9)
    efficiency34.Draw("E")
    efficiency34.GetYaxis().SetTitle("Efficiency")

    histgen45 = rdf_gen_cut_45.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", 7, -1, 1),"cos_theta_beam")
    histreco45 = rdf_reco_cut_45.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", 7, -1, 1),"cos_theta_beam")
    histreco45.SetDefaultSumw2()

    efficiency45 = histreco45.Clone()
    efficiency45.Sumw2()
    efficiency45.Divide(histreco45.GetPtr(), histgen45.GetPtr(), 1, 1, "B")
    efficiency45.Draw("Same E")
    efficiency45.SetLineColor(8)

    histgen56 = rdf_gen_cut_56.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", 7, -1, 1),"cos_theta_beam")
    histreco56 = rdf_reco_cut_56.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", 7, -1, 1),"cos_theta_beam")
    histreco56.SetDefaultSumw2()

    efficiency56 = histreco45.Clone()
    efficiency56.Sumw2()
    efficiency56.Divide(histreco56.GetPtr(), histgen56.GetPtr(), 1, 1, "B")
    efficiency56.Draw("Same E")
    efficiency56.SetLineColor(6)

    histgen36 = rdf_gen_cut_36.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", 7, -1, 1),"cos_theta_beam")
    histreco36 = rdf_reco_cut_36.Histo1D(("cosThetaHistogram", "Detector Efficiencies for different p_{T} ranges; cos(#theta_{beam})", 7, -1, 1),"cos_theta_beam")
    histreco36.SetDefaultSumw2()

    efficiency36 = histreco36.Clone()
    efficiency36.Sumw2()
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
    c4.SaveAs("Output.pdf)")