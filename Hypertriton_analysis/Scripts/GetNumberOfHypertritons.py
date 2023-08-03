import ROOT
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json


def CreatePlots(i, output="Output.pdf"):
    c = ROOT.TCanvas()
    inv_mass = rdf.Filter(f"{bins[i]} < cos_theta_beam and cos_theta_beam < {bins[i+1]} ")\
    .Histo1D(("InvariantMassHistogram", "3<p_{T}<6,  "+f"{str(bins[i])[:5]}"+"< cos(#theta_{beam})<"+f"{str(bins[i+1])[:5]}; m[GeV]", 100, 2.96, 3.04),"m")
    
    inv_mass.SetFillColor(ROOT.kBlue-10)
    inv_mass_copy = inv_mass.Clone()

    inv_mass.Draw()

    fit_func = ROOT.TF1("fit_func", "gausn(0) + [3]*expo(4)")
    fit_func.SetParLimits(1, 2.98, 2.998)
    fit_func.SetParLimits(2, 0.0001, 0.02)
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
        

    # plots cos(theta*) distribution
    c1 = ROOT.TCanvas()
    cos_theta = rdf.Histo1D(("CosThetaHistogram", "cos(theta*)_wrt_beam", 100, -1, 1),"cos_theta_beam")
    cos_theta.Draw()
    c1.Draw()
    c1.SaveAs("Output.pdf(")


    # create plots with signal fit
    for i in range(7):
        CreatePlots(i, "Output.pdf")

    # plot number of Hypertritons
    with open("./DetectorEfficiency_3_pt_6.json", 'r') as f:
        detector_efficiency = np.array(json.load(f))

    c_Eff = ROOT.TCanvas()
    centered_bins = (bins[:-1] + bins[1:]) / 2
    # hist = ROOT.TH1D("hist", "Number of Hypertritons per cos(#theta_{beam}) bin", 7, -1, 1)

    # hist.Draw()

    gr = ROOT.TGraphErrors(7, centered_bins, np.array(number_of_Hypertritons), np.array([1/7]*7), np.array(errorbars))
    gr.SetFillColor(ROOT.kOrange+6)
    gr.SetTitle("Number of Hypertritons for 3<p_{T}<6")#, "cos(#theta_{beam})", "Counts")
    gr.GetXaxis().SetTitle("cos(#theta_{beam})")
    gr.GetYaxis().SetTitle("Counts")
    gr.GetYaxis().SetRangeUser(0,400)
    # gr.SetMarkerColor(4)
    # gr.SetMarkerStyle(21)


    ROOT.gStyle.SetBarWidth(1.17)
    ROOT.gStyle.SetTitleFontSize(0.045)

    gr.Draw("ABZ")

    c_Eff.Draw()
    c_Eff.SaveAs("Output.pdf")

    #----------------------------------
    # plot efficiencies
    pt_min = 3
    pt_max = 6

    # load the files that have been created with HypertritonRestframeBoost.py
    rdf_gen = ROOT.RDF.FromCSV("generatedHypertritons_cut.csv")
    rdf_reco = ROOT.RDF.FromCSV("reconstrucedHypertritons_cut.csv")
    print(rdf_reco.GetColumnNames())

    pt_distr = ROOT.TCanvas()
    # NOTE: here the and must be used otherwise it will be interpreted wrong!!
    rdf_gen_cut = rdf_gen.Filter(f"{pt_min} < pt and pt < {pt_max}")
    rdf_reco_cut = rdf_reco.Filter(f"{pt_min} < pt and pt < {pt_max} ")


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
    c_Eff_2.SaveAs("Output.pdf")

    #--------------------------------------

    c_Eff_3 = ROOT.TCanvas()


    gr = ROOT.TGraphErrors(7, centered_bins, np.array(number_of_Hypertritons)/detector_efficiency, np.array([1/7]*7), np.array(errorbars)/detector_efficiency)
    gr.SetFillColor(ROOT.kOrange+6)
    gr.SetTitle("Number of Hypertritons(corrected) for 3<p_{T}<6")#, "cos(#theta_{beam})", "Counts")
    gr.GetXaxis().SetTitle("cos(#theta_{beam})")
    gr.GetYaxis().SetTitle("Counts")
    gr.GetYaxis().SetRangeUser(0,400)
    # gr.SetMarkerColor(4)
    # gr.SetMarkerStyle(21)


    ROOT.gStyle.SetBarWidth(1.17)
    ROOT.gStyle.SetTitleFontSize(0.045)

    gr.Draw("ABZ")

    c_Eff_3.Draw()
    c_Eff_3.SaveAs("Output.pdf)")