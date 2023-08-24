import numpy as np
import ROOT
import json


def ProcedureForPtRange(pt_min, pt_max, rdf, bkgdf):
    bins = 5000
    BDT_bins = 100

    # filter for pt range
    bkgdf_filtered = bkgdf.Filter(f"{pt_min} < pt and pt < {pt_max}")
    rdf_filtered = rdf.Filter(f"{pt_min} < pt and pt < {pt_max} and -15 < model_output and 15 > model_output")
    
    # get minimum and maximum BDT Score
    BDT_min = int(np.ceil(rdf_filtered.Min("model_output").GetValue()))
    BDT_max = int(np.floor(rdf_filtered.Max("model_output").GetValue()))
    print(f"{BDT_min} < BDTScore < {BDT_max}")
    
    # specifiy BDT cuts
    BDTCuts = np.linspace(BDT_min, BDT_max, BDT_bins+1)
    
    # get Efficiency
    hist = rdf_filtered.Histo1D(("BDTscore", "BDT Score", bins, BDT_min, BDT_max),"model_output")
    hist_cum = hist.GetCumulative()
    hist_cum.Scale(1/hist_cum.GetMaximum())
    hist_cum.Scale(-1)
    
    for bin in range(1, bins+1):
        hist_cum.AddBinContent(bin)

    # create empty histogram for Significance
    Signhisto = ROOT.TH1D("Signhisto", "Significance", BDT_bins, BDT_min, BDT_max)

    #loop through BDT cuts
    for i in range(BDT_bins):

        BDTCut = BDTCuts[i]
        print(BDTCut)

        # Obtain Efficiency from histogram
        BDTEfficiency = hist_cum.GetBinContent(hist_cum.FindBin(BDTCut))

        # calculate Number of Hypertritron Candidates
        Nhyp = Nhyp_factor * BDTEfficiency

        # Obtain Number of Background Candidates
        # Filtering: pt range, m range (determined by fit), ct
        Nbkg = bkgdf_filtered.Filter(f"model_output > {BDTCut}").Count().GetValue()

        # calculate Significance
        Sign = Nhyp/np.sqrt(Nhyp+Nbkg)

        # Fill Histogram
        Signhisto.SetBinContent(i + 1, Sign)
    
    
    
    Effhisto = hist_cum.Rebin(int(bins/BDT_bins))
    Effhisto.Scale(1/Effhisto.GetMaximum())

    Signeffhisto = Effhisto * Signhisto
    
    # get Significance and Best BDT cut
    return Signhisto.GetBinContent(Signeffhisto.GetMaximumBin()), Signeffhisto.GetXaxis().GetBinCenter(Signeffhisto.GetMaximumBin())
    

    
if __name__ == "__main__":
    
    # load Dataframe
    rdf = ROOT.RDataFrame("df", [
        "DataframesForBDTEfficiency/Dataframe_2.0_pt_2.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_2.5_pt_3.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_3.0_pt_3.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_3.5_pt_4.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_4.0_pt_4.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_4.5_pt_5.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_5.0_pt_5.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_5.5_pt_6.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_6.0_pt_6.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_6.5_pt_7.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_7.0_pt_7.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_7.5_pt_8.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_8.0_pt_8.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/Dataframe_8.5_pt_9.0_BDTEfficiency.root" ,
    ])
    bkgdf = ROOT.RDataFrame("df", [
        "DataframesForBDTEfficiency/DataframeBkg_2.0_pt_2.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_2.5_pt_3.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_3.0_pt_3.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_3.5_pt_4.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_4.0_pt_4.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_4.5_pt_5.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_5.0_pt_5.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_5.5_pt_6.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_6.0_pt_6.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_6.5_pt_7.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_7.0_pt_7.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_7.5_pt_8.0_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_8.0_pt_8.5_BDTEfficiency.root" ,
        "DataframesForBDTEfficiency/DataframeBkg_8.5_pt_9.0_BDTEfficiency.root" ,
    ])

    # Calculate Nhyp factor outside the loop
    Nhyp_factor = 2 * 278971416 * 2.6e-5 * 0.25

    #Filter Background
    bkgdf = bkgdf.Filter(" m > 2.987048438549096 and m < 2.99568866477378 and 1 < ct and ct < 35")
    rdf = rdf.Filter(f" 1 < ct and ct < 35")


    # pt ranges:
    max_int_size = 7#GeV/c
    min_int_size = 1.6#GeV/c #before: 1.6
    step_size = 0.2#GeV/c
    maximumPt=9#GeV/c
    minimumPt=2#GeV/c

    pt_ranges=[]
    significances=[]

    for i in range(int((max_int_size-min_int_size)/step_size)+1):
        int_size = max_int_size - i*step_size
        print("pt interval size", int_size)
        for j in range(i+1):
            pt_min = minimumPt+step_size*j
            pt_max = pt_min + int_size
            print(f"{pt_min } < pt < {pt_max}")
            pt_ranges.append((pt_min, pt_max))
            significances.append(ProcedureForPtRange(pt_min, pt_max, rdf, bkgdf))
            
    with open(f"ptSignificances{step_size}GeV.json", 'w') as f:
        json.dump(significances, f, indent=2) 

    with open(f"ptIntervals{step_size}GeV.json", 'w') as f:
        json.dump(pt_ranges, f, indent=2) 

    Best_pT = pt_ranges[significances.index(max(significances))]
    print("best pt interval:", Best_pT)