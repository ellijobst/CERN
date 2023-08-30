import numpy as np
import ROOT
import logging
import json

def procedure(rdf_filtered, bkgdf_filtered, gendf_filtered):
    
    nreco = rdf_filtered.Count().GetValue()*1/0.99
    ngen = gendf_filtered.Count().GetValue()

    RecoEfficiency = nreco/ngen

    binnumber = 85
    # Define Arrays
    eff = np.linspace(0.1, 0.95, binnumber)
    result = [GetSignificance(x, rdf_filtered, bkgdf_filtered, RecoEfficiency) for x in eff]
    sign = [x[0] for x in result]
    nhyp = [x[1] for x in result]
    nbkg = [x[2] for x in result]
    
    signxeff = eff*sign
    
    return eff, sign, signxeff, nhyp, nbkg

def GetSignificance(x, rdf_filtered, bkgdf_filtered, recoeff):
    print(x)
    
    # rdf_new = rdf_filtered.Filter(f"BDTEfficiency < {x}")
    bkgdf_new = bkgdf_filtered.Filter(f"BDTEfficiency < {x}")
    
    Nbkg = bkgdf_new.Count().GetValue()
    Nhyp = 2*278971416*2.6e-5*0.25*x*recoeff
    
    return (Nhyp/np.sqrt(Nhyp+Nbkg), Nhyp, Nbkg)



if __name__ == "__main__":
    logging.basicConfig(filename='ptOptimization-28082023.log', level=logging.DEBUG)
    logging.info("Creating Dataframes")
    
    # ROOT.EnableImplicitMT(8)
    
    # Load Dataframes
    file_path_rdf =  "../Output/DataframesForptOptimization/DataframeNew_*_pt_*_BDTEfficiency.root"
    rdf = ROOT.RDataFrame("df", file_path_rdf)
    
    file_path_bkgdf =  "../Output/DataframesForptOptimization/DataframeBkgNew_*_pt_*_BDTEfficiency.root"
    bkgdf = ROOT.RDataFrame("df", file_path_bkgdf)
    
    gendf = ROOT.RDataFrame("df", "../Output/DataframesForptOptimization/GeneratedHypertritons.root")
 
    # Define Variables for Fit
    m = ROOT.RooRealVar("m", "m [GeV]", 2.97, 3.01)
    m0 = ROOT.RooRealVar("m0", "m0", 2.992,  2.990, 2.994)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.002,  0.001, 0.003)
    alphal = ROOT.RooRealVar("alphal", "alpha L", 2.05,1, 2.1)
    nl = ROOT.RooRealVar("nl", "n L", 4.7,  3.9, 6.5)
    alphar = ROOT.RooRealVar("alphar", "alpha right", 2.08, 1, 2.2)
    nr = ROOT.RooRealVar("nr", "n right", 3.7, 2.9 , 6.5)


    # Define FitFunction
    fit_func = ROOT.RooCrystalBall("fit_func", "Crystal Ball PDF", m, m0, sigma, alphal, nl, alphar, nr)


    # Define Factor for Number of Hypertriton Candidates
    Nhyp_factor = 2 * 278971416 * 2.6e-5 * 0.25

    # Define pt ranges:
    max_int_size = 7#GeV/c
    min_int_size = 1.6#GeV/c #before: 1.6
    step_size = 0.2#GeV/c
    maximumPt=9#GeV/c
    minimumPt=2#GeV/c
    
    # Define Output Dictionary
    Significances = {}
    
    # Loop through pt ranges
    for i in range(int((max_int_size-min_int_size)/step_size)+1):
        
        int_size = max_int_size - i*step_size
        logging.warning(f"pt interval size: {int_size}")
        
        for j in range(i+1):
            
            pt_min = minimumPt+step_size*j
            pt_max = pt_min + int_size
            logging.warning(f"{pt_min } < pt < {pt_max}")
            
            # Cut on pt
            rdf_filtered = rdf.Filter(f"{pt_min} < pt  and pt <= {pt_max}")
            bkgdf_filtered = bkgdf.Filter(f"{pt_min} < pt  and pt <= {pt_max}")
            gendf_filtered = gendf.Filter(f"{pt_min} < pt  and pt <= {pt_max}")
            
            # Create Histogram for Fit
            inv_mass = rdf_filtered.Histo1D(("InvariantMassHistogram", "hist; m[GeV]", 80, 2.97, 3.01),"m")
            dh = ROOT.RooDataHist("dh", "dataset with m", m, inv_mass.GetPtr())
            
            # Fit
            fit_func.fitTo(dh)

            m_min = m0.getVal() - 3*sigma.getVal()
            m_max = m0.getVal() + 3*sigma.getVal()

            # Filter Dataframes
            bkgdf_filtered = bkgdf_filtered.Filter(f"m > {m_min} and m < {m_max} and 1 < ct and ct < 35")
            rdf_filtered = rdf_filtered.Filter(f" 1 < ct and ct < 35")
            
            # Calculate Efficiency, Significance, etc.
            eff, sign, signxeff, nhyp, nbkg =  procedure(rdf_filtered, bkgdf_filtered, gendf_filtered)
            
            # Append the Significance to the Dictionary
            BestBDTSign = sign[signxeff.tolist().index(max(signxeff.tolist()))]
            logging.info("BestBDTSign: {BestBDTSign}")
            Significances[f"{pt_min}<pt<{pt_max}"] = BestBDTSign

            inv_mass = None
            rdf_filered = None
            bkgdf_filtered = None
            gendf_filtered = None
            
            
    with open('ptSignificancesNew.json', 'w') as f:
        # write the dictionary to the file in JSON format
        json.dump(Significances, f)