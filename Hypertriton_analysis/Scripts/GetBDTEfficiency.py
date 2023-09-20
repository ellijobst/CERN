import ROOT

cppcode = """
#include <ROOT/RDataFrame.hxx>
#include <TH1.h>

float GetAll(const float &BDT_score, TH1D &hist_cum) {
    int bins = 5000;

    // Find and return the value at the specified BDT score bin
    int bin = hist_cum.FindBin(BDT_score);
    return hist_cum.GetBinContent(bin);
}
struct GetBDTEfficiency {
    TH1D *fHist1D;
    GetBDTEfficiency(TH1D *h) : fHist1D(h) {}

    float operator()(const float &model_output) {
       return GetAll(model_output, *fHist1D);
    }
};
"""

ROOT.gInterpreter.ProcessLine(cppcode)

for pt2 in range(4,18,1):
    
    # Loop through pt in steps of 0.5GeV
    pt=pt2/2
    print(pt)
    
    
    # Define input file names 
    # NOTE: these Files had the model trained and then applied to it for 0.5GeV pt ranges
#     filenamebkg = f"DataframesForBDTEfficiency/DataframeBkg_{pt}_pt_{pt+0.5}.root"
    filename = f"DataframeMCEPangle_3.0_pt_6.4_3.root"
    filenamedata = f"DataframeDataEPangle_3.0_pt_6.4_3.root"
    
    # Create Dataframes
    rdf = ROOT.RDataFrame("df", filename)
#     rdfbkg = ROOT.RDataFrame("df", filenamebkg)
    rdfdata = ROOT.RDataFrame("df", filenamedata)
    
    # Set Number of bins, and get Number of generated Hypertriton Candidates for this pt range
#     Ngenpt = ngen[pt2-4]
    bins = 5000
    
    # define Histogram that converts BDTScore into BDTEfficiency with MCs
    hist = rdf.Histo1D(("BDTScoreHist", "BDT Score", bins, -15, 15),"model_output")
    hist_cum = hist.GetCumulative()
    maximum = hist_cum.GetMaximum()
    hist_cum.Scale(1/(maximum))
    hist_cum.Scale(-1)
    
    for bin in range(1, bins+1):
        hist_cum.AddBinContent(bin)
        
#     hist_cum.Scale(maximum/Ngenpt)
#     h = hist_cum

    # Define BDTEfficiencies for the Dataframes
    # rdf = rdf.Define("BDTEfficiency", ROOT.GetBDTEfficiency(hist_cum), ["model_output"])
    rdfdata = rdfdata.Define("BDTEfficiency", ROOT.GetBDTEfficiency(hist_cum), ["model_output"])
    # rdfbkg = rdfbkg.Define("BDTEfficiency", ROOT.GetBDTEfficiency(hist_cum), ["model_output"])
    
    # Cut on BDTEfficiency and save modified DataFrames to .root file
    treeName = "df"

    # fileName = f"DataframeMCEPangle_3.0_pt_6.4_BDTEfficiency_3.root"
#     fileNamebkg = f"DataframesForBDTEfficiency/DataframeBkgNew_{pt}_pt_{pt+0.5}_BDTEfficiency.root"
    fileNamedata = f"DataframeDataEPangle_3.0_pt_6.4_BDTEfficiency_3.root"
    
    # rdf = rdf.Filter(f"BDTEfficiency < 0.99").Snapshot(treeName, fileName)
    rdfdata = rdfdata.Filter(f"BDTEfficiency < 0.99").Snapshot(treeName, fileNamedata)
#     rdfbkg = rdfbkg.Filter(f"BDTEfficiency < 0.99").Snapshot(treeName, fileNamebkg, {"pt", "Matter",  "centrality", "ct", "m", "model_output", "BDTEfficiency"})

    print("Root df saved.")
    break
    