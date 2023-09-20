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

def AddBDTEfficiency():
    ROOT.gInterpreter.ProcessLine(cppcode)
    
    # Create Dataframes
    rdfMC = ROOT.RDataFrame("df", InfilenameMC)
    rdfData = ROOT.RDataFrame("df", InfilenameData)
    rdfBkg = ROOT.RDataFrame("df", InfilenameBkg)


    # Set Number of bins for BDTEfficiency Histogram
    bins = 5000

    # define Histogram that converts BDTScore into BDTEfficiency with MCs
    hist = rdfMC.Histo1D(("BDTScoreHist", "BDT Score", bins, -15, 15),"model_output")
    hist_cum = hist.GetCumulative()
    maximum = hist_cum.GetMaximum()
    hist_cum.Scale(1/(maximum))
    hist_cum.Scale(-1)

    for bin in range(1, bins+1):
        hist_cum.AddBinContent(bin)


    # Add BDTEfficiencies to the Dataframes
    rdfMC = rdfMC.Define("BDTEfficiency", ROOT.GetBDTEfficiency(hist_cum), ["model_output"])
    rdfData = rdfData.Define("BDTEfficiency", ROOT.GetBDTEfficiency(hist_cum), ["model_output"])
    rdfBkg = rdfBkg.Define("BDTEfficiency", ROOT.GetBDTEfficiency(hist_cum), ["model_output"])

    # Cut on BDTEfficiency and save modified DataFrames to .root file
    treeName = "df"

    # save dataframes to root file and cut on BDTEfficiency
    rdfMC = rdfMC.Filter(f"BDTEfficiency < 0.99").Snapshot(treeName, OutfileNameMC)
    rdfData = rdfData.Filter(f"BDTEfficiency < 0.99").Snapshot(treeName, OutfileNameData)
    rdfBkg = rdfBkg.Filter(f"BDTEfficiency < 0.99").Snapshot(treeName, OutfileNameBkg)

    print("Root df saved.")

if __name__ == "__main__":
    # Define input file names 
    InfilenameMC = f"DataframeMCEPangle_3.0_pt_6.4_3.root"
    InfilenameData = f"DataframeDataEPangle_3.0_pt_6.4_3.root"
    InfilenameBkg = "Backgroundfile.root"

    # Define output file names
    OutfileNameMC = f"DataframeMCEPangle_3.0_pt_6.4_BDTEfficiency_3.root"
    OutfileNameData = f"DataframeDataEPangle_3.0_pt_6.4_BDTEfficiency_3.root"
    OutfileNameBkg = f"DataframeBkgEPangle_3.0_pt_6.4_BDTEfficiency_3.root"

    AddBDTEfficiency()

