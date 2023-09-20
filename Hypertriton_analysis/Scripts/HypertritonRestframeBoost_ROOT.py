import ROOT

def HypertritonRestFrameBoost(df_gen, fileName):

    # specify masses of decay daughters
    mHe3 = 2.80839 #GeV
    mPi = 0.139570 #GeV

    # define lorentzvectors
    df_gen = df_gen.Define("LorentzvectorLABHe3", f"ROOT::Math::PxPyPzMVector(pxHe3, pyHe3, pzHe3, {mHe3})")
    df_gen = df_gen.Define("LorentzvectorLABPi", f"ROOT::Math::PxPyPzMVector(pxPi, pyPi, pzPi, {mPi})")

    # get Hypertriton data
    df_gen = df_gen.Define("pxHyp", "pxPi+pxHe3")
    df_gen = df_gen.Define("pyHyp", "pyPi+pyHe3")
    df_gen = df_gen.Define("pzHyp", "pzPi+pzHe3")

    # get energies
    df_gen = df_gen.Define("EHe3", f"sqrt(pxHe3*pxHe3+pyHe3*pyHe3+pzHe3*pzHe3+{mHe3}*{mHe3})")
    df_gen = df_gen.Define("EPi", f"sqrt(pxPi*pxPi+pyPi*pyPi+pzPi*pzPi+{mPi}*{mPi})")
    df_gen = df_gen.Define("EHyp", "EPi + EHe3")

    # get mass hyp
    df_gen = df_gen.Define("mHyp", "sqrt(EHyp*EHyp - pxHyp*pxHyp - pyHyp*pyHyp - pzHyp*pzHyp)")

    # get transverse momentum
    df_gen = df_gen.Define("pTHe3", "LorentzvectorLABHe3.Pt()")
    df_gen = df_gen.Define("pTPi", "LorentzvectorLABPi.Pt()")


    # boost in CMS of hyp
    df_gen = df_gen.Define("LorentzvectorLABHyp", f"ROOT::Math::PxPyPzMVector(pxHyp, pyHyp, pzHyp, mHyp)")

    df_gen = df_gen.Define("LorentzvectorCMSHyp", "BoostToCMHyp(LorentzvectorLABHyp, LorentzvectorLABHyp)")
    df_gen = df_gen.Define("LorentzvectorCMSHe3", "BoostToCMHyp(LorentzvectorLABHyp, LorentzvectorLABHe3)")
    df_gen = df_gen.Define("LorentzvectorCMSPi", "BoostToCMHyp(LorentzvectorLABHyp, LorentzvectorLABPi)")


    # calculate cos theta* wrt beam
    df_gen = df_gen.Define("CosThetaWrtBeam", "GetCosThetaBeam(LorentzvectorCMSHe3)")

    # calculate cos theta* wrt to pHyp
    df_gen = df_gen.Define("CosThetaWrtpHyp", "GetCosTheta(LorentzvectorCMSHe3, LorentzvectorLABHyp)")

    # calculate average event plane angle
    df_gen = df_gen.Define("EPangle", "(EPangleV0A + EPangleV0C)/2")

    treeName = "df"
    fileName = fileName
    df_gen.Snapshot(treeName, fileName, {"Phi", "EPangle", "CosThetaWrtBeam", "BDTEfficiency", "CosThetaWrtpHyp", "centrality", "ct", "Matter", "pt", "Rapidity", "m", "PseudoRapidityHe3"})

if __name__ == "__main__":
    ROOT.gInterpreter.ProcessLine('#include "FunctionsForHypertritonRestFrameBoost.h"')

    file_path_rdf =  "DataframeData_*_pt_*_BDTEfficiency.root"
    rdf = ROOT.RDataFrame("df", file_path_rdf)

    HypertritonRestFrameBoost(rdf, "DataframeDatawCosTheta")


    file_path_rdfMC =  "DataframeMC_*_pt_*_BDTEfficiency.root"
    rdf = ROOT.RDataFrame("df", file_path_rdfMC)
    HypertritonRestFrameBoost(rdf, "DataframeMCwCosTheta")




