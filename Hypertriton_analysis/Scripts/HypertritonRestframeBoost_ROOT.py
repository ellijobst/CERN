import ROOT


cppcode = """
#include "Math/Vector4D.h"
#include "Math/LorentzVector.h"
#include "Math/Boost.h"
#include "Math/Vector3D.h"

// Perform Lorentz boost to center-of-mass frame
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> BoostToCMHyp(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& vHyp, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& v){
    // Get the beta vector for the boost to the center-of-mass frame
    ROOT::Math::DisplacementVector3D boostVector = vHyp.BoostToCM();

    // Create a ROOT::Math::Boost object using the beta vector
    ROOT::Math::Boost boost(boostVector);

    // Perform the Lorentz boost using the created Boost object
    return boost(v);
};

// Get Cos(theta*) w.r.t. Beam Axis
float GetCosThetaBeam(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& v){
    float px = v.X();
    float py = v.Y();
    float pz = v.Z();
    return pz/sqrt(px*px+py*py+pz*pz);
}

// Get Cos(theta*) w.r.t. an arbitrary Axis
float GetCosTheta(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& v, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& w){
    ROOT::Math::XYZVector V(v.X(), v.Y(), v.Z());
    ROOT::Math::XYZVector W(w.X(), w.Y(), w.Z());
    double s = V.Dot(W);
    return s/sqrt(V.Mag2()*W.Mag2());
}
"""


# get dataframes and apply preselection
# df_gen = ROOT.RDataFrame("GenTable", "./Data/SignalTable_B_20g7.root").Filter("rapidity < 0.5")
# df_reco = ROOT.RDataFrame("SignalTable", "./Data/SignalTable_B_20g7.root").Filter("Rapidity < 0.5 and PseudoRapidityHe3 < 0.8 and PseudoRapidityPion < 0.8 ")

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


    # calculate cos theta* (only for helium)
    df_gen = df_gen.Define("CosThetaWrtBeam", "GetCosThetaBeam(LorentzvectorCMSHe3)")


    df_gen = df_gen.Define("CosThetaWrtpHyp", "GetCosTheta(LorentzvectorCMSHe3, LorentzvectorLABHyp)")

    treeName = "df"
    fileName = fileName
    df_gen.Snapshot(treeName, fileName, {"CosThetaWrtBeam", "BDTEfficiency", "CosThetaWrtpHyp", "centrality", "ct", "Matter", "pt", "Rapidity", "m", "PseudoRapidityHe3"})

if __name__ == "__main__":
    ROOT.gInterpreter.ProcessLine('#include "MyFunctions.h"')

    file_path_rdf =  "DataframeData_*_pt_*_BDTEfficiency.root"
    rdf = ROOT.RDataFrame("df", file_path_rdf)

    HypertritonRestFrameBoost(rdf, "DataframeDatawCosTheta")


    file_path_rdfMC =  "DataframeMC_*_pt_*_BDTEfficiency.root"
    rdf = ROOT.RDataFrame("df", file_path_rdfMC)
    HypertritonRestFrameBoost(rdf, "DataframeMCwCosTheta")




