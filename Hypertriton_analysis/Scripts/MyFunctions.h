#include "Math/Vector4D.h"
#include "Math/LorentzVector.h"
#include "Math/Boost.h"
#include "Math/Vector3D.h"

// Perform Lorentz boost to center-of-mass frame
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> BoostToCMHyp(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& vHyp, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& v){
    // Get the beta vector for the boost to the center-of-mass frame
    ROOT::Math::XYZVector boostVector = vHyp.BoostToCM();

    // Create a ROOT::Math::Boost object using the beta vector
    ROOT::Math::Boost boost(boostVector);

    // Perform the Lorentz boost using the created Boost object
    return boost(v);
}

float GetCosThetaBeam(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& v){
    float px = v.X();
    float py = v.Y();
    float pz = v.Z();
    return pz/sqrt(px*px+py*py+pz*pz);
}


float GetCosTheta(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& v, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& w){
    ROOT::Math::XYZVector V(v.X(), v.Y(), v.Z());
    ROOT::Math::XYZVector W(w.X(), w.Y(), w.Z());
    double s = V.Dot(W);
    return s/sqrt(V.Mag2()*W.Mag2());
}